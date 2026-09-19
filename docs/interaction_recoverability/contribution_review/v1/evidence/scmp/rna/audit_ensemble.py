"""Audit the ensemble machinery against exact enumeration on tiny windows.

Everything here is checked against ground truth computed by enumerating the whole
family, which is the only way to know that a sampler is unbiased and a bracket is sound
rather than merely self-consistent.

Checks:
  * marginals normalise, and for the validated grammar match enumeration exactly
  * every sample is a valid member of the family
  * the sampled E[g] agrees with the exact E[g] within its recorded standard error
  * the deterministic mean bracket from the affine envelopes contains the exact E[g]
  * yhat = b + a^T mu + gamma E[g] equals the exact ensemble mean of f
  * no probability floor deletes anything

Run:  python -m scmp.rna.audit_ensemble --config configs/rna_tiny.yaml
"""
from __future__ import annotations

import argparse, hashlib, itertools, json, math, os, random, sys, time

import numpy as np
import torch
import yaml

from ..model import CertifiedModel, backbone_adjacency
from ..oracles import RNANonCrossingOracle
from ..oracles.bruteforce import enumerate_rna
from .ensemble import (UNVALIDATED, VALIDATED, DeclaredGibbsPrior, ViennaPrior,
                       ensemble_bracket, expected_residual, mean_bracket,
                       predict_ensemble)
from .features import build_covariates, per_node_covariates


class Audit:
    def __init__(self):
        self.checks = 0
        self.failures = []

    def check(self, ok, label, **kw):
        self.checks += 1
        if not ok:
            self.failures.append({"check": label, **kw})
        return ok


def exact_distribution(prior, members):
    """Exact probability of every member under the declared Gibbs prior."""
    coef = prior.coefficients()
    logw = np.array([sum(coef.get(e, 0.0) for e in m) for m in members])
    w = np.exp(logw - logw.max())
    return w / w.sum()


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m scmp.rna.audit_ensemble")
    ap.add_argument("--config", default="configs/rna_tiny.yaml")
    args = ap.parse_args(argv)
    cfg = yaml.safe_load(open(args.config))
    fam = cfg["family"]
    a = Audit()
    rng = random.Random(cfg["audit"]["seed"])
    rows = []

    print(f"RNA ensemble audit: {cfg['audit']['name']}")
    print(f"probability floor  : {fam['probability_floor']} "
          f"(declared={fam['declare_floor_explicitly']}) -- no hidden deletion")
    print(f"chemistry vocab    : {cfg['covariates']['chemistry_vocabulary']}")
    print(f"assay vocab        : {cfg['covariates']['assay_vocabulary']}\n")

    for w in cfg["windows"]:
        seq, n = w["sequence"], w["length"]
        o = RNANonCrossingOracle(seq, fam["min_loop"], fam["canonical_only"])
        ground = o.ground_set()
        members = [tuple(sorted(s)) for s in
                   enumerate_rna(seq, fam["min_loop"], fam["canonical_only"])]

        # no hidden deletion: the ground set must be every canonical pair
        expect = [(i, j) for i in range(n) for j in range(i + fam["min_loop"] + 1, n)
                  if (seq[i], seq[j]) in
                  {("A", "U"), ("U", "A"), ("G", "C"), ("C", "G"), ("G", "U"), ("U", "G")}]
        a.check(list(ground) == expect, "no_probability_floor_deletion",
                window=w["name"], ground=len(ground), expected=len(expect))

        X_flat = build_covariates("GCGCAAAAGCGCAAAAGCGCA"[:21], "unmodified",
                                  "reporter_yfp",
                                  cfg["covariates"]["chemistry_vocabulary"],
                                  cfg["covariates"]["assay_vocabulary"])
        Xn = torch.tensor(per_node_covariates(X_flat, n), dtype=torch.float64)
        B = backbone_adjacency(n)

        for pspec in cfg["priors"]:
            if pspec["kind"] == "declared_gibbs":
                theta = {e: pspec["theta_scale"] * ((k % 5) - 2) / 2.0
                         for k, e in enumerate(ground)}
                prior = DeclaredGibbsPrior(o, theta)
                want_status = VALIDATED
            else:
                prior = ViennaPrior(seq, o)
                want_status = UNVALIDATED

            M = prior.marginals()
            a.check(M.status == want_status, "marginal_status_declared",
                    window=w["name"], prior=prior.name, got=M.status)
            a.check(M.normalisation_ok, "marginals_in_unit_interval",
                    window=w["name"], prior=prior.name)

            if pspec["kind"] == "declared_gibbs":
                p = exact_distribution(prior, members)
                mu_exact = np.array([sum(pi for pi, m in zip(p, members) if e in m)
                                     for e in ground])
                err = float(np.max(np.abs(M.mu - mu_exact)))
                a.check(err < 1e-9, "exact_marginals_match_enumeration",
                        window=w["name"], max_abs_err=err)
            else:
                mu_exact = None

            for gamma in cfg["model"]["gamma"]:
                for seed in cfg["model"]["seeds"]:
                    m = CertifiedModel(Xn.shape[1], ground, n,
                                       d_hid=cfg["model"]["d_hid"], gamma=gamma,
                                       seed=seed)
                    g = torch.Generator().manual_seed(seed + 7)
                    with torch.no_grad():
                        for pp in m.parameters():
                            pp.add_(torch.randn(pp.shape, generator=g, dtype=pp.dtype))

                    # exact E[g] and exact ensemble mean of f, by enumeration
                    if pspec["kind"] == "declared_gibbs":
                        with torch.no_grad():
                            gs = np.array([float(m.g(Xn, m.adjacency(m.indicator(s), B)))
                                           for s in members]) if m.residual_active \
                                else np.zeros(len(members))
                            fs = np.array([float(m(Xn, m.indicator(s), B))
                                           for s in members])
                        eg_exact = float(p @ gs)
                        ef_exact = float(p @ fs)
                    else:
                        eg_exact = ef_exact = None

                    br = mean_bracket(m, Xn, B, M.mu)
                    for ns in cfg["sampling"]["n_samples"]:
                        E = expected_residual(m, Xn, B, prior, ns, rng)
                        a.check(E.n_valid == E.n_samples or not m.residual_active,
                                "every_sample_is_a_valid_member", window=w["name"],
                                prior=prior.name, valid=E.n_valid, n=E.n_samples)
                        if eg_exact is not None and m.residual_active:
                            a.check(br[0] - 1e-9 <= eg_exact <= br[1] + 1e-9,
                                    "deterministic_bracket_contains_exact_mean",
                                    window=w["name"], bracket=br, exact=eg_exact)
                            tol = 6 * (E.stderr if E.stderr == E.stderr else 0) + 1e-9
                            a.check(abs(E.mean - eg_exact) <= tol,
                                    "sampled_mean_within_recorded_error",
                                    window=w["name"], n=ns, got=E.mean,
                                    exact=eg_exact, stderr=E.stderr)
                            yhat = predict_ensemble(m, Xn, M.mu, eg_exact)
                            a.check(abs(yhat - ef_exact) < 1e-9,
                                    "yhat_equals_exact_ensemble_mean_of_f",
                                    window=w["name"], yhat=yhat, exact=ef_exact)
                        rows.append({"window": w["name"], "prior": prior.name,
                                     "grammar": M.grammar, "status": M.status,
                                     "gamma": gamma, "seed": seed, "n_samples": ns,
                                     "Eg_sampled": E.mean, "stderr": E.stderr,
                                     "Eg_exact": eg_exact, "bracket": br,
                                     "bracket_width": br[1] - br[0],
                                     "n_valid": E.n_valid})
        print(f"  {w['name']:>5}  n={n:3d}  |E|={len(ground):3d}  |F|={len(members):5d}"
              f"  checks so far {a.checks}")

    print(f"\nchecks: {a.checks}   failures: {len(a.failures)}")
    for f in a.failures[:10]:
        print(f"  FAIL {f}")

    ex = [r for r in rows if r["Eg_exact"] is not None and r["stderr"] == r["stderr"]]
    if ex:
        errs = [abs(r["Eg_sampled"] - r["Eg_exact"]) for r in ex]
        ses = [r["stderr"] for r in ex]
        widths = [r["bracket_width"] for r in ex]
        print(f"\nsampled vs exact E[g] : median |err| {np.median(errs):.4f}, "
              f"median recorded SE {np.median(ses):.4f}")
        print(f"deterministic bracket : median width {np.median(widths):.4f} "
              f"(no sampling error)")
    print("\nthermodynamic scope: only the native target window is folded; guide "
          "chemistry enters X and is never handed to the energy model")

    out = cfg["audit"]["out_dir"]
    os.makedirs(out, exist_ok=True)
    path = os.path.join(out, f"ensemble-{time.strftime('%Y%m%d-%H%M%S')}.json")
    payload = json.dumps({"config": args.config, "checks": a.checks,
                          "failures": a.failures, "rows": rows,
                          "numerical_status": "float_unverified"},
                         indent=2, default=str)
    open(path, "w").write(payload)
    open(path + ".sha256", "w").write(hashlib.sha256(payload.encode()).hexdigest() + "\n")
    print(f"written: {path}")
    return 1 if a.failures else 0


if __name__ == "__main__":
    sys.exit(main())
