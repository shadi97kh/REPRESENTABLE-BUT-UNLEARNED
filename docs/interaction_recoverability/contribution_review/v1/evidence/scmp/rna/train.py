"""RNA on-target efficacy fit against ENSEMBLE readouts. Dry run only until authorised.

The objective is yhat = b + a^T mu + gamma * E[g(X, z)], compared with the measured
ensemble readout. No latent structure is ever paired with the readout.

REFUSES to run without --dry-run unless `pilot.authorized_budget` is set. The dry run
also refuses to proceed silently past the construct question: if the dataset's assayed
context is a reporter and substitution has not been explicitly permitted, the plan is
reported as BLOCKED rather than quietly folding native windows.
"""
from __future__ import annotations

import argparse, json, os, sys, time
import yaml

from ..oracles import RNANonCrossingOracle
from .features import ContextSubstitutionRefused


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m scmp.rna.train")
    ap.add_argument("--config", default="configs/rna_pilot.yaml")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args(argv)
    cfg = yaml.safe_load(open(a.config))
    p = cfg["pilot"]
    ws = cfg["window_sensitivity"]
    blockers = []

    print(f"RNA train: {a.config}")
    print(f"task            : {cfg['task']}")
    print(f"objective       : {cfg['fit']['objective']}  "
          f"yhat = b + a^T mu + gamma E[g]")
    print(f"per-structure labels: {cfg['fit']['label_per_latent_structure']} "
          f"(must be false)\n")

    if cfg["fit"]["label_per_latent_structure"]:
        blockers.append("config would label latent structures with the ensemble readout")

    ds = cfg["dataset"]
    print(f"dataset         : {ds['name']}")
    print(f"assayed context : {ds['context_source']}")
    print(f"substitution    : allowed={ds['allow_context_substitution']}")
    if ds["context_source"] == "reporter_construct" and \
            not ds["allow_context_substitution"]:
        blockers.append(
            f"{ds['name']}: assayed context is a reporter construct and native "
            f"substitution is not permitted, so no target window can be built for it")

    print(f"\nwindow sensitivity: lengths {ws['lengths']} nt, min_loop "
          f"{ws['min_loop']}, canonical_only {ws['canonical_only']}, "
          f"floor {ws['probability_floor']}")
    if ws["probability_floor"] != 0.0:
        blockers.append("probability floor is non-zero: pairs would be deleted")
    print(f"isoform policy    : {ws['isoform_policy']}  "
          f"fold_exon_unions={ws['fold_exon_unions']}")
    if ws["fold_exon_unions"]:
        blockers.append("config would fold exon unions")

    pr = cfg["prior"]
    print(f"\nprior             : {pr['kind']}")
    print(f"exact marginals only for a validated grammar: "
          f"{pr['exact_marginals_require_validated_grammar']}")

    # cost of the plan, MEASURED on one representative window per length with the
    # sampler that would actually be used
    import random as _random, time as _time
    from .ensemble import DeclaredGibbsPrior
    n_draws = cfg["fit"]["n_samples"]
    print(f"\n{'length':>7} {'|E|':>7} {'lattice':>10} {'sample s':>10}  note")
    est = []
    for L in ws["lengths"]:
        seq = ("GCGCAAAAGCGCAUAUAUAU" * 20)[:L]
        o = RNANonCrossingOracle(seq, ws["min_loop"], ws["canonical_only"])
        k = len(o.ground_set())
        prior = DeclaredGibbsPrior(o, {e: 0.0 for e in o.ground_set()})
        t = _time.time()
        draws = prior.sample(n_draws, _random.Random(0))
        dt = _time.time() - t
        valid = all(o.feasibility(forced=s).feasible for s in draws[:8])
        est.append({"length": L, "candidates": k, "sample_seconds": dt,
                    "n_draws": n_draws, "spot_checked_valid": valid})
        note = "tractable (stochastic backtrace)" if dt < 5.0 else "too slow"
        print(f"{L:7d} {k:7d} {'10^%.0f' % (k * 0.30103):>10} {dt:10.2f}  {note}")
        if dt >= 5.0:
            blockers.append(f"length {L}: {n_draws} draws take {dt:.1f}s")
        if not valid:
            blockers.append(f"length {L}: sampler produced an infeasible draw")

    print(f"\nsplit unit       : {cfg['splits']['unit']}")
    print(f"seeds            : {cfg['fit']['seeds']}   "
          f"n_samples {cfg['fit']['n_samples']}")
    print(f"authorised budget: {p['authorized_budget']}")

    if not a.dry_run and p["authorized_budget"] is None:
        print("\nREFUSING TO RUN: no authorised compute budget.")
        return 3

    print(f"\nblockers ({len(blockers)}):")
    for b in blockers:
        print(f"  - {b}")

    out = p.get("out_dir", "runs/rna_pilot")
    os.makedirs(out, exist_ok=True)
    path = os.path.join(out, f"dryrun-{time.strftime('%Y%m%d-%H%M%S')}.json")
    json.dump({"config": cfg, "window_cost": est, "blockers": blockers,
               "dry_run": a.dry_run}, open(path, "w"), indent=2)
    print(f"\nplan written: {path}")
    if a.dry_run:
        print("DRY RUN ONLY. Nothing was launched, no window was folded, no fit was run.")
    return 1 if blockers else 0


if __name__ == "__main__":
    sys.exit(main())
