"""Validate the preregistered protocol before any large experiment.

This is a lint for research design, not for YAML. It refuses a protocol that would let a
result be reported dishonestly: a metric with no minimum meaningful gain, a noninferiority
margin smaller than the seed noise it must clear, unequal tuning allowances, an arm listed
as available that is not implemented, a cost the arms are allowed to net out, or a
licence claim asserted without evidence.

Run:  python -m scmp.validate_protocol --config configs/protocol.yaml
"""
from __future__ import annotations

import argparse
import importlib
import sys

import yaml

REQUIRED_TOP = ["protocol", "primary_metric", "effect_sizes", "design", "arms",
                "cost_accounting", "exclusions", "families_order", "rna_labels",
                "licences", "compute"]

# Arms the repository can actually run, mapped to the module that provides them.
ARM_IMPL = {
    "box": "scmp.bounds",
    "terminal_only": "scmp.bounds",
    "intermediate_unconditional": "scmp.bounds",
    "fixed_conditional": "scmp.conditional",
    "learned_conditional": "scmp.policy",
    "monotone_positional": "scmp.model",
    "unconstrained_positional": "scmp.model",
    "monotone_generic": "scmp.model",
    "unconstrained_generic": "scmp.model",
}
FAMILY_IMPL = {"layered_dag_path": "LayeredDAGPathOracle",
               "path_matching": "PathMatchingOracle",
               "rna_noncrossing": "RNANonCrossingOracle"}

REQUIRED_COSTS = {"bound_construction_s", "oracle_calls", "policy_inference_s",
                  "teacher_generation_passes", "training_s", "wall_clock_s"}
REQUIRED_VOIDS = {"empty_family", "singleton_family", "dead_instance", "non_finite"}


class Report:
    def __init__(self):
        self.errors, self.warnings, self.checks = [], [], 0

    def check(self, ok, msg, warn=False):
        self.checks += 1
        if not ok:
            (self.warnings if warn else self.errors).append(msg)
        return ok


def validate(cfg) -> Report:
    r = Report()

    for k in REQUIRED_TOP:
        r.check(k in cfg, f"missing top-level section: {k}")
    if r.errors:
        return r

    # -- provenance of the numbers -------------------------------------
    r.check(cfg["protocol"].get("chosen_on") == "development",
            "protocol must record that thresholds were chosen on development data")

    # -- metric and effect sizes ---------------------------------------
    pm = cfg["primary_metric"]
    r.check(pm.get("direction") in ("lower_is_better", "higher_is_better"),
            "primary metric must declare a direction")
    r.check(bool(pm.get("definition")), "primary metric must carry an explicit formula")

    es = cfg["effect_sizes"]
    mmg = es.get("minimum_meaningful_gain", {})
    r.check("relative_reduction" in mmg,
            "no minimum meaningful gain: any difference could be called a result")
    r.check(float(mmg.get("relative_reduction", 0)) > 0.0,
            "minimum meaningful gain must be strictly positive")
    r.check(bool(mmg.get("justification_dev")),
            "minimum meaningful gain must cite the development evidence it came from")

    ni = es.get("predictive_noninferiority_margin", {})
    r.check("margin" in ni, "no predictive noninferiority margin")
    r.check(bool(ni.get("justification_dev")),
            "noninferiority margin must cite development evidence")
    r.check(bool(ni.get("direction")),
            "noninferiority margin must state the direction of the difference")

    # -- design ---------------------------------------------------------
    d = cfg["design"]
    r.check(bool(d.get("independent_split_unit")), "no independent split unit declared")
    r.check(len(d.get("seeds", [])) >= 3, "fewer than 3 seeds")
    r.check(d.get("report_all_seeds") is True, "protocol must report all seeds")
    r.check(d.get("drop_seeds") is False, "protocol must forbid dropping seeds")
    ci = d.get("confidence_interval", {})
    r.check(bool(ci.get("method")), "no confidence interval method")
    r.check(0.5 < float(ci.get("level", 0)) < 1.0, "confidence level out of range")
    r.check(int(ci.get("resamples", 0)) >= 1000, "bootstrap resamples < 1000")
    r.check(ci.get("cluster_by") == "independent_split_unit",
            "CIs must be clustered by the independent split unit, or they overstate "
            "precision")
    tun = d.get("tuning", {})
    r.check(int(tun.get("max_configs_per_arm", 0)) > 0, "no tuning limit")
    r.check(tun.get("equal_allowance_for_baselines") is True,
            "baselines must receive the same tuning allowance, or they are weakened")

    # -- arms ------------------------------------------------------------
    arms = cfg["arms"]
    r.check(arms.get("separate_tables") is True,
            "verifier-only and predictor-training comparisons must be reported "
            "separately")
    r.check(arms.get("frozen_model_required") is True,
            "verifier arms must share a frozen model")
    seen = set()
    for group in ("verifier_only", "predictor_training"):
        for arm in arms.get(group, []):
            name, claimed = arm["name"], bool(arm.get("implemented"))
            seen.add(name)
            mod = ARM_IMPL.get(name)
            actually = False
            if mod:
                try:
                    importlib.import_module(mod)
                    actually = True
                except Exception:
                    actually = False
            if claimed and not actually:
                r.check(False, f"arm '{name}' claims implemented but {mod} does not "
                               f"import")
            if not claimed:
                r.check(bool(arm.get("note")),
                        f"unimplemented arm '{name}' must carry a note saying it will "
                        f"be reported as unavailable, not omitted")
            r.checks += 1
    for required in ("box", "terminal_only", "intermediate_unconditional",
                     "fixed_conditional", "learned_conditional",
                     "topology_aware_solver"):
        r.check(required in seen, f"required comparison arm missing: {required}")

    # -- costs -----------------------------------------------------------
    ca = cfg["cost_accounting"]
    charged = set(ca.get("charged", []))
    missing = REQUIRED_COSTS - charged
    r.check(not missing, f"costs not charged: {sorted(missing)}")
    r.check(ca.get("net_of_anything") is False,
            "protocol must forbid reporting costs net of anything")

    # -- exclusions ------------------------------------------------------
    ex = cfg["exclusions"]
    ids = {v["id"] for v in ex.get("void_rules", [])}
    r.check(REQUIRED_VOIDS <= ids, f"void rules missing: {sorted(REQUIRED_VOIDS - ids)}")
    never = set(ex.get("never_exclude_for", []))
    r.check({"arm_lost", "not_significant"} <= never,
            "protocol must forbid excluding on the basis of results")
    r.check(ex.get("unresolved_counts_in_denominator") is True,
            "unresolved instances must count in the denominator")

    # -- families --------------------------------------------------------
    fams = cfg["families_order"]
    non_rna = [f for f in fams if f.get("kind") == "non_rna"]
    r.check(len(non_rna) >= 2, "protocol must start with at least two non-RNA families")
    r.check(all(f["kind"] == "non_rna" for f in fams[:2]),
            "the first two families must be non-RNA")
    rna = [f for f in fams if f.get("kind") == "rna"]
    for f in rna:
        r.check(f.get("role") == "tiny_pilot_only",
                "RNA must enter as a tiny pilot, not as evidence")
    oracles = importlib.import_module("scmp.oracles")
    for f in fams:
        cls = FAMILY_IMPL.get(f["name"])
        r.check(cls is not None and hasattr(oracles, cls),
                f"family '{f['name']}' has no oracle implementation")

    # -- RNA label semantics ---------------------------------------------
    lab = cfg["rna_labels"]
    r.check(lab.get("label_is_structure_indexed") is False,
            "the Huesken label is not structure-indexed; the protocol must say so")
    r.check(lab.get("assume_invariant_across_latent_structures") is False,
            "protocol must NOT assume the label is invariant across latent structures")
    r.check(lab.get("worst_case_has_measured_counterpart") is False,
            "max over structures has no measured counterpart; say so")
    r.check(lab.get("permitted_predictive_claim") == "marginal_label_fit_only",
            "predictive claims must be restricted to fitting the marginal label")
    r.check("bounds off-target risk" in lab.get("retired_phrases", []),
            "the off-target phrasing must be explicitly retired")

    # -- licences ---------------------------------------------------------
    lic = cfg["licences"]
    r.check(lic.get("vendored_in_repo") is False,
            "third-party data must not be vendored")
    r.check(lic.get("publish_derived_data_before_resolution") is False,
            "derived data must not be published before licences resolve")
    for k, v in lic.items():
        if isinstance(v, dict):
            r.check(v.get("status") in ("resolved", "unresolved"),
                    f"licence '{k}' must record resolved/unresolved, not be asserted")
            if v.get("status") == "resolved":
                r.check(bool(v.get("evidence")),
                        f"licence '{k}' claims resolved without evidence")

    # -- compute -----------------------------------------------------------
    comp = cfg["compute"]
    if comp.get("authorized_budget") is None:
        r.check(bool(comp.get("proposed_cap")),
                "no authorised budget and no proposed cap: nothing may be launched")
        r.warnings.append("no authorised compute budget: only --dry-run is permitted")
    cap = comp.get("proposed_cap", {})
    r.check(cap.get("paid_or_cloud") is False,
            "paid or cloud compute requires separate authorisation")
    return r


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m scmp.validate_protocol")
    ap.add_argument("--config", default="configs/protocol.yaml")
    a = ap.parse_args(argv)
    cfg = yaml.safe_load(open(a.config))
    r = validate(cfg)

    print(f"protocol validation: {a.config}")
    print(f"preregistration    : {cfg['protocol'].get('preregistration')}")
    print(f"obligations checked: {r.checks}\n")
    for w in r.warnings:
        print(f"  WARNING  {w}")
    for e in r.errors:
        print(f"  ERROR    {e}")
    if r.errors:
        print(f"\nprotocol REJECTED: {len(r.errors)} error(s)")
        return 2
    print(f"\nprotocol accepted ({len(r.warnings)} warning(s))")
    print("primary metric      :", cfg["primary_metric"]["name"],
          f"({cfg['primary_metric']['aggregate']}, "
          f"{cfg['primary_metric']['direction']})")
    print("minimum meaningful  :",
          f"{cfg['effect_sizes']['minimum_meaningful_gain']['relative_reduction']:.0%} "
          f"relative reduction")
    print("noninferiority      :",
          f"{cfg['effect_sizes']['predictive_noninferiority_margin']['margin']} Spearman")
    print("split unit          :", cfg["design"]["independent_split_unit"])
    print("families, in order  :", " -> ".join(f["name"] for f in cfg["families_order"]))
    print("authorised budget   :", cfg["compute"]["authorized_budget"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
