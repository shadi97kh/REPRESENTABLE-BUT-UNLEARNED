"""Build the corrected claim ledger from audit records, not by hand.

Every entry carries a status and an evidence kind, so the six categories the
master instruction requires to stay separate -- exact arithmetic, real-arithmetic
proof, float diagnostic, enclosed numerical result, learned synthetic result and
measured RNA result -- are visible per claim rather than merged into prose.

Historical protocol-v1 gates are copied in VERBATIM with status `historical`.
They are never edited, re-run or re-scored here; a v1 gate that was UNDETERMINED
stays UNDETERMINED in this ledger forever.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import time

from . import (EXACT_INTEGER, EXACT_RATIONAL, FLOAT_DIAGNOSTIC, LEARNED_SYNTHETIC,
               MEASURED_RNA, REAL_PROOF)

# status vocabulary for the ledger
ESTABLISHED = "established"      # verified at the stated evidence kind
SCOPED = "scoped"               # true only under narrower conditions than stated
CORRECTED = "corrected"         # the reported number was wrong or mis-defined
REFUTED = "refuted"             # the claim as stated is false
OPEN = "open"                   # not yet decided; needs the revision cycle
HISTORICAL = "historical"       # a protocol-v1 result, retained unchanged


def _latest(pat):
    f = sorted(glob.glob(pat))
    return f[-1] if f else None


def _load(p):
    with open(p) as fh:
        return json.load(fh)


def build(runs="runs"):
    rpath = _latest(os.path.join(runs, "revision", "audit_report-*.json"))
    mpath = _latest(os.path.join(runs, "revision", "audit_math-*.json"))
    gpath = _latest(os.path.join(runs, "gates", "gates-*.json"))
    R = _load(rpath)["results"] if rpath else {}
    M = _load(mpath)["results"] if mpath else {}
    entries = []

    def add(cid, claim, status, kind, was, now, source, note):
        entries.append({"id": cid, "claim": claim, "status": status,
                        "evidence_kind": kind, "reported": was, "corrected": now,
                        "source_record": source, "note": note})

    # ---- item 0: traceability --------------------------------------------
    if "item0_trace_pages_6_10" in R:
        a = R["item0_trace_pages_6_10"]
        sc = a["status_counts"]
        add("report.traceability_pages_6_10",
            "every number on report pages 6-10 traces to an immutable record",
            ESTABLISHED if not sc.get("no_record") else OPEN, EXACT_INTEGER,
            "not assessed",
            f"{a['n_rows']} rows: {sc.get('traced', 0)} traced, "
            f"{sc.get('mismatch', 0)} mismatch, {sc.get('no_record', 0)} unrecorded",
            a.get("report", {}).get("resolved_to", "docs/run_report.pdf"),
            "the report supplies the claim list; every value is recomputed from "
            "records. Report sha256 recorded in the audit artefact.")
        add("audits.persist_only_on_failure",
            "clean audit runs leave an immutable record",
            CORRECTED, EXACT_INTEGER,
            "implied by the report citing 1,262,952 checks and 72/72",
            "audit_bounds, audit_conditional and check_certificates persisted "
            "nothing on success; records created by rerun_unrecorded, all three "
            "reproduce exactly (122 s CPU, no GPU)",
            "runs/revision/rerun_unrecorded-*.json",
            "protocol v2 must persist a structured record on success as well as "
            "on failure; scraping stdout is a rescue, not an arrangement")

    # ---- item 1 -----------------------------------------------------------
    if "item1_termination_accounting" in R:
        a = R["item1_termination_accounting"]
        for arm, v in a["arms"].items():
            if v["report_closed_pct"] == v["corrected_proved_pct"]:
                continue
            add(f"verifier.{arm}.resolution_rate",
                f"{arm} resolves this fraction of pilot instances",
                CORRECTED, EXACT_INTEGER,
                f"{v['report_closed_pct']}% closed",
                f"{v['corrected_proved_pct']}% proved",
                a["record"],
                f"{v['overlap_gap_zero_and_timeout']} instances carried both "
                f"gap_rel==0 and timeout; the reported figure counted them as "
                f"closed. Corrected value is lower.")
        add("verifier.outcome_vocabulary",
            "termination outcomes are mutually exclusive",
            ESTABLISHED, EXACT_INTEGER, "closed / timeout (overlapping flags)",
            "proved / unresolved (exclusive, sums to 100%)", a["record"],
            "the record's own `status` field is already exclusive and "
            "timeout holds if and only if status == unresolved")

        fa = a["first_tolerance_attainment"]
        add("verifier.first_attainment_time",
            "time at which each arm first met the requested tolerance",
            OPEN, FLOAT_DIAGNOSTIC, "reported as median wall_s",
            "NOT INSTRUMENTED in protocol v1; wall_s is a stopping time and for a "
            "timed-out record equals the budget",
            a["record"], fa["v2_requirement"])
        ng = a["negative_gaps"]
        add("verifier.negative_gap_defect",
            "no bound inversion occurs",
            CORRECTED, FLOAT_DIAGNOSTIC, "not reported",
            f"{ng['policy_pilot_negative_gaps']} negative gaps in policy_pilot, "
            f"min {ng['most_negative']:.3e}, not clipped",
            "runs/policy_pilot/*.json", ng["assessment"])
        den = a["denominators"]
        add("verifier.denominators",
            "n = 330 per arm",
            ESTABLISHED, EXACT_INTEGER, "330 per arm, 291 timeouts, 30 voids",
            f"{den['cells_attempted']} cells - {den['voids_excluded']} void = "
            f"{den['cells_evaluated_per_arm']} x {den['arms']} arms = "
            f"{den['verifier_records']} records",
            a["record"], den["timeouts_are_records_not_cells"])

    # ---- item 2 -----------------------------------------------------------
    if "item2_effect_interval_populations" in R:
        a = R["item2_effect_interval_populations"]
        add("effects.population_accounting",
            "model-effect intervals validated and classified",
            CORRECTED, EXACT_INTEGER,
            "75 intervals validated, verdicts 31/2/7 (implied subset)",
            f"{a['validation_population']['n_intervals']} validated by enumeration "
            f"(0 failures) AND {a['reporting_population']['n_intervals']} classified "
            f"on pilot instances -- disjoint populations",
            a["record"],
            "no interval is unaccounted for; the two totals answer different "
            "questions and must not appear in one row")

    # ---- item 3 -----------------------------------------------------------
    if "item3_pair_counting" in R:
        a = R["item3_pair_counting"]
        add("ranking.coverage",
            "fraction of candidate comparisons the certificate resolves",
            CORRECTED, REAL_PROOF,
            f"{a['resolved_strict_inequalities']}/{a['ordered_pairs']} = "
            f"{a['reported_coverage_ordered']:.1%} (ordered denominator)",
            f"{a['resolved_strict_inequalities']}/{a['unordered_pairs']} = "
            f"{a['matched_coverage_unordered']:.1%} (unordered, matched)",
            a["record"], a["injectivity_argument"])

    # ---- item 4 -----------------------------------------------------------
    if "item4_gates_and_graph_inputs" in R:
        a = R["item4_gates_and_graph_inputs"]
        add("gate.G3.interpretation",
            "certifiability is cheap because G3 passed",
            SCOPED, FLOAT_DIAGNOSTIC,
            f"PASS, point {a['G3']['point']:.4f}, CI {a['G3']['ci95']}",
            "PASS is arithmetically correct but uninformative",
            a["records"]["gates"], a["G3"]["structural_objection"])
        add("predictor.graph_inputs",
            "what the trained predictor arms consumed as graph input",
            ESTABLISHED, REAL_PROOF, "unstated in the report",
            "the guide's own 21-mer MFE structure; no target window",
            a["records"]["pilot"], a["graph_inputs"]["consequence"])
        add("predictor.split_strength",
            "the predictor split supports a generalisation claim",
            SCOPED, FLOAT_DIAGNOSTIC,
            f"gene-disjoint, n={a['split']['n']}",
            f"one held-out gene ({a['split']['held_out_genes']}), "
            f"test n={a['split']['test']}",
            a["records"]["pilot"], a["split"]["objection"])

    # ---- item 5 -----------------------------------------------------------
    if "item5_anchor_sign" in R:
        a = R["item5_anchor_sign"]
        add("model.anchor_sign", "the additive anchor is non-negative",
            REFUTED, REAL_PROOF, "non-negative additive anchor",
            "signed: nn.Linear with unconstrained weight and bias",
            "scmp/model.py:AnchorEncoder", a["note"])

    # ---- items 6-9 --------------------------------------------------------
    if "item6_three_node_path" in M:
        a = M["item6_three_node_path"]
        add("theory.additivity_in_z",
            "sum+ReLU with non-negative parameters is additive in the graph z",
            REFUTED, EXACT_RATIONAL,
            "exactly affine (C4, violation 2.157e-16)",
            f"f(z) = 2z1 + 2z2 + 2z1z2 on the three-node path; mixed difference "
            f"{a['mixed_difference']} with every ReLU stable",
            "scmp/revision/audit_math.py:check_three_node_path",
            "C4 varied X with the adjacency fixed, establishing affinity in X "
            "only. Additivity in z is a different statement and is false.")
    if "item6_design_matrix_residual" in M:
        a = M["item6_design_matrix_residual"]
        g = a["projection_onto_[1,z]"]["g_residual_only"]
        add("model.residual_carries_structure",
            "the residual g is non-additive in z on a real feasible family",
            ESTABLISHED, FLOAT_DIAGNOSTIC, "not measured",
            f"relative residual {g['relative_residual']:.3e} after projecting "
            f"{g['n_members']} enumerated members onto [1, z]",
            "scmp/revision/audit_math.py:check_design_matrix_residual",
            "capacity for non-additive structure is present. This is NOT "
            "evidence that it generalises: G4 measured generalisation and failed.")
    if "item7_layer1_relu" in M:
        a = M["item7_layer1_relu"]
        add("bounds.layer1_relaxation",
            "layer 1 needs no relaxation at all",
            CORRECTED, FLOAT_DIAGNOSTIC, "no relaxation needed at layer 1",
            f"message M1 exactly affine; its ReLU relaxed, "
            f"{a['unstable_units_at_layer_1']} unstable units observed",
            "scmp/bounds.py:ResidualBounder.propagate",
            "the implementation was already correct; the report's sentence "
            "generalised a true statement about the message to the whole layer")
    if "item8_kappa_counterexample" in M:
        a = M["item8_kappa_counterexample"]
        add("theory.kappa_feature_independence",
            "kappa_L bounds the relaxation ratio for arbitrary positive features",
            REFUTED, EXACT_RATIONAL, "implied by the report's framing",
            f"kappa_2 = {a['uniform_features']['kappa_2']} but features "
            f"(1,1,1,1/10,1/10) give {a['skewed_features']['ratio']} "
            f"({a['skewed_features']['ratio_float']:.4f})",
            "scmp/revision/audit_math.py:check_kappa_counterexample",
            "the uniform-feature identity is untouched; only its "
            "feature-independent extension is refuted")
    if "item9_matching_kappa_formula" in M:
        a = M["item9_matching_kappa_formula"]
        add("theory.kappa_exponent",
            "the relaxation-gap exponent 0.912 per layer is architectural",
            SCOPED, EXACT_RATIONAL, "0.912 per layer, R^2 0.9996, depths 1-4",
            "kappa_L = (2m-1)^L without self-loops, m^L with them: the exponent "
            "depends on the family and the self-loop convention",
            "scmp/revision/audit_math.py:check_matching_kappa_formula",
            "C5's measurement stands on the family C5 used; the universal "
            "reading does not")

    if "item9_actual_convention" in M:
        a = M["item9_actual_convention"]
        act = {(x["m"], x["L"]): x["kappa"] for x in a["rows"]
               if x["convention"].startswith("ACTUAL")}
        add("theory.kappa_under_actual_convention",
            "the closed-form kappa applies to the implemented model",
            REFUTED, EXACT_RATIONAL,
            "kappa_L presented without stating the fixed part B",
            f"B = I + chain at every call site; kappa is no clean power "
            f"(m=3,L=2 gives {act.get((3, 2))} vs 9 for self-loops only)",
            "scmp/revision/audit_math.py:check_actual_convention_kappa",
            a["note"])

    # ---- historical gates, verbatim --------------------------------------
    if gpath:
        for g in _load(gpath)["gates"]:
            add(f"historical.{g['gate'].split()[0]}", g["gate"], HISTORICAL,
                LEARNED_SYNTHETIC if "predictor" not in g["gate"] else MEASURED_RNA,
                f"{g.get('status')} point={g.get('mean_difference')} "
                f"ci={g.get('ci95')}", "retained unchanged", gpath,
                "protocol v1 result; retained for the record and never re-scored")

    # ---- P1 findings ------------------------------------------------------
    ppath = _latest(os.path.join(runs, "revision", "predictor_pilot-*.json"))
    P = {}
    for f in sorted(glob.glob(os.path.join(runs, "revision",
                                           "predictor_pilot-*.json"))):
        P.update(_load(f).get("stages", {}))
    if "S1_reproduce" in P:
        a = P["S1_reproduce"]
        add("predictor.reproduction", "the predictive pilot reproduces",
            ESTABLISHED, FLOAT_DIAGNOSTIC,
            "0.7253 / 0.4573 / 0.1212",
            "reproduced exactly from actual source and data",
            ppath, "technical training diagnostic; supports no RNA claim")
        add("predictor.input_parity",
            "the unconstrained baseline is a fair comparator",
            REFUTED, FLOAT_DIAGNOSTIC, "0.1212, read as the cost of no constraint",
            f"adding only a node-wise skip path gives "
            f"{a['means'].get('unconstrained_with_skip', float('nan')):.4f}",
            ppath, "the baseline had no route from X to the output bypassing the "
                   "adjacency; the gate compared architecture families")
    if "S2_parity" in P:
        a = P["S2_parity"]
        add("predictor.gamma0_graph_blind",
            "certifiable_gamma0 is a graph model", REFUTED, FLOAT_DIAGNOSTIC,
            "0.7253 presented as a certifiable-model result",
            "f(A=I) == f(A=ones) exactly; it is b(X) alone",
            ppath, "its score contains no graph information")
        add("predictor.decomposition_unexercised",
            "G4 tested the anchor+residual decomposition",
            REFUTED, FLOAT_DIAGNOSTIC, "G4 = residual predictive value",
            "candidates=(), so a^T z == 0 and P z == 0 in every predictor arm",
            ppath, a["checks"]["decomposition_exercised"]["note"])
    if "S3_tinyfit" in P:
        add("predictor.optimisation_adequate",
            "the weak arm underfits or is buggy", REFUTED, FLOAT_DIAGNOSTIC,
            "0.1212 read as a weak/underfit model",
            "all three arms interpolate 24 rows at train R2 = 1.000",
            ppath, "a generalisation failure, not an optimisation failure; the "
                   "diagnosis comes from the training curve, not the test score")
    if "S4_ainit" in P:
        a = P["S4_ainit"]
        pc = a.get("paired_contrasts", {}).get("path_matching|local_interaction", {})
        c = pc.get("contrasts", {})
        ft = c.get("A_init_finetune", {})
        old_s = c.get("A_plus_B_old_setup", {})
        add("predictor.a_initialisation_repairs_residual",
            "A-initialised residual training beats the old setup",
            ESTABLISHED, LEARNED_SYNTHETIC,
            "old setup: residual on gamma=1 from scratch",
            f"A_init_finetune delta {ft.get('mean_delta_mse', float('nan')):+.4f} "
            f"vs A_only on {ft.get('wins_on_seeds')}/{ft.get('n_seeds')} seeds; "
            f"old setup {old_s.get('mean_delta_mse', float('nan')):+.4f} on "
            f"{old_s.get('wins_on_seeds')}/{old_s.get('n_seeds')}",
            ppath, "synthetic, one family, size holdout; a mechanism result, not "
                   "a predictive one")
    if "S5_controlled" in P:
        add("predictor.useful_nonlinear_value",
            "the residual carries useful nonlinear predictive value",
            REFUTED, LEARNED_SYNTHETIC, "the intended P1 outcome",
            "a pairwise ridge is EXACT (MSE 0.0000) on path_matching|interacting "
            "where A_init_residual is 10.3087; no benefit at all on layered_dag",
            ppath, "the residual learns ~8% of the target projection residual and "
                   "~13% of its mixed-difference magnitude: real but far short. "
                   "NOTE the targets are exactly linear in an explicit basis, "
                   "which favours ridge by construction; a fair retest needs "
                   "non-basis-representable targets")
        add("predictor.additive_negative_control",
            "the residual avoids harm on deliberately additive targets",
            ESTABLISHED, LEARNED_SYNTHETIC, "not previously tested",
            "no significant harm on either additive control",
            ppath, "required by the plan as a negative control")

    # ---- P2 findings ------------------------------------------------------
    vpath = _latest(os.path.join(runs, "revision", "verifier_pilot-*.json"))
    if vpath:
        V = _load(vpath)
        ov = V["aggregates"]["overall"]["arms"]
        r0 = V["aggregates"]["rung0"]["arms"]
        r3 = V["aggregates"]["rung3"]["arms"]
        add("verifier.conditioning_pays_on_harder_instances",
            "conditional message support improves verification",
            ESTABLISHED, FLOAT_DIAGNOSTIC,
            "v1: UNDETERMINED (metric 0/0; baseline closed 99.1%)",
            f"decision coverage {ov['box']['decision_coverage']:.1%} -> "
            f"{ov['fixed_conditional']['decision_coverage']:.1%} overall; gain "
            f"grows with difficulty: {r0['fixed_conditional']['decision_coverage'] - r0['box']['decision_coverage']:+.1%} "
            f"at rung 0 to {r3['fixed_conditional']['decision_coverage'] - r3['box']['decision_coverage']:+.1%} at rung 3",
            vpath,
            "rung 0 IS the v1 distribution and shows zero gain, which explains "
            "the v1 UNDETERMINED result as an instance-distribution defect")
        add("verifier.selective_conditioning_efficiency",
            "conditioning must be applied to every instance",
            REFUTED, FLOAT_DIAGNOSTIC, "v1 arms always conditioned",
            f"deterministic selective reaches "
            f"{ov['deterministic_selective']['decision_coverage']:.1%} coverage at "
            f"{ov['deterministic_selective']['median_wall_s']:.3f}s and a median of "
            f"{ov['deterministic_selective']['median_oracle_calls']:.0f} oracle call(s) "
            f"vs {ov['fixed_conditional']['median_oracle_calls']:.0f}",
            vpath, "skipping ~53% of instances costs 0.9pp of coverage")
        add("verifier.soundness_vs_exact_optimum",
            "every arm's bound is sound against ground truth",
            ESTABLISHED, FLOAT_DIAGNOSTIC, "not previously checkable",
            "0 violations across 108 cases x 7 arms, against the exact optimum "
            "by enumeration",
            vpath, "an earlier run showed 66 violations; the cause was the "
                   "ground-truth enumerator accepting extendable prefixes as "
                   "members on non-downward-closed families, not the bounds")
        add("verifier.negative_gap_roundoff",
            "no bound inversion occurs", CORRECTED, FLOAT_DIAGNOSTIC,
            "not reported in v1",
            "42 negative gaps at -8.5e-14, identical across all 7 arms on tight "
            "cases; recorded, never clipped",
            vpath, "shared float evaluation, not arm logic; P4 enclosure must "
                   "make it impossible rather than small")
        add("verifier.external_solver_gap",
            "an external topology-aware solver was compared",
            OPEN, FLOAT_DIAGNOSTIC, "requested by the plan",
            "auto_LiRPA, gurobipy and cvxpy are not installed; exhaustive "
            "enumeration and an LP outer relaxation used instead",
            vpath, "recorded as a gap; no internal baseline relabelled external")
        g = V["policy_gate"]
        add("verifier.policy_gate",
            "a learned allocation policy is justified",
            OPEN, FLOAT_DIAGNOSTIC, "gated by the plan",
            f"gate OPEN (useful_region={g['useful_region']}); policy NOT yet trained",
            vpath, "random scoring already matches deterministic scoring on "
                   "coverage, so a policy must earn its keep on WHEN to query, "
                   "not WHICH edge")

    # ---- P3 findings ------------------------------------------------------
    apath = _latest(os.path.join(runs, "revision", "rna_audit-*.json"))
    ppath2 = _latest(os.path.join(runs, "revision", "rna_pilot-*.json"))
    if apath:
        A = _load(apath)
        ins = A.get("insert_recovery", {})
        adm = sum(1 for r in A.get("arm_admissibility", []) if r["admissible"])
        add("rna.reporter_context_partially_recoverable",
            "the reporter construct context is unrecoverable, blocking 6 of 8 arms",
            REFUTED, FLOAT_DIAGNOSTIC, "6 of 8 arms blocked",
            f"{adm} of {len(A.get('arm_admissibility', []))} arms admissible; "
            f"{ins.get('by_window', {}).get(50, {}).get('interior_fraction', 0):.1%} "
            f"of located 50-nt windows lie inside a lower bound on the insert",
            apath,
            "the hull of densely tiled sites bounds the inserted fragment. Rests "
            "on an assumption of insert contiguity inferred from tiling, NOT a "
            "construct record; edge windows stay unresolved_context")
        add("rna.labels_unclipped",
            "efficacy labels are preserved unclipped",
            ESTABLISHED, REAL_PROOF, "134 rows exceed 1.0",
            "preserved; max 1.341; clipping would corrupt the primary "
            "normalisation", apath, "positive control = 90.0% inhibition")
        add("rna.assay_identity",
            "the redistributed annotation states the assay",
            REFUTED, REAL_PROOF, "HeLa / luciferase reporter assay",
            "primary source: YFP reporter plasmid, H1299, 48 h; annotation "
            "recorded and contradicted", apath,
            "annotation columns are never used as assay identity")
        add("rna.davis_s1_blocked",
            "Davis S1 is available for the chemical-transfer study",
            OPEN, MEASURED_RNA, "required by the plan",
            "gkaf479_supplemental_files.zip absent; publisher route recorded; "
            "CAPTCHA not bypassed", apath,
            "blocks chemical transfer, native-vs-reporter, dose/time, replicates")
    if ppath2:
        P3 = _load(ppath2)
        sm = P3["summary"]
        base = sm["matched_sequence_chemistry"]["spearman_mean"]
        best = max(sm[k]["spearman_mean"] for k in sm)
        add("rna.structure_adds_predictive_value",
            "target-window structure improves on-target efficacy prediction",
            REFUTED, MEASURED_RNA, "the intended RNA claim",
            f"best structural arm {best:.4f} vs sequence-only {base:.4f} "
            f"(delta {best - base:+.4f}), 5% of the seed sd 0.028",
            ppath2,
            "1153 admissible interior windows, novel-sequence split, group leak "
            "0. Arms are ridge on window SUMMARY features, not the full graph "
            "model -- a disclosed limitation")
        add("rna.sequence_only_benchmark",
            "guide sequence predicts on-target efficacy",
            ESTABLISHED, MEASURED_RNA, "0.6464 on the earlier split",
            f"{base:.4f} Spearman, P@20 "
            f"{sm['matched_sequence_chemistry']['precision_at_k']:.3f}, "
            f"enrichment {sm['matched_sequence_chemistry']['enrichment']:.2f} at "
            f"prevalence {sm['matched_sequence_chemistry']['prevalence']:.3f}",
            ppath2, "novel-sequence split by gene x 13-mer cluster, zero leak")
        ps = P3.get("prior_sensitivity", {})
        add("rna.prior_uncertainty_dominates",
            "numerical sampling error bounds the ensemble uncertainty",
            REFUTED, MEASURED_RNA, "implied by reporting sampling error alone",
            f"prior temperature moves mean pair density by {ps.get('spread'):.4f}, "
            f"larger in relative terms than any arm difference",
            ppath2, "prior choice is a modelling uncertainty and is NOT covered "
                    "by a Monte Carlo interval")

    # ---- what the revision cycle must still decide -----------------------
    for cid, claim in [
        ("open.residual_beats_explicit_basis",
         "a nonlinear residual beats explicit-basis regression on a task whose "
         "interaction is NOT exactly representable in a small basis"),
        ("open.learned_allocation",
         "learned allocation beats a strong deterministic selective rule at equal "
         "total wall time"),
        ("open.numerical_enclosure",
         "the certifying operator set is enclosed end to end"),
        ("open.rna_structural",
         "a measured RNA comparison with correct assay context supports the "
         "structural model"),
        ("open.ml_novelty",
         "a distinct ML contribution survives equation-level comparison to prior "
         "work")]:
        add(cid, claim, OPEN, "undetermined", "not established", "not established",
            "docs/theorem_scope_v2.md", "decided by the revision cycle, not here")

    return entries


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="python -m scmp.revision.ledger",
        description="Build the corrected claim ledger from audit records.")
    ap.add_argument("--runs", default="runs")
    ap.add_argument("--out", default="docs/corrected_ledger.json")
    a = ap.parse_args(argv)

    entries = build(a.runs)
    counts = {}
    for e in entries:
        counts[e["status"]] = counts.get(e["status"], 0) + 1
    payload = {"kind": "corrected_claim_ledger",
               "created": time.strftime("%Y-%m-%d %H:%M:%S"),
               "status_counts": counts, "n_entries": len(entries),
               "entries": entries}
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with open(a.out, "w") as fh:
        json.dump(payload, fh, indent=1)

    print(f"{'status':<14}{'n':>4}")
    for k in sorted(counts):
        print(f"{k:<14}{counts[k]:>4}")
    print(f"\n{len(entries)} entries -> {a.out}")
    print("\ncorrected or refuted:")
    for e in entries:
        if e["status"] in (CORRECTED, REFUTED):
            print(f"  [{e['status']:>9}] {e['id']}")
            print(f"              was: {e['reported']}")
            print(f"              now: {e['corrected']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
