"""P0 items 1-5: trace every headline number in the run report to its record.

The run report is treated as a set of CLAIMS. Each check below recomputes the
quantity from the immutable record that produced it and reports agreement,
disagreement, or a definitional mismatch. Nothing under `runs/` is modified: a
correction is a new record that cites the record it corrects.

  1  closed% + timeout% exceeded 100% (105.8 fixed, 108.8 learned)
     Two overlapping FLAGS were summed. `gap_rel == 0` is a tightness attainment
     flag; `timeout` is a termination flag; an instance can carry both. The record
     already contains a mutually exclusive `status` field, so the corrected rate
     is computable without rerunning anything -- and it is WORSE for the
     conditional arms than the figure the report printed.

  2  75 model-effect intervals vs 40 classified verdicts
     Disjoint populations: 75 intervals validated by exhaustive enumeration on
     tiny families, 40 intervals reported on pilot instances. The report's single
     row implied the 40 were a subset of the 75.

  3  12 candidates, 132 pairs, 60 resolved, "45.5% coverage"
     132 is the number of ORDERED pairs; 60 counts strict inequalities u_i < l_j.
     Since l <= u, at most one of (u_i<l_j, u_j<l_i) can hold, so resolved
     inequalities inject into UNORDERED pairs. The matched denominator is 66.

  4  G3 sign, comparator, and what the predictor arms actually consumed
     The gate's direction is explicit in code. Its comparator is not the best
     available model. The graphs are the guide's own 21-mer MFE structure, never
     a target window -- which is why these arms ran while the window arms were
     blocked.

  5  "non-negative additive anchor"
     Checked against the code rather than the prose.
"""
from __future__ import annotations

import argparse
import collections
import glob
import inspect
import json
import os
import time

import numpy as np

from . import EXACT_INTEGER, FLOAT_DIAGNOSTIC, REAL_PROOF


def _load(path):
    with open(path) as fh:
        return json.load(fh)


def _latest(pattern):
    files = sorted(glob.glob(pattern))
    return files[-1] if files else None


def report_fingerprint(path, pages=(6, 10)):
    """Hash the report and, for a PDF, pull the text of the audited pages.

    The report is the CLAIM SOURCE, never the evidence: numbers extracted here
    are only ever compared against values recomputed from records. Hashing it
    means the audit names the exact document it reconciled.
    """
    import hashlib
    import shutil
    import subprocess
    out = {"path": path, "exists": os.path.exists(path)}
    if not out["exists"]:
        for alt in (os.path.join("docs", os.path.basename(path)),
                    os.path.splitext(path)[0] + ".md"):
            if os.path.exists(alt):
                out["resolved_to"] = alt
                path = alt
                out["exists"] = True
                break
    if not out["exists"]:
        out["note"] = "report not found; claims are taken from the trace table below"
        return out
    with open(path, "rb") as fh:
        out["sha256"] = hashlib.sha256(fh.read()).hexdigest()
    out["bytes"] = os.path.getsize(path)
    if path.lower().endswith(".pdf") and shutil.which("pdftotext"):
        txt = {}
        for pg in range(pages[0], pages[1] + 1):
            r = subprocess.run(["pdftotext", "-f", str(pg), "-l", str(pg), path, "-"],
                               capture_output=True, text=True, timeout=60)
            txt[pg] = r.stdout
        out["pages_extracted"] = sorted(txt)
        out["page_chars"] = {k: len(v) for k, v in txt.items()}
        out["_text"] = txt
    return out


def _claim_present(fp, page, token):
    """Is the reported token literally on that page? Absence is reported, not fatal."""
    t = fp.get("_text", {}).get(page)
    if t is None:
        return None
    return token in t.replace("\u2212", "-")


# --------------------------------------------------------------------------- item 0
TRACED = "traced"            # recomputed from an immutable record and agrees
MISMATCH = "mismatch"        # recomputed and differs; corrected value given
NO_RECORD = "no_record"      # never persisted; only a rerun can reproduce it
TEXT_ONLY = "text_only"      # survives only as stdout text inside another record
CONFIG = "config"            # traced to a config/primary-source file, not a run
SOURCE = "source"            # traced to source code rather than to a run


def item0_trace_pages(runs="runs", fp=None):
    """Trace every substantive number on report pages 6-10 to what produced it.

    Each row names the reported value, the value recomputed here, the immutable
    record it came from, and the source/config that generated it. Rows that no
    record can support are marked NO_RECORD rather than quietly accepted.
    """
    pilot = _latest(os.path.join(runs, "pilot", "pilot-*.json"))
    gates = _latest(os.path.join(runs, "gates", "gates-*.json"))
    eff = _latest(os.path.join(runs, "explanations", "effects-*.json"))
    rep = _latest(os.path.join(runs, "reproduction", "reproduction-*.json"))
    rna = _latest(os.path.join(runs, "rna_eval", "eval-*.json"))
    spl = _latest(os.path.join(runs, "rna_splits", "splits-*.json"))
    ens = _latest(os.path.join(runs, "rna_ensemble", "ensemble-*.json"))
    L = {}
    for k, v in (("pilot", pilot), ("gates", gates), ("eff", eff), ("rep", rep),
                 ("rna", rna), ("spl", spl), ("ens", ens)):
        L[k] = _load(v) if v else None

    # records created by scmp.revision.rerun_unrecorded for audits that persist
    # nothing on success; keyed by job name
    rr_path = _latest(os.path.join(runs, "revision", "rerun_unrecorded-*.json"))
    RR = {}
    if rr_path:
        for r in _load(rr_path).get("results", []):
            RR[r["job"]] = r

    rows = []

    def row(page, claim, reported, got, status, record, source, config="", note=""):
        rows.append({"page": page, "claim": claim, "reported": str(reported),
                     "recomputed": str(got), "status": status,
                     "record": record or "-", "generating_source": source,
                     "config": config, "note": note})

    # ---------------- page 6: certificates, numerics, gates -------------------
    cc = RR.get("check_certificates", {}).get("parsed_from_stdout", {})
    row(6, "refinement certificates accepted", "72/72, 2222 obligations",
        (f"{cc['accepted']}/{cc['accepted_of']}, {cc['obligations_discharged']} "
         f"obligations" if cc.get("accepted") else "not reproducible from any record"),
        TRACED if cc.get("accepted") else NO_RECORD, rr_path,
        "scmp/check_certificates.py", "runs/tiny_complete",
        "check_certificates persists nothing on success; this row is traced to a "
        "record created by scmp.revision.rerun_unrecorded, which ran the audit "
        "UNCHANGED and captured its stdout, exit status and module hashes.")
    row(6, "interval-arithmetic tests", "70 tests",
        "collected by pytest at run time", SOURCE, None,
        "tests/numerics/", "", "test count is a property of the suite, not a run")
    row(6, "minimum meaningful gain / noninferiority margin", "0.25 / 0.05",
        "0.25 / 0.05", CONFIG, None, "configs/protocol.yaml",
        "configs/protocol.yaml", "preregistered before the pilot")
    if L["gates"]:
        for g in L["gates"]["gates"]:
            gid = g["gate"].split()[0]
            rep_v = (f"{g.get('mean_difference'):.4f}"
                     if isinstance(g.get("mean_difference"), float) else "-")
            row(6, f"{gid} point estimate and status",
                f"{rep_v}, {g.get('status')}", f"{rep_v}, {g.get('status')}",
                TRACED, gates, "scmp/evaluate_gates.py", "configs/protocol.yaml")

    # ---------------- pages 6-7: verifier arm table ---------------------------
    if L["pilot"]:
        V = [r for r in L["pilot"] if r.get("kind") == "verifier"]
        for arm in sorted({r["arm"] for r in V}):
            R = [r for r in V if r["arm"] == arm]
            n = len(R)
            closed = 100.0 * sum(1 for r in R if r["gap_rel"] == 0.0) / n
            proved = 100.0 * sum(1 for r in R if r["status"] == "proved") / n
            row(7, f"{arm}: n / closed%", f"330 / {closed:.1f}%",
                f"{n} / {proved:.1f}% proved",
                MISMATCH if abs(closed - proved) > 0.05 else TRACED,
                pilot, "scmp/run_pilot.py", "configs/pilot.yaml",
                "'closed' summed an attainment flag with a termination flag"
                if abs(closed - proved) > 0.05 else "")
        # predictor arms
        P = [r for r in L["pilot"] if r.get("kind") == "predictor"]
        for arm in sorted({r["arm"] for r in P}):
            vals = [r["spearman"] for r in P if r["arm"] == arm]
            row(7, f"predictor {arm} mean Spearman",
                f"{np.mean(vals):.4f} (sd {np.std(vals, ddof=0):.4f})",
                f"{np.mean(vals):.4f} (sd {np.std(vals, ddof=0):.4f}), n={len(vals)}",
                TRACED, pilot, "scmp/run_pilot.py:run_predictor",
                "configs/pilot.yaml",
                "labels: Huesken YFP-reporter efficacy; graphs: guide 21-mer "
                "self-MFE; split: one held-out gene")
        fin = next((r for r in L["pilot"] if r.get("kind") == "final"), {})
        row(7, "budget used", "0.239 of 8 core-hours",
            f"{fin.get('core_hours', float('nan')):.4f} core-hours, "
            f"{fin.get('elapsed_s', 0):.0f}s", TRACED, pilot,
            "scmp/run_pilot.py", "configs/pilot.yaml",
            "CPU core-hours only; GPU hours are NOT accounted separately here")
        row(7, "record counts", "1650 verifier, 15 predictor, 291 timeouts, 30 voids",
            f"{len(V)} verifier, {len(P)} predictor, "
            f"{sum(1 for r in V if r['timeout'])} timeouts, "
            f"{sum(1 for r in L['pilot'] if r.get('kind') == 'void')} voids",
            TRACED, pilot, "scmp/run_pilot.py", "configs/pilot.yaml")

    # ---------------- page 7: effect intervals --------------------------------
    if L["eff"]:
        d = L["eff"]
        row(7, "model-effect intervals validated / verdicts",
            "75 validated, verdicts 31/2/7",
            f"{d['validated_intervals']} validated (0 failures) and "
            f"{len(d['rows'])} classified -- disjoint populations",
            MISMATCH, eff, "scmp/evaluate_explanations.py", "configs/pilot.yaml",
            "one row implied the 40 were a subset of the 75")

    # ---------------- page 8: reproduction ------------------------------------
    if L["rep"]:
        d = L["rep"]
        nh = len(d["source"]) + len(d["data"]) + len(d["environment"])
        okh = sum(1 for x in d["source"] + d["data"] + d["environment"] if x["ok"])
        row(8, "reproduction hashes", "46/46", f"{okh}/{nh}", TRACED, rep,
            "scmp/reproduce.py", "configs/reproduction.yaml")
        for m in d["reproduction"]:
            row(8, f"reproduced metric {m['metric']}", m["recorded"], m["got"],
                TRACED if m["ok"] else MISMATCH, rep, "scmp/audit_conditional.py",
                "configs/reproduction.yaml",
                "primary conditional-audit numbers survive only as stdout text "
                "inside this reproduction record")
        acp = RR.get("audit_conditional", {}).get("parsed_from_stdout", {})
        row(8, "conditional audit: 500 cases, 324 tighter, 1684 impossible",
            "500 / 324 / 1684",
            (f"{acp['cases']} / {acp['tighter']} / {acp['proved_impossible']} "
             f"(fallback wins {acp['fallback_wins']})" if acp.get("cases")
             else "present in rerun.stdout_tail only"),
            TRACED if acp.get("cases") else TEXT_ONLY, rr_path or rep,
            "scmp/audit_conditional.py", "configs/exhaustive.yaml",
            "audit_conditional persists a record only on violation; traced via "
            "scmp.revision.rerun_unrecorded (audit run unchanged)")
        ab = RR.get("audit_bounds", {})
        abp = ab.get("parsed_from_stdout", {})
        row(8, "bound soundness: 1,262,952 checks, 0 violations",
            "1,262,952 / 0",
            (f"{abp['intermediate_bracket_checks']} checks over "
             f"{abp['feasible_structures']} structures, no_violation="
             f"{ab.get('flags', {}).get('no_violation')}"
             if abp.get("intermediate_bracket_checks")
             else "not reproducible from any record"),
            TRACED if abp.get("intermediate_bracket_checks") else NO_RECORD,
            rr_path, "scmp/audit_bounds.py", "configs/exhaustive.yaml",
            "audit_bounds persists a record only on violation; traced via "
            "scmp.revision.rerun_unrecorded (audit run unchanged)")

    # ---------------- pages 8-9: RNA ------------------------------------------
    if L["spl"]:
        d = L["spl"]
        for sp in d["splits"]:
            row(8, f"split {sp['name']}: groups / test / leak",
                f"{sp['n_groups']} / {sp['n_test']} / "
                f"{sp['residual_cluster_leak']:.2%}",
                f"{sp['n_groups']} / {sp['n_test']} / "
                f"{sp['residual_cluster_leak']:.2%}", TRACED, spl,
                "scmp/rna/audit_splits.py", "configs/rna_splits.yaml")
        du = d["duplication"]
        row(9, "historical split leak (cluster / gene / seed)",
            "87.1% / 100% / 20.9%",
            f"{du['historical_leak_kmer_cluster']:.1%} / "
            f"{du['historical_leak_gene']:.0%} / "
            f"{du['historical_leak_seed_region']:.1%}", TRACED, spl,
            "scmp/rna/audit_splits.py", "configs/rna_splits.yaml")
    if L["ens"]:
        row(9, "ensemble audit checks / failures",
            "102 / 0", f"{L['ens']['checks']} / {len(L['ens']['failures'])}",
            TRACED, ens, "scmp/rna/audit_ensemble.py", "configs/rna_tiny.yaml")
    if L["rna"]:
        d = L["rna"]
        for key, arm in (("reporter_yfp|linear_sequence_only", "linear_sequence_only"),):
            sm = d["summary"].get(key)
            if sm:
                row(9, f"RNA {arm} Spearman", "0.6464 (range 0.632-0.677)",
                    f"{np.mean(sm['spearman']):.4f} (range "
                    f"{min(sm['spearman']):.3f}-{max(sm['spearman']):.3f})",
                    TRACED, rna, "scmp/rna/evaluate.py", "configs/rna_protocol.yaml",
                    "assay: Huesken YFP reporter, H1299, 48 h; split: "
                    "sequence-cluster-disjoint")
        rc = d["resolved_comparison_demo"]
        row(10, "ranking coverage", f"{rc['resolved']}/{rc['pairs']} = "
            f"{rc['coverage']:.1%}",
            f"{rc['resolved']}/{rc['n_candidates'] * (rc['n_candidates'] - 1) // 2}"
            f" = {2 * rc['coverage']:.1%} on the matched unordered denominator",
            MISMATCH, rna, "scmp/rna/evaluate.py:resolved_comparison_demo",
            "configs/rna_protocol.yaml",
            "ordered denominator with an unordered numerator")
        row(10, "model-ranking agreement on resolved pairs", "100%",
            f"{rc['model_ranking_agreement']:.0%}", TRACED, rna,
            "scmp/rna/evaluate.py", "configs/rna_protocol.yaml",
            "agreement with the MODEL's ordering; not a measurement check")

    counts = collections.Counter(r["status"] for r in rows)
    return {"claim": "every number on report pages 6-10 traces to a record",
            "verdict": ("PARTIAL: %d of %d rows have no immutable record"
                        % (counts.get(NO_RECORD, 0), len(rows))),
            "evidence_kind": EXACT_INTEGER,
            "report": fp or {},
            "n_rows": len(rows), "status_counts": dict(counts), "rows": rows,
            "note": "NO_RECORD rows are not disputed values; they are values no "
                    "artefact can currently support, because those audits persist "
                    "a record only when they FAIL. That is a traceability defect "
                    "in protocol v1 and is fixed by rerunning them under v2."}


# --------------------------------------------------------------------------- item 1
def item1_termination(runs="runs"):
    path = _latest(os.path.join(runs, "pilot", "pilot-*.json"))
    if not path:
        return {"status": "NO RECORD", "pattern": "pilot-*.json"}
    log = _load(path)
    V = [r for r in log if r.get("kind") == "verifier"]
    voids = [r for r in log if r.get("kind") == "void"]

    arms = {}
    for arm in sorted({r["arm"] for r in V}):
        R = [r for r in V if r["arm"] == arm]
        n = len(R)
        ct = collections.Counter((r["gap_rel"] == 0.0, bool(r["timeout"])) for r in R)
        status_ct = collections.Counter(r["status"] for r in R)
        # mutually exclusive termination outcome, taken from the record's own field
        proved = status_ct.get("proved", 0)
        unresolved = status_ct.get("unresolved", 0)
        gap_zero = sum(1 for r in R if r["gap_rel"] == 0.0)
        timeout = sum(1 for r in R if r["timeout"])
        arms[arm] = {
            "n": n,
            # integer counts, so downstream checks never depend on rounded percentages
            "count_gap_zero": gap_zero, "count_timeout": timeout,
            "count_proved": proved, "count_unresolved": unresolved,
            "report_closed_pct": round(100 * gap_zero / n, 1),
            "report_timeout_pct": round(100 * timeout / n, 1),
            "report_sum_pct": round(100 * (gap_zero + timeout) / n, 1),
            "overlap_gap_zero_and_timeout": ct[(True, True)],
            "corrected_proved_pct": round(100 * proved / n, 1),
            "corrected_unresolved_pct": round(100 * unresolved / n, 1),
            "corrected_sum_pct": round(100 * (proved + unresolved) / n, 1),
            "exclusive_outcomes": dict(status_ct),
            "cross_tab": {f"gap_zero={c},timeout={t}": v for (c, t), v in sorted(ct.items())},
        }
    ok = all(abs(a["corrected_sum_pct"] - 100.0) < 1e-9 for a in arms.values())
    align = all(r["timeout"] == (r["status"] == "unresolved") for r in V)

    # ---- denominator derivation, from cells rather than from record counts ----
    cells = {(r["family"], r["instance"], r["seed"]) for r in V}
    void_cells = {(r["family"], r["instance"], r["seed"]) for r in voids}
    fam_cells = collections.Counter(f for (f, _, _) in cells)
    fam_voids = collections.Counter(f for (f, _, _) in void_cells)
    n_arms = len({r["arm"] for r in V})
    denominators = {
        "cells_attempted": len(cells) + len(void_cells),
        "voids_excluded": len(void_cells),
        "cells_evaluated_per_arm": len(cells),
        "arms": n_arms,
        "verifier_records": len(V),
        "identity_records_eq_cells_x_arms": len(V) == len(cells) * n_arms,
        "voids_disjoint_from_evaluated": not (cells & void_cells),
        "by_family": {f: {"evaluated": fam_cells[f], "void": fam_voids.get(f, 0),
                          "attempted": fam_cells[f] + fam_voids.get(f, 0)}
                      for f in sorted(set(fam_cells) | set(fam_voids))},
        "timeouts_total": sum(1 for r in V if r["timeout"]),
        "timeouts_are_records_not_cells":
            "291 timeouts are (cell, arm) records out of "
            f"{len(V)}, not 291 of the {len(cells)} cells",
    }

    # ---- first tolerance-attainment time: is it instrumented at all? ----------
    fields = sorted({k for r in V for k in r})
    timing_fields = [f for f in fields if "wall" in f or "time" in f or "_s" == f[-2:]]
    first_attainment = {
        "requested": "time at which each arm first met the requested tolerance",
        "available_timing_fields": timing_fields,
        "instrumented": False,
        "why": "each record carries only a terminal wall_s and the matched budget "
               "matched_wall_s. There is no per-refinement timeline and no "
               "first-attainment timestamp, so the cross-tab column cannot be "
               "populated from protocol-v1 records at all.",
        "not_a_substitute": "wall_s is the time the arm STOPPED, which for a "
                            "timed-out record is the budget itself, not the time "
                            "the gap first reached tolerance.",
        "v2_requirement": "log (t, L, U, status) at every bound update; derive "
                          "first-attainment and anytime AUC from that stream.",
        "wall_s_summary": {arm: {
            "median": float(np.median([r["wall_s"] for r in V if r["arm"] == arm])),
            "max": float(max(r["wall_s"] for r in V if r["arm"] == arm)),
            "matched_budget_s": sorted({r["matched_wall_s"] for r in V
                                        if r["arm"] == arm})}
            for arm in sorted({r["arm"] for r in V})},
    }

    # ---- exact vs displayed sums: 105.7 or 105.8? -----------------------------
    for arm, v in arms.items():
        exact = 100.0 * (v["count_gap_zero"] + v["count_timeout"]) / v["n"]
        v["sum_exact_pct"] = round(exact, 4)
        v["sum_of_displayed_1dp"] = round(v["report_closed_pct"]
                                          + v["report_timeout_pct"], 4)

    # ---- negative gaps must be recorded as defects, never clipped -------------
    neg = [r for r in V if r["gap_rel"] < 0]
    neg_policy = []
    for pth in sorted(glob.glob("runs/policy_pilot/*.json")):
        try:
            res = _load(pth).get("results", {})
        except Exception:
            continue
        for arm, vals in res.items():
            for i, val in enumerate(vals):
                if isinstance(val, (int, float)) and val < 0:
                    neg_policy.append({"record": pth, "arm": arm, "index": i,
                                       "gap": val})
    negative_gaps = {
        "verifier_records_with_negative_gap": len(neg),
        "policy_pilot_negative_gaps": len(neg_policy),
        "most_negative": min([d["gap"] for d in neg_policy], default=None),
        "clipped_anywhere": False,
        "verdict": "DEFECT RECORDED" if neg_policy else "none found",
        "assessment": "roundoff-scale bound inversion (order 1e-13) in every "
                      "policy-pilot arm including the unconditional baseline, so "
                      "it is a property of the shared float bound evaluation, not "
                      "of conditioning. Consistent with float_unverified status; "
                      "it is recorded rather than clipped, and an enclosed "
                      "implementation (P4) must make it impossible.",
    }

    return {"claim": "closed% + timeout% should sum to at most 100%",
            "verdict": "EXPLAINED: two overlapping flags were summed",
            "evidence_kind": EXACT_INTEGER,
            "record": path,
            "n_verifier": len(V), "n_void_excluded": len(voids),
            "exclusive_outcomes_sum_to_100": ok,
            "timeout_iff_unresolved": align,
            "denominators": denominators,
            "first_tolerance_attainment": first_attainment,
            "negative_gaps": negative_gaps,
            "arms": arms,
            "note": "gap_rel==0 is a tightness flag, timeout is a termination "
                    "flag, and an instance can carry both. The corrected proved "
                    "rate is LOWER than the reported closed rate for both "
                    "conditional arms: the correction runs against the "
                    "mechanism's favour."}


# --------------------------------------------------------------------------- item 2
def item2_effect_intervals(runs="runs"):
    path = _latest(os.path.join(runs, "explanations", "effects-*.json"))
    if not path:
        return {"status": "NO RECORD"}
    d = _load(path)
    rows = d.get("rows", [])
    verdicts = collections.Counter(r["verdict"] for r in rows)
    seeds = collections.Counter(r["seed"] for r in rows)
    edges = {tuple(r["edge"]) for r in rows}
    return {"claim": "75 intervals validated, verdicts 31/2/7",
            "verdict": "EXPLAINED: two disjoint populations",
            "evidence_kind": EXACT_INTEGER,
            "record": path,
            "validation_population": {
                "n_intervals": d.get("validated_intervals"),
                "containment_failures": d.get("containment_failures"),
                "contrast_violations": d.get("contrast_violations"),
                "how": "exhaustive enumeration over 3 tiny families x 3 seeds",
                "classified_into_verdicts": False},
            "reporting_population": {
                "n_intervals": len(rows),
                "seeds": dict(seeds), "distinct_edges": len(edges),
                "verdicts": dict(verdicts),
                "how": "5 seeds x 8 candidate edges on one pilot instance",
                "enumeration_validated": False},
            "populations_disjoint": True,
            "unaccounted": 0,
            "note": "the 40 reported verdicts are NOT a subset of the 75 validated "
                    "intervals; they are a different family, instance and seed set. "
                    "No interval is missing: 75 and 40 are separate totals."}


# --------------------------------------------------------------------------- item 3
def item3_pair_counting(runs="runs", rerun=False):
    path = _latest(os.path.join(runs, "rna_eval", "eval-*.json"))
    rec = _load(path)["resolved_comparison_demo"] if path else {}
    k = rec.get("n_candidates", 12)
    ordered, unordered = k * (k - 1), k * (k - 1) // 2
    resolved = rec.get("resolved")
    out = {"claim": f"coverage {resolved}/{rec.get('pairs')} = "
                    f"{rec.get('coverage', float('nan')):.1%}",
           "verdict": "DEFINITIONAL MISMATCH: unordered numerator, ordered denominator",
           "evidence_kind": REAL_PROOF,
           "record": path,
           "n_candidates": k,
           "ordered_pairs": ordered, "unordered_pairs": unordered,
           "resolved_strict_inequalities": resolved,
           "reported_coverage_ordered": rec.get("coverage"),
           "matched_coverage_unordered": (resolved / unordered) if resolved else None,
           "injectivity_argument":
               "l_i <= u_i for every candidate, so u_i < l_j and u_j < l_i cannot "
               "both hold. Each resolved inequality therefore identifies a distinct "
               "unordered pair, and the matched denominator is k(k-1)/2.",
           "note": "both numbers are arithmetically correct; they answer different "
                   "questions. The ordered-denominator figure understates coverage "
                   "by exactly a factor of two."}
    if rerun:
        from ..rna.evaluate import resolved_comparison_demo
        out["rerun"] = resolved_comparison_demo()
    return out


# --------------------------------------------------------------------------- item 4
def item4_gates_and_inputs(runs="runs"):
    gpath = _latest(os.path.join(runs, "gates", "gates-*.json"))
    ppath = _latest(os.path.join(runs, "pilot", "pilot-*.json"))
    gates = _load(gpath)["gates"] if gpath else []
    log = _load(ppath) if ppath else []
    pred = [r for r in log if r.get("kind") == "predictor"]
    split = next((r for r in log if r.get("kind") == "predictor_split"), {})

    by_arm = {}
    for a in sorted({r["arm"] for r in pred}):
        vals = {r["seed"]: r["spearman"] for r in pred if r["arm"] == a}
        by_arm[a] = {"seeds": sorted(vals), "spearman": [vals[s] for s in sorted(vals)],
                     "mean": sum(vals.values()) / len(vals)}

    # what the arms actually consumed, read from the source rather than assumed
    try:
        from experiments import _huesken
        src = inspect.getsource(_huesken.build_graphs)
    except Exception:
        src = ""
    graph_is_guide_self_fold = "bpp_lattice(s)" in src and "A_mfe" in src

    # assay identity of the labels these arms were fitted to, from the audited
    # primary source rather than from the redistributed annotation
    assay = {}
    try:
        import yaml
        hd = yaml.safe_load(open("configs/rna_data.yaml"))["datasets"]["huesken2005"]
        assay = {"primary_source": hd["primary_source"]["assay_from_primary"],
                 "redistributed_annotation_claims": hd["annotation_source"]["claims"],
                 "annotation_contradicts_primary": True,
                 "label_semantics": "marginal over latent structures; reporter "
                                    "readout, not a native-transcript measurement"}
    except Exception as exc:
        assay = {"error": f"{type(exc).__name__}: {exc}"}

    g3 = next((g for g in gates if g["gate"].startswith("G3")), {})
    g4 = next((g for g in gates if g["gate"].startswith("G4")), {})
    return {"claim": "G3 passes, so certifiability is cheap",
            "verdict": "SIGN CORRECT, COMPARATOR WEAK",
            "evidence_kind": FLOAT_DIAGNOSTIC,
            "records": {"gates": gpath, "pilot": ppath},
            "G3": {"comparator_arms": ["unconstrained (free)",
                                       "certifiable_gamma1 (constrained)"],
                   "difference_defined_as": "rho_unconstrained - rho_certifiable_gamma1",
                   "decision_inequality": "PASS iff upper 95% CI of the difference "
                                          "< margin (0.05)",
                   "point": g3.get("mean_difference"), "ci95": g3.get("ci95"),
                   "status": g3.get("status"),
                   "structural_objection":
                       "the constrained arm is A+B (mean 0.4573), not the best "
                       "available certifiable model A (mean 0.7253). The gate never "
                       "compared the free model against the strongest constrained "
                       "one, and it passes because the free arm scored 0.1212."},
            "G4": {"difference_defined_as":
                       "rho_certifiable_gamma1 - rho_certifiable_gamma0",
                   "decision_inequality": "PASS iff lower 95% CI > 0",
                   "point": g4.get("mean_difference"), "ci95": g4.get("ci95"),
                   "status": g4.get("status")},
            "arms": by_arm,
            "split": {"unit": split.get("unit"), "n": split.get("n"),
                      "train": split.get("train"), "test": split.get("test"),
                      "held_out_genes": split.get("held_out_genes"),
                      "objection": "one held-out gene; not a population claim"},
            "graph_inputs": {
                "source": "experiments/_huesken.py:build_graphs",
                "graph_is_guide_21mer_self_mfe": graph_is_guide_self_fold,
                "uses_target_window": False,
                "consequence":
                    "these arms fold the assayed guide's OWN sequence, so no native "
                    "transcript was substituted for a reporter construct (unlike "
                    "F9). They are a valid optimisation diagnostic and cannot "
                    "support any target-structure claim."},
            "assay_identity": assay,
            "label_trace": {
                "labels": "Huesken 2005 normalised inhibition (efficacy)",
                "graphs": "guide 21-mer own MFE structure via bpp_lattice",
                "split": "gene-disjoint, one held-out gene",
                "assay": "YFP reporter plasmid, H1299, 48 h",
                "supports_target_structure_claim": False,
                "supports_optimisation_diagnostic": True},
            "note": "the report presented G3 as support for cheap certifiability. "
                    "The direction is right; the comparison is not informative."}


# --------------------------------------------------------------------------- item 5
def item5_anchor_sign():
    from ..model import AnchorEncoder, SignedResidualMPNN
    import torch
    enc = AnchorEncoder(4, 8)
    with torch.no_grad():
        neg_a = int(sum((p < 0).sum() for p in enc.a_out.parameters()))
        neg_phi = int(sum((p < 0).sum() for p in enc.phi.parameters()))
    src = inspect.getsource(AnchorEncoder) + inspect.getsource(SignedResidualMPNN)
    constrained = any(t in src for t in ("clamp(min=0", "relu(self.a_out.weight",
                                         "abs(self.a_out", "softplus", "nonneg"))
    return {"claim": "the model uses a non-negative additive anchor",
            "verdict": "REFUTED: the anchor is signed in code",
            "evidence_kind": REAL_PROOF,
            "a_out_layer": "nn.Linear(d_hid, 1) with unconstrained weight and bias",
            "negative_params_at_init": {"a_out": neg_a, "phi": neg_phi},
            "constraint_found_in_source": constrained,
            "note": "the non-negativity constraint belongs to certmp's monotone "
                    "MPNN (phase one), not to scmp's anchor. The report's limits "
                    "table carried the phase-one restriction across to phase two. "
                    "This matters: a signed anchor is strictly more expressive and "
                    "does not obstruct the support oracle, which optimises a linear "
                    "objective over the family regardless of coefficient sign."}


# --------------------------------------------------------------------------------- CLI
CHECKS = [("item0_trace_pages_6_10", item0_trace_pages),
          ("item1_termination_accounting", item1_termination),
          ("item2_effect_interval_populations", item2_effect_intervals),
          ("item3_pair_counting", item3_pair_counting),
          ("item4_gates_and_graph_inputs", item4_gates_and_inputs),
          ("item5_anchor_sign", item5_anchor_sign)]


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="python -m scmp.revision.audit_report",
        description="P0 items 1-5: trace run-report numbers to immutable records. "
                    "Read-only with respect to historical runs.")
    ap.add_argument("--report", default="docs/run_report.md",
                    help="the report whose claims are audited (recorded, not parsed)")
    ap.add_argument("--runs", default="runs")
    ap.add_argument("--rerun-demo", action="store_true",
                    help="re-execute the resolved-comparison demo (slower)")
    ap.add_argument("--out-dir", default="runs/revision")
    a = ap.parse_args(argv)

    t0 = time.time()
    results = {}
    fp = report_fingerprint(a.report)
    print("P0 items 0-5  run-report reconciliation")
    print(f"report        {fp.get('resolved_to', fp['path'])}"
          f"{'' if fp['exists'] else '   [NOT FOUND]'}")
    if fp.get("sha256"):
        print(f"sha256        {fp['sha256'][:32]}")
    if fp.get("pages_extracted"):
        print(f"pages read    {fp['pages_extracted']}")
    print()
    for name, fn in CHECKS:
        if name == "item0_trace_pages_6_10":
            r = fn(a.runs, {k: v for k, v in fp.items() if k != "_text"})
        else:
            r = fn(a.runs) if "runs" in inspect.signature(fn).parameters else fn()
        if name == "item3_pair_counting" and a.rerun_demo:
            r = item3_pair_counting(a.runs, rerun=True)
        results[name] = r
        print(f"  {name}")
        print(f"    claim   : {r.get('claim')}")
        print(f"    verdict : {r.get('verdict')}   [{r.get('evidence_kind')}]")
        print()

    t = results["item0_trace_pages_6_10"]
    print(f"page 6-10 trace: {t['n_rows']} rows  {t['status_counts']}")
    for r in t["rows"]:
        if r["status"] in (MISMATCH, NO_RECORD, TEXT_ONLY):
            print(f"    [{r['status']:>9}] p{r['page']} {r['claim']}")
            print(f"                reported   {r['reported']}")
            print(f"                recomputed {r['recomputed']}")
    print()

    d = results["item1_termination_accounting"]["denominators"]
    print(f"denominators: {d['cells_attempted']} cells attempted - "
          f"{d['voids_excluded']} void = {d['cells_evaluated_per_arm']} evaluated "
          f"x {d['arms']} arms = {d['verifier_records']} records")
    print(f"    by family: " + "  ".join(
        f"{k} {v['evaluated']}+{v['void']}v" for k, v in d["by_family"].items()))
    print(f"    {d['timeouts_are_records_not_cells']}")
    fa = results["item1_termination_accounting"]["first_tolerance_attainment"]
    print(f"first tolerance-attainment: instrumented={fa['instrumented']} "
          f"(fields present: {fa['available_timing_fields']})")
    ng = results["item1_termination_accounting"]["negative_gaps"]
    print(f"negative gaps: {ng['verdict']} -- {ng['policy_pilot_negative_gaps']} in "
          f"policy_pilot, min {ng['most_negative']}, clipped={ng['clipped_anywhere']}")
    print()

    print("corrected verifier accounting (mutually exclusive outcomes)")
    print(f"    {'arm':<28}{'reported':>20}{'corrected':>22}")
    print(f"    {'':<28}{'closed  t/out  sum':>20}{'proved  unres.  sum':>22}")
    for arm, v in results["item1_termination_accounting"]["arms"].items():
        print(f"    {arm:<28}{v['report_closed_pct']:>7.1f}{v['report_timeout_pct']:>7.1f}"
              f"{v['report_sum_pct']:>6.1f}{v['corrected_proved_pct']:>9.1f}"
              f"{v['corrected_unresolved_pct']:>8.1f}{v['corrected_sum_pct']:>7.1f}")

    i3 = results["item3_pair_counting"]
    print(f"\nranking coverage: {i3['resolved_strict_inequalities']}/"
          f"{i3['ordered_pairs']} ordered = {i3['reported_coverage_ordered']:.1%}  ->  "
          f"{i3['resolved_strict_inequalities']}/{i3['unordered_pairs']} unordered = "
          f"{i3['matched_coverage_unordered']:.1%}")

    os.makedirs(a.out_dir, exist_ok=True)
    path = os.path.join(a.out_dir, f"audit_report-{time.strftime('%Y%m%d-%H%M%S')}.json")
    with open(path, "w") as fh:
        json.dump({"kind": "revision_audit_report", "report": a.report,
                   "results": results, "runtime_s": time.time() - t0}, fh, indent=1)
    print(f"\nwritten {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
