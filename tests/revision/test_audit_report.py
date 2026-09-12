"""The P0 record-level reconciliations, as regression tests.

These read the immutable pilot/gate/effect records. They assert the ARITHMETIC
relationships that the corrections rest on, so a future change to the accounting
cannot quietly reintroduce the overlapping-flag bug.
"""
import json
import os

import pytest

from scmp.revision.audit_report import (item1_termination, item2_effect_intervals,
                                        item3_pair_counting, item4_gates_and_inputs,
                                        item5_anchor_sign)

pytestmark = pytest.mark.skipif(not os.path.isdir("runs/pilot"),
                                reason="immutable pilot records not present")


# ------------------------------------------------------------ item 1
def test_exclusive_outcomes_sum_to_exactly_100_percent():
    r = item1_termination()
    assert r["exclusive_outcomes_sum_to_100"] is True
    for arm, v in r["arms"].items():
        assert set(v["exclusive_outcomes"]) <= {"proved", "unresolved"}
        assert v["corrected_proved_pct"] + v["corrected_unresolved_pct"] == \
            pytest.approx(100.0, abs=0.05), arm


def test_timeout_holds_exactly_when_unresolved():
    assert item1_termination()["timeout_iff_unresolved"] is True


def test_overlap_explains_the_excess_and_correction_is_unfavourable():
    """The corrected rate must be <= the reported one wherever flags overlapped."""
    r = item1_termination()
    overlapping = {a: v for a, v in r["arms"].items()
                   if v["overlap_gap_zero_and_timeout"] > 0}
    assert overlapping, "expected the conditional arms to overlap"
    for arm, v in overlapping.items():
        assert v["report_sum_pct"] > 100.0, arm
        assert v["corrected_proved_pct"] < v["report_closed_pct"], arm
        # exact integer identity: the excess IS the double-counted instances
        assert v["count_gap_zero"] + v["count_timeout"] - v["n"] == \
            v["overlap_gap_zero_and_timeout"], arm
        # and the corrected proved count excludes exactly those
        assert v["count_proved"] == v["count_gap_zero"] - \
            v["overlap_gap_zero_and_timeout"], arm


# ------------------------------------------------------------ item 2
def test_effect_interval_populations_are_disjoint_and_accounted():
    r = item2_effect_intervals()
    assert r["validation_population"]["containment_failures"] == 0
    assert r["reporting_population"]["n_intervals"] == \
        sum(r["reporting_population"]["verdicts"].values())
    assert r["validation_population"]["classified_into_verdicts"] is False
    assert r["unaccounted"] == 0


# ------------------------------------------------------------ item 3
def test_pair_counting_denominators():
    r = item3_pair_counting()
    k = r["n_candidates"]
    assert r["ordered_pairs"] == k * (k - 1)
    assert r["unordered_pairs"] == k * (k - 1) // 2
    # resolved inequalities inject into unordered pairs, so they cannot exceed them
    assert r["resolved_strict_inequalities"] <= r["unordered_pairs"]
    assert r["matched_coverage_unordered"] == pytest.approx(
        2 * r["reported_coverage_ordered"])


# ------------------------------------------------------------ item 4
def test_g3_direction_and_comparator_are_explicit():
    r = item4_gates_and_inputs()
    assert r["G3"]["difference_defined_as"] == \
        "rho_unconstrained - rho_certifiable_gamma1"
    means = {a: v["mean"] for a, v in r["arms"].items()}
    # the gate's constrained arm is not the best certifiable arm available
    assert means["certifiable_gamma0"] > means["certifiable_gamma1"]
    assert means["certifiable_gamma1"] > means["unconstrained"]


def test_predictor_graphs_are_guide_self_folds_not_target_windows():
    r = item4_gates_and_inputs()
    assert r["graph_inputs"]["graph_is_guide_21mer_self_mfe"] is True
    assert r["graph_inputs"]["uses_target_window"] is False


def test_split_is_a_single_held_out_gene():
    r = item4_gates_and_inputs()
    assert len(r["split"]["held_out_genes"]) == 1


# ------------------------------------------------------------ item 5
def test_anchor_is_signed_not_nonnegative():
    r = item5_anchor_sign()
    assert r["constraint_found_in_source"] is False
    assert r["negative_params_at_init"]["a_out"] > 0
    assert r["verdict"].startswith("REFUTED")


# ------------------------------------------------------------ ledger
def test_ledger_retains_historical_gates_unchanged():
    from scmp.revision.ledger import HISTORICAL, build
    entries = build()
    hist = [e for e in entries if e["status"] == HISTORICAL]
    assert len(hist) == 4                       # G1..G4
    assert any("UNDETERMINED" in e["reported"] for e in hist)
    for e in hist:
        assert e["corrected"] == "retained unchanged"
