"""The P0 mathematical corrections, as regression tests.

These pin exact rational values. If a future refactor changes the graph encoding
or the walk-score convention, these fail rather than silently reverting to the
claims the run report made.
"""
from fractions import Fraction

import pytest

from scmp.revision import EXACT_RATIONAL
from scmp.revision.audit_math import (adjacency, check_bounded_feature_repair,
                                      check_design_matrix_residual,
                                      check_kappa_counterexample,
                                      check_layer1_relu,
                                      check_matching_kappa_formula,
                                      check_three_node_path, matchings,
                                      walk_score)


# ------------------------------------------------------------------ enumeration
def test_matchings_of_triangle_excludes_conflicting_pairs():
    """Grown incrementally, so this must still be exactly the conflict-free set."""
    got = set(matchings([(0, 1), (1, 2), (0, 2)]))
    assert got == {(), ((0, 1),), ((1, 2),), ((0, 2),)}


def test_matchings_count_matches_known_value_for_K4():
    # K4 has 10 matchings: 1 empty + 6 single + 3 perfect
    assert len(matchings([(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)])) == 10


def test_walk_score_is_sum_of_squared_degrees_at_L2():
    """1^T A^2 1 = ||A1||^2 = sum_i d_i^2, in exact arithmetic."""
    A = adjacency(5, [(0, 1), (1, 2), (3, 4)])
    assert walk_score(A, [Fraction(1)] * 5, 2) == Fraction(1 + 4 + 1 + 1 + 1)


# ------------------------------------------------------- item 6: additivity in z
def test_three_node_path_is_not_additive_in_z():
    r = check_three_node_path()
    assert r["observed"] == ["0", "2", "2", "6"]
    assert r["mixed_difference"] == "2"          # f = 2z1 + 2z2 + 2z1z2
    assert r["all_relus_stable"] is True         # nonlinearity is NOT from clipping
    assert r["evidence_kind"] == EXACT_RATIONAL
    assert r["verdict"] == "REFUTED"


def test_real_model_residual_is_not_additive_in_z():
    """The actual torch model on an enumerated feasible family, not a toy."""
    r = check_design_matrix_residual(seed=0)
    g = r["projection_onto_[1,z]"]["g_residual_only"]
    assert g["n_members"] > 1 and g["n_edges"] > 1
    assert g["additive_in_z"] is False
    assert g["relative_residual"] > 1e-6


# --------------------------------------------------------- item 7: layer-1 ReLU
def test_layer1_relu_is_relaxed_and_can_be_unstable():
    r = check_layer1_relu(seed=0)
    assert r["message_M1_exactly_affine"] is True
    assert r["relu_relaxation_applied_at_H1"] is True   # the code does relax it
    assert r["unstable_units_at_layer_1"] >= 1          # and it is needed
    # ReLU(z1+z2-1) is affine-preactivated yet non-additive
    assert r["relu_demo_mixed_difference"] == 1


# ------------------------------------------------- item 8: kappa counterexample
def test_kappa_is_not_feature_independent():
    r = check_kappa_counterexample()
    assert r["uniform_features"]["kappa_2"] == "2"
    assert r["skewed_features"]["ratio"] == "31/11"
    assert Fraction(31, 11) > 2
    assert r["exceeds_kappa"] is True
    assert r["evidence_kind"] == EXACT_RATIONAL


def test_bounded_feature_repair_holds():
    r = check_bounded_feature_repair()
    assert r["holds"] is True
    # the repair is loose: it must not be mistaken for the tight ratio
    assert Fraction(r["bound"]) > Fraction(r["observed_ratio"])


# -------------------------------------------------- item 9: matching kappa form
@pytest.mark.parametrize("m,L", [(1, 1), (2, 1), (2, 2), (3, 2), (3, 3)])
def test_matching_kappa_closed_forms(m, L):
    """(2m-1)^L without self-loops, m^L with a self-loop at every vertex."""
    r = check_matching_kappa_formula(max_m=m, layers=(L,))
    rows = [x for x in r["rows"] if x["m"] == m and x["L"] == L]
    assert rows, "closed form not evaluated"
    for row in rows:
        assert row["agrees"], row
        want = m ** L if row["self_loops"] else (2 * m - 1) ** L
        assert Fraction(row["kappa"]) == want


def test_exponent_depends_on_family_not_only_depth():
    """Two conventions on the same family give different exponents in L."""
    r = check_matching_kappa_formula(max_m=3, layers=(2,))
    by = {(x["m"], x["self_loops"]): Fraction(x["kappa"]) for x in r["rows"]}
    assert by[(3, False)] == 25 and by[(3, True)] == 9      # 5^2 vs 3^2
    assert by[(2, False)] == 9 and by[(2, True)] == 4       # 3^2 vs 2^2
    assert r["verdict"] == "REFUTED as universal"


# --------------------------- item 9b: the repository's actual convention
def test_actual_convention_is_self_loops_plus_chain():
    from scmp.revision.audit_math import check_actual_convention_kappa
    r = check_actual_convention_kappa(max_m=2, layers=(1, 2))
    assert r["backbone_adjacency_defaults"] == {"self_loops": True, "chain": True}


def test_actual_convention_kappa_matches_neither_closed_form():
    """B = I + chain is in every member, so the candidate-only formulas do not apply."""
    from fractions import Fraction

    from scmp.revision.audit_math import check_actual_convention_kappa
    r = check_actual_convention_kappa(max_m=3, layers=(1, 2, 3))
    act = {(x["m"], x["L"]): Fraction(x["kappa"]) for x in r["rows"]
           if x["convention"].startswith("ACTUAL")}
    assert act[(3, 2)] == Fraction(177, 41)
    assert act[(2, 2)] == Fraction(61, 25)
    # strictly between 1 and the self-loop-only value: the shared core compresses it
    for (m, L), k in act.items():
        if m > 1:
            assert 1 < k < Fraction(m) ** L, (m, L, k)


def test_actual_convention_is_not_a_clean_power():
    from fractions import Fraction

    from scmp.revision.audit_math import check_actual_convention_kappa
    r = check_actual_convention_kappa(max_m=3, layers=(1, 2))
    act = {(x["m"], x["L"]): Fraction(x["kappa"]) for x in r["rows"]
           if x["convention"].startswith("ACTUAL")}
    assert act[(3, 2)] != act[(3, 1)] ** 2      # kappa_2 != kappa_1^2
