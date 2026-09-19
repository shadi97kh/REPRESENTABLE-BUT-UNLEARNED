"""Randomised audit of the support oracles against independent enumeration.

Every quantity the oracles expose is checked against `scmp.oracles.bruteforce`, which
never calls the recurrence. Exact-integer quantities must match exactly; floating-point
quantities are compared with a tolerance and are reported separately, because a float
agreement is a diagnostic and never evidence of exactness.

Run:
    python -m scmp.audit_oracles --max-n 10 --cases 200 --seed 7027
"""
from __future__ import annotations

import argparse
import math
import random
import sys
import time

from .oracles import (EXACT_INTEGER, FLOAT_DIAGNOSTIC, INFEASIBLE, NEG_INF,
                      LayeredDAGPathOracle, RNANonCrossingOracle, derivation_edges,
                      score_from_derivation)
from .oracles.bruteforce import (brute_count, brute_log_partition, brute_marginals,
                                 brute_support, enumerate_dag_paths, enumerate_rna,
                                 rna_is_valid)

# Non-crossing partial matchings on n points are counted by the Motzkin numbers
# (OEIS A001006). A closed-form check independent of both the DP and the enumerator.
MOTZKIN = [1, 1, 2, 4, 9, 21, 51, 127, 323, 835, 2188, 5798, 15511, 41835, 113634]

FLOAT_TOL = 1e-9


class Audit:
    def __init__(self):
        self.exact_checks = 0
        self.float_checks = 0
        self.failures = []

    def exact(self, ok, label, got=None, want=None):
        self.exact_checks += 1
        if not ok:
            self.failures.append(("exact", label, got, want))

    def approx(self, got, want, label, tol=FLOAT_TOL):
        self.float_checks += 1
        if want == NEG_INF or got == NEG_INF:
            ok = (got == want)
        else:
            ok = abs(got - want) <= tol * max(1.0, abs(want))
        if not ok:
            self.failures.append(("float", label, got, want))


def random_sequence(rng, n):
    return "".join(rng.choice("ACGU") for _ in range(n))


def audit_rna_case(a: Audit, rng, n):
    min_loop = rng.choice([0, 1, 2, 3])
    canon = rng.choice([True, False])
    seq = random_sequence(rng, n)
    o = RNANonCrossingOracle(seq, min_loop=min_loop, canonical_only=canon)
    members = enumerate_rna(seq, min_loop, canon)
    tag = f"rna n={n} min_loop={min_loop} canon={canon}"

    c = o.count()
    a.exact(c.numerical_status == EXACT_INTEGER, f"{tag}: count status")
    a.exact(c.value == brute_count(members), f"{tag}: count", c.value, len(members))

    ground = o.ground_set()
    coef = {e: rng.randint(-9, 9) for e in ground}
    s = o.support(coef)
    exp, _ = brute_support(members, coef)
    a.exact(s.numerical_status == EXACT_INTEGER, f"{tag}: support status")
    a.exact(s.value == exp, f"{tag}: support", s.value, exp)

    if s.termination_status != INFEASIBLE:
        a.exact(rna_is_valid(s.witness, seq, min_loop, canon), f"{tag}: witness valid")
        a.exact(sum(coef.get(e, 0) for e in s.witness) == s.value,
                f"{tag}: witness scores as reported")
        prods = s.stats["productions"]
        edges = derivation_edges(prods)
        a.exact(len(edges) == len(set(edges)), f"{tag}: P injective on derivation")
        a.exact(tuple(sorted(edges)) == s.witness, f"{tag}: P edge multiset == witness")
        a.exact(score_from_derivation(prods, coef) == s.value,
                f"{tag}: score reconstructs through P")

    fcoef = {e: rng.uniform(-2, 2) for e in ground}
    a.approx(o.log_partition(fcoef).value, brute_log_partition(members, fcoef),
             f"{tag}: log_partition")
    if members and ground:
        got = o.marginals(fcoef).value
        want = brute_marginals(members, fcoef, ground)
        for e in ground:
            a.approx(got[e], want[e], f"{tag}: marginal {e}")

    # masks: a mixture of consistent and deliberately inconsistent
    for _ in range(3):
        k = rng.randint(0, min(2, len(ground)))
        forced = tuple(rng.sample(list(ground), k)) if ground else ()
        j = rng.randint(0, min(2, len(ground)))
        forbidden = tuple(rng.sample(list(ground), j)) if ground else ()
        dp = o.count(forced, forbidden)
        bf = len(enumerate_rna(seq, min_loop, canon, forced, forbidden))
        a.exact(dp.value == bf, f"{tag}: masked count {forced}/{forbidden}", dp.value, bf)
        feas = o.feasibility(forced, forbidden)
        a.exact(bool(feas) == (bf > 0), f"{tag}: feasibility agrees with enumeration")

    if canon is False and min_loop == 0 and n < len(MOTZKIN):
        a.exact(c.value == MOTZKIN[n], f"{tag}: Motzkin M_{n}", c.value, MOTZKIN[n])


def audit_dag_case(a: Audit, rng, n_layers, width):
    sizes = [1] + [rng.randint(1, width) for _ in range(n_layers - 1)]
    edges = [(l, u, v)
             for l in range(len(sizes) - 1)
             for u in range(sizes[l]) for v in range(sizes[l + 1])
             if rng.random() < 0.7]
    o = LayeredDAGPathOracle(sizes, edges)
    paths = enumerate_dag_paths(sizes, edges)
    tag = f"dag layers={sizes}"

    c = o.count()
    a.exact(c.value == len(paths), f"{tag}: count", c.value, len(paths))

    ground = o.ground_set()
    coef = {e: rng.randint(-9, 9) for e in ground}
    s = o.support(coef)
    exp, _ = brute_support(paths, coef)
    a.exact(s.value == exp, f"{tag}: support", s.value, exp)
    if s.termination_status != INFEASIBLE:
        a.exact(s.witness in set(paths), f"{tag}: witness is a real path")
        a.exact(derivation_edges(s.stats["productions"]) == list(s.witness),
                f"{tag}: P identity")

    fcoef = {e: rng.uniform(-2, 2) for e in ground}
    a.approx(o.log_partition(fcoef).value, brute_log_partition(paths, fcoef),
             f"{tag}: log_partition")
    if paths and ground:
        got = o.marginals(fcoef).value
        want = brute_marginals(paths, fcoef, ground)
        for e in ground:
            a.approx(got[e], want[e], f"{tag}: marginal {e}")

    for _ in range(2):
        k = rng.randint(0, min(2, len(ground)))
        forced = tuple(rng.sample(list(ground), k)) if ground else ()
        dp = o.count(forced)
        bf = len(enumerate_dag_paths(sizes, edges, forced))
        a.exact(dp.value == bf, f"{tag}: masked count {forced}", dp.value, bf)


def complexity_table(max_n, seed):
    """Measured chart values and materialised productions, not asserted ones."""
    rng = random.Random(seed)
    rows = []
    for n in range(1, max_n + 1):
        seq = random_sequence(rng, n)
        o = RNANonCrossingOracle(seq, min_loop=0, canonical_only=False)
        st = o.count().stats
        rows.append((n, st["chart_values"], st["materialized_productions"],
                     st["ground_set_size"]))
    return rows


def fit_exponent(rows, col):
    pts = [(math.log(n), math.log(v)) for n, *rest in rows
           for v in [rest[col]] if n >= 3 and v > 0]
    if len(pts) < 2:
        return float("nan")
    mx = sum(p[0] for p in pts) / len(pts)
    my = sum(p[1] for p in pts) / len(pts)
    num = sum((x - mx) * (y - my) for x, y in pts)
    den = sum((x - mx) ** 2 for x, _ in pts)
    return num / den if den else float("nan")


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m scmp.audit_oracles")
    ap.add_argument("--max-n", type=int, default=10)
    ap.add_argument("--cases", type=int, default=200)
    ap.add_argument("--seed", type=int, default=7027)
    args = ap.parse_args(argv)

    rng = random.Random(args.seed)
    a = Audit()
    t0 = time.time()

    print(f"scmp oracle audit  max-n={args.max_n} cases={args.cases} seed={args.seed}")
    print("independent enumerator: scmp.oracles.bruteforce (never calls the recurrence)\n")

    n_rna = n_dag = 0
    for i in range(args.cases):
        if i % 4 == 3:
            audit_dag_case(a, rng, rng.randint(2, 5), rng.randint(1, 3))
            n_dag += 1
        else:
            audit_rna_case(a, rng, rng.randint(0, args.max_n))
            n_rna += 1

    print(f"cases            : {n_rna} RNA, {n_dag} layered-DAG")
    print(f"exact  (integer) : {a.exact_checks} checks")
    print(f"float  (diagnostic): {a.float_checks} checks, tol {FLOAT_TOL:g}")

    print("\nMotzkin cross-check (OEIS A001006), min_loop=0 canonical_only=False:")
    ok_motz = True
    for n in range(0, min(args.max_n, len(MOTZKIN) - 1) + 1):
        o = RNANonCrossingOracle("A" * n, min_loop=0, canonical_only=False)
        got = o.count().value
        mark = "ok" if got == MOTZKIN[n] else "MISMATCH"
        ok_motz &= (got == MOTZKIN[n])
        print(f"  n={n:2d}  count={got:8d}  M_{n}={MOTZKIN[n]:8d}  {mark}")
    a.exact(ok_motz, "Motzkin sequence")

    print("\ncomplexity, measured on the densest family (min_loop=0, canonical_only=False):")
    rows = complexity_table(args.max_n, args.seed)
    print(f"  {'n':>3} {'chart values':>13} {'productions':>12} {'|E|':>6} "
          f"{'chart/n^2':>10} {'prod/n^3':>10}")
    for n, ch, pr, gs in rows:
        c2 = ch / n**2 if n else float("nan")
        p3 = pr / n**3 if n else float("nan")
        print(f"  {n:3d} {ch:13d} {pr:12d} {gs:6d} {c2:10.3f} {p3:10.3f}")
    # The closed forms are exact, so check them rather than trusting a log-log slope:
    # one chart cell per interval [i..j] including the empty ones, and one production
    # per (i, j, k) with i <= k <= j, plus one "unpaired" production per non-empty cell.
    closed_ok = True
    for n, ch, pr, _ in rows:
        want_ch = (n + 1) * (n + 2) // 2
        want_pr = (n + 2) * (n + 1) * n // 6
        if ch != want_ch or pr != want_pr:
            closed_ok = False
            print(f"  n={n}: chart {ch} vs {want_ch}, productions {pr} vs {want_pr}")
    a.exact(closed_ok, "closed-form chart/production counts")
    print(f"  closed form: chart values = (n+1)(n+2)/2 = Theta(n^2), "
          f"productions = C(n+2,3) = Theta(n^3)  -> {'verified exactly' if closed_ok else 'MISMATCH'}")
    e_chart = fit_exponent(rows, 0)
    e_prod = fit_exponent(rows, 1)
    print(f"  log-log fit over n=3..{args.max_n}: chart n^{e_chart:.2f}, "
          f"productions n^{e_prod:.2f} — biased low by lower-order terms at this range; "
          f"the closed forms above are the claim, not the fit")

    dt = time.time() - t0
    print(f"\nruntime: {dt:.2f}s")
    if a.failures:
        print(f"\nFAILURES: {len(a.failures)}")
        for kind, label, got, want in a.failures[:20]:
            print(f"  [{kind}] {label}: got={got!r} want={want!r}")
        if len(a.failures) > 20:
            print(f"  ... and {len(a.failures) - 20} more")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
