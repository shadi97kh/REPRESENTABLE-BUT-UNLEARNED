"""Support-oracle interface shared by every graph family.

A family F is a set of subsets of a finite ground set E (the "edges"). An oracle answers,
under a forced/forbidden mask:

    support(c)        max_{x in F} <c, x>, with a witness
    count()           |F|, exact integer
    log_partition(c)  log sum_{x in F} exp(<c, x>)
    marginals(c)      per-edge Pr[e in x] under the Gibbs measure exp(<c,x>)/Z
    feasibility()     is F non-empty under the mask, with a reason if not
    condition()       a new oracle with the mask baked in

Every result carries the certificate fields required of this project: assumptions,
numerical status, family hash, bound and its direction, witness, termination status.

NUMERICAL STATUS is not decoration. `exact_integer` means the value was produced by
Python integer arithmetic with no rounding anywhere. `float_diagnostic` means the value
is a floating-point quantity reported for inspection and must not be used to support an
exactness claim. Callers are expected to branch on it.
"""
from __future__ import annotations

import hashlib
import json
import math
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Any, Iterable, Sequence

NEG_INF = float("-inf")

EXACT_INTEGER = "exact_integer"
FLOAT_DIAGNOSTIC = "float_diagnostic"
FLOAT_UNVERIFIED = "float_unverified"

PROVED = "proved"
INFEASIBLE = "infeasible"
EMPTY_GROUND_SET = "empty_ground_set"


# --------------------------------------------------------------------------- semirings
class Semiring:
    """(zero, one, add, mul). `select` is used only by argmax-style semirings."""
    name = "abstract"
    zero: Any = None
    one: Any = None
    tracks_argmax = False

    def add(self, a, b): raise NotImplementedError
    def mul(self, a, b): raise NotImplementedError
    def weight(self, coefficient): raise NotImplementedError


class CountingSemiring(Semiring):
    """Exact structure counting over Python ints. Coefficients are ignored."""
    name = "counting"
    zero = 0
    one = 1

    def add(self, a, b): return a + b
    def mul(self, a, b): return a * b
    def weight(self, coefficient): return 1


class MaxPlusSemiring(Semiring):
    """Maximisation. Exact whenever the coefficients are integers."""
    name = "max_plus"
    zero = NEG_INF
    one = 0
    tracks_argmax = True

    def add(self, a, b): return a if a >= b else b
    def mul(self, a, b):
        if a == NEG_INF or b == NEG_INF: return NEG_INF
        return a + b
    def weight(self, coefficient): return coefficient


class LogSumExpSemiring(Semiring):
    """log-sum-exp. Floating point by construction: never exact."""
    name = "log_sum_exp"
    zero = NEG_INF
    one = 0.0

    def add(self, a, b):
        if a == NEG_INF: return b
        if b == NEG_INF: return a
        hi, lo = (a, b) if a >= b else (b, a)
        return hi + math.log1p(math.exp(lo - hi))

    def mul(self, a, b):
        if a == NEG_INF or b == NEG_INF: return NEG_INF
        return a + b

    def weight(self, coefficient): return float(coefficient)


COUNTING = CountingSemiring()
MAX_PLUS = MaxPlusSemiring()
LOG_SUM_EXP = LogSumExpSemiring()


# ----------------------------------------------------------------------------- results
@dataclass(frozen=True)
class Feasibility:
    feasible: bool
    reason: str
    termination_status: str
    family_hash: str
    assumptions: tuple = ()

    def __bool__(self): return self.feasible


@dataclass(frozen=True)
class OracleResult:
    """One answer plus everything needed to audit it."""
    value: Any
    witness: Any
    numerical_status: str
    termination_status: str
    family_hash: str
    bound_direction: str            # "exact" | "upper" | "lower" | "none"
    assumptions: tuple = ()
    stats: dict = field(default_factory=dict)

    @property
    def is_exact(self) -> bool:
        return self.numerical_status == EXACT_INTEGER


def family_hash(payload: dict) -> str:
    """Stable hash of everything that defines the family, mask included."""
    blob = json.dumps(payload, sort_keys=True, default=str).encode()
    return hashlib.sha256(blob).hexdigest()[:16]


def normalise_mask(forced: Iterable | None, forbidden: Iterable | None):
    f = frozenset(tuple(e) for e in (forced or ()))
    b = frozenset(tuple(e) for e in (forbidden or ()))
    return f, b


# ------------------------------------------------------------------------------ oracle
class SupportOracle(ABC):
    """A conditionable support oracle over a finite family of edge subsets."""

    family_id: str = "abstract"

    # -- family description -------------------------------------------------
    @abstractmethod
    def ground_set(self) -> tuple:
        """Every edge that could appear in some member of the family."""

    @abstractmethod
    def describe(self) -> dict:
        """JSON-able description; feeds `family_hash`."""

    def assumptions(self) -> tuple:
        return ()

    # -- core queries -------------------------------------------------------
    @abstractmethod
    def feasibility(self, forced=(), forbidden=()) -> Feasibility: ...

    @abstractmethod
    def support(self, coefficients, forced=(), forbidden=()) -> OracleResult:
        """max <c, x> over the conditioned family, with a witness member."""

    @abstractmethod
    def count(self, forced=(), forbidden=()) -> OracleResult:
        """|F| under the mask. Exact integer."""

    @abstractmethod
    def log_partition(self, coefficients, forced=(), forbidden=()) -> OracleResult:
        """log sum exp <c, x>. Float diagnostic."""

    def marginals(self, coefficients, forced=(), forbidden=()) -> OracleResult:
        """Pr[e in x] for every e in the ground set, under exp(<c,x>)/Z.

        Computed as Z(forced + {e}) / Z, which is exact in the algebra and needs only
        `log_partition`. Cost is |E| + 1 partition evaluations; an inside-outside pass
        would be one, and is the obvious optimisation once |E| grows.
        """
        forced, forbidden = normalise_mask(forced, forbidden)
        base = self.log_partition(coefficients, forced, forbidden)
        h = base.family_hash
        if base.value == NEG_INF:
            return OracleResult({}, None, FLOAT_DIAGNOSTIC, INFEASIBLE, h, "none",
                                self.assumptions(), {"partition_calls": 1})
        out, calls = {}, 1
        for e in self.ground_set():
            if e in forbidden:
                out[e] = 0.0
                continue
            if e in forced:
                out[e] = 1.0
                continue
            num = self.log_partition(coefficients, forced | {e}, forbidden)
            calls += 1
            out[e] = 0.0 if num.value == NEG_INF else math.exp(num.value - base.value)
        return OracleResult(out, None, FLOAT_DIAGNOSTIC, PROVED, h, "exact",
                            self.assumptions(), {"partition_calls": calls})

    # -- conditioning -------------------------------------------------------
    def condition(self, forced=(), forbidden=()) -> "ConditionedOracle":
        return ConditionedOracle(self, *normalise_mask(forced, forbidden))

    # -- instrumentation ----------------------------------------------------
    def stats(self) -> dict:
        return {}


class ConditionedOracle(SupportOracle):
    """`base` with a mask permanently applied. Masks compose by union."""

    def __init__(self, base: SupportOracle, forced: frozenset, forbidden: frozenset):
        self.base = base
        self.forced = forced
        self.forbidden = forbidden
        self.family_id = f"{base.family_id}|conditioned"

    def _merge(self, forced, forbidden):
        f, b = normalise_mask(forced, forbidden)
        return self.forced | f, self.forbidden | b

    def ground_set(self): return self.base.ground_set()
    def assumptions(self): return self.base.assumptions()

    def describe(self):
        d = dict(self.base.describe())
        d["forced"] = sorted(self.forced)
        d["forbidden"] = sorted(self.forbidden)
        return d

    def feasibility(self, forced=(), forbidden=()):
        return self.base.feasibility(*self._merge(forced, forbidden))

    def support(self, coefficients, forced=(), forbidden=()):
        return self.base.support(coefficients, *self._merge(forced, forbidden))

    def count(self, forced=(), forbidden=()):
        return self.base.count(*self._merge(forced, forbidden))

    def log_partition(self, coefficients, forced=(), forbidden=()):
        return self.base.log_partition(coefficients, *self._merge(forced, forbidden))

    def marginals(self, coefficients, forced=(), forbidden=()):
        return self.base.marginals(coefficients, *self._merge(forced, forbidden))

    def condition(self, forced=(), forbidden=()):
        f, b = self._merge(forced, forbidden)
        return ConditionedOracle(self.base, f, b)

    def stats(self): return self.base.stats()


# ------------------------------------------------------------------------- derivations
@dataclass(frozen=True)
class Production:
    """One rule application in a derivation.

    `edge` is the production-to-edge map P: the edge this production introduces, or
    None for productions that introduce none. A derivation's edge multiset must equal
    the structure it derives, and the map must be injective on edge-introducing
    productions — that is what makes the grammar unambiguous.
    """
    kind: str
    span: tuple
    edge: tuple | None = None


def derivation_edges(productions: Sequence[Production]) -> list:
    return [p.edge for p in productions if p.edge is not None]


def score_from_derivation(productions: Sequence[Production], coefficients: dict):
    """Recompute <c, x> from the derivation alone, via P. Used to validate P."""
    return sum(coefficients.get(e, 0) for e in derivation_edges(productions))
