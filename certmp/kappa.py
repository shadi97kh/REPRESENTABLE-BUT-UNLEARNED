"""kappa_L: the worst-case relaxation gap as a purely combinatorial quantity.

For a certified depth-L sum-aggregation MPNN (affine on the non-negative orthant, see
the affinity result), with uniform features and no bias, f(A) is proportional to the
length-L walk count W_L(A) = 1' A^L 1. Hence for a downward-closed feasible family F with
lattice top TOP,

    kappa_L(F) = W_L(TOP) / max_{G in F} W_L(G)

equals the relaxation gap exactly, and biases only add lower-order terms, so kappa_L is a
TIGHT UPPER BOUND on the gap for the whole certified class. It needs no network, no
weights and no training: it is a ratio of walk counts.

Everything here is additive. Nothing in this module imports or modifies models.py
internals; models are used only through their public .forward(X, A).
"""
import itertools
import numpy as np

def adjacency(n, edges, self_loops=True):
    A = np.eye(n) if self_loops else np.zeros((n, n))
    for (i, j) in edges:
        A[i, j] = A[j, i] = 1.0
    return A

def walk_count(A, L):
    """W_L(A) = 1' A^L 1, the number of length-L walks (self-loops included).

    Computed as L matrix-VECTOR products rather than a matrix power: 1'A^L 1 = 1'(A^L 1),
    which is O(L n^2) instead of O(L n^3) and identical to floating-point associativity.
    This matters because c8 evaluates it once per sampled structure.
    """
    A = np.asarray(A, float)
    v = np.ones(A.shape[0])
    for _ in range(L):
        v = A @ v
    return float(v.sum())

def kappa_L(n, top_edges, feasible_edge_sets, L):
    """kappa_L = W_L(top) / max over the feasible family. Returns the ratio and the
    argmax walk count so callers can report both."""
    wt = walk_count(adjacency(n, top_edges), L)
    ws = [walk_count(adjacency(n, g), L) for g in feasible_edge_sets]
    best = max(ws)
    return wt / best, wt, best

def kappa_from_adjacency(A_top, A_feasible, L):
    """Same, when structures are already adjacency matrices (e.g. sampled RNA structures)."""
    wt = walk_count(A_top, L)
    ws = [walk_count(A, L) for A in A_feasible]
    best = max(ws)
    return wt / best, wt, best

def measured_gap(model, X, A_top, A_feasible):
    """Empirical gap f(top) / max_{G in F} f(G) for any model exposing .forward(X, A)."""
    top = model.forward(X, A_top)
    best = max(model.forward(X, A) for A in A_feasible)
    return top / best, top, best

def enumerate_family(n, candidate_edges, predicate, max_edges=None):
    """Exhaustive enumeration of a downward-closed family over a candidate edge set.
    predicate(edges, n) -> bool. Exact, so max over F is real rather than sampled.
    Guard: refuses beyond 2^20 subsets."""
    k = len(candidate_edges)
    if k > 20:
        raise ValueError(f"|candidate|={k}: 2^{k} subsets is too many to enumerate")
    out = []
    top = max_edges if max_edges is not None else k
    for r in range(top + 1):
        for s in itertools.combinations(candidate_edges, r):
            e = list(s)
            if predicate(e, n):
                out.append(e)
    return out

# ---- standard downward-closed families, ordered by exclusion strength ----
def _degrees(edges):
    d = {}
    for (i, j) in edges:
        d[i] = d.get(i, 0) + 1
        d[j] = d.get(j, 0) + 1
    return d

def fam_matching(edges, n):  return all(v <= 1 for v in _degrees(edges).values())
def fam_deg2(edges, n):      return all(v <= 2 for v in _degrees(edges).values())
def fam_deg3(edges, n):      return all(v <= 3 for v in _degrees(edges).values())
def fam_join_closed(edges, n): return True          # all subsets: no exclusion at all

def fam_forest(edges, n):
    par = list(range(n))
    def find(a):
        while par[a] != a:
            par[a] = par[par[a]]; a = par[a]
        return a
    for (i, j) in edges:
        ri, rj = find(i), find(j)
        if ri == rj: return False
        par[ri] = rj
    return True

FAMILIES = [
    ("matching (deg<=1)",      fam_matching),
    ("deg<=2",                 fam_deg2),
    ("deg<=3",                 fam_deg3),
    ("forest (acyclic)",       fam_forest),
    ("all subsets (join-closed)", fam_join_closed),
]

def structure_to_edges(dotbracket):
    """Dot-bracket secondary structure -> 0-indexed base-pair edge list."""
    stack, edges = [], []
    for i, c in enumerate(dotbracket):
        if c == "(": stack.append(i)
        elif c == ")" and stack: edges.append((stack.pop(), i))
    return edges
