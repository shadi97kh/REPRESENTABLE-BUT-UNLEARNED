"""Saturate a secondary structure to a MAXIMAL valid one, and reason about the downward
closure such structures generate.

Why this matters. The edge lattice is the full power set of base pairs, which contains
enormous numbers of edge subsets that are not valid secondary structures: crossing pairs,
bases with two partners, hairpins too short. That is where f7 and f9 measured their slack.

A valid structure is downward closed under the monotone model in the following sense: if
S is a subset of T's pair set, monotonicity gives model(S) <= model(T). So a set of
MAXIMAL valid structures {T_i} bounds every structure contained in any of them, and the
bound is max_i model(T_i). Every pair in that union is realisable, so the bound cannot be
attained at an impossible configuration the way the lattice bound can.

Validity here is ViennaRNA's: canonical pairs only, at least MIN_LOOP unpaired bases inside
a hairpin, each base in at most one pair, and no crossing pairs.
"""
CANONICAL = {("A", "U"), ("U", "A"), ("G", "C"), ("C", "G"), ("G", "U"), ("U", "G")}
MIN_LOOP = 3


def loop_regions(pairs, n):
    """Group unpaired positions by their innermost enclosing pair.

    Two positions can only pair with each other if they sit in the same loop; a pair
    spanning two different loops would cross an existing pair.
    """
    partner = {}
    for (i, j) in pairs:
        partner[i] = j; partner[j] = i
    regions, stack = {}, []
    for i in range(n):
        if i in partner:
            if partner[i] > i: stack.append(i)          # opening
            else: stack.pop()                            # closing
        else:
            regions.setdefault(stack[-1] if stack else -1, []).append(i)
    return regions


def _fill(region, seq, out):
    """Greedy nested fill of one loop region, preferring the shortest span so that each
    added pair blocks as little of the region as possible."""
    if len(region) < 2: return
    best = None
    for a in range(len(region)):
        i = region[a]
        for b in range(a + 1, len(region)):
            j = region[b]
            if j - i <= MIN_LOOP: continue
            if (seq[i], seq[j]) in CANONICAL:
                if best is None or (j - i) < best[2]:
                    best = (i, j, j - i, a, b)
                break                                    # nearest partner for this i
    if best is None: return
    i, j, _, a, b = best
    out.append((i, j))
    _fill(region[a + 1:b], seq, out)                     # strictly inside the new pair
    _fill(region[:a] + region[b + 1:], seq, out)         # outside it, same loop


def saturate(pairs, seq, n, max_passes=8):
    """Add pairs until none can be added. Returns a frozenset: a MAXIMAL valid structure
    containing `pairs`. Deterministic, so identical inputs give identical outputs."""
    cur = set(pairs)
    for _ in range(max_passes):
        added = []
        for _, region in sorted(loop_regions(cur, n).items()):
            _fill(region, seq, added)
        if not added: break
        cur |= set(added)
    return frozenset(cur)


def is_valid(pairs, seq, n):
    seen = set()
    ps = sorted(pairs)
    for (i, j) in ps:
        if not (0 <= i < j < n): return False
        if j - i <= MIN_LOOP: return False
        if (seq[i], seq[j]) not in CANONICAL: return False
        if i in seen or j in seen: return False
        seen.add(i); seen.add(j)
    for x in range(len(ps)):
        i, j = ps[x]
        for y in range(x + 1, len(ps)):
            a, b = ps[y]
            if i < a < j < b or a < i < b < j: return False   # crossing
    return True


def is_maximal(pairs, seq, n):
    """No further pair can be added without breaking validity."""
    for _, region in sorted(loop_regions(set(pairs), n).items()):
        add = []
        _fill(region, seq, add)
        if add: return False
    return True
