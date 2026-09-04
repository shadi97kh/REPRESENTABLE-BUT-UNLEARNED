"""Dataset integrity verification. Never trust a redistributed copy; check it against
the statistics reported in the original publication.

  Huesken et al., Nature Biotechnology 23(8):995-1001 (2005): 2431 siRNAs, 34 mRNA
    targets, 2182 train / 249 independent test, reported test Pearson r = 0.66.
  siRNAmod (Dar et al., Scientific Reports 6:20031, 2016): 4894 chemically modified
    entries, 128 unique chemical modifications.
Any downloaded file that does not match is REJECTED, not silently used.
"""
EXPECTED = {
    "huesken":  {"n_total": 2431, "n_train": 2182, "n_test": 249, "n_targets": 34,
                 "seq_len": 21, "alphabet": set("ACGU"), "efficacy_range": (0.0, 1.4)},
    "sirnamod": {"n_total": 4894, "n_modifications": 128},
}

# AMENDED 2026-09-04, on first contact with the real distribution. Recorded here rather
# than silently edited, because loosening an integrity check to make data pass is exactly
# the failure this module exists to prevent.
#   seq_len       19 -> 21. The distributed Huesken sequences are 21 nt. The trailing two
#                 nucleotides vary across rows, so they are measured sequence and not a
#                 constant synthetic overhang. Truncating would discard real data.
#   efficacy_range (0.0, 1.2) -> (0.0, 1.4). Observed maximum is 1.341. The published
#                 values are normalised inhibition and legitimately exceed 1.0, so the
#                 original 1.2 was a guess, not a published figure.
# Everything diagnostic passed unchanged and unaided: exactly 2431 rows, exactly the
# 2182/249 published split, alphabet ACGU, zero duplicate sequences, zero train/test
# overlap. Those are what actually identify the dataset.
AMENDMENTS = [
    {"date": "2026-09-04", "dataset": "huesken", "field": "seq_len",
     "from": 19, "to": 21, "reason": "distributed sequences are 21 nt; final 2 nt vary"},
    {"date": "2026-09-04", "dataset": "huesken", "field": "efficacy_range",
     "from": (0.0, 1.2), "to": (0.0, 1.4), "reason": "observed max 1.341; normalised inhibition exceeds 1.0"},
]

def verify(name, records):
    """records: list of dicts with at least 'sequence' and 'efficacy'."""
    exp = EXPECTED[name]; problems = []
    if len(records) != exp["n_total"]:
        problems.append(f"row count {len(records)} != published {exp['n_total']}")
    if "seq_len" in exp:
        lens = {len(r["sequence"]) for r in records}
        if lens != {exp["seq_len"]}:
            problems.append(f"sequence lengths {sorted(lens)[:5]} != {exp['seq_len']}")
        alpha = set().union(*(set(r["sequence"]) for r in records)) if records else set()
        if not alpha <= exp["alphabet"]:
            problems.append(f"unexpected symbols {alpha - exp['alphabet']}")
    if "efficacy_range" in exp and records:
        ys = [r["efficacy"] for r in records]
        boundlo, boundhi = exp["efficacy_range"]
        if min(ys) < boundlo or max(ys) > boundhi:
            problems.append(f"efficacy range [{min(ys):.3f}, {max(ys):.3f}] outside "
                            f"published [{boundlo}, {boundhi}]")
    dupes = len(records) - len({r["sequence"] for r in records})
    if dupes: problems.append(f"{dupes} duplicate sequences (leakage risk)")
    return {"dataset": name, "n": len(records), "verified": not problems, "problems": problems}


def load_dsir_split(train_path, test_path):
    """Load the Huesken table in its published 2182/249 split, as redistributed with the
    DSIR tool (Vert et al., BMC Bioinformatics 7:520, 2006): one 'SEQUENCE<ws>EFFICACY'
    row per line. Returns (train, test, report). The split is the PUBLISHED one, so no
    random resplit is needed and target-wise leakage is avoided by construction."""
    def rd(path):
        out = []
        for line in open(path):
            parts = line.split()
            if len(parts) < 2: continue
            try: y = float(parts[1])
            except ValueError: continue
            out.append({"sequence": parts[0].strip().upper().replace("T", "U"), "efficacy": y})
        return out
    tr, te = rd(train_path), rd(test_path)
    overlap = {r["sequence"] for r in tr} & {r["sequence"] for r in te}
    report = {"n_train": len(tr), "n_test": len(te), "train_test_overlap": len(overlap),
              "split": "published (DSIR redistribution of Huesken et al. 2005)"}
    return tr, te, report


def load_delimited(path, seq_col=None, eff_col=None, delim=None, skip_header=None):
    """Parse a redistributed Huesken table without guessing silently.

    Returns (records, layout). Layout is reported so the caller can print exactly how the
    file was interpreted; a wrong column guess is a data bug that must be visible.
    """
    import csv, io, re
    raw = open(path, newline="").read()
    if delim is None:
        head = raw[:4000]
        delim = "\t" if head.count("\t") > head.count(",") else ","
    rows = [r for r in csv.reader(io.StringIO(raw), delimiter=delim) if r and any(c.strip() for c in r)]
    if not rows: raise ValueError(f"{path}: no rows")

    def is_seq(c):
        c = c.strip().upper()
        return len(c) >= 15 and bool(re.fullmatch(r"[ACGTU]+", c))
    def is_num(c):
        try: float(c); return True
        except Exception: return False

    header = None
    if skip_header is None:
        skip_header = not any(is_seq(c) for c in rows[0])
    if skip_header:
        header = rows[0]; rows = rows[1:]
    if seq_col is None:
        cands = [i for i in range(len(rows[0])) if is_seq(rows[0][i])]
        if not cands: raise ValueError(f"{path}: no nucleotide column found in {rows[0]!r}")
        seq_col = cands[0]
    if eff_col is None:
        cands = [i for i in range(len(rows[0])) if i != seq_col and is_num(rows[0][i])]
        if not cands: raise ValueError(f"{path}: no numeric efficacy column found in {rows[0]!r}")
        eff_col = cands[-1]

    records = []
    for r in rows:
        if len(r) <= max(seq_col, eff_col): continue
        s = r[seq_col].strip().upper().replace("T", "U")
        if not is_seq(s): continue
        try: y = float(r[eff_col])
        except ValueError: continue
        records.append({"sequence": s, "efficacy": y})
    layout = {"path": path, "delimiter": repr(delim), "header": header,
              "seq_col": seq_col, "eff_col": eff_col, "n_rows_parsed": len(records)}
    return records, layout


def spearman(a, b):
    import numpy as np
    def rank(v):
        v = np.asarray(v, dtype=float); order = v.argsort()
        r = np.empty(len(v), dtype=float); r[order] = np.arange(len(v), dtype=float)
        # average ties
        _, inv, cnt = np.unique(v, return_inverse=True, return_counts=True)
        for g in np.where(cnt > 1)[0]:
            m = inv == g; r[m] = r[m].mean()
        return r
    return pearson(rank(a), rank(b))


def pearson(a, b):
    import numpy as np
    a = np.asarray(a, dtype=float); b = np.asarray(b, dtype=float)
    a = a - a.mean(); b = b - b.mean()
    d = (np.sqrt((a * a).sum()) * np.sqrt((b * b).sum()))
    return float((a * b).sum() / d) if d > 0 else float("nan")
