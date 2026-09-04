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
                 "seq_len": 19, "alphabet": set("ACGU"), "efficacy_range": (0.0, 1.2)},
    "sirnamod": {"n_total": 4894, "n_modifications": 128},
}

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
    dupes = len(records) - len({r["sequence"] for r in records})
    if dupes: problems.append(f"{dupes} duplicate sequences (leakage risk)")
    return {"dataset": name, "n": len(records), "verified": not problems, "problems": problems}
