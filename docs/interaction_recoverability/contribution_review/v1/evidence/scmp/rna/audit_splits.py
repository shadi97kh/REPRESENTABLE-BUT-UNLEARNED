"""Cross-dataset duplication audit, then held-out split construction.

ORDER MATTERS. Duplication is audited first. A gene-disjoint split that still shares
near-identical guide sequences across the boundary leaks, and the leak is invisible in
the split statistics. So clusters are built before splits, and every proposed split is
reported with the leakage that remains after it.

SPLIT UNITS: gene, sequence cluster, study, patent family. A unit that does not exist in
the available data is reported as INFEASIBLE, not silently skipped -- a single-study
corpus cannot support a study-held-out claim, and saying so is the point of the audit.

Run:  python -m scmp.rna.audit_splits --config configs/rna_splits.yaml
"""
from __future__ import annotations

import argparse, csv, hashlib, json, os, random, sys, time
from collections import Counter, defaultdict

import yaml


def load_huesken():
    tr = [l.split() for l in open("data/TrainAll2182.txt") if l.strip()]
    te = [l.split() for l in open("data/TestAll249.txt") if l.strip()]
    rows = []
    for split, part in (("train", tr), ("test", te)):
        for r in part:
            rows.append({"seq": r[0].upper().replace("T", "U"), "y": float(r[1]),
                         "historical_split": split, "study": "huesken2005",
                         "patent_family": None})
    gene = {}
    ann = "data/Huesken_2431_annotated.tsv"
    if os.path.exists(ann):
        with open(ann) as fh:
            rd = csv.reader(fh, delimiter="\t")
            hdr = next(rd)
            gi = hdr.index("Gene")
            si = next(k for k, h in enumerate(hdr) if "21mer" in h.lower())
            ai = hdr.index("Accession_number") if "Accession_number" in hdr else None
            for row in rd:
                if len(row) > max(gi, si):
                    key = row[si].strip().upper().replace("T", "U")
                    gene[key] = (row[gi].strip(),
                                 row[ai].strip() if ai is not None else None)
    for r in rows:
        g = gene.get(r["seq"], (None, None))
        r["gene"], r["accession"] = g
    return rows


def kmer_clusters(seqs, k, ):
    """Single-linkage clusters over a shared k-mer. Union-find, no pairwise scan."""
    parent = list(range(len(seqs)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x

    def union(x, y):
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[rx] = ry

    index = defaultdict(list)
    for i, s in enumerate(seqs):
        for p in range(len(s) - k + 1):
            index[s[p:p + k]].append(i)
    for _, members in index.items():
        for m in members[1:]:
            union(members[0], m)
    return [find(i) for i in range(len(seqs))]


def leakage(rows, assign, unit_key):
    """Fraction of test rows sharing a cluster with a train row, after the split."""
    tr_units = {r[unit_key] for r, s in zip(rows, assign) if s == "train"}
    te = [r for r, s in zip(rows, assign) if s == "test"]
    if not te:
        return float("nan"), 0
    leaked = sum(1 for r in te if r[unit_key] in tr_units)
    return leaked / len(te), leaked


def group_split(rows, key, frac, seed):
    groups = defaultdict(list)
    for i, r in enumerate(rows):
        groups[r[key]].append(i)
    keys = sorted(groups, key=lambda g: (-len(groups[g]), str(g)))
    rng = random.Random(seed)
    rng.shuffle(keys)
    target = frac * len(rows)
    test_keys, n = set(), 0
    for g in keys:
        if n >= target:
            break
        test_keys.add(g); n += len(groups[g])
    assign = ["test" if rows[i][key] in test_keys else "train"
              for i in range(len(rows))]
    return assign, test_keys


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m scmp.rna.audit_splits")
    ap.add_argument("--config", default="configs/rna_splits.yaml")
    args = ap.parse_args(argv)
    cfg = yaml.safe_load(open(args.config))
    dup = cfg["duplication"]
    report = {"config": args.config, "duplication": {}, "splits": [],
              "cross_dataset": {}, "findings": []}

    def note(sev, msg, **kw):
        report["findings"].append({"severity": sev, "finding": msg, **kw})

    print(f"RNA split audit: {cfg['audit']['name']}")
    print(f"available: {cfg['datasets_available']}   "
          f"pending: {cfg['datasets_pending']}\n")

    rows = load_huesken()
    seqs = [r["seq"] for r in rows]

    # ---------------- duplication, before any split ----------------------
    print("duplication audit (runs before splits)")
    exact = Counter(seqs)
    n_exact_dup = sum(c - 1 for c in exact.values() if c > 1)
    print(f"  exact duplicate guides      : {n_exact_dup}")

    cl = kmer_clusters(seqs, dup["kmer_cluster"]["k"])
    for r, c in zip(rows, cl):
        r["kmer_cluster"] = c
    sizes = Counter(cl)
    multi = sum(1 for v in sizes.values() if v > 1)
    print(f"  {dup['kmer_cluster']['k']}-mer clusters            : {len(sizes)} "
          f"({multi} with >1 member, largest {max(sizes.values())})")

    s0, s1 = dup["seed_region"]["start"], dup["seed_region"]["end"]
    for r in rows:
        r["seed_region"] = r["seq"][s0:s1]
    seed_groups = Counter(r["seed_region"] for r in rows)
    print(f"  distinct seed regions (2-8) : {len(seed_groups)}, largest group "
          f"{max(seed_groups.values())}")

    genes = Counter(r["gene"] for r in rows)
    print(f"  genes                       : {len(genes)}")
    report["duplication"] = {
        "exact_duplicates": n_exact_dup, "n_kmer_clusters": len(sizes),
        "kmer_clusters_multi": multi, "largest_kmer_cluster": max(sizes.values()),
        "n_seed_regions": len(seed_groups),
        "largest_seed_group": max(seed_groups.values()), "n_genes": len(genes)}

    # leakage in the historical split, measured not assumed
    hist = [r["historical_split"] for r in rows]
    for unit in ("kmer_cluster", "gene", "seed_region"):
        frac, n = leakage(rows, hist, unit)
        print(f"  historical split leakage by {unit:<13}: {frac:.1%} of test rows "
              f"({n})")
        report["duplication"][f"historical_leak_{unit}"] = frac
        if frac > 0.5:
            note("WARN", f"historical split leaks heavily by {unit}",
                 fraction=round(frac, 4))

    # ---------------- cross-dataset duplication --------------------------
    print("\ncross-dataset duplication")
    for d in cfg["datasets_pending"]:
        print(f"  {d:<12}: PENDING -- not acquired, so no cross-dataset check is "
              f"possible")
        report["cross_dataset"][d] = "pending_acquisition"
        note("PENDING", f"cross-dataset duplication against {d} not yet checkable")
    note("BLOCKING", "no multi-dataset split may be built until the pending datasets "
                     "are acquired and de-duplicated against huesken2005")

    # ---------------- proposed splits ------------------------------------
    print(f"\nproposed held-out splits (target test fraction "
          f"{cfg['target_test_fraction']:.0%})")
    print(f"  {'split':>18} {'unit':>16} {'groups':>7} {'test rows':>10} "
          f"{'residual leak':>14}  status")
    for spec in cfg["splits"]:
        unit = {"gene": "gene", "sequence_cluster": "kmer_cluster",
                "study": "study", "patent_family": "patent_family"}[spec["name"]]
        present = [r[unit] for r in rows if r[unit] is not None]
        n_groups = len(set(present))
        if n_groups == 0:
            status, ntest, lk = "INFEASIBLE (unit absent)", 0, float("nan")
            note("INFEASIBLE", f"{spec['name']} split: unit '{unit}' absent from the "
                               f"available data")
        elif n_groups < 2:
            status, ntest, lk = f"INFEASIBLE (only {n_groups} group)", 0, float("nan")
            note("INFEASIBLE",
                 f"{spec['name']} split: a single {unit} cannot be held out; "
                 f"the corpus cannot support a {spec['name']}-generalisation claim",
                 n_groups=n_groups)
        else:
            assign, tk = group_split(rows, unit, cfg["target_test_fraction"],
                                     cfg["seed"])
            ntest = sum(1 for a in assign if a == "test")
            lk, _ = leakage(rows, assign, "kmer_cluster")
            status = "ok" if lk < 0.01 else f"leaks {lk:.1%} by sequence cluster"
            report["splits"].append({"name": spec["name"], "unit": unit,
                                     "n_groups": n_groups, "n_test": ntest,
                                     "residual_cluster_leak": lk,
                                     "test_groups": sorted(map(str, tk))[:10]})
            if lk >= 0.01:
                note("WARN", f"{spec['name']} split leaves sequence-cluster leakage",
                     fraction=round(lk, 4))
        print(f"  {spec['name']:>18} {unit:>16} {n_groups:7d} {ntest:10d} "
              f"{lk:14.2%}  {status}")

    print(f"\nfindings")
    for f in report["findings"]:
        print(f"  [{f['severity']:<10}] {f['finding']}")
        for k, v in f.items():
            if k not in ("severity", "finding"):
                print(f"               {k}: {v}")

    out = cfg["audit"]["out_dir"]
    os.makedirs(out, exist_ok=True)
    path = os.path.join(out, f"splits-{time.strftime('%Y%m%d-%H%M%S')}.json")
    payload = json.dumps(report, indent=2, default=str)
    with open(path, "w") as fh:
        fh.write(payload)
    with open(path + ".sha256", "w") as fh:
        fh.write(hashlib.sha256(payload.encode()).hexdigest() + "\n")
    print(f"\nwritten: {path}")
    blocking = [f for f in report["findings"] if f["severity"] == "BLOCKING"]
    return 1 if blocking else 0


if __name__ == "__main__":
    sys.exit(main())
