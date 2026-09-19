"""Data audit and manifest builder for the on-target siRNA efficacy study.

TASK FRAMING. On-target efficacy and ranking. No off-target or safety quantity is
computed, claimed, or implied anywhere in this package.

WHAT THIS CATCHES. Two failures that are easy to commit silently and expensive to find
later:

  1. ASSAY / CONSTRUCT IDENTITY. A redistributed annotation can disagree with the primary
     paper about what was actually measured. The audit compares them and reports the
     conflict rather than picking one.

  2. SILENT CONSTRUCT SUBSTITUTION. If the assay used a reporter plasmid, the sequence
     context around the target site is the reporter's, not the native transcript's.
     Substituting native flanking sequence produces a model conditioned on RNA that was
     never in the experiment. The audit refuses to mark such a context as resolved.

Run:  python -m scmp.rna.audit_data --config configs/rna_data.yaml
"""
from __future__ import annotations

import argparse, hashlib, json, os, sys, time
from collections import Counter

import yaml


def sha256_file(path, n=16):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()[:n]


class Audit:
    def __init__(self):
        self.findings = []
        self.checks = 0

    def note(self, severity, dataset, msg, **extra):
        self.checks += 1
        self.findings.append({"severity": severity, "dataset": dataset,
                              "finding": msg, **extra})

    def ok(self):
        self.checks += 1

    def blocking(self):
        return [f for f in self.findings if f["severity"] == "BLOCKING"]


# ---------------------------------------------------------------------- huesken
def audit_huesken(spec, a):
    man = {"dataset": "huesken2005", "role": spec["role"],
           "task": "on_target_efficacy", "version": "DSIR redistribution",
           "source_url": "https://biodev.cea.fr/DSIR/data/",
           "primary": spec["primary_source"], "hashes": {}, "fields": {},
           "exclusions": [], "provenance_conflicts": [],
           "context_resolution": {}}

    for path, want in spec["files"].items():
        if not os.path.exists(path):
            a.note("BLOCKING", "huesken2005", f"missing input {path}")
            continue
        got = sha256_file(path)
        man["hashes"][path] = got
        if got != want:
            a.note("BLOCKING", "huesken2005",
                   f"hash moved for {path}", want=want, got=got)
        else:
            a.ok()

    # -- records, guide length, unclipped labels, historical split --------
    tr = [l.split() for l in open("data/TrainAll2182.txt") if l.strip()]
    te = [l.split() for l in open("data/TestAll249.txt") if l.strip()]
    seqs = [r[0].upper().replace("T", "U") for r in tr + te]
    ys = [float(r[1]) for r in tr + te]
    exp = spec["expect"]

    man["n_records"] = len(seqs)
    man["n_unique_sequences"] = len(set(seqs))
    a.ok() if len(seqs) == exp["n_records"] else a.note(
        "BLOCKING", "huesken2005", "record count changed",
        want=exp["n_records"], got=len(seqs))
    a.ok() if len(set(seqs)) == len(seqs) else a.note(
        "WARN", "huesken2005", "duplicate guide sequences present",
        n_dup=len(seqs) - len(set(seqs)))

    lens = Counter(len(s) for s in seqs)
    man["guide_lengths"] = dict(lens)
    if set(lens) != {exp["guide_length"]}:
        a.note("BLOCKING", "huesken2005", "guide length changed",
               want=exp["guide_length"], got=sorted(lens))
    else:
        a.ok()

    man["label_min"], man["label_max"] = min(ys), max(ys)
    man["label_unclipped"] = max(ys) > 1.0
    if not man["label_unclipped"]:
        a.note("BLOCKING", "huesken2005",
               "labels appear clipped at 1.0; the published scale exceeds it")
    elif max(ys) < exp["label_max_at_least"]:
        a.note("WARN", "huesken2005", "label maximum lower than expected",
               got=max(ys), want_at_least=exp["label_max_at_least"])
    else:
        a.ok()

    split = [len(tr), len(te)]
    man["historical_split"] = split
    a.ok() if split == exp["historical_split"] else a.note(
        "BLOCKING", "huesken2005", "historical split changed",
        want=exp["historical_split"], got=split)

    # -- mapping coverage -------------------------------------------------
    ann = "data/Huesken_2431_annotated.tsv"
    gene_of = {}
    if os.path.exists(ann):
        import csv
        with open(ann) as fh:
            rd = csv.reader(fh, delimiter="\t")
            hdr = next(rd)
            gi = hdr.index("Gene") if "Gene" in hdr else None
            si = next((k for k, h in enumerate(hdr) if "21mer" in h.lower()), None)
            for row in rd:
                if gi is not None and si is not None and len(row) > max(gi, si):
                    gene_of[row[si].strip().upper().replace("T", "U")] = row[gi].strip()
    mapped = sum(1 for s in seqs if s in gene_of)
    man["mapping_coverage"] = {"annotated": mapped, "total": len(seqs),
                               "fraction": mapped / max(len(seqs), 1),
                               "n_genes": len(set(gene_of.values()))}
    a.ok() if mapped == len(seqs) else a.note(
        "WARN", "huesken2005", "annotation does not cover every guide",
        mapped=mapped, total=len(seqs))

    # -- assay / construct identity, primary vs redistributed annotation ---
    prim = spec["primary_source"]["assay_from_primary"]
    claim = spec["annotation_source"]["claims"]
    man["assay"] = {"construct": prim["construct"], "n_constructs": prim["n_constructs"],
                    "cell_line": prim["cell_line"], "timepoint_h": prim["timepoint_h"],
                    "readout": prim["readout"],
                    "normalisation": prim["normalisation"],
                    "evidence": prim["evidence"],
                    "primary_access": spec["primary_source"]["access"]}
    if claim.get("cell_line") != prim["cell_line"]:
        a.note("BLOCKING", "huesken2005",
               "assay identity conflict: cell line",
               primary=prim["cell_line"],
               redistributed_annotation=claim.get("cell_line"),
               action="use the primary source; the annotation's assay fields are unsafe")
        man["provenance_conflicts"].append(
            {"field": "cell_line", "primary": prim["cell_line"],
             "annotation": claim.get("cell_line")})
    if "luciferase" in str(claim.get("assay", "")).lower() and \
            "yfp" in prim["readout"].lower():
        a.note("BLOCKING", "huesken2005",
               "assay identity conflict: reporter protein",
               primary=prim["readout"], redistributed_annotation=claim.get("assay"))
        man["provenance_conflicts"].append(
            {"field": "reporter", "primary": prim["readout"],
             "annotation": claim.get("assay")})

    # -- context resolution ----------------------------------------------
    reporter = "reporter" in prim["construct"].lower()
    man["context_resolution"] = {
        "assayed_context": "reporter plasmid 3' UTR" if reporter else "native transcript",
        "native_accession": None, "species": "human (target 3' UTR fragments)",
        "strand": "sense of the reporter transcript", "isoform": None,
        "status": "UNRESOLVED" if reporter else "resolved",
        "reason": ("the assay used a reporter construct, so the flanking context is the "
                   "reporter's, not the native transcript's; the construct sequence is "
                   "not in hand") if reporter else "",
        "substitution_permitted": False}
    if reporter:
        a.note("BLOCKING", "huesken2005",
               "assayed context is a reporter construct; native transcript flanking "
               "sequence MUST NOT be substituted for it",
               affects="any model conditioned on sequence context around the target site",
               prior_use="scmp/experiments f9/f10 used GENCODE native windows")

    man["fields"] = {"sequence": "Antisense_21mer (guide, 21 nt)",
                     "chemistry": "unmodified RNA (no modification field present)",
                     "target": "Gene symbol via siRNAEfficacyDB join",
                     "assay": "YFP reporter, H1299, 48 h (primary)",
                     "condition": "not recorded per-row; single assay condition"}
    man["exclusions"] = [
        {"rule": "rows failing published-statistics integrity check", "n": 0},
        {"rule": "guides whose gene annotation is missing",
         "n": len(seqs) - mapped},
    ]
    man["licence"] = {"status": "unresolved",
                      "note": "DSIR redistribution terms not established"}
    return man


# ------------------------------------------------------------------ davis 2025
def audit_davis(spec, a):
    man = {"dataset": "davis2025", "role": spec["role"], "status": spec["status"],
           "task": "on_target_efficacy", "citation": spec["citation"],
           "doi": spec["doi"], "efficacy_table": spec["efficacy_table"],
           "fields": {"sequence": "siRNA sequence",
                      "chemistry": "chemical scaffold: " +
                                   ", ".join(spec["chemical_scaffolds"]),
                      "target": ", ".join(spec["target_genes"]),
                      "assay": ", ".join(spec["assay_types"]),
                      "condition": "percent target expression"},
           "geo": {"accession": spec["geo_accession"],
                   "content": spec["geo_content"],
                   "is_efficacy_labels": spec["geo_is_efficacy_labels"]},
           "hashes": {}, "exclusions": [], "context_resolution": {}}
    a.ok()
    if spec["geo_is_efficacy_labels"]:
        a.note("BLOCKING", "davis2025",
               "config asserts the GEO accession contains efficacy labels; it does not")
    else:
        a.note("INFO", "davis2025",
               f"{spec['geo_accession']} is {spec['geo_content']} and must never be "
               f"used as knockdown labels")
    a.note("PENDING", "davis2025",
           f"{spec['efficacy_table']} not yet downloaded; acquire before this becomes "
           f"the main therapeutic dataset", doi=spec["doi"])
    man["context_resolution"] = {
        "assayed_context": "BOTH: native (QuantiGene) and reporter (luciferase) arms",
        "status": "RESOLVABLE",
        "reason": "the study labels each measurement with its assay type, so native and "
                  "reporter rows can be separated rather than pooled",
        "substitution_permitted": False,
        "required_action": "split by assay_type; native rows may use the native "
                           "transcript context, reporter rows may not"}
    man["licence"] = {"status": "unresolved",
                      "note": "NAR is open access; supplementary terms not yet checked"}
    man["strengths"] = ["fully chemically modified siRNAs (therapeutically relevant)",
                        "explicit chemistry field",
                        "assay type recorded per measurement",
                        "four therapeutic targets"]
    return man


# ----------------------------------------------------------------- cmsirnadb
def audit_cmsirnadb(spec, a):
    man = {"dataset": "cmsirnadb", "role": spec["role"], "status": spec["status"],
           "citation": spec["citation"], "doi": spec["doi"], "url": spec["url"],
           "reported": spec["reported"], "unresolved": spec["unresolved"],
           "hashes": {}, "fields": {}, "exclusions": []}
    a.note("PENDING", "cmsirnadb", "optional dataset; not acquired")
    a.note("WARN", "cmsirnadb",
           f"{spec['reported']['sequences']} sequences drawn from "
           f"{spec['reported']['patents']} patents: entries within a patent family are "
           f"NOT independent and require a patent-family split")
    for u in spec["unresolved"]:
        a.note("WARN", "cmsirnadb", f"unresolved before use: {u}")
    man["licence"] = {"status": "unresolved"}
    man["blocking_before_use"] = [
        "records vs unique sequences must be distinguished",
        "chemical variants of the same guide must not be treated as independent rows",
        "condition fields must be recovered per record",
        "patent-family grouping must be recovered for splitting",
    ]
    return man


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m scmp.rna.audit_data")
    ap.add_argument("--config", default="configs/rna_data.yaml")
    args = ap.parse_args(argv)
    cfg = yaml.safe_load(open(args.config))
    a = Audit()

    print(f"RNA data audit: {cfg['audit']['name']}")
    print(f"task           : {cfg['audit']['task']}  "
          f"(no off-target or safety quantity is computed)\n")

    manifests = {}
    ds = cfg["datasets"]
    manifests["huesken2005"] = audit_huesken(ds["huesken2005"], a)
    manifests["davis2025"] = audit_davis(ds["davis2025"], a)
    manifests["cmsirnadb"] = audit_cmsirnadb(ds["cmsirnadb"], a)

    h = manifests["huesken2005"]
    print("huesken2005 (historical benchmark, preserved)")
    print(f"  records / unique      : {h['n_records']} / {h['n_unique_sequences']}")
    print(f"  guide lengths         : {h['guide_lengths']}")
    print(f"  labels (unclipped)    : [{h['label_min']:.3f}, {h['label_max']:.3f}]  "
          f"unclipped={h['label_unclipped']}")
    print(f"  historical split      : {h['historical_split']}")
    print(f"  mapping coverage      : {h['mapping_coverage']['annotated']}/"
          f"{h['mapping_coverage']['total']} "
          f"({h['mapping_coverage']['fraction']:.1%}), "
          f"{h['mapping_coverage']['n_genes']} genes")
    print(f"  assay (primary source): {h['assay']['readout']}, "
          f"{h['assay']['cell_line']}, {h['assay']['timepoint_h']} h, "
          f"{h['assay']['n_constructs']} constructs")
    print(f"  context               : {h['context_resolution']['assayed_context']} "
          f"-> {h['context_resolution']['status']}")

    print("\ndavis2025 (candidate main therapeutic dataset)")
    d = manifests["davis2025"]
    print(f"  status                : {d['status']}")
    print(f"  efficacy labels       : {d['efficacy_table']}")
    print(f"  GEO {d['geo']['accession']:<12}      : {d['geo']['content']} "
          f"-- efficacy labels: {d['geo']['is_efficacy_labels']}")
    print(f"  context               : {d['context_resolution']['status']} "
          f"({d['context_resolution']['assayed_context']})")

    print("\ncmsirnadb (optional)")
    c = manifests["cmsirnadb"]
    print(f"  reported              : {c['reported']}")
    print(f"  unresolved            : {len(c['unresolved'])} items incl. "
          f"records-vs-unique and patent-family grouping")

    print(f"\nfindings ({a.checks} checks)")
    for sev in ("BLOCKING", "WARN", "PENDING", "INFO"):
        rows = [f for f in a.findings if f["severity"] == sev]
        for f in rows:
            print(f"  [{sev:<8}] {f['dataset']:<12} {f['finding']}")
            for k, v in f.items():
                if k not in ("severity", "dataset", "finding"):
                    print(f"             {k}: {v}")

    out = cfg["audit"]["out_dir"]
    os.makedirs(out, exist_ok=True)
    path = os.path.join(out, f"manifest-{time.strftime('%Y%m%d-%H%M%S')}.json")
    payload = json.dumps({"config": args.config, "task": cfg["audit"]["task"],
                          "manifests": manifests, "findings": a.findings,
                          "checks": a.checks}, indent=2, default=str)
    with open(path, "w") as fh:
        fh.write(payload)
    with open(path + ".sha256", "w") as fh:
        fh.write(hashlib.sha256(payload.encode()).hexdigest() + "\n")
    print(f"\nmanifests written: {path}")

    blocking = a.blocking()
    if blocking:
        print(f"\n{len(blocking)} BLOCKING finding(s): the dataset may be used for the "
              f"historical benchmark, but a context-conditioned model may not be "
              f"trained on substituted native sequence.")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
