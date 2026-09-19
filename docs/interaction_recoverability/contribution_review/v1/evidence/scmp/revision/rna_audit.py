"""P3: audit the locally available Huesken and Davis resources.

Three things this module refuses to do, because each was a real failure mode:

  * it never clips efficacy above 1.0. The primary normalisation puts the
    positive control at 90.0% inhibition and the least active siRNA at 0%, so
    values above 100% are attainable. 134 rows (5.5%) exceed it.
  * it never treats the redistributed annotation's assay columns as assay
    identity. They say HeLa/luciferase; the primary source says YFP reporter in
    H1299. The columns are recorded, and contradicted.
  * it never substitutes a plausible accession or invents a construct. Where the
    reporter context is unresolved it says so and the row is excluded.

WHAT IT DOES INVESTIGATE. Before declaring reporter context unrecoverable, it
maps each guide's reverse complement into GENCODE isoforms and takes, per gene,
the convex hull of the tiled site positions. Dense contiguous tiling is what an
inserted fragment produces, so that hull is a LOWER BOUND on the insert. A window
lying wholly inside the hull has native flanking sequence that is inside the
insert under one stated assumption -- that the insert is contiguous native
sequence. That is an inference from tiling, not a construct record, and it is
labelled as such wherever it is used.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time

import numpy as np

from . import FLOAT_DIAGNOSTIC, REAL_PROOF
from .provenance import sha256_file

RC = str.maketrans("ACGU", "UGCA")


def revcomp(s):
    return s.translate(RC)[::-1]


# ------------------------------------------------------------------ huesken
def audit_huesken(cfg):
    hc = cfg["datasets"]["huesken2005"]
    files = {k: {"path": v, "exists": os.path.exists(v),
                 "sha256": sha256_file(v) if os.path.exists(v) else None,
                 "bytes": os.path.getsize(v) if os.path.exists(v) else None}
             for k, v in hc["local_files"].items()}
    if not files["annotated"]["exists"]:
        return {"status": "ABSENT", "files": files}

    rows = list(csv.DictReader(open(hc["local_files"]["annotated"]), delimiter="\t"))
    pct = np.array([float(r["%Inhibition"]) for r in rows])
    eff = pct / 100.0
    ann_cell = sorted({r["Cell"] for r in rows})
    ann_tech = sorted({r["Technology"] for r in rows})
    prim = hc["assay_primary_source"]

    sys.path.insert(0, ".")
    from experiments._huesken import gene_map, load_verified
    recs, rep, ver = load_verified()
    gm = gene_map()
    genes = [gm.get(r["sequence"], "?") for r in recs]
    reff = np.array([r["efficacy"] for r in recs])

    return {
        "status": "PRESENT", "files": files,
        "n_annotated_rows": len(rows),
        "n_split_records": len(recs), "published_split": rep,
        "verification": {"verified": bool(ver["verified"]),
                         "problems": ver.get("problems", [])},
        "label": {
            "field": hc["label"]["field"], "scale": hc["label"]["scale"],
            "min": float(eff.min()), "max": float(eff.max()),
            "mean": float(eff.mean()),
            "n_above_1": int((eff > 1.0).sum()),
            "frac_above_1": float((eff > 1.0).mean()),
            "clipped": False,
            "split_records_match_annotated_range":
                bool(abs(reff.max() - eff.max()) < 1e-9),
            "note": "values above 1.0 are legitimate under the primary "
                    "normalisation and are preserved unclipped"},
        "assay_identity": {
            "primary_source": prim,
            "redistributed_annotation": {"cell_line": ann_cell,
                                         "technology": ann_tech},
            "agrees_on_cell_line": bool(ann_cell == [prim["cell_line"]]),
            "agrees_on_readout": False,
            "resolution": "the primary source governs; the annotation columns "
                          "are recorded and contradicted, never used as identity"},
        "chemistry": {"variants_present": 1, "description": hc["chemistry"],
                      "chemical_transfer_feasible": False},
        "genes": {"n_distinct": len({g for g in genes if g != "?"}),
                  "unmapped_rows": int(sum(1 for g in genes if g == "?"))},
        "evidence_kind": REAL_PROOF,
    }


# -------------------------------------------------------------------- davis
def audit_davis(cfg):
    dc = cfg["datasets"]["davis2025"]
    req = dc["required_file"]
    target = req["place_at"]
    present = os.path.exists(target)
    out = {"doi": dc["doi"], "licence": dc["licence"],
           "required_file": req["name"], "expected_path": target,
           "present": present,
           "do_not_substitute": dc["do_not_substitute"]}
    if not present:
        out["status"] = "MISSING"
        out["what_is_needed"] = {
            "file": req["name"],
            "contains": req["contains"],
            "route": req["publisher_route"],
            "place_at": target,
            "why_not_automated": "the PMC mirror serves a reCAPTCHA challenge. "
                                 "It was not bypassed and no challenge HTML was "
                                 "parsed. A human download is required.",
            "licence_permits_redistribution": True}
        out["blocked_capabilities"] = [
            "chemical scaffold covariates", "chemical-transfer evaluation",
            "native vs reporter assay separation", "dose/time covariates",
            "replicate structure", "cross-dataset de-duplication"]
    else:
        out["status"] = "PRESENT"
        out["sha256"] = sha256_file(target)
        out["bytes"] = os.path.getsize(target)
        out["next"] = ("schema audit: map guide/site, chemical scaffold, assay, "
                       "dose/time, replicate, target, accession, species, "
                       "outcome; separate reporter and native rows and their "
                       "incompatible label scales")
    return out


# ------------------------------------------------- insert / flank recovery
def audit_insert_recovery(cfg, max_genes=None):
    """Bound the inserted fragment from the tiled sites, before declaring loss."""
    ic = cfg["insert_recovery"]
    sys.path.insert(0, ".")
    from experiments._huesken import gene_map, load_verified
    from experiments.f9_target_context import GENCODE, load_transcripts

    if not os.path.exists(GENCODE):
        return {"status": "GENCODE_ABSENT", "path": GENCODE}
    recs, _, _ = load_verified()
    gm = gene_map()
    genes = sorted({gm.get(r["sequence"], "?") for r in recs} - {"?"})
    tx = load_transcripts(GENCODE, genes)

    spans, located, isoform_ambiguous = {}, 0, 0
    for r in recs:
        g = gm.get(r["sequence"], "?")
        if g not in tx:
            continue
        site = revcomp(r["sequence"].upper().replace("T", "U"))
        hits = []
        for k, iso in enumerate(tx[g]):
            p = iso.find(site)
            if p < 0:
                p = iso.find(site[2:])          # tolerate the 2-nt 3' overhang
            if p >= 0:
                hits.append((k, p, len(iso)))
        if not hits:
            continue
        located += 1
        if len({h[0] for h in hits}) > 1:
            isoform_ambiguous += 1              # recorded, never silently folded
        k, p, L = hits[0]
        d = spans.setdefault(g, {"pos": [], "len": L, "isoform_index": k,
                                 "n_isoforms": len(tx[g])})
        d["pos"].append((p, p + len(site)))

    per_gene, by_window = {}, {}
    for g, d in spans.items():
        lo = min(a for a, _ in d["pos"])
        hi = max(b for _, b in d["pos"])
        per_gene[g] = {"n_sites": len(d["pos"]), "transcript_len": d["len"],
                       "hull_start": lo, "hull_end": hi, "hull_nt": hi - lo,
                       "transcript_fraction": (hi - lo) / max(d["len"], 1),
                       "n_isoforms": d["n_isoforms"]}
    for W in ic["windows"]:
        inte = edge = 0
        for g, d in spans.items():
            lo = min(a for a, _ in d["pos"])
            hi = max(b for _, b in d["pos"])
            for (a, b) in d["pos"]:
                c = (a + b) // 2
                if c - W // 2 >= lo and c + W // 2 <= hi:
                    inte += 1
                else:
                    edge += 1
        by_window[W] = {"interior": inte, "edge": edge,
                        "interior_fraction": inte / max(inte + edge, 1),
                        "interior_status": ic["status_for_interior"],
                        "edge_status": ic["status_for_edge"]}

    return {"status": "INVESTIGATED",
            "method": ic["method"], "assumption": ic["assumption"],
            "genes_with_transcripts": len(tx), "genes_with_sites": len(spans),
            "sites_located": located, "sites_total": len(recs),
            "located_fraction": located / max(len(recs), 1),
            "isoform_ambiguous_sites": isoform_ambiguous,
            "per_gene": per_gene, "by_window": by_window,
            "verdict": ("reporter context is PARTIALLY recoverable: interior "
                        "windows lie inside a lower bound on the insert"),
            "caveats": [
                "the hull is a LOWER bound on the insert; the true insert may "
                "extend further, so edge windows are unresolved rather than "
                "proven to be vector sequence",
                "contiguity of the insert is inferred from dense tiling, not "
                "read from a construct record",
                "vector backbone sequence is not recovered at all",
                f"only {located}/{len(recs)} sites are locatable in GENCODE",
                "isoform-ambiguous sites are counted and never folded into a "
                "union"],
            "evidence_kind": FLOAT_DIAGNOSTIC}


# ------------------------------------------------------------------- splits
def audit_splits(cfg, insert):
    """Novel-sequence grouping, and an explicit verdict on chemical transfer."""
    sys.path.insert(0, ".")
    from experiments._huesken import gene_map, load_verified
    recs, _, _ = load_verified()
    gm = gene_map()

    def cluster(seq, k=13):
        s = seq.upper().replace("T", "U")
        return min(s[i:i + k] for i in range(max(len(s) - k + 1, 1)))

    groups = {}
    for r in recs:
        g = gm.get(r["sequence"], "?")
        key = (g, cluster(r["sequence"]))     # site proxy + sequence cluster
        groups.setdefault(key, []).append(r)

    sizes = np.array([len(v) for v in groups.values()])
    gene_counts = {}
    for r in recs:
        g = gm.get(r["sequence"], "?")
        gene_counts[g] = gene_counts.get(g, 0) + 1

    sc = cfg["splits"]
    return {
        "novel_sequence": {
            "group_by": sc["novel_sequence"]["group_by"],
            "n_groups": len(groups),
            "n_rows": len(recs),
            "median_group_size": float(np.median(sizes)),
            "max_group_size": int(sizes.max()),
            "singleton_groups": int((sizes == 1).sum()),
            "feasible": True,
            "note": sc["novel_sequence"]["note"]},
        "chemical_transfer": {
            "status": sc["chemical_transfer"]["status"],
            "reason": "Huesken contains a single chemistry (unmodified). A "
                      "transfer evaluation needs at least two scaffolds on the "
                      "same site, which only Davis S1 can supply.",
            "feasible": False},
        "target_wise": {
            "n_target_genes": len(gene_counts),
            "rows_per_gene": dict(sorted(gene_counts.items(),
                                         key=lambda kv: -kv[1])[:10]),
            "min_rows_per_gene": int(min(gene_counts.values())),
            "note": "results must be shown per target gene; a row-wise bootstrap "
                    "over this few genes does not license a population claim"},
        "forbidden": sc["forbidden"],
        "evidence_kind": REAL_PROOF}


# --------------------------------------------------- arm admissibility table
def arm_admissibility(huesken, davis, insert, splits):
    """Admissibility decided by RESOLVED INPUTS, never by dataset name."""
    have_seq = huesken.get("status") == "PRESENT"
    have_chem = huesken.get("chemistry", {}).get("chemical_transfer_feasible", False)
    interior = insert.get("by_window", {}).get(50, {}).get("interior", 0) \
        if insert.get("status") == "INVESTIGATED" else 0
    have_window_inferred = interior > 0
    have_window_recorded = davis.get("status") == "PRESENT"

    def row(arm, needs, ok, why):
        return {"arm": arm, "requires": needs, "admissible": bool(ok), "note": why}

    return [
        row("matched_sequence_chemistry", ["guide sequence"], have_seq,
            "guide sequence is fully resolved; chemistry is constant here"),
        row("target_sequence", ["target site sequence"], have_seq,
            "the site is the reverse complement of the guide, so it carries the "
            "same information; a null result against sequence-only is expected "
            "by construction and is not evidence about target context"),
        row("mfe_gnn", ["target window", "folding"], have_window_inferred,
            "admissible ONLY on interior windows under the stated insert "
            "assumption; results must be labelled inferred_insert_interior"),
        row("sampled_ensemble_gnn", ["target window", "prior", "sampler"],
            have_window_inferred,
            "same restriction; prior sensitivity and sampling error must be "
            "reported alongside"),
        row("anchor_A", ["target window"], have_window_inferred,
            "same restriction"),
        row("repaired_A_plus_B", ["target window", "A-initialised residual"],
            have_window_inferred,
            "same restriction; P1 showed the residual is not predictively "
            "useful against explicit-basis regression"),
        row("chemical_transfer", [">=2 chemical scaffolds on one site"],
            have_chem, "blocked: Huesken has a single chemistry"),
        row("native_vs_reporter", ["assay-typed rows"], have_window_recorded,
            "blocked: needs Davis S1 assay-type column"),
    ]


# --------------------------------------------------------------------- CLI
def main(argv=None):
    import yaml
    ap = argparse.ArgumentParser(
        prog="python -m scmp.revision.rna_audit",
        description="Audit local Huesken and Davis resources. Investigates "
                    "insert recovery before declaring reporter context lost.")
    ap.add_argument("--config", default="configs/revision_rna_data.yaml")
    ap.add_argument("--skip-insert", action="store_true",
                    help="skip the GENCODE mapping (slow)")
    ap.add_argument("--out-dir", default="runs/revision")
    a = ap.parse_args(argv)
    cfg = yaml.safe_load(open(a.config))
    t0 = time.time()

    hue = audit_huesken(cfg)
    dav = audit_davis(cfg)
    ins = {"status": "SKIPPED"} if a.skip_insert else audit_insert_recovery(cfg)
    spl = audit_splits(cfg, ins)
    arms = arm_admissibility(hue, dav, ins, spl)

    print("HUESKEN 2005")
    print(f"  rows {hue.get('n_annotated_rows')} annotated / "
          f"{hue.get('n_split_records')} in the published split; "
          f"verified={hue.get('verification', {}).get('verified')}")
    lb = hue.get("label", {})
    print(f"  label  range [{lb.get('min'):.4f}, {lb.get('max'):.4f}]  "
          f"{lb.get('n_above_1')} rows above 1.0 ({lb.get('frac_above_1'):.2%})  "
          f"clipped={lb.get('clipped')}")
    ai = hue.get("assay_identity", {})
    print(f"  assay  PRIMARY {ai['primary_source']['construct']}, "
          f"{ai['primary_source']['cell_line']}, "
          f"{ai['primary_source']['timepoint_h']}h")
    print(f"         ANNOTATION claims {ai['redistributed_annotation']['cell_line']} / "
          f"{ai['redistributed_annotation']['technology']} -> CONTRADICTED")
    print(f"  genes  {hue.get('genes', {}).get('n_distinct')}   "
          f"chemistry variants {hue.get('chemistry', {}).get('variants_present')}")

    print(f"\nDAVIS 2025   status {dav['status']}")
    if dav["status"] == "MISSING":
        w = dav["what_is_needed"]
        print(f"  required   {w['file']}  ({w['contains']})")
        print(f"  route      {w['route']}")
        print(f"  place at   {w['place_at']}")
        print(f"  not automated: {w['why_not_automated']}")
        print(f"  blocks     {', '.join(dav['blocked_capabilities'])}")

    if ins.get("status") == "INVESTIGATED":
        print(f"\nINSERT RECOVERY  ({ins['sites_located']}/{ins['sites_total']} "
              f"sites located, {ins['located_fraction']:.1%}; "
              f"{ins['genes_with_sites']} genes)")
        print(f"  {'window':>7}{'interior':>10}{'edge':>7}{'interior%':>11}")
        for W, v in ins["by_window"].items():
            print(f"  {W:>7}{v['interior']:>10}{v['edge']:>7}"
                  f"{v['interior_fraction']:>10.1%}")
        print(f"  verdict: {ins['verdict']}")
        print(f"  isoform-ambiguous sites recorded: "
              f"{ins['isoform_ambiguous_sites']} (never folded into a union)")

    print(f"\nSPLITS")
    ns = spl["novel_sequence"]
    print(f"  novel-sequence groups {ns['n_groups']} over {ns['n_rows']} rows "
          f"(median size {ns['median_group_size']:.0f}, "
          f"{ns['singleton_groups']} singletons)  feasible={ns['feasible']}")
    print(f"  chemical transfer: {spl['chemical_transfer']['status']} -- "
          f"{spl['chemical_transfer']['reason']}")
    print(f"  target genes {spl['target_wise']['n_target_genes']}, "
          f"min rows/gene {spl['target_wise']['min_rows_per_gene']}")

    print(f"\nARM ADMISSIBILITY (by resolved inputs, not dataset name)")
    for r in arms:
        print(f"  {'YES' if r['admissible'] else ' NO':>3}  {r['arm']:<28}"
              f"needs {', '.join(r['requires'])}")

    res = {"kind": "revision_rna_audit", "config": a.config,
           "config_sha256": sha256_file(a.config),
           "huesken": hue, "davis": dav, "insert_recovery": ins,
           "splits": spl, "arm_admissibility": arms,
           "metrics_predeclared": cfg["metrics"],
           "f9_quarantine": cfg["f9_native_windows"],
           "runtime_s": round(time.time() - t0, 2)}
    os.makedirs(a.out_dir, exist_ok=True)
    path = os.path.join(a.out_dir,
                        f"rna_audit-{time.strftime('%Y%m%d-%H%M%S')}.json")
    with open(path, "w") as fh:
        json.dump(res, fh, indent=1, default=str)
    with open(path + ".sha256", "w") as fh:
        fh.write(sha256_file(path) + "  " + os.path.basename(path) + "\n")
    print(f"\nwritten {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
