"""Shared Huesken loading and graph construction, so every experiment that touches the
dataset uses the same verified path and the same published split."""
import os, sys
import numpy as np
from certmp import data as cdata
from certmp.ensemble import bpp_lattice
from certmp.reach import build_A

TRAIN_PATH = os.environ.get("HUESKEN_TRAIN", "data/TrainAll2182.txt")
TEST_PATH  = os.environ.get("HUESKEN_TEST",  "data/TestAll249.txt")
ANNOT_PATH = os.environ.get("HUESKEN_ANNOT", "data/Huesken_2431_annotated.tsv")
VAL_SEED   = 20260904


def pairs_of(structure):
    stack, pairs = [], set()
    for i, c in enumerate(structure):
        if c == "(": stack.append(i)
        elif c == ")": pairs.add((stack.pop(), i))
    return pairs


def load_verified():
    """Load and REJECT unless the published statistics match. Exits non-zero if absent."""
    for p in (TRAIN_PATH, TEST_PATH):
        if not os.path.exists(p):
            print(f"DATASET ABSENT: {p}. Run `make data`. Nothing is synthesised.")
            sys.exit(2)
    tr, te, rep = cdata.load_dsir_split(TRAIN_PATH, TEST_PATH)
    records = tr + te
    ver = cdata.verify("huesken", records)
    if not ver["verified"]:
        print(f"REJECTED: {ver['problems']}")
        sys.exit(3)
    return records, rep, ver


def gene_map():
    """21-mer -> gene symbol, from the siRNAEfficacyDB join. Empty dict if absent."""
    if not os.path.exists(ANNOT_PATH): return {}
    import csv
    out = {}
    with open(ANNOT_PATH) as fh:
        rd = csv.reader(fh, delimiter="\t"); hdr = next(rd)
        gi = next((i for i, h in enumerate(hdr) if h.lower() == "gene"), None)
        si = next((i for i, h in enumerate(hdr) if "21mer" in h.lower()), None)
        if gi is None or si is None: return {}
        for row in rd:
            if len(row) > max(gi, si):
                out[row[si].strip().upper().replace("T", "U")] = row[gi].strip()
    return out


def build_graphs(seqs, feature_fn):
    """MFE-conditioned graph per sequence, which is what current pipelines use."""
    out = []
    for s in seqs:
        lat = bpp_lattice(s)
        n = lat["n"]
        backbone = [(i, i + 1) for i in range(n - 1)]
        A_mfe = build_A(n, backbone + sorted(pairs_of(lat["mfe_structure"])), [], [])
        out.append(dict(seq=s, n=n, X=feature_fn(s), A_mfe=A_mfe,
                        mfe_structure=lat["mfe_structure"]))
    return out


def published_split(n_total, n_train, val_size=200, seed=VAL_SEED):
    idx_te = np.arange(n_train, n_total)
    shuf = np.random.RandomState(seed).permutation(np.arange(n_train))
    return shuf[val_size:], shuf[:val_size], idx_te
