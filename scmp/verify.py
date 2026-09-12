"""Run refinement over a config and write one replayable certificate per case.

Run:  python -m scmp.verify --config configs/tiny_complete.yaml
"""
from __future__ import annotations

import argparse, itertools, json, os, sys, time
import torch, yaml

from .check_certificates import rebuild
from .oracles.bruteforce import enumerate_rna
from .refine import refine


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m scmp.verify")
    ap.add_argument("--config", default="configs/tiny_complete.yaml")
    a = ap.parse_args(argv)
    cfg = yaml.safe_load(open(a.config))
    v = cfg["verify"]
    out = v["out_dir"]
    os.makedirs(out, exist_ok=True)
    for old in os.listdir(out):
        if old.endswith(".json"):
            os.remove(os.path.join(out, old))

    t0 = time.time()
    grid = list(itertools.product(cfg["models"]["d_hid"], cfg["models"]["gamma"],
                                 cfg["models"]["init_seed"],
                                 cfg["models"]["perturb_scale"],
                                 cfg["models"]["perturb_seed"]))
    n_proved = n_unresolved = 0
    written = 0
    print(f"verify  config={a.config}  out={out}")
    for fam in cfg["families"]:
        members = enumerate_rna(fam["seq"], fam["min_loop"], fam["canonical_only"])
        for d_hid, gamma, iseed, scale, pseed in grid:
            for tag, nodes in (("full", v["max_nodes"]),
                               ("tight", cfg["stress"]["tight_node_budget"])):
                meta = {"seq": fam["seq"], "min_loop": fam["min_loop"],
                        "canonical_only": fam["canonical_only"], "d_hid": d_hid,
                        "gamma": gamma, "init_seed": iseed,
                        "perturb_scale": scale, "perturb_seed": pseed,
                        "max_nodes": nodes, "tag": tag}
                o, m, X, B = rebuild(meta)
                r = refine(m, X, o, B, max_nodes=nodes,
                           time_limit=v["time_limit"], tol=v["tol"],
                           fairness_period=v["fairness_period"],
                           reference_mode=v["reference_mode"])
                n_proved += r.status == "proved"
                n_unresolved += r.status == "unresolved"
                name = (f"{fam['seq']}_h{d_hid}_g{gamma}_s{iseed}"
                        f"_p{scale}_{tag}.json")
                with open(os.path.join(out, name), "w") as fh:
                    json.dump({"meta": meta, "certificate": r.certificate()}, fh,
                              indent=1)
                written += 1
        print(f"  {fam['seq']:14s} |F|={len(members):4d}  "
              f"cases={len(grid)*2}")
    print(f"\ncertificates written : {written} -> {out}")
    print(f"proved               : {n_proved}")
    print(f"unresolved (budget)  : {n_unresolved}")
    print(f"numerical status     : inherited from the bounds (see docs/bounds_proof.md)")
    print(f"runtime              : {time.time()-t0:.1f}s")
    print(f"\nnow check them independently:")
    print(f"  python -m scmp.check_certificates --run {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
