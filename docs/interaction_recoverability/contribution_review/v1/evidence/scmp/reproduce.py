"""Reproduce the winning development configuration under verified pins.

WHAT THIS DOES AND DOES NOT DO. It verifies that the *current* interpreter matches the
recorded pins and that every tracked source file and data input still hashes to what the
environment manifest recorded, then reruns the winning configuration in a clean
subprocess and compares the numbers against the recorded ones. It does NOT provision a
fresh virtualenv or container; that is a real limitation and is reported as one, not
glossed.

Output is an IMMUTABLE result file: the payload is written, then hashed, and the hash is
stored beside it. `export_results` refuses any file whose content hash does not match.

Run:  python -m scmp.reproduce --config configs/reproduction.yaml
"""
from __future__ import annotations

import argparse, hashlib, importlib.metadata as md, json, os, re, subprocess, sys, time
import yaml


def sha256_file(path, n=16):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()[:n]


def check_environment(cfg, report):
    env = cfg["environment"]
    import numpy, torch, RNA
    got = {"python": ".".join(sys.version.split()[0].split(".")[:2]),
           "numpy": numpy.__version__, "torch": torch.__version__.split("+")[0]}
    try:
        got["ViennaRNA"] = md.version("ViennaRNA")
    except Exception:
        got["ViennaRNA"] = "unknown"
    for k, want in env["require"].items():
        ok = got.get(k, "").startswith(want)
        report["environment"].append({"item": k, "want": want, "got": got.get(k),
                                      "ok": ok})
    report["environment"].append(
        {"item": "RNA.__version__", "want": env["rna_module_version"],
         "got": RNA.__version__, "ok": RNA.__version__ == env["rna_module_version"],
         "note": "module attribute differs from the distribution version by design"})
    lock = env["lock"]
    report["environment"].append({"item": "lockfile", "want": "present", "got": lock,
                                  "ok": os.path.exists(lock)})


def check_source_hashes(cfg, report):
    """Every hash the environment manifest recorded must still hold."""
    manifest = cfg["environment"]["manifest"]
    text = open(manifest).read()
    rows = re.findall(r"^\| `([^`]+)` \| `([0-9a-f]{12})` \|$", text, re.M)
    for path, want in rows:
        if not os.path.exists(path):
            report["source"].append({"path": path, "want": want, "got": None,
                                     "ok": False})
            continue
        got = sha256_file(path, 12)
        report["source"].append({"path": path, "want": want, "got": got,
                                 "ok": got == want})


def check_data_hashes(cfg, report):
    for path, want in (cfg.get("data_hashes") or {}).items():
        if not os.path.exists(path):
            report["data"].append({"path": path, "want": want, "got": None, "ok": False,
                                   "note": "gitignored third-party input; fetch first"})
            continue
        got = sha256_file(path, 16)
        report["data"].append({"path": path, "want": want, "got": got,
                               "ok": got == want})


NUM = r"([-+]?\d+\.?\d*)"


def rerun(cfg, report):
    wc = cfg["winning_config"]
    env = dict(os.environ)
    env["PYTHONPATH"] = "."          # clean import path, nothing inherited
    env.pop("CUDA_VISIBLE_DEVICES", None)
    t0 = time.time()
    proc = subprocess.run([str(x) for x in wc["command"]], capture_output=True,
                          text=True, env=env, cwd=".")
    out = proc.stdout
    report["rerun"] = {"returncode": proc.returncode, "seconds": time.time() - t0,
                       "command": [str(x) for x in wc["command"]],
                       "stdout_tail": out[-2000:]}
    got = {}
    m = re.search(r"cases audited\s*:\s*(\d+)", out)
    if m: got["cases"] = int(m.group(1))
    m = re.search(r"tighter than baseline\s*:\s*(\d+)", out)
    if m: got["tighter_than_baseline"] = int(m.group(1))
    m = re.search(r"edges proved impossible\s*:\s*(\d+)", out)
    if m: got["proved_impossible"] = int(m.group(1))
    for name, label in (("baseline_median_gap", "baseline"),
                        ("conditional_median_gap", "largest_gap")):
        m = re.search(rf"^\s*{label}\s+{NUM}\s+{NUM}", out, re.M)
        if m: got[name] = float(m.group(1))
    tol = float(wc["tolerance"])
    for k, want in wc["recorded"].items():
        g = got.get(k)
        ok = g is not None and (abs(g - want) <= max(tol, tol * abs(want))
                                if isinstance(want, float) else g == want)
        report["reproduction"].append({"metric": k, "recorded": want, "got": g,
                                       "ok": bool(ok)})


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m scmp.reproduce")
    ap.add_argument("--config", default="configs/reproduction.yaml")
    a = ap.parse_args(argv)
    cfg = yaml.safe_load(open(a.config))
    rc = cfg["reproduction"]
    report = {"config": a.config, "name": rc["name"], "environment": [], "source": [],
              "data": [], "reproduction": [], "rerun": {},
              "limitation": "verifies the current interpreter against pins; does not "
                            "provision a fresh virtualenv or container",
              "started": time.strftime("%Y-%m-%d %H:%M:%S")}

    print(f"reproduction: {rc['name']}  config={a.config}\n")
    check_environment(cfg, report)
    check_source_hashes(cfg, report)
    check_data_hashes(cfg, report)

    def summarise(key, label):
        rows = report[key]
        bad = [r for r in rows if not r["ok"]]
        print(f"{label:>18}: {len(rows)-len(bad)}/{len(rows)} verified")
        for r in bad[:6]:
            print(f"                    MISMATCH {r.get('path', r.get('item'))}: "
                  f"want {r['want']} got {r['got']}")
        return bad

    bad_env = summarise("environment", "environment")
    bad_src = summarise("source", "source hashes")
    bad_dat = summarise("data", "data hashes")

    strict = rc.get("strict_hashes", True)
    if strict and (bad_src or bad_dat):
        report["status"] = "aborted_hash_mismatch"
        print("\nABORTING: a source or data hash moved. Reproduction of a different "
              "tree is not a reproduction.")
    else:
        print(f"\nrerunning winning config: "
              f"{' '.join(str(x) for x in cfg['winning_config']['command'])}")
        rerun(cfg, report)
        rows = report["reproduction"]
        bad = [r for r in rows if not r["ok"]]
        print(f"\n{'metric':>26} {'recorded':>12} {'reproduced':>12}  ok")
        for r in rows:
            print(f"{r['metric']:>26} {str(r['recorded']):>12} {str(r['got']):>12}  "
                  f"{'yes' if r['ok'] else 'NO'}")
        report["status"] = "reproduced" if not bad else "diverged"
        print(f"\nstatus: {report['status']}  "
              f"({len(rows)-len(bad)}/{len(rows)} metrics matched, "
              f"{report['rerun']['seconds']:.1f}s)")

    out = rc["out_dir"]
    os.makedirs(out, exist_ok=True)
    path = os.path.join(out, f"reproduction-{time.strftime('%Y%m%d-%H%M%S')}.json")
    payload = json.dumps(report, indent=2, default=str)
    with open(path, "w") as fh:
        fh.write(payload)
    digest = hashlib.sha256(payload.encode()).hexdigest()
    with open(path + ".sha256", "w") as fh:
        fh.write(digest + "\n")
    print(f"\nimmutable result: {path}")
    print(f"content sha256  : {digest[:32]}...")
    return 0 if report["status"] == "reproduced" else 1


if __name__ == "__main__":
    sys.exit(main())
