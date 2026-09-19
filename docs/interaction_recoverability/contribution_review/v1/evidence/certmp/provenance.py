import json, os, sys, time, hashlib, platform, subprocess

def _git(*a):
    try: return subprocess.check_output(["git",*a], stderr=subprocess.DEVNULL).decode().strip()
    except Exception: return "unavailable"

def new_run(name, config):
    ts = time.strftime("%Y%m%d-%H%M%S"); d = os.path.join("runs", ts); os.makedirs(d, exist_ok=True)
    vers = {}
    for mod in ("numpy", "RNA", "torch"):
        try: vers[mod] = __import__(mod).__version__
        except Exception: vers[mod] = "absent"
    json.dump({"experiment": name, "timestamp": ts, "git_sha": _git("rev-parse","HEAD"),
               "git_dirty": _git("status","--porcelain") != "", "python": sys.version.split()[0],
               "platform": platform.platform(), "versions": vers, "config": config,
               "config_sha256": hashlib.sha256(json.dumps(config, sort_keys=True, default=str).encode()).hexdigest()[:16]},
              open(os.path.join(d,"provenance.json"),"w"), indent=2)
    return d

def record(rundir, name, payload):
    json.dump(payload, open(os.path.join(rundir, f"{name}.json"),"w"), indent=2, default=str)
    print(f"[provenance] {rundir}/{name}.json")
