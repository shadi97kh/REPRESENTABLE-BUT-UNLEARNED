"""Static artifact verification and exact completed-child command journal.

This does not import any numerical or learner implementation, replay a check,
or update a prior artifact. Run under resource_gate.py.
"""
import ast
import hashlib
import json
from pathlib import Path
import re
import subprocess

ROOT = Path(__file__).resolve().parents[4]
HERE = Path(__file__).resolve().parent
RUN = ROOT / "runs/interaction_recoverability/learned-target-v1-20260912-164000"


def main():
    hashes = json.loads((HERE / "preservation.json").read_text())
    changed = [rel for rel, sha in hashes.items()
               if not (ROOT / rel).is_file()
               or hashlib.sha256((ROOT / rel).read_bytes()).hexdigest() != sha]
    assert not changed, changed
    for name in ["checks.py", "checks_attempt_job33.py", "resource_gate.py", "finalize.py"]:
        ast.parse((HERE / name).read_text(), filename=name)
    missing = []
    for name in ["proof.md", "decision.md"]:
        body = (HERE / name).read_text()
        assert not [c for c in body if ord(c) < 32 and c not in "\n\t"]
        assert body.count(r"\[") == body.count(r"\]")
        for line in body.splitlines():
            assert line.count("$") % 2 == 0, (name, line)
        for target in re.findall(r"\]\(([^)]+)\)", body):
            if "://" in target or target == "commands_and_resources.md":
                continue
            if not (HERE / target.split("#")[0]).is_file():
                missing.append([name, target])
    assert not missing, missing
    checked = json.loads((HERE / "checks.json").read_text())
    assert checked["status"] == "passed"
    assert checked["numerical_workers"] == checked["numerical_threads"] == 1
    assert checked["gpu_s"] == 0
    rows = [json.loads(line) for line in (RUN / "ledger.jsonl").read_text().splitlines()]
    marker = next(r for r in rows if r.get("task") == "local-five-law-v1")
    completed = {r["job_id"]: r for r in rows
                 if r["event"] == "end" and r["job_id"] > marker["prior_through_job"]}
    journal = []
    for row in rows:
        if row["event"] == "start" and row["job_id"] in completed:
            journal.append({"start": row, "end": completed[row["job_id"]]})
    with (HERE / "command_records.json").open("x") as handle:
        json.dump(journal, handle, indent=2)
    tracked = subprocess.run(["git", "diff", "--stat"], cwd=ROOT,
                             text=True, capture_output=True, check=True).stdout
    status = subprocess.run(["git", "status", "--short"], cwd=ROOT,
                            text=True, capture_output=True, check=True).stdout
    output = {
        "status": "passed",
        "preserved_file_count": len(hashes),
        "changed_preserved_files": changed,
        "tracked_diff_stat": tracked,
        "git_status_short": status,
        "checks_reexecuted": False,
        "journal_completed_jobs": sorted(completed),
        "note": "Finalization's own completed status is in the current ledger and subsequent resource audit. The final resource note is authored after that audit.",
        "artifact_sha256": {
            name: hashlib.sha256((HERE / name).read_bytes()).hexdigest()
            for name in ["proof.md", "decision.md", "checks.py", "checks.json",
                         "checks_attempt_job33.py", "failure_job33.json",
                         "resource_gate.py", "finalize.py", "command_records.json"]
        }
    }
    with (HERE / "final_checks.json").open("x") as handle:
        json.dump(output, handle, indent=2)
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
