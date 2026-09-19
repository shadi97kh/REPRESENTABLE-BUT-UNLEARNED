"""Static preservation/syntax/link review. Does not import or run the estimator."""
import ast
import hashlib
import json
from pathlib import Path
import re
import subprocess

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
RUN = ROOT / "runs/interaction_recoverability/learned-target-v1-20260912-164000"


def main():
    protected = json.loads((HERE/"preservation.json").read_text())
    changed = [p for p,h in protected.items() if not (ROOT/p).is_file()
               or hashlib.sha256((ROOT/p).read_bytes()).hexdigest() != h]
    assert not changed, changed
    for name in ["estimator.py","checks.py","resource_gate.py","finalize.py"]:
        ast.parse((HERE/name).read_text(), filename=name)
    source = ast.parse((HERE/"estimator.py").read_text())
    flag = [node.value for node in source.body if isinstance(node,ast.Assign)
            and any(isinstance(t,ast.Name) and t.id=="SAMPLED_EXECUTION_ENABLED" for t in node.targets)]
    assert len(flag)==1 and isinstance(flag[0],ast.Constant) and flag[0].value is False
    missing = []
    for name in ["theorem.md","estimator.md","novelty.md","decision.md"]:
        body=(HERE/name).read_text()
        assert not [c for c in body if ord(c)<32 and c not in "\n\t"]
        assert body.count(r"\[")==body.count(r"\]")
        for target in re.findall(r"\]\(([^)]+)\)",body):
            if "://" in target or target=="commands_and_resources.md":
                continue
            if not (HERE/target.split("#")[0]).is_file():
                missing.append([name,target])
    assert not missing, missing
    result=json.loads((HERE/"checks.json").read_text())
    assert result["status"]=="passed" and result["gpu_s"]==0
    rows=[json.loads(x) for x in (RUN/"ledger.jsonl").read_text().splitlines()]
    marker=next(x for x in rows if x.get("task")=="local-attainment-v1")
    completed={x["job_id"]:x for x in rows if x["event"]=="end" and x["job_id"]>marker["prior_through_job"]}
    journal=[dict(start=x,end=completed[x["job_id"]]) for x in rows
             if x["event"]=="start" and x["job_id"] in completed]
    with (HERE/"command_records.json").open("x") as f:
        json.dump(journal,f,indent=2)
    tracked=subprocess.run(["git","diff","--stat"],cwd=ROOT,text=True,capture_output=True,check=True).stdout
    names=["theorem.md","estimator.md","novelty.md","decision.md","estimator.py",
           "checks.py","checks.json","resource_gate.py","finalize.py","command_records.json"]
    output=dict(status="passed",preserved_file_count=len(protected),changed_preserved_files=changed,
                tracked_diff_stat=tracked,sampled_execution_enabled=False,numerical_checks_rerun=False,
                journal_completed_jobs=sorted(completed),
                artifact_sha256={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in names},
                note="The completed status of this static job is in the live ledger. The final resource audit/note are authored afterward under the existing overhead debit.")
    with (HERE/"final_checks.json").open("x") as f:
        json.dump(output,f,indent=2)
    print(json.dumps(output,indent=2))


if __name__=="__main__":
    main()
