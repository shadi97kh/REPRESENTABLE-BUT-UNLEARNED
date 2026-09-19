#!/usr/bin/env bash
set -eu

# EXECUTED command journal, through completed job 27
# The existing outputs are exclusive and intentionally refuse replay.
# This file contains the tested operational forms, not authorization to rerun fits.
# All discovery/read/image/report argv and statuses: command_records.json; ledger.jsonl is authoritative.
# Wrapper prefixes below are reconstructed from recorded phase/timeout; child argv are exact.

# Initialization ran successfully before job 1:
python -m interaction_recoverability.learned_target.budget --run-root runs/interaction_recoverability/learned-target-v1-20260912-164000 --init

# Job 3: EXECUTED; exit 0; snapshot resource_snapshot-3.json
python -m interaction_recoverability.learned_target.budget --run-root runs/interaction_recoverability/learned-target-v1-20260912-164000 --phase compute --job-seconds 180.0 -- python -m pytest -q tests/interaction_recoverability_learned_target

# Job 5: EXECUTED; exit 0; snapshot resource_snapshot-5.json
python -m interaction_recoverability.learned_target.budget --run-root runs/interaction_recoverability/learned-target-v1-20260912-164000 --phase compute --job-seconds 600.0 -- python -m interaction_recoverability.learned_target --help

# Job 6: EXECUTED; exit 0; snapshot resource_snapshot-6.json
python -m interaction_recoverability.learned_target.budget --run-root runs/interaction_recoverability/learned-target-v1-20260912-164000 --phase compute --job-seconds 180.0 -- python -m interaction_recoverability.learned_target validate --run-root runs/interaction_recoverability/learned-target-v1-20260912-164000

# Job 7: EXECUTED; exit 0; snapshot resource_snapshot-7.json
python -m interaction_recoverability.learned_target.budget --run-root runs/interaction_recoverability/learned-target-v1-20260912-164000 --phase compute --job-seconds 1200.0 -- python -m interaction_recoverability.learned_target development --run-root runs/interaction_recoverability/learned-target-v1-20260912-164000

# Job 10: EXECUTED; exit 124; snapshot resource_snapshot-10.json
python -m interaction_recoverability.learned_target.budget --run-root runs/interaction_recoverability/learned-target-v1-20260912-164000 --phase compute --job-seconds 1.0 -- python -c 'import subprocess,sys,time; from pathlib import Path; p=subprocess.Popen([sys.executable,"-c","import time; time.sleep(60)"],start_new_session=True); Path("runs/interaction_recoverability/learned-target-v1-20260912-164000/watchdog_child_pid.txt").write_text(str(p.pid)); time.sleep(60)'

# Job 12: EXECUTED; exit 0; snapshot resource_snapshot-12.json
python -m interaction_recoverability.learned_target.budget --run-root runs/interaction_recoverability/learned-target-v1-20260912-164000 --phase compute --job-seconds 60.0 -- python -m interaction_recoverability.learned_target freeze --run-root runs/interaction_recoverability/learned-target-v1-20260912-164000

# Job 13: EXECUTED; exit 0; snapshot resource_snapshot-13.json
python -m interaction_recoverability.learned_target.budget --run-root runs/interaction_recoverability/learned-target-v1-20260912-164000 --phase compute --job-seconds 120.0 -- python -m interaction_recoverability.learned_target generate --run-root runs/interaction_recoverability/learned-target-v1-20260912-164000

# Job 14: EXECUTED; exit 0; snapshot resource_snapshot-14.json
python -m interaction_recoverability.learned_target.budget --run-root runs/interaction_recoverability/learned-target-v1-20260912-164000 --phase compute --job-seconds 1200.0 -- python -m interaction_recoverability.learned_target train --run-root runs/interaction_recoverability/learned-target-v1-20260912-164000 --start 0 --stop 12

# Job 17: EXECUTED; exit 0; snapshot resource_snapshot-17.json
python -m interaction_recoverability.learned_target.budget --run-root runs/interaction_recoverability/learned-target-v1-20260912-164000 --phase compute --job-seconds 1200.0 -- python -m interaction_recoverability.learned_target train --run-root runs/interaction_recoverability/learned-target-v1-20260912-164000 --start 12 --stop 24

# Job 18: EXECUTED; exit 0; snapshot resource_snapshot-18.json
python -m interaction_recoverability.learned_target.budget --run-root runs/interaction_recoverability/learned-target-v1-20260912-164000 --phase report --job-seconds 120.0 -- python -m interaction_recoverability.learned_target evaluate --run-root runs/interaction_recoverability/learned-target-v1-20260912-164000

# Job 19: EXECUTED; exit 0; snapshot resource_snapshot-19.json
python -m interaction_recoverability.learned_target.budget --run-root runs/interaction_recoverability/learned-target-v1-20260912-164000 --phase compute --job-seconds 120.0 -- python -m interaction_recoverability.learned_target plot --run-root runs/interaction_recoverability/learned-target-v1-20260912-164000

# Job 25: EXECUTED; exit 0; snapshot resource_snapshot-25.json
python -m interaction_recoverability.learned_target.budget --run-root runs/interaction_recoverability/learned-target-v1-20260912-164000 --phase report --job-seconds 120.0 -- python docs/interaction_recoverability/learned_target_v1/compile_report.py --run-root runs/interaction_recoverability/learned-target-v1-20260912-164000 report

# Job 26: EXECUTED; exit 0; snapshot resource_snapshot-26.json
python -m interaction_recoverability.learned_target.budget --run-root runs/interaction_recoverability/learned-target-v1-20260912-164000 --phase report --job-seconds 60.0 -- python -m interaction_recoverability.learned_target verify --run-root runs/interaction_recoverability/learned-target-v1-20260912-164000

# Report finalization itself runs after this journal cutoff; its exact command and
# eventual exit status are in ledger.jsonl and the final resource snapshot.
# No future optional fits are listed or authorized here.
