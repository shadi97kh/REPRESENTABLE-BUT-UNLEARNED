#!/usr/bin/env bash
# Exact admitted commands for the completed continuation, from /home/shadi/iclr2027.
# Historical record, not an automatic replay script: registration and output files refuse overwrite.
# No sample fitting is authorized by these commands.
PYTHONDONTWRITEBYTECODE=1 python docs/interaction_recoverability/direct_target/v1/resource_gate.py register
PYTHONDONTWRITEBYTECODE=1 python docs/interaction_recoverability/direct_target/v1/resource_gate.py --seconds 60 run -- python -B docs/interaction_recoverability/direct_target/v1/checks.py
PYTHONDONTWRITEBYTECODE=1 python docs/interaction_recoverability/direct_target/v1/resource_gate.py --seconds 30 run -- cat interaction_recoverability/unknown_profiles/core.py docs/interaction_recoverability/learned_target_v1/novelty.md docs/interaction_recoverability/learned_target_v1/source_manifest.json
PYTHONDONTWRITEBYTECODE=1 python docs/interaction_recoverability/direct_target/v1/resource_gate.py --seconds 30 run -- python -B docs/interaction_recoverability/direct_target/v1/final_checks.py
PYTHONDONTWRITEBYTECODE=1 python docs/interaction_recoverability/direct_target/v1/resource_gate.py verify
PYTHONDONTWRITEBYTECODE=1 python -B docs/interaction_recoverability/direct_target/v1/finish_report.py
