# Static discovery/read commands actually issued by the primary agent before registration.
# These do not constitute a replay instruction or fitting authorization.
pwd; git status --short; for p in /AGENTS.md /home/AGENTS.md /home/shadi/AGENTS.md /home/shadi/iclr2027/AGENTS.md; do if [ -f "$p" ]; then cat "$p"; fi; done; cat Codex_Direct_Target_Theorem_Command.md
rg --files -g AGENTS.md -g '!runs/**' -g '!data/**'; cat docs/interaction_recoverability/unknown_profiles/model.md docs/interaction_recoverability/unknown_profiles/proofs.md docs/interaction_recoverability/unknown_profiles/rate_audit.md
cat interaction_recoverability/learned_target/budget.py; tail -3 runs/interaction_recoverability/learned-target-v1-20260912-164000/ledger.jsonl; cat runs/interaction_recoverability/learned-target-v1-20260912-164000/resource_snapshot-28.json; rg --files docs/interaction_recoverability/direct_target
cat docs/interaction_recoverability/learned_target_v1/method.md docs/interaction_recoverability/learned_target_v1/theorem.md docs/interaction_recoverability/learned_target_v1/review.md docs/interaction_recoverability/learned_target_v1/adversarial_review_20260912.md; rg --files interaction_recoverability/unknown_profiles; rg --files runs/interaction_recoverability | rg 'sources/|trabs|donoho|quantile|statistic|optimal'
# Subsequent admitted shell commands and finalization are in commands.sh.
# File-authoring calls, primary web lookups, and static agent reads are described in commands_and_resources.md.
