# Exact static document edit after independent proof review; covered by overhead debit.
# Does not change any historical proof or execute numerical calculations.
python -B - <<'PY'
from pathlib import Path
p=Path('docs/interaction_recoverability/direct_target/v1/modulus.md')
s=p.read_text()
assert 'C0=65.6233873491624' in s
p.write_text(s.replace('C0=65.6233873491624', 'C0=49.3625+20.38 sqrt(2/pi) (approximately 65.6233873491624)'))
PY
