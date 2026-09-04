PY ?= python3
export PYTHONPATH := .
.PHONY: setup test f1 f2 f3 f3b f4 f5 f6 all freeze
setup:  ; $(PY) -m pip install -r requirements.txt
test:   ; $(PY) tests/test_sanity.py
f1:     ; $(PY) -m experiments.f1_extremality
f2:     ; $(PY) -m experiments.f2_lattice_size
f3:     ; $(PY) -m experiments.f3_certificate
f3b:    ; $(PY) -m experiments.f3b_relaxation_slack
f4:     ; $(PY) -m experiments.f4_topk_stability
f5:     ; $(PY) -m experiments.f5_stress_extremality
f6:     ; $(PY) -m experiments.f6_trained
all: test f1 f2 f3 f3b f4 f5
freeze: ; git add -A && git commit -m "freeze: prereg + scaffold" && git tag -f prereg
