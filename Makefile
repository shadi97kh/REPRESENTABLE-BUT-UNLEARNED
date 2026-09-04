PY ?= python3
.PHONY: setup test f1 f2 f3 f4 all freeze
setup:  ; $(PY) -m pip install -r requirements.txt
test:   ; $(PY) tests/test_sanity.py
f1:     ; $(PY) -m experiments.f1_extremality
f2:     ; $(PY) -m experiments.f2_lattice_size
f3:     ; $(PY) -m experiments.f3_certificate
f4:     ; $(PY) -m experiments.f4_topk_stability
all: test f1 f2 f3 f4
freeze: ; git add -A && git commit -m "freeze: prereg + scaffold" && git tag -f prereg
