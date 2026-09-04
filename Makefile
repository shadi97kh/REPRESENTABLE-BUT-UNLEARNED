PY ?= python3
export PYTHONPATH := .
.PHONY: setup data gencode test f1 f2 f3 f3b f4 f5 f6 f7 f8 f9 f10 f10b c3 c4 c5 c6 all theory freeze
setup:  ; $(PY) -m pip install -r requirements.txt

# Huesken et al. 2005 in its published 2182/249 split, as redistributed with DSIR
# (Vert et al., BMC Bioinformatics 7:520, 2006). Not vendored: it is third-party data.
# certmp.data.verify rejects the files if they do not match the published statistics.
gencode:
	mkdir -p data/gencode
	curl -fsSL -o data/gencode/gencode.v47.pc_transcripts.fa.gz \
	  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.pc_transcripts.fa.gz

data:
	mkdir -p data
	curl -fsSL -o data/TrainAll2182.txt https://biodev.cea.fr/DSIR/data/TrainAll2182.txt
	curl -fsSL -o data/TestAll249.txt   https://biodev.cea.fr/DSIR/data/TestAll249.txt
	$(PY) -c "from certmp import data as d; \
	  tr,te,rep=d.load_dsir_split('data/TrainAll2182.txt','data/TestAll249.txt'); \
	  v=d.verify('huesken',tr+te); print(rep); print(v); \
	  raise SystemExit(0 if v['verified'] else 1)"

test:   ; $(PY) tests/test_sanity.py
f1:     ; $(PY) -m experiments.f1_extremality
f2:     ; $(PY) -m experiments.f2_lattice_size
f3:     ; $(PY) -m experiments.f3_certificate
f3b:    ; $(PY) -m experiments.f3b_relaxation_slack
f4:     ; $(PY) -m experiments.f4_topk_stability
f5:     ; $(PY) -m experiments.f5_stress_extremality
f6:     ; $(PY) -m experiments.f6_trained
f7:     ; $(PY) -m experiments.f7_soundness
f8:     ; $(PY) -m experiments.f8_baseline
f9:     ; $(PY) -m experiments.f9_target_context
f10:    ; $(PY) -m experiments.f10_maximal
f10b:   ; $(PY) -m experiments.f10b_convergence
c3:     ; $(PY) -m experiments.c3_aggregator_characterization
c4:     ; $(PY) -m experiments.c4_affinity
c5:     ; $(PY) -m experiments.c5_slack_scaling
c6:     ; $(PY) -m experiments.c6_closure
theory: c3 c4 c5 c6
all: test f1 f2 f3b f5 f6 f7 f8 f9 f10 f10b theory   # f3, f4 retired
freeze: ; git add -A && git commit -m "freeze: prereg + scaffold" && git tag -f prereg
