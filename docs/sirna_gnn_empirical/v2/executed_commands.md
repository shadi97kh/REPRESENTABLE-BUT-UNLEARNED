# Exact executed stage commands

Commands below are verbatim runner-expanded commands; full logs and failed attempts remain in the run. No prepared-only fit is presented as executed.

```sh
/opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v2/audit.py
```

Exit 0; wall 16.451 s; child CPU 16.572 s.

```sh
/opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v2/plan.py
```

Exit 0; wall 5.375 s; child CPU 5.54 s.

```sh
/opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v2/factorial.py
```

Exit -15; wall 57.929 s; child CPU 57.795 s.

```sh
/opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v2/factorial.py
```

Exit 0; wall 62.24 s; child CPU 73.879 s.

```sh
/opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v2/grouped.py
```

Exit 0; wall 260.513 s; child CPU 302.492 s.

```sh
/opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v2/b2.py
```

Exit 0; wall 56.626 s; child CPU 57.382 s.

```sh
/opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v2/evaluate.py
```

Exit 0; wall 12.789 s; child CPU 12.836 s.

```sh
/opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v2/external.py
```

Exit 0; wall 6.277 s; child CPU 6.476 s.

```sh
/opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v2/figures.py
```

Exit 0; wall 10.034 s; child CPU 10.11 s.

```sh
/opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v2/source_checks.py
```

Exit 0; wall 2.319 s; child CPU 2.423 s.

```sh
/opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v2/paper.py
```

Exit 0; wall 35.627 s; child CPU 35.688 s.

```sh
/opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v2/integrity.py
```

Exit 0; wall 5.726 s; child CPU 5.856 s.

```sh
/opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v2/external_uncertainty.py
```

Exit 0; wall 2.118 s; child CPU 2.229 s.

```sh
/opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v2/paper.py
```

Exit 0; wall 35.779 s; child CPU 35.849 s.

```sh
/opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v2/integrity.py
```

Exit 0; wall 5.723 s; child CPU 5.836 s.

```sh
/opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v2/paper.py
```

Exit 0; wall 35.828 s; child CPU 35.931 s.

```sh
/opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v2/integrity.py
```

Exit 0; wall 5.725 s; child CPU 5.828 s.

Additional direct administrative commands included `source_checks.py`, `prepare_paper.py`, `paper_tables.py`, `tectonic --keep-logs --keep-intermediates .../v3/main.tex`, `pdftotext`, `pdftoppm`, source HTTP acquisition and preservation hashing. Source checks were subsequently recorded as their own runner stage. Their earlier direct administrative cost is not falsely included in the stage ledger.

Resumable entry point: `bash sirna_gnn_empirical/v2/run_all.sh`. In a scientific-only archive, pass `audit source_checks plan factorial grouped b2 external evaluate figures`; paper/integrity stages additionally require the manuscript/full original workspace. New-machine fitting uses the same scripts with a fresh run path and the preserved data; the included completed manifests verify current outputs without training again.
