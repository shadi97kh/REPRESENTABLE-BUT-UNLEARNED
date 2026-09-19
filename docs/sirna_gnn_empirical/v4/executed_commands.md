# Executed commands and accounting

The scientific entry point was invoked as:

```bash
bash sirna_gnn_empirical/v4/run_all.sh adjudicate plan fit
```

Every following child command is copied from its actual receipt, with exit status; a receipt is not a prepared command.

```bash
/opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v4/adjudicate.py
```

Exit 0; wall 2.425s; child CPU 2.489s; log `logs/adjudicate-1789371169567026054.log`.

```bash
/opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v4/plan.py
```

Exit 0; wall 1.020s; child CPU 1.094s; log `logs/plan-1789371172003173150.log`.

```bash
/opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v4/train_focused.py
```

Exit 0; wall 4210.790s; child CPU 4423.589s; log `logs/fit-1789371173026278885.log`.

```bash
/opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v4/evaluate.py
```

Exit 0; wall 18.765s; child CPU 19.283s; log `logs/evaluate-1789375415218757896.log`.

```bash
/opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v4/visibility.py
```

Exit 0; wall 8.472s; child CPU 8.617s; log `logs/visibility-1789375434586597322.log`.

```bash
/opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v4/scientific_audit.py
```

Exit 0; wall 12.137s; child CPU 12.209s; log `logs/audit-1789375443496018637.log`.

```bash
/opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v4/paper.py
```

Exit 0; wall 1.538s; child CPU 1.672s; log `logs/paper-1789375456068481086.log`.

```bash
/opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v4/paper.py
```

Exit 0; wall 1.541s; child CPU 1.690s; log `logs/paper-1789375548681280925.log`.

```bash
/opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v4/paper.py
```

Exit 0; wall 1.504s; child CPU 1.652s; log `logs/paper-1789375644207368899.log`.

```bash
/opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v4/paper.py
```

Exit 0; wall 1.513s; child CPU 1.662s; log `logs/paper-1789375804469021774.log`.

```bash
/opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v4/paper.py
```

Exit 0; wall 1.512s; child CPU 1.661s; log `logs/paper-1789375892335534330.log`.

```bash
/opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v4/validate_delivery.py
```

Exit 0; wall 12.250s; child CPU 12.313s; log `logs/validate-1789375941427447754.log`.

```bash
/opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v4/paper.py
```

Exit 0; wall 1.874s; child CPU 2.022s; log `logs/paper-1789376157593086472.log`.

```bash
/opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v4/paper.py
```

Exit 0; wall 1.829s; child CPU 1.978s; log `logs/paper-1789376237662822876.log`.

```bash
/opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v4/validate_delivery.py
```

Exit 0; wall 12.242s; child CPU 12.318s; log `logs/validate-1789376266028372878.log`.

```bash
/opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v4/package.py
```

Exit 0; wall 141.693s; child CPU 141.562s; log `logs/package-1789376481353977586.log`.

```bash
/opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v4/review_checks.py
```

Exit 0; wall 3.269s; child CPU 3.254s; log `logs/review_checks-1789371445923651442.log`.

```bash
/opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v4/source_checks.py
```

Exit 0; wall 8.497s; child CPU 6.821s; log `logs/source_checks-1789371675785658695.log`.

```bash
/opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v4/ranking.py
```

Exit 0; wall 2.704s; child CPU 2.615s; log `logs/ranking-1789371857429116180.log`.

```bash
/opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v4/ranking.py
```

Exit 0; wall 3.027s; child CPU 2.875s; log `logs/ranking-1789372188439550922.log`.

Additional executed direct preparation/authoring commands include `python sirna_gnn_empirical/v4/write_dataset.py`, script syntax checks, `prepare_minimal.py` (one-time initial source creation), iterative paper/build commands, rendering and source-ZIP checks. Their direct interactive costs are not all in the scientific stage receipt and are not assigned invented zero cost. Initial unsuccessful parsing/acquisition and ranking-member-count attempts are retained in the run attempt directories. Additional recorded typesetting attempts are verbatim in `authoring_commands.jsonl`; generation costs overlap the paper-stage receipts and are not summed twice. Final external build commands and their measured costs are verbatim in the manuscript artifacts `validation_commands.json`.

The default resumable entry point is `bash sirna_gnn_empirical/v4/run_all.sh`. It verifies completed stage hashes and skips them; this is distinct from executing a fresh fit. Do not delete completion records to reinterpret an old result as a new experiment. No historical allowance wrapper or ledger is used.
