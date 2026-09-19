# Executed commands

Resumable entry point actually launched: `bash sirna_gnn_empirical/v3/run_all.sh activity`; the queued continuation actually launched pair, b3_evaluate, evaluate and analysis_supplement after activity. The following table is generated from original command records, including failures. Completed-stage verification is not new fitting.

| stage | command | exit_code | wall_s | child_cpu_s | log |
|---|---|---|---|---|---|
| legacy_analysis | /opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v3/legacy_analysis.py | 0 | 75.337861 | 75.421102 | logs/legacy_analysis-1789349804327633453.log |
| protocol | /opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v3/protocol.py | 0 | 1.015337 | 1.089518 | logs/protocol-1789350127831779743.log |
| checks | /opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v3/checks.py | 1 | 3.497328 | 3.643814 | logs/checks-1789350128851738451.log |
| checks | /opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v3/checks.py | 0 | 5.667620 | 5.832451 | logs/checks-1789350200654156312.log |
| activity | /opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v3/activity.py | 0 | 11903.680519 | 12067.060765 | logs/activity-1789350461798247368.log |
| pair | /opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v3/pair.py | 0 | 833.725028 | 837.582085 | logs/pair-1789362379476827760.log |
| b3_evaluate | /opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v3/b3_evaluate.py | 0 | 17.261607 | 17.959085 | logs/b3_evaluate-1789363218313257027.log |
| evaluate | /opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v3/evaluate.py | 0 | 14.582723 | 14.657842 | logs/evaluate-1789363237214631679.log |
| analysis_supplement | /opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v3/analysis_supplement.py | 0 | 42.174583 | 43.055538 | logs/analysis_supplement-1789363253680780935.log |
| source_checks | /opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v3/source_checks.py | 0 | 0.321578 | 0.395792 | logs/source_checks-1789363367538887702.log |
| parameter_audit | /opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v3/parameter_audit.py | 0 | 114.508992 | 123.160566 | logs/parameter_audit-1789363369163029241.log |
| figures | /opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v3/figures.py | 0 | 10.022731 | 10.095333 | logs/figures-1789365055868958927.log |
| paper | /opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v3/paper.py | 0 | 2.171644 | 2.245599 | logs/paper-1789365067278453280.log |
| paper | /opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v3/paper.py | 0 | 2.161546 | 2.235508 | logs/paper-1789365192310315669.log |
| deliver | /opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v3/deliver.py | 1 | 27.070208 | 27.554989 | logs/deliver-1789365791919862422.log |

## Executed manuscript authoring attempts

| command | exit_code | wall_s | child_cpu_s | log |
|---|---|---|---|---|
| /opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v3/build_audit.py | 1 | 0.375629 | 0.375345 | authoring_attempts/build_audit-1789363535000947559.log |
| /opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v3/build_audit.py | 1 | 0.382728 | 0.382409 | authoring_attempts/build_audit-1789363588053753594.log |
| /opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v3/build_audit.py | 1 | 0.378114 | 0.377831 | authoring_attempts/build_audit-1789363604683167957.log |
| /opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v3/build_audit.py | 0 | 47.332222 | 47.173477 | authoring_attempts/build_audit-1789363671290426057.log |
| /opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v3/build_audit.py | 0 | 48.086434 | 48.230772 | authoring_attempts/build_audit-1789364094530040240.log |
| /opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v3/build_audit.py | 0 | 48.179628 | 48.326516 | authoring_attempts/build_audit-1789364380531131164.log |
| /opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v3/build_audit.py | 0 | 48.278040 | 48.419068 | authoring_attempts/build_audit-1789364679947892345.log |
| /opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v3/figures.py | 0 | 10.078257 | 10.119570 | authoring_attempts/figures-1789363522792258941.log |
| /opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v3/figures.py | 0 | 10.095548 | 10.168037 | authoring_attempts/figures-1789364082263622832.log |
| /opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v3/figures.py | 1 | 0.022931 | 0.022620 | authoring_attempts/figures-1789364349365597535.log |
| /opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v3/figures.py | 0 | 10.030507 | 10.095273 | authoring_attempts/figures-1789364368306818153.log |
| /opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v3/figures.py | 0 | 10.019811 | 10.092531 | authoring_attempts/figures-1789364667776694975.log |
| /opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v3/paper.py | 0 | 2.129364 | 2.201176 | authoring_attempts/paper-1789363532871059968.log |
| /opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v3/paper.py | 0 | 2.114839 | 2.188780 | authoring_attempts/paper-1789363602567667867.log |
| /opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v3/paper.py | 0 | 2.152062 | 2.226044 | authoring_attempts/paper-1789363669137833216.log |
| /opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v3/paper.py | 0 | 2.169788 | 2.243727 | authoring_attempts/paper-1789364092359745819.log |
| /opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v3/paper.py | 0 | 2.192765 | 2.266684 | authoring_attempts/paper-1789364378337852224.log |
| /opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v3/paper.py | 0 | 2.150352 | 2.224345 | authoring_attempts/paper-1789364677797026714.log |

These commands perform plotting, table generation or document compilation/rendering, with zero additional model fits. Nested build logs are the same work, not an additive resource estimate. The isolated CairoSVG retry is recorded in administrative_preflight/workflow_row_convention_attempts.json.

Source attempt `sources/acquire_official_pdb.record.json`:

```json
{
  "command": [
    "/opt/miniforge3/bin/python",
    "-m",
    "gdown",
    "https://drive.google.com/uc?id=1F7cNJXMNPSjFb0UvDkDRHTkt9Tt4EGWe",
    "-O",
    "/home/shadi/iclr2027/runs/sirna_gnn_empirical/v3-20260914T013314Z/sources/ensirna_mod_pdb.zip",
    "--no-cookies"
  ],
  "cwd": "/home/shadi/iclr2027",
  "exit_code": 1,
  "wall_s": 0.587351321009919,
  "child_cpu_s": 0.23764799999999986,
  "max_child_rss_kib": 151936,
  "log": "acquire_official_pdb.log"
}
```

Source attempt `sources/install_dependencies.record.json`:

```json
{
  "command": [
    "/opt/miniforge3/bin/python",
    "-m",
    "pip",
    "install",
    "--target",
    "/home/shadi/iclr2027/sirna_gnn_empirical/v3/_deps",
    "gdown",
    "tensorboard"
  ],
  "cwd": "/home/shadi/iclr2027",
  "exit_code": 0,
  "wall_s": 6.366850656922907,
  "child_cpu_s": 4.854154,
  "max_child_rss_kib": 151936,
  "log": "install_dependencies.log"
}
```

Source attempt `sources/official_preprocess_help.record.json`:

```json
{
  "command": [
    "/opt/miniforge3/bin/python",
    "-m",
    "data.get_pdb",
    "--help"
  ],
  "cwd": "/home/shadi/iclr2027/handoffs/codex_revision_v1/ENsiRNA_official/ENsiRNA-mod",
  "exit_code": 0,
  "wall_s": 1.6809359630569816,
  "child_cpu_s": 1.7299979999999997,
  "max_child_rss_kib": 747892,
  "log": "official_preprocess_help.log"
}
```

Source attempt `sources/rosetta_probe.record.json`:

```json
{
  "command": [
    "/app/rosetta/rosetta.binary.linux.release-371/main/source/bin/rna_denovo.static.linuxgccrelease",
    "-help"
  ],
  "error": "[Errno 2] No such file or directory: '/app/rosetta/rosetta.binary.linux.release-371/main/source/bin/rna_denovo.static.linuxgccrelease'",
  "exit_code": null,
  "reason": "Required geometry generator is not installed at official configured path"
}
```

Primary acquisition/Bramsen parsing, source and citation inspection, authoring, workflow export, builds, rendering, copying and hashing also occurred. Their source scripts and retained evidence identify the operations; not every interactive administrative call has an exact process-time meter. No prepared published fit is listed as executed. Later independent compilation and packaging commands have separate verification records. No historical wrapper, ledger write, new wet-lab work, backend development, repository push or external contact occurred.
