"""Render paper figures from immutable saved results; no model imports or fitting.

The default destination is paper/figures. Existing destinations are refused.
Archived pilot PNG/SVG files are copied byte for byte. Their PDF companions
embed those PNGs; use the original SVGs for vector publication.
"""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[1]
PILOT = "runs/interaction_recoverability/learned-target-v1-20260912-164000"
P4 = "runs/revision/p4_gpu-20260911-171200-784773-"
LOCAL = "docs/interaction_recoverability/local_five_law/v1/checks.json"
ATTAIN = "docs/interaction_recoverability/local_attainment/v1/checks.json"
INSTABILITY = "runs/interaction_recoverability/prep-20260912-082154/numerics-v2/diagnostics.json"
SLACK = "runs/20260904-034559/c5_slack_scaling.json"
KAPPA = "runs/20260904-043204/c7_kappa_characterization.json"
COVERAGE = "runs/20260904-020253/f7_soundness.json"
THEOREM = "paper/research/docs/interaction_recoverability/local_attainment/v1/theorem.md"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "paper/figures")
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyBboxPatch
    from PIL import Image, ImageDraw
    from reportlab.pdfgen.canvas import Canvas
    plt.rcParams.update({
        "font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.alpha": .18, "svg.fonttype": "none",
        "svg.hashsalt": "iclr2027-saved-results-v1", "pdf.fonttype": 42,
        "savefig.dpi": 200, "font.family": "DejaVu Sans",
    })
    colors = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#6F6F6F", "#E69F00"]
    records = []

    def js(path):
        return json.loads((ROOT / path).read_text())

    def csvs(path):
        with (ROOT / path).open(newline="") as f:
            return list(csv.DictReader(f))

    def record(name, title, caption, inputs, **extra):
        records.append(dict(
            id=name, title=title, caption=caption,
            sources=[dict(path=p, sha256=sha(ROOT / p)) for p in inputs],
            outputs={ext: dict(path=str((out / (name+"."+ext)).relative_to(ROOT))
                              if out.is_relative_to(ROOT) else str(out / (name+"."+ext)),
                              sha256=sha(out / (name+"."+ext))) for ext in ["png", "svg", "pdf"]},
            **extra,
        ))

    def save(fig, name, title, caption, inputs, **extra):
        fig.savefig(out / (name+".png"), bbox_inches="tight")
        fig.savefig(out / (name+".svg"), bbox_inches="tight", metadata={"Date": None})
        fig.savefig(out / (name+".pdf"), bbox_inches="tight",
                    metadata={"CreationDate": None, "ModDate": None})
        plt.close(fig)
        record(name, title, caption, inputs, pdf_kind="vector", **extra)

    # Static diagram: no simulated values or unobserved truth supplied to a learner.
    fig, ax = plt.subplots(figsize=(11, 5.7))
    ax.set(xlim=(0, 11), ylim=(0, 6)); ax.axis("off")
    def box(x, y, w, h, text, color):
        ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.09",
                                   facecolor=color, edgecolor="#334155", linewidth=.8))
        ax.text(x+w/2,y+h/2,text,ha="center",va="center",fontsize=10)
    box(.15,3.9,3.05,1.4,"Five independent source groups\nReference + four single actions\nn per group; total N = 5n","#DBEAFE")
    box(3.85,3.9,3.05,1.4,"Empirical quantile curves\nKnown latent law and warp\nTwo calibrated unit anchors","#DBEAFE")
    box(7.55,3.9,3.05,1.4,"Contracting central identity\nFinite depth; normalization\nNo profile derivatives required","#DCFCE7")
    box(7.55,1.45,3.05,1.4,"Source-only calibration\nAll four coefficients unknown\nDeclared local search region","#DCFCE7")
    box(3.85,1.45,3.05,1.4,"Signed target integration\nOutward lattice + tail control\nShared-reference covariance","#FEF3C7")
    box(.15,1.45,3.05,1.4,"Target contrast estimate\nPointwise local expansion\nJoint outcomes not observed","#FEF3C7")
    for a,b in [((3.3,4.6),(3.75,4.6)),((7,4.6),(7.45,4.6)),
                ((9.08,3.8),(9.08,2.95)),((7.45,2.15),(7,2.15)),
                ((3.75,2.15),(3.3,2.15))]:
        ax.annotate("",xy=b,xytext=a,arrowprops=dict(arrowstyle="->",lw=1.6))
    ax.text(5.5,.55,"Mathematical construction with certified-arithmetic specification.\n"
            "Reference backend is uncertified; sampled-data execution remains disabled.",
            ha="center",va="center",fontsize=11,color="#7C2D12")
    save(fig,"fig01_five_law_construction","Five-law direct-target construction",
         "Static schematic of the saved local attainment construction. Both profiles and all four coefficients are unknown; the latent law, warp, private unit anchors and declared neighborhood are assumed known. This is a mathematical estimator specification, not an executed data analysis. Efficiency, global rates and biological applicability are not established.",
         [THEOREM])

    archived = [
        ("fig02_pilot_interaction_error","figure1_interaction_error","Learned pilot: interaction error",
         "Archived figure, unchanged. All three replicates per cell are shown. The frozen pilot has 24 datasets and ten methods. This descriptive comparison does not establish a statistical rate or significance."),
        ("fig03_pilot_ablation_cost","figure2_ablation_cost","Learned pilot: ablations and recorded costs",
         "Archived figure, unchanged. Adaptive aggregate MSE improved by 3.5569% against the strongest complete control, below the frozen 10% screen. Shared fitted work is charged in the saved component costs; the later audit qualifies these as component accounting, not an independently optimized equal-runtime frontier."),
        ("fig04_pilot_witness_nonlinearity","figure3_witness_nonlinearity","Learned pilot: witnesses and nonlinear remainder",
         "Archived figure, unchanged. Finite direction witnesses are lower diagnostics, not upper uncertainty bounds over the full profile class. The observed finite-path nonlinear remainder is retained. Near-diagonal witnesses do not demonstrate the proposed mechanism."),
    ]
    for name, old, title, caption in archived:
        inputs=[f"{PILOT}/evaluation/{old}.{ext}" for ext in ["png","svg"]]
        for source in inputs:
            shutil.copyfile(ROOT/source, out/(name+Path(source).suffix))
        with Image.open(out/(name+".png")) as im:
            w,h=im.size
        width=720; height=width*h/w
        pdf=Canvas(str(out/(name+".pdf")),pagesize=(width,height),invariant=1)
        pdf.drawImage(str(out/(name+".png")),0,0,width=width,height=height)
        pdf.showPage(); pdf.save()
        record(name,title,caption,inputs+[f"{PILOT}/evaluation/raw_results.csv",
               f"{PILOT}/evaluation/summary.json",f"{PILOT}/evaluation/figure_data.json"],
               pdf_kind="raster archival copy; original SVG is vector",
               archived_png_svg_unchanged=True)

    d=js(LOCAL)
    fig,axs=plt.subplots(1,2,figsize=(10,3.8),layout="constrained")
    for i,row in enumerate(d["derivative_checks"]):
        points=sorted((int(k),v) for k,v in row["neumann_max_errors"].items())
        axs[0].plot(*zip(*points),marker=["o","s"][i],color=colors[i],label=row["kind"])
    axs[0].set(yscale="log",xlabel="Saved contraction depth",ylabel="Maximum reconstruction error",
               title="A  Central-series diagnostics",xticks=[32,128,512])
    axs[0].legend(fontsize=8)
    rows=d["profile_functional_checks"]
    for i,kind in enumerate(["trig","gaussian_holdout"]):
        rr=rows[2*i:2*i+2]
        axs[1].plot([r["order"] for r in rr],[r["absolute_error"] for r in rr],
                    marker=["o","s"][i],color=colors[i],label=kind)
    axs[1].set(yscale="log",xlabel="Saved quadrature order",ylabel="Absolute functional residual",
               title="B  Quadrature diagnostics",xticks=[64,128])
    save(fig,"fig05_local_analytic_checks","Local five-law analytic checks",
         "Saved deterministic trigonometric and Gaussian-holdout fixtures. Depths and quadrature orders are the original check settings. These are reconstruction and floating-point refinement diagnostics, not empirical error rates, variance estimates, or certified integration enclosures. The earlier failed absolute-threshold check is retained in the artifact archive.",
         [LOCAL],plotted_rows={"derivative_checks":d["derivative_checks"],"profile_functional_checks":rows})

    d=js(ATTAIN); rows=d["domain_checks"]; x=list(range(len(rows)))
    fig,ax=plt.subplots(figsize=(7.8,4.4),layout="constrained")
    ax.plot(x,[float(r["lower_effective_rank"]) for r in rows],"o-",color=colors[0],label="Lower effective rank")
    ax.plot(x,[float(r["required_rank"]) for r in rows],"s--",color=colors[1],label="Required rank (sufficient condition)")
    ax.set(yscale="log",xticks=x,xticklabels=[r"$10^{12}$",r"$10^{32}$",r"$10^{64}$"],
           xlabel="Per-source count n in saved schedule diagnostic",ylabel="Rank",
           title="The conservative construction has extreme sufficient requirements")
    for i,r in enumerate(rows):
        ax.annotate("passes" if r["domain_condition"] else "fails",
                    (i,float(r["lower_effective_rank"])),xytext=(0,10),textcoords="offset points",
                    ha="center",fontsize=9)
    ax.margins(y=.2);ax.legend(loc="lower right",fontsize=9)
    save(fig,"fig06_attainment_practicality","Attainment: sufficient-condition limitation",
         "The three saved public-count schedule checks are plotted without rerunning the estimator. This conservative sufficient quantile-domain condition fails at n=10^12 and 10^32 and passes at n=10^64. It is not a necessary sample-size bound or a demonstrated learning curve. The saved precision requirements are 18,438, 31,725 and 52,986 bits; the reference implementation remains uncertified and sample-disabled.",
         [ATTAIN],plotted_rows=rows)

    d=js(INSTABILITY);rows=d["rows"]
    fig,ax=plt.subplots(figsize=(7.5,4),layout="constrained")
    for key,label,color in [("axis2_w1","Observed-axis Wasserstein distance",colors[0]),
                            ("primary_interaction_gap","Target interaction gap",colors[1])]:
        ax.plot([r["delta"] for r in rows],
                [r[key]["value"] if isinstance(r[key],dict) else r[key] for r in rows],
                "o-",label=label,color=color)
    ax.set(xscale="log",yscale="log",xlabel="Saved separation parameter delta",
           ylabel="Distance / absolute contrast gap",title="Historical instability construction")
    ax.legend(fontsize=9)
    save(fig,"supp01_instability","Instability under weak separation",
         "Saved analytic construction values show shrinking observed-axis Wasserstein distance while a target contrast gap persists. Wasserstein proximity alone is not statistical indistinguishability. The complete-experiment Hellinger argument is separate in the proof; floating-point quadrature values are not certified enclosures. This construction is distinct from the calibrated local five-law neighborhood.",
         [INSTABILITY])

    p=P4+"figure1_prediction_ranking.csv";rows=csvs(p)
    fig,axs=plt.subplots(1,2,figsize=(11,4),layout="constrained",sharey=True)
    for i,r in enumerate(rows):
        center=float(r["delta_spearman"]);lo=float(r["delta_ci95_low"]);hi=float(r["delta_ci95_high"])
        axs[0].errorbar(center,i,xerr=[[center-lo],[hi-center]],fmt="o",color=colors[i],capsize=3)
        axs[1].plot(float(r["macro_target_ndcg_at_5"]),i,"o",color=colors[i])
    axs[0].axvline(0,color="black",lw=.8,linestyle="--")
    axs[0].set(yticks=range(len(rows)),yticklabels=[r["arm"] for r in rows],
               xlabel="Delta Spearman vs guide ridge (saved 95% CI)",title="A  Prediction comparison")
    axs[0].invert_yaxis()
    axs[1].set(xlabel="Macro target NDCG@5",title="B  Within-target ranking")
    save(fig,"supp02_p4_prediction","P4 prediction and ranking: negative outcome",
         "All six saved arms, 262 held-out rows, and previously computed paired confidence intervals. No resampling was performed for this export. The model capability gate failed; these RNA comparisons are exploratory in an inferred contiguous reporter context. Actual A+B Spearman is 0.2508 versus 0.6475 for guide ridge. Prediction and ranking do not certify biological ordering.",
         [p],plotted_rows=rows)

    p=P4+"figure2_coverage_total_compute_v2.csv"
    rows=[r for r in csvs(p) if int(r["amortization_instances"])==100]
    strata=sorted({r["stratum"] for r in rows});arms=list(dict.fromkeys(r["arm"] for r in rows))
    fig,axs=plt.subplots(1,2,figsize=(10,4.5),layout="constrained")
    markers=["o","s","^","D","x"]
    for ax,stratum in zip(axs,strata):
        for i,arm in enumerate(arms):
            rr=sorted([r for r in rows if r["stratum"]==stratum and r["arm"]==arm],key=lambda r:float(r["total_budget_s"]))
            ax.plot([float(r["total_budget_s"]) for r in rr],[float(r["decision_coverage"]) for r in rr],
                    color=colors[i],marker=markers[i],markersize=10-i,markerfacecolor="none",
                    lw=1,linestyle=["-","--",":","-.",":"][i],label=arm,zorder=10-i)
        ax.set(xscale="log",ylim=(-.05,1.12),xlabel="Total seconds; 100-instance amortization",
               ylabel="Decision coverage",title=stratum.replace("_"," "))
    axs[1].legend(loc="lower right",fontsize=8)
    save(fig,"supp03_p4_coverage_cost","P4 coverage at matched total cost",
         "All five arms and both saved strata at the primary 100-instance amortization. Overlapping curves are retained. The original family saturates at coverage 1 by 0.05 seconds; in the wider family exhaustive enumeration reaches 1 while the other arms reach 0.75. Missing gap values remain missing in the CSV. Numerical status is float_unverified and the capability gate failed, so this is a diagnostic, not learned superiority.",
         [p],plotted_rows=rows)

    p=P4+"target_summary.csv";rows=[r for r in csvs(p) if int(r["n"])>=5]
    fig,ax=plt.subplots(figsize=(8.8,5.7),layout="constrained")
    for i,r in enumerate(rows):
        a=float(r["reference_selected_activity"]);b=float(r["selected_activity"])
        ax.plot([a,b],[i,i],color="#AAAAAA",lw=1)
        ax.plot(a,i,"o",color=colors[0],label="Guide ridge" if i==0 else None)
        ax.plot(b,i,"s",color=colors[1],label="Actual A+B" if i==0 else None)
    ax.set(yticks=range(len(rows)),yticklabels=[f'{r["gene"]} (n={r["n"]})' for r in rows],
           xlabel="Mean measured activity of selected top five",
           title="Saved target pools: shortlist utility")
    ax.invert_yaxis();ax.legend(fontsize=9)
    save(fig,"supp04_p4_target_utility","P4 target-level shortlist utility",
         "All 16 saved held-out target pools with at least five candidates. Paired points are the saved selected-activity summaries for actual A+B and guide ridge. The 262-row candidate table is retained for traceability. These are observed assay selection summaries under the recorded reporter-context assumptions, not treatment recommendations, causal interaction estimates or biological certificates.",
         [p,P4+"figure3_target_candidates.csv"],plotted_rows=rows)

    d=js(SLACK);k=js(KAPPA)
    fig,axs=plt.subplots(1,2,figsize=(10,4),layout="constrained")
    rr=d["part_a_unbounded"]
    axs[0].loglog([r["m"] for r in rr],[r["slack"] for r in rr],"o-",color=colors[0])
    axs[0].set(xlabel="Matching-family size m",ylabel="Saved relaxation slack",
               title="A  Slack can grow without bound")
    for i,family in enumerate(dict.fromkeys(r["family"] for r in k["rows"])):
        rr=[r for r in k["rows"] if r["family"]==family]
        axs[1].scatter([r["kappa"] for r in rr],[r["gap_nobias"] for r in rr],
                       s=28,label=family,color=colors[i%len(colors)])
    values=[r["kappa"] for r in k["rows"]]
    axs[1].plot([min(values),max(values)],[min(values),max(values)],"--",color="#555555",lw=.8)
    axs[1].set(xscale="log",yscale="log",xlabel="Walk-count ratio kappa",
               ylabel="Observed no-bias gap",title="B  Exact characterization fixtures")
    axs[1].legend(fontsize=7)
    save(fig,"supp05_structured_relaxation","Historical structured-family relaxation",
         "Historical C5/C7 saved diagnostics. Left: non-union-closed matching-family relaxation slack grows from 8.96 to 2290.30 across the recorded sizes. Right: the no-bias fixture gap agrees with its walk-count ratio to floating-point precision. These conditional algebraic results do not establish predictive utility or biological validity.",
         [SLACK,KAPPA])

    d=js(COVERAGE);rows=d["rows"];x=list(range(len(rows)))
    labels=["all edges" if not r["canonical_only"] else f'canonical\nfloor {r["floor"]:g}' for r in rows]
    fig,axs=plt.subplots(2,1,figsize=(10,6),sharex=True,layout="constrained")
    axs[0].plot(x,[r["median_coverage"] for r in rows],"o-",label="Median coverage",color=colors[0])
    axs[0].plot(x,[r["min_coverage"] for r in rows],"s--",label="Minimum coverage",color=colors[1])
    axs[0].set(ylabel="Sampled structure coverage",ylim=(-.02,1.05),title="Probability truncation trades coverage for tightness")
    axs[0].legend(fontsize=9)
    axs[1].plot(x,[r["median_tightness"] for r in rows],"o-",color=colors[2])
    axs[1].set(yscale="log",ylabel="Median bound tightness ratio",xticks=x,xticklabels=labels)
    axs[1].tick_params(axis="x",labelsize=8)
    save(fig,"supp06_ensemble_coverage","Historical ensemble coverage limitation",
         "All nine saved F7 configurations, 20,000 samples per configuration. Coverage is measured against sampled structures, not proved ensemble probability. At positive truncation floors the lattice can exclude valid structures; the saved violation counts are retained in the input JSON. Zero observed violations cannot prove soundness. These records explain why lattice exactness must not be described as universal Boltzmann-ensemble coverage.",
         [COVERAGE],plotted_rows=rows)

    manifest=dict(version=1,scope="Rendering and copying existing artifacts only; no new numerical experiments",
                  figure_count=len(records),figures=records,
                  renderer_sha256=sha(Path(__file__)),numerical_workers=1,numerical_threads=1,
                  affinity=sorted(os.sched_getaffinity(0)),gpu_s=0)
    (out/"manifest.json").write_text(json.dumps(manifest,indent=2,allow_nan=False)+"\n")
    # A visual overview is also derived solely from the exported figures.
    tw,th=580,350;sheet=Image.new("RGB",(tw*3,th*4),"white");draw=ImageDraw.Draw(sheet)
    for i,r in enumerate(records):
        with Image.open(out/(r["id"]+".png")) as im:
            im=im.convert("RGB");im.thumbnail((tw-20,th-35))
            x=(i%3)*tw+(tw-im.width)//2;y=(i//3)*th+28
            sheet.paste(im,(x,y))
        draw.text(((i%3)*tw+10,(i//3)*th+8),r["id"],fill="black")
    sheet.save(out/"overview.png")
    print(json.dumps(dict(figures=len(records),output=str(out),new_fits=0,new_experiments=0)))


if __name__=="__main__":
    main()
