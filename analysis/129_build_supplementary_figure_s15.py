import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(os.environ.get("ICCA_VM_PROJECT_ROOT", Path(__file__).resolve().parents[1]))
RESULTS = ROOT / "results" / "116_secact_robustness_benchmarks"
FIGURES = ROOT / "figures" / "manuscript_v1"
SUPPLEMENT = FIGURES / "supplementary"
SUPPLEMENT.mkdir(parents=True, exist_ok=True)


def panel_label(axis, label):
    axis.text(-0.12, 1.08, label, transform=axis.transAxes, fontsize=12, fontweight="bold", va="top")


def main():
    fdr_summary = pd.read_csv(RESULTS / "secact_global_fdr_summary.tsv", sep="\t")
    benchmark = pd.read_csv(RESULTS / "published_signature_benchmark.tsv", sep="\t")
    anchor = pd.read_csv(RESULTS / "leave_one_anchor_out.tsv", sep="\t")
    bulk = pd.read_csv(RESULTS / "bulk_leave_one_cohort_out.tsv", sep="\t")
    secact_loo = pd.read_csv(RESULTS / "secact_leave_one_cohort_out.tsv", sep="\t")
    figure, axes = plt.subplots(2, 2, figsize=(11, 8.2))

    axis = axes[0, 0]
    labels = ["Nominal\nboth", "BH-FDR\nboth", "Fisher\ncombined"]
    values = [int(fdr_summary.loc[fdr_summary.metric == "nominal_positive_both", "value"].iloc[0]), int(fdr_summary.loc[fdr_summary.metric == "global_FDR_positive_both", "value"].iloc[0]), int(fdr_summary.loc[fdr_summary.metric == "combined_Fisher_FDR_positive", "value"].iloc[0])]
    bars = axis.bar(labels, values, color=["#9BBBD4", "#3B82A0", "#C97B63"], width=0.62)
    axis.set_ylim(0, 340)
    axis.set_ylabel("SecAct groups")
    axis.set_title("Global multiple-testing audit", fontsize=10, pad=12)
    axis.spines[["top", "right"]].set_visible(False)
    for bar, value in zip(bars, values):
        axis.text(bar.get_x() + bar.get_width() / 2, value + 8, str(value), ha="center", va="bottom", fontsize=9)
    axis.text(0.02, 0.03, "657 nonredundant groups; Fisher is sensitivity only", transform=axis.transAxes, fontsize=7, color="#5D6670")
    panel_label(axis, "A")

    axis = axes[0, 1]
    display_names = {"six_anchor_program": "Six-anchor program", "published_hallmark_hypoxia": "Hallmark hypoxia", "published_hallmark_emt": "Hallmark EMT", "published_hallmark_angiogenesis": "Hallmark angiogenesis", "published_hallmark_inflammatory_response": "Hallmark inflammatory response", "published_hallmark_tgf_beta_signaling": "Hallmark TGF-β"}
    plot_data = benchmark.copy()
    plot_data["display"] = plot_data["signature"].map(display_names)
    plot_data = plot_data.dropna(subset=["display"])
    matrix = plot_data.pivot(index="display", columns="cohort", values="spearman_with_anchor_expression_score")
    ordered = [display_names[key] for key in ["six_anchor_program", "published_hallmark_hypoxia", "published_hallmark_emt", "published_hallmark_angiogenesis", "published_hallmark_inflammatory_response", "published_hallmark_tgf_beta_signaling"]]
    matrix = matrix.reindex(ordered)
    image = axis.imshow(matrix.to_numpy(float), cmap="coolwarm", vmin=-1, vmax=1, aspect="auto")
    axis.set_xticks(np.arange(len(matrix.columns)), ["GSE125449", "GSE138709"], rotation=25, ha="right")
    axis.set_yticks(np.arange(len(matrix.index)), matrix.index)
    axis.tick_params(axis="y", labelsize=7.5)
    for row_index in range(matrix.shape[0]):
        for column_index in range(matrix.shape[1]):
            axis.text(column_index, row_index, f"{matrix.iloc[row_index, column_index]:.2f}", ha="center", va="center", fontsize=8)
    axis.set_title("Published-program benchmark", fontsize=10, pad=12)
    axis.set_xlabel("Cohort")
    figure.colorbar(image, ax=axis, fraction=0.046, pad=0.04, label="Spearman ρ")
    panel_label(axis, "B")

    axis = axes[1, 0]
    for cohort, color in [("GSE138709", "#2B7A9B"), ("GSE125449", "#C97B63")]:
        subset = anchor[anchor.cohort == cohort]
        axis.plot(subset["omitted_anchor"], subset["spearman_with_full_anchor_score"], marker="o", linewidth=1.6, label=cohort, color=color)
    axis.axhline(0.8, color="#9AA3AD", linestyle="--", linewidth=0.8)
    axis.set_ylim(0.55, 1.03)
    axis.set_ylabel("ρ with full six-anchor score")
    axis.set_title("Leave-one-anchor-out stability", fontsize=10, pad=12)
    axis.tick_params(axis="x", rotation=35)
    axis.legend(frameon=False, fontsize=8, loc="lower left")
    axis.spines[["top", "right"]].set_visible(False)
    panel_label(axis, "C")

    axis = axes[1, 1]
    bulk_summary = bulk.groupby("omitted_cohort", as_index=False)["direction_consistent_remaining"].sum()
    omitted_order = ["GSE32225", "GSE107943", "GSE26566"]
    bulk_summary["omitted_cohort"] = pd.Categorical(bulk_summary["omitted_cohort"], omitted_order, ordered=True)
    bulk_summary = bulk_summary.sort_values("omitted_cohort")
    bars = axis.bar(bulk_summary["omitted_cohort"].astype(str), bulk_summary["direction_consistent_remaining"], color="#5E9C89", width=0.58)
    axis.set_ylim(0, 7.8)
    axis.set_ylabel("Direction-consistent groups / 7")
    axis.set_title("Leave-one-cohort-out directionality", fontsize=10, pad=12)
    axis.tick_params(axis="x", rotation=25)
    axis.spines[["top", "right"]].set_visible(False)
    for bar in bars:
        axis.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.2, f"{int(bar.get_height())}/7", ha="center", fontsize=9)
    secact_text = secact_loo.groupby("omitted_cohort")["global_FDR_positive_retained"].sum().to_dict()
    axis.text(0.02, 0.03, f"SecAct one-cohort FDR+ after exclusion: {secact_text.get('GSE138709', 0)} or {secact_text.get('GSE125449', 0)} groups; descriptive only", transform=axis.transAxes, fontsize=7, color="#5D6670")
    panel_label(axis, "D")

    figure.suptitle("Supplementary Fig. S15 | Robustness and calibration of the VM-associated molecular coordinate", fontsize=13, fontweight="bold", x=0.03, ha="left")
    figure.tight_layout(rect=[0, 0, 1, 0.96])
    figure.savefig(RESULTS / "FigS15_robustness_calibration_composite.png", dpi=400)
    figure.savefig(RESULTS / "FigS15_robustness_calibration_composite.pdf")
    figure.savefig(SUPPLEMENT / "FigS15_robustness_calibration_composite.png", dpi=400)
    figure.savefig(SUPPLEMENT / "FigS15_robustness_calibration_composite.pdf")
    plt.close(figure)


if __name__ == "__main__":
    main()
