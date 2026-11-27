import csv
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def read_groups(path):
    groups = {
        "AD_only": {"percent": [], "run": []},
        "all_three": {"percent": [], "run": []},
    }
    with open(path, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            group = row["detection_pattern"]
            if group not in groups:
                continue
            groups[group]["percent"].append(float(row["hydrophobic_percent"]))
            groups[group]["run"].append(int(row["max_hydrophobic_run_aa"]))

    for group in groups:
        if not groups[group]["percent"]:
            raise ValueError(f"{group}에 해당하는 서열이 없습니다.")
    return groups


def save_plots(groups, folder):
    folder.mkdir(parents=True, exist_ok=True)
    ad = groups["AD_only"]
    common = groups["all_three"]
    ad_label = f"AD only (n={len(ad['percent'])})"
    common_label = f"All three (n={len(common['percent'])})"

    # 기존 그래프 1. 소수성 비율 상자그림
    fig, ax = plt.subplots(figsize=(7.4, 5))
    ax.boxplot([ad["percent"], common["percent"]])
    ax.set_xticks([1, 2], [ad_label, common_label])
    ax.set_ylabel("Hydrophobic residues (%)")
    ax.set_title("Hydrophobic proportion by detection group")
    fig.tight_layout()
    fig.savefig(folder / "hydrophobic_percent_boxplot.png", dpi=200)
    plt.close(fig)

    # 기존 그래프 2. 소수성 비율 vs 최대 연속 길이 산점도
    fig, ax = plt.subplots(figsize=(7.4, 5))
    ax.scatter(common["percent"], common["run"], label=common_label,
               color="steelblue", alpha=0.4)
    ax.scatter(ad["percent"], ad["run"], label=ad_label,
               color="darkorange", marker="^", alpha=0.5)
    ax.set_xlabel("Hydrophobic residues (%)")
    ax.set_ylabel("Longest hydrophobic run (aa)")
    ax.set_title("Hydrophobic proportion and longest run")
    ax.legend()
    fig.tight_layout()
    fig.savefig(folder / "hydrophobic_run_scatter.png", dpi=200)
    plt.close(fig)

    # 추가 그래프 3. 최대 연속 소수성 길이별 펩타이드 비율 막대그래프
    # 서열 수가 다르므로 각 집단 안에서 비율(%)로 변환하여 공정하게 비교한다.
    max_run = max(max(ad["run"]), max(common["run"]))
    lengths = list(range(max_run + 1))
    ad_counts = Counter(ad["run"])
    common_counts = Counter(common["run"])
    ad_pct = [100 * ad_counts[i] / len(ad["run"]) for i in lengths]
    common_pct = [100 * common_counts[i] / len(common["run"]) for i in lengths]

    fig, ax = plt.subplots(figsize=(9, 5.5))
    x = np.arange(len(lengths))
    width = 0.36
    bars1 = ax.bar(x - width/2, ad_pct, width, label=ad_label,
                   color="darkorange", edgecolor="white")
    bars2 = ax.bar(x + width/2, common_pct, width, label=common_label,
                   color="steelblue", edgecolor="white")
    ax.set_xticks(x, [str(i) for i in lengths])
    ax.set_xlabel("Longest consecutive hydrophobic run (aa)")
    ax.set_ylabel("Peptides within each group (%)")
    ax.set_title("Distribution of longest hydrophobic runs")
    ax.set_ylim(0, max(ad_pct + common_pct) * 1.19)
    ax.spines[['top','right']].set_visible(False)
    ax.yaxis.grid(True, alpha=0.17)
    ax.set_axisbelow(True)
    ax.legend(frameon=False)
    for bars in (bars1,bars2):
        for b in bars:
            if b.get_height() >= 3:
                ax.annotate(f"{b.get_height():.1f}%", (b.get_x()+b.get_width()/2, b.get_height()),
                            xytext=(0,3), textcoords="offset points", ha="center", fontsize=9)
    fig.tight_layout()
    fig.savefig(folder / "hydrophobic_run_barplot.png", dpi=220)
    plt.close(fig)

    print("AD_only 연속 길이별 개수:", sorted(ad_counts.items()))
    print("all_three 연속 길이별 개수:", sorted(common_counts.items()))
    print("그래프 3개 저장 완료:", folder)


def save_similar_hydrophobicity_plot(csv_path, folder):
    """Compare average longest runs within 10pp hydrophobicity bands.

    Export a summary CSV and a fourth bar plot.
    Each data point is a UNIQUE PEPTIDE SEQUENCE, not a patient.
    """
    from statistics import mean, median
    rows = []
    with open(csv_path, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if row["detection_pattern"] not in ("AD_only", "all_three"):
                continue
            rows.append({
                "group": row["detection_pattern"],
                "percent": float(row["hydrophobic_percent"]),
                "run": int(row["max_hydrophobic_run_aa"]),
                "length": int(row["length_aa"]),
            })
    bands = list(range(0, 70, 10))
    result = []
    for low in bands:
        for group in ("AD_only", "all_three"):
            subset = [r for r in rows if r["group"] == group and low <= r["percent"] < low+10]
            n = len(subset)
            result.append({
                "bin": f"{low}-{low+10}%",
                "lower_pct": low,
                "upper_pct_exclusive": low+10,
                "group": group,
                "n_peptides": n,
                "mean_longest_run_aa": round(mean(r["run"] for r in subset), 5) if n else "",
                "median_longest_run_aa": median(r["run"] for r in subset) if n else "",
                "pct_run_ge_2": round(100*sum(r["run"] >= 2 for r in subset)/n, 3) if n else "",
                "mean_length_aa": round(mean(r["length"] for r in subset), 3) if n else "",
            })
    folder.mkdir(parents=True, exist_ok=True)
    with open(folder / "similar_percent_run_summary.csv", "w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(result[0]))
        writer.writeheader()
        writer.writerows(result)

    fig, ax = plt.subplots(figsize=(10.3, 6.2))
    xs = np.arange(len(bands))
    width = 0.38
    for group, delta, color, label in [
        ("AD_only", -width/2, "darkorange", "AD only"),
        ("all_three", width/2, "steelblue", "Detected in all three"),
    ]:
        stats = [next(r for r in result if r["lower_pct"] == low and r["group"] == group) for low in bands]
        values = [float(r["mean_longest_run_aa"]) if r["n_peptides"] else 0.0 for r in stats]
        bars = ax.bar(xs+delta, values, width, label=label, color=color, alpha=0.88, edgecolor="white")
        for bar, st in zip(bars, stats):
            if st["n_peptides"] < 10:
                bar.set_hatch("///")
                bar.set_edgecolor("white")
            ax.annotate(f'n={st["n_peptides"]}', (bar.get_x()+bar.get_width()/2, bar.get_height()),
                        xytext=(0, 4), textcoords="offset points", ha="center", fontsize=9)
    ax.set_xticks(xs, [f"{low}-{low+10}%" for low in bands])
    ax.set_xlabel("Hydrophobic residue percentage (10-percentage-point bins)")
    ax.set_ylabel("Mean longest hydrophobic run (aa)")
    ax.set_title("Longest hydrophobic run within similar hydrophobicity ranges")
    ax.yaxis.grid(True, alpha=0.25)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_ylim(0, 4.25)
    ax.legend(frameon=False, loc="upper left")
    fig.text(0.12, 0.017, "Hatched bars: n < 10 per group/bin; each observation is a peptide, not a patient.", fontsize=9)
    fig.tight_layout(rect=[0, 0.04, 1, 1])
    fig.savefig(folder / "similar_percent_run_barplot.png", dpi=220)
    plt.close(fig)
    print("비슷한 소수성 비율 구간별 연속 길이 비교 완료")


if __name__ == "__main__":
    folder = Path(__file__).resolve().parent / "results"
    data_path = folder / "peptide_summary.csv"
    groups = read_groups(data_path)
    save_plots(groups, folder)
    save_similar_hydrophobicity_plot(data_path, folder)
