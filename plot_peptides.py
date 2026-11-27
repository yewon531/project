import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # 서버에서 화면을 띄우지 않고 PNG 파일로 저장한다.
import matplotlib.pyplot as plt


def read_groups(path):
    # 그룹마다 소수성 비율과 최대 연속 길이를 같은 순서로 담는다.
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

            percent = float(row["hydrophobic_percent"])
            run = int(row["max_hydrophobic_run_aa"])

            groups[group]["percent"].append(percent)
            groups[group]["run"].append(run)

    for group in groups:
        if len(groups[group]["percent"]) == 0:
            raise ValueError(f"{group}에 해당하는 서열이 없습니다.")

    return groups


def save_plots(groups, folder):
    ad = groups["AD_only"]
    common = groups["all_three"]

    ad_label = f"AD only (n={len(ad['percent'])})"
    common_label = f"All three (n={len(common['percent'])})"

    # 1. 상자그림: 두 그룹의 소수성 비율 분포를 비교한다.
    fig, ax = plt.subplots()
    ax.boxplot([ad["percent"], common["percent"]])
    ax.set_xticks([1, 2])
    ax.set_xticklabels([ad_label, common_label])
    ax.set_ylabel("Hydrophobic residues (%)")
    ax.set_title("Hydrophobic proportion by detection group")
    fig.tight_layout()
    fig.savefig(folder / "hydrophobic_percent_boxplot.png", dpi=200)
    plt.close(fig)

    # 2. 산점도: 점 하나는 고유 펩타이드 서열 하나다.
    fig, ax = plt.subplots()
    ax.scatter(common["percent"], common["run"],
               label=common_label, color="steelblue", alpha=0.4)
    ax.scatter(ad["percent"], ad["run"],
               label=ad_label, color="darkorange", marker="^", alpha=0.5)
    ax.set_xlabel("Hydrophobic residues (%)")
    ax.set_ylabel("Longest hydrophobic run (aa)")
    ax.set_title("Hydrophobic proportion and longest run")
    ax.legend()
    fig.tight_layout()
    fig.savefig(folder / "hydrophobic_run_scatter.png", dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    folder = Path(__file__).resolve().parent / "results"
    groups = read_groups(folder / "peptide_summary.csv")
    save_plots(groups, folder)
    print(f"그래프 2개를 저장했습니다: {folder}")
