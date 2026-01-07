import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def find_latest(pattern: str) -> str:
    files = glob.glob(pattern)
    if not files:
        raise FileNotFoundError(f"No files found for pattern: {pattern}")
    return max(files, key=os.path.getmtime)


def load_bo(results_dir: str):
    csv_path = find_latest(os.path.join(results_dir, "bo_single_run_*.csv"))
    print(f"[INFO] Using BO CSV: {csv_path}")
    df = pd.read_csv(csv_path)

    df["elapsed_time_s"] -= df["elapsed_time_s"].min()

    return df, csv_path


def load_ga(results_dir: str):
    csv_path = find_latest(os.path.join(results_dir, "ga_single_run_*.csv"))
    print(f"[INFO] Using GA CSV: {csv_path}")
    df = pd.read_csv(csv_path)

    df["elapsed_time_s"] -= df["elapsed_time_s"].min()

    return df, csv_path


def set_paper_style():
    plt.rcParams.update({
        "font.size": 10,
        "axes.labelsize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
        "lines.linewidth": 1.4,
        "lines.markersize": 4,
        "figure.dpi": 120,
    })

def compute_axis_limits(series_bo, series_ga, padding=0.1):
    vmin = min(series_bo.min(), series_ga.min())
    vmax = max(series_bo.max(), series_ga.max())

    # 防止上下限相同的情况
    if abs(vmax - vmin) < 1e-6:
        vmax = vmin + 1e-3

    # 按 padding 扩展，通常为 10%
    span = vmax - vmin
    vmin -= padding * span
    vmax += padding * span

    return vmin, vmax


def plot_three_metrics_vs_iter(df_bo, df_ga, save_dir):
    fig, axes = plt.subplots(3, 1, figsize=(4, 6), sharex=True)

    # ========== 自动计算纵轴范围 ==========
    y_min, y_max = compute_axis_limits(df_bo["y_best"], df_ga["y_best"])
    lcc_min, lcc_max = compute_axis_limits(df_bo["lcc_best"], df_ga["lcc_best_norm"])
    qmax_min, qmax_max = compute_axis_limits(df_bo["qmax_best"], df_ga["qmax_best_norm"])

    # ========== Plot Objective ==========
    ax = axes[0]
    ax.plot(df_bo["iter"], df_bo["y_best"], marker="o", label="BO")
    ax.plot(df_ga["gen"], df_ga["y_best"], marker="s", linestyle="--", label="GA")
    ax.set_ylabel("Best objective")
    ax.set_ylim(y_min, y_max)
    ax.legend(frameon=False)

    # ========== Plot LCC ==========
    ax = axes[1]
    ax.plot(df_bo["iter"], df_bo["lcc_best"], marker="o", label="BO")
    ax.plot(df_ga["gen"], df_ga["lcc_best_norm"], marker="s", linestyle="--", label="GA")
    ax.set_ylabel("Best LCC")
    ax.set_ylim(lcc_min, lcc_max)

    # ========== Plot Qmax ==========
    ax = axes[2]
    ax.plot(df_bo["iter"], df_bo["qmax_best"], marker="o", label="BO")
    ax.plot(df_ga["gen"], df_ga["qmax_best_norm"], marker="s", linestyle="--", label="GA")
    ax.set_ylabel("Best Qmax")
    ax.set_xlabel("Iteration")
    ax.set_ylim(qmax_min, qmax_max)

    for ax in axes:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    plt.tight_layout()

    out_png = os.path.join(save_dir, "combined_BO_GA_metrics_vs_iteration.png")
    out_pdf = os.path.join(save_dir, "combined_BO_GA_metrics_vs_iteration.pdf")
    fig.savefig(out_png, dpi=300)
    fig.savefig(out_pdf)

    print(f"[SAVE] {out_png}")
    print(f"[SAVE] {out_pdf}")


def plot_three_metrics_vs_time(df_bo, df_ga, save_dir):
    fig, axes = plt.subplots(3, 1, figsize=(4, 6), sharex=True)

    # ========== 自动计算纵轴范围 ==========
    y_min, y_max = compute_axis_limits(df_bo["y_best"], df_ga["y_best"])
    lcc_min, lcc_max = compute_axis_limits(df_bo["lcc_best"], df_ga["lcc_best_norm"])
    qmax_min, qmax_max = compute_axis_limits(df_bo["qmax_best"], df_ga["qmax_best_norm"])

    # ========== Objective ==========
    ax = axes[0]
    ax.plot(df_bo["elapsed_time_s"], df_bo["y_best"], marker="o", label="BO")
    ax.plot(df_ga["elapsed_time_s"], df_ga["y_best"], marker="s", linestyle="--", label="GA")
    ax.set_ylabel("Best objective")
    ax.set_ylim(y_min, y_max)
    ax.legend(frameon=False)

    # ========== LCC ==========
    ax = axes[1]
    ax.plot(df_bo["elapsed_time_s"], df_bo["lcc_best"], marker="o", label="BO")
    ax.plot(df_ga["elapsed_time_s"], df_ga["lcc_best_norm"], marker="s", linestyle="--", label="GA")
    ax.set_ylabel("Best LCC")
    ax.set_ylim(lcc_min, lcc_max)

    # ========== Qmax ==========
    ax = axes[2]
    ax.plot(df_bo["elapsed_time_s"], df_bo["qmax_best"], marker="o", label="BO")
    ax.plot(df_ga["elapsed_time_s"], df_ga["qmax_best_norm"], marker="s", linestyle="--", label="GA")
    ax.set_ylabel("Best Qmax")
    ax.set_xlabel("Time (s)")
    ax.set_ylim(qmax_min, qmax_max)

    for ax in axes:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    plt.tight_layout()

    out_png = os.path.join(save_dir, "combined_BO_GA_metrics_vs_time.png")
    out_pdf = os.path.join(save_dir, "combined_BO_GA_metrics_vs_time.pdf")
    fig.savefig(out_png, dpi=300)
    fig.savefig(out_pdf)

    print(f"[SAVE] {out_png}")
    print(f"[SAVE] {out_pdf}")
    
def plot_bo_metrics_vs_iteration(df_bo, save_dir):
    fig, axes = plt.subplots(3, 1, figsize=(4, 6), sharex=True)

    # Objective
    axes[0].plot(df_bo["iter"], df_bo["y_best"], marker="o")
    axes[0].set_ylabel("Objective")
    axes[0].spines["top"].set_visible(False)
    axes[0].spines["right"].set_visible(False)

    # LCC
    axes[1].plot(df_bo["iter"], df_bo["lcc_best"], marker="o")
    axes[1].set_ylabel("LCC")
    axes[1].spines["top"].set_visible(False)
    axes[1].spines["right"].set_visible(False)

    # Qmax
    axes[2].plot(df_bo["iter"], df_bo["qmax_best"], marker="o")
    axes[2].set_ylabel("Qmax")
    axes[2].set_xlabel("Iteration")
    axes[2].spines["top"].set_visible(False)
    axes[2].spines["right"].set_visible(False)

    plt.tight_layout()

    out_png = os.path.join(save_dir, "BO_metrics_vs_iteration.png")
    out_pdf = os.path.join(save_dir, "BO_metrics_vs_iteration.pdf")
    fig.savefig(out_png, dpi=300)
    fig.savefig(out_pdf)
    print(f"[SAVE] {out_png}")
    print(f"[SAVE] {out_pdf}")

def plot_ga_metrics_vs_iteration(df_ga, save_dir):
    fig, axes = plt.subplots(3, 1, figsize=(4, 6), sharex=True)

    # Objective
    axes[0].plot(df_ga["gen"], df_ga["y_best"], marker="s")
    axes[0].set_ylabel("Objective")
    axes[0].spines["top"].set_visible(False)
    axes[0].spines["right"].set_visible(False)

    # LCC
    axes[1].plot(df_ga["gen"], df_ga["lcc_best_norm"], marker="s")
    axes[1].set_ylabel("LCC")
    axes[1].spines["top"].set_visible(False)
    axes[1].spines["right"].set_visible(False)

    # Qmax
    axes[2].plot(df_ga["gen"], df_ga["qmax_best_norm"], marker="s")
    axes[2].set_ylabel("Qmax")
    axes[2].set_xlabel("Generation")
    axes[2].spines["top"].set_visible(False)
    axes[2].spines["right"].set_visible(False)

    plt.tight_layout()

    out_png = os.path.join(save_dir, "GA_metrics_vs_iteration.png")
    out_pdf = os.path.join(save_dir, "GA_metrics_vs_iteration.pdf")
    fig.savefig(out_png, dpi=300)
    fig.savefig(out_pdf)
    print(f"[SAVE] {out_png}")
    print(f"[SAVE] {out_pdf}")


def main():

    results_dir = os.environ.get("RESULTS_DIR", "results")
    os.makedirs(results_dir, exist_ok=True)

    df_bo, _ = load_bo(results_dir)
    df_ga, _ = load_ga(results_dir)

    set_paper_style()

    plot_bo_metrics_vs_iteration(df_bo, results_dir)
    plot_ga_metrics_vs_iteration(df_ga, results_dir)
    plot_three_metrics_vs_iter(df_bo, df_ga, results_dir)
    plot_three_metrics_vs_time(df_bo, df_ga, results_dir)

    print("\n[Done] All comparison plots generated.\n")


if __name__ == "__main__":
    main()
