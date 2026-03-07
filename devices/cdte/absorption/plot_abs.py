from pathlib import Path
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator


def parse_abs_file(path: Path):
    lambdas = []
    alphas = []
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        s = line.strip()
        if not s or s.startswith("/"):
            continue
        parts = s.replace(",", " ").split()
        if len(parts) < 2:
            continue
        try:
            lam = float(parts[0])
            alpha = float(parts[1])
        except ValueError:
            continue
        lambdas.append(lam)
        alphas.append(alpha)
    return lambdas, alphas


def main():
    parser = argparse.ArgumentParser(description="Plot alpha vs lambda from .abs files.")
    parser.add_argument(
        "--root",
        default=".",
        help="Root folder to search for .abs files (default: current directory).",
    )
    parser.add_argument(
        "--output",
        default="alpha_vs_lambda_all_abs.png",
        help="Output image path (default: alpha_vs_lambda_all_abs.png).",
    )
    parser.add_argument(
        "--files",
        nargs="+",
        default=None,
        help="Specific .abs files to plot (relative to --root or absolute paths).",
    )
    parser.add_argument(
        "--title",
        default=None,
        help="Optional figure title. Default is no title for publication-style figures.",
    )
    parser.add_argument(
        "--logy",
        action="store_true",
        help="Use log scale on y-axis.",
    )
    args = parser.parse_args()

    root = Path(args.root).resolve()
    if args.files:
        abs_files = []
        for value in args.files:
            p = Path(value)
            if not p.is_absolute():
                p = root / p
            abs_files.append(p.resolve())
        abs_files = sorted(abs_files)
    else:
        abs_files = sorted(root.rglob("*.abs"))
    if not abs_files:
        raise SystemExit(f"No .abs files found under: {root}")

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "font.size": 10,
            "axes.labelsize": 10,
            "axes.titlesize": 10,
            "legend.fontsize": 9,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "axes.linewidth": 1.0,
            "mathtext.default": "regular",
        }
    )

    fig, ax = plt.subplots(figsize=(5.2, 3.8))
    colors = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d"]
    plotted = 0

    for i, abs_file in enumerate(abs_files):
        lambdas, alphas = parse_abs_file(abs_file)
        if not lambdas:
            continue
        label = abs_file.stem if abs_file.parent == root else f"{abs_file.parent.name}/{abs_file.stem}"
        ax.plot(lambdas, alphas, linewidth=2.0, color=colors[i % len(colors)], label=label)
        plotted += 1

    if plotted == 0:
        raise SystemExit("No numeric data could be parsed from any .abs file.")

    if args.logy:
        ax.set_yscale("log")

    ax.set_xlabel(r"Wavelength, $\lambda$ (nm)")
    ax.set_ylabel(r"Absorption coefficient, $\alpha$ (m$^{-1}$)")
    if args.title:
        ax.set_title(args.title, pad=8)
    ax.grid(axis="y", linestyle=":", linewidth=0.8, color="0.75")
    ax.tick_params(which="both", direction="in", top=True, right=True, length=4, width=1.0)
    ax.tick_params(which="minor", length=2.5, width=0.8)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    if not args.logy:
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))

    ax.legend(frameon=False, loc="best", handlelength=2.0)
    fig.tight_layout()

    output = Path(args.output).resolve()
    fig.savefig(output, dpi=600)
    fig.savefig(output.with_suffix(".pdf"))
    print(output)


if __name__ == "__main__":
    main()
