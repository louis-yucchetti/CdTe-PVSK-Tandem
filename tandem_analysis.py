#!/usr/bin/env python3
"""Build and plot a 2T PVSK/CdTe tandem from SCAPS IV/QE exports."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Tuple

import matplotlib.pyplot as plt
import numpy as np


def apply_publication_style() -> None:
    """Apply a clean journal-style matplotlib theme."""
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "dejavuserif",
            "font.size": 11,
            "axes.labelsize": 12,
            "axes.titlesize": 12,
            "axes.linewidth": 1.0,
            "legend.fontsize": 10,
            "legend.frameon": False,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.major.size": 5,
            "ytick.major.size": 5,
            "xtick.minor.size": 2.5,
            "ytick.minor.size": 2.5,
            "xtick.major.width": 1.0,
            "ytick.major.width": 1.0,
            "xtick.minor.width": 0.8,
            "ytick.minor.width": 0.8,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.03,
        }
    )


def style_axis(ax: plt.Axes) -> None:
    """Apply consistent axis styling for publication figures."""
    ax.minorticks_on()
    ax.tick_params(which="both", top=True, right=True)
    ax.grid(True, which="major", linestyle=":", linewidth=0.6, alpha=0.35)


def read_scaps_two_columns(path: Path) -> Tuple[np.ndarray, np.ndarray]:
    """Read first two numeric columns from SCAPS text output."""
    x_vals = []
    y_vals = []

    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            try:
                x_val = float(parts[0])
                y_val = float(parts[1])
            except ValueError:
                continue
            x_vals.append(x_val)
            y_vals.append(y_val)

    if not x_vals:
        raise ValueError(f"No numeric data found in {path}")

    return np.asarray(x_vals, dtype=float), np.asarray(y_vals, dtype=float)


def sort_unique(x: np.ndarray, y: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Sort by x and keep first occurrence when duplicates are present."""
    order = np.argsort(x)
    x_sorted = x[order]
    y_sorted = y[order]
    x_unique, unique_idx = np.unique(x_sorted, return_index=True)
    return x_unique, y_sorted[unique_idx]


def interpolate_y_at_x(x: np.ndarray, y: np.ndarray, x_target: float) -> float:
    """Interpolate y(x_target) if x_target is inside bounds, else NaN."""
    if x_target < np.min(x) or x_target > np.max(x):
        return float("nan")
    return float(np.interp(x_target, x, y))


def interpolate_x_at_y(x: np.ndarray, y: np.ndarray, y_target: float) -> float:
    """Interpolate x(y_target) after sorting by y."""
    y_sorted, x_sorted = sort_unique(y, x)
    if y_target < np.min(y_sorted) or y_target > np.max(y_sorted):
        return float("nan")
    return float(np.interp(y_target, y_sorted, x_sorted))


def extract_pv_metrics(
    voltage: np.ndarray, current: np.ndarray, pin_mw_cm2: float = 100.0
) -> Dict[str, float]:
    """Extract Voc, Jsc, FF, eta, Vmpp, Jmpp from a JV curve."""
    v, j = sort_unique(voltage, current)

    jsc_signed = interpolate_y_at_x(v, j, 0.0)
    jsc = abs(jsc_signed) if np.isfinite(jsc_signed) else float("nan")
    voc = interpolate_x_at_y(v, j, 0.0)

    power_gen = -v * j  # mW/cm2 for V in V and J in mA/cm2
    idx_mpp = int(np.argmax(power_gen))
    pmax = float(power_gen[idx_mpp])
    vmpp = float(v[idx_mpp])
    jmpp = abs(float(j[idx_mpp]))

    ff = float("nan")
    if np.isfinite(voc) and np.isfinite(jsc) and voc > 0 and jsc > 0:
        ff = 100.0 * pmax / (voc * jsc)

    eta = 100.0 * pmax / pin_mw_cm2

    return {
        "voc": voc,
        "jsc": jsc,
        "ff": ff,
        "eta": eta,
        "vmpp": vmpp,
        "jmpp": jmpp,
        "pmax": pmax,
    }


def voltage_from_current(
    voltage: np.ndarray, current: np.ndarray, current_grid: np.ndarray
) -> np.ndarray:
    """Interpolate V(J) on a target current grid."""
    j_sorted, v_sorted = sort_unique(current, voltage)
    if current_grid.min() < j_sorted.min() or current_grid.max() > j_sorted.max():
        raise ValueError(
            "Current grid extends outside available JV data range. "
            "Provide broader JV data (including reverse bias) for strict 2T construction."
        )
    return np.interp(current_grid, j_sorted, v_sorted)


def construct_tandem_2t(
    v_top: np.ndarray,
    j_top: np.ndarray,
    v_bottom: np.ndarray,
    j_bottom: np.ndarray,
    pin_mw_cm2: float = 100.0,
    n_points: int = 3000,
) -> Dict[str, object]:
    """Construct full 2T tandem JV by current-domain series connection."""
    j_min_common = max(np.min(j_top), np.min(j_bottom))
    j_max_common = min(np.max(j_top), np.max(j_bottom))
    if j_min_common >= j_max_common:
        raise ValueError("No overlapping current range between top and bottom JV curves.")

    current_grid = np.linspace(j_min_common, j_max_common, n_points)
    v_top_at_j = voltage_from_current(v_top, j_top, current_grid)
    v_bottom_at_j = voltage_from_current(v_bottom, j_bottom, current_grid)
    v_tandem = v_top_at_j + v_bottom_at_j

    metrics = extract_pv_metrics(v_tandem, current_grid, pin_mw_cm2=pin_mw_cm2)
    idx_vmin = int(np.argmin(v_tandem))
    status: Dict[str, float | str] = {
        "jsc_source": "tandem_curve",
        "voc_source": "tandem_curve",
        "v_tandem_min": float(v_tandem[idx_vmin]),
        "j_at_v_tandem_min": float(current_grid[idx_vmin]),
        "v_top_at_v_tandem_min": float(v_top_at_j[idx_vmin]),
        "v_bottom_at_v_tandem_min": float(v_bottom_at_j[idx_vmin]),
    }

    if not np.isfinite(metrics["jsc"]):
        jsc_top = abs(interpolate_y_at_x(v_top, j_top, 0.0))
        jsc_bottom = abs(interpolate_y_at_x(v_bottom, j_bottom, 0.0))
        if np.isfinite(jsc_top) and np.isfinite(jsc_bottom):
            metrics["jsc"] = min(jsc_top, jsc_bottom)
            status["jsc_source"] = "min_subcell_jsc"

    if not np.isfinite(metrics["voc"]):
        voc_top = interpolate_x_at_y(v_top, j_top, 0.0)
        voc_bottom = interpolate_x_at_y(v_bottom, j_bottom, 0.0)
        if np.isfinite(voc_top) and np.isfinite(voc_bottom):
            metrics["voc"] = voc_top + voc_bottom
            status["voc_source"] = "sum_subcell_voc"

    if (
        not np.isfinite(metrics["ff"])
        and np.isfinite(metrics["voc"])
        and np.isfinite(metrics["jsc"])
        and metrics["voc"] > 0
        and metrics["jsc"] > 0
    ):
        metrics["ff"] = 100.0 * metrics["pmax"] / (metrics["voc"] * metrics["jsc"])

    return {"voltage": v_tandem, "current": current_grid, "metrics": metrics, "status": status}


def construct_tandem_4t_equivalent(
    metrics_top: Dict[str, float], metrics_bottom: Dict[str, float], pin_mw_cm2: float = 100.0
) -> Dict[str, float]:
    """Build equivalent 4T metrics from independent subcell power addition."""
    pmax = metrics_top["pmax"] + metrics_bottom["pmax"]
    voc = metrics_top["voc"] + metrics_bottom["voc"]
    jsc = metrics_top["jsc"] + metrics_bottom["jsc"]
    vmpp = metrics_top["vmpp"] + metrics_bottom["vmpp"]
    jmpp = pmax / vmpp if np.isfinite(vmpp) and vmpp != 0 else float("nan")
    ff = (
        100.0 * pmax / (voc * jsc)
        if np.isfinite(voc) and np.isfinite(jsc) and voc > 0 and jsc > 0
        else float("nan")
    )
    eta = 100.0 * pmax / pin_mw_cm2

    return {
        "voc": voc,
        "jsc": jsc,
        "ff": ff,
        "eta": eta,
        "vmpp": vmpp,
        "jmpp": jmpp,
        "pmax": pmax,
    }


def construct_tandem_optical_eqe(
    wavelength_top: np.ndarray,
    eqe_top_percent: np.ndarray,
    wavelength_bottom: np.ndarray,
    eqe_bottom_percent: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Build interpolated top/bottom EQE and optical tandem EQE."""
    wl_top, eqe_top = sort_unique(wavelength_top, np.clip(eqe_top_percent, 0.0, 100.0))
    wl_bottom, eqe_bottom = sort_unique(
        wavelength_bottom, np.clip(eqe_bottom_percent, 0.0, 100.0)
    )

    wl_common = np.union1d(wl_top, wl_bottom)
    eqe_top_interp = np.interp(wl_common, wl_top, eqe_top)
    eqe_bottom_interp = np.interp(wl_common, wl_bottom, eqe_bottom)
    # Optical stack response: top + transmitted bottom (clip to 100%).
    eqe_tandem_optical = np.clip(eqe_top_interp + eqe_bottom_interp, 0.0, 100.0)

    return wl_common, eqe_top_interp, eqe_bottom_interp, eqe_tandem_optical


def save_metrics_csv(
    out_path: Path,
    metrics_cdte_unfiltered: Dict[str, float],
    metrics_cdte_filtered: Dict[str, float],
    metrics_top: Dict[str, float],
    metrics_tandem_2t: Dict[str, float],
    metrics_tandem_4t: Dict[str, float],
) -> None:
    header = (
        "cell,model,Voc_V,Jsc_mAcm2,FF_percent,eta_percent,"
        "Vmpp_V,Jmpp_mAcm2,Pmax_mWcm2"
    )
    lines = [header]

    rows = [
        ("CdTe_unfiltered", "single_junction", metrics_cdte_unfiltered),
        ("CdTe_filtered", "single_junction_filtered", metrics_cdte_filtered),
        ("PVSK_top", "single_junction", metrics_top),
        ("Tandem_2T", "series_current_matched", metrics_tandem_2t),
        ("Tandem_4T", "equivalent_power_addition", metrics_tandem_4t),
    ]
    for name, model, row in rows:
        lines.append(
            ",".join(
                [
                    name,
                    model,
                    f"{row['voc']:.6f}",
                    f"{row['jsc']:.6f}",
                    f"{row['ff']:.6f}",
                    f"{row['eta']:.6f}",
                    f"{row['vmpp']:.6f}",
                    f"{row['jmpp']:.6f}",
                    f"{row['pmax']:.6f}",
                ]
            )
        )

    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot JV and EQE curves and extract PV parameters for PVSK/CdTe tandem."
    )
    parser.add_argument("--pvsk-iv", default="PVSK_louis.iv", help="Top-cell JV file")
    parser.add_argument(
        "--cdte-iv", default="CdTe_filtered_louis.iv", help="Bottom-cell JV file"
    )
    parser.add_argument(
        "--cdte-unfiltered-iv",
        default="CdTe_louis.iv",
        help="Unfiltered CdTe JV file",
    )
    parser.add_argument("--pvsk-qe", default="PVSK_louis.qe", help="Top-cell EQE file")
    parser.add_argument(
        "--cdte-qe", default="CdTe_filtered_louis.qe", help="Bottom-cell EQE file"
    )
    parser.add_argument(
        "--cdte-unfiltered-qe",
        default="CdTe_louis.qe",
        help="Unfiltered CdTe EQE file for comparison plotting",
    )
    parser.add_argument(
        "--pin-mw-cm2",
        type=float,
        default=100.0,
        help="Incident power density in mW/cm^2 (default: 100)",
    )
    parser.add_argument(
        "--outdir",
        default=".",
        help="Output directory for plots and csv files",
    )
    parser.add_argument(
        "--eqe-min-nm",
        type=float,
        default=250.0,
        help="Minimum wavelength (nm) shown in EQE plot (default: 250)",
    )
    parser.add_argument(
        "--eqe-max-nm",
        type=float,
        default=900.0,
        help="Maximum wavelength (nm) shown in EQE plot (default: 900)",
    )
    args = parser.parse_args()

    if args.eqe_min_nm >= args.eqe_max_nm:
        raise ValueError("--eqe-min-nm must be strictly smaller than --eqe-max-nm.")

    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    apply_publication_style()

    pvsk_iv_path = Path(args.pvsk_iv)
    cdte_filtered_iv_path = Path(args.cdte_iv)
    cdte_unfiltered_iv_path = Path(args.cdte_unfiltered_iv)
    pvsk_qe_path = Path(args.pvsk_qe)
    cdte_qe_path = Path(args.cdte_qe)
    cdte_unfiltered_qe_path = Path(args.cdte_unfiltered_qe)

    v_top, j_top = read_scaps_two_columns(pvsk_iv_path)
    v_cdte_filtered, j_cdte_filtered = read_scaps_two_columns(cdte_filtered_iv_path)
    v_cdte_unfiltered, j_cdte_unfiltered = read_scaps_two_columns(cdte_unfiltered_iv_path)
    wl_top, eqe_top = read_scaps_two_columns(pvsk_qe_path)
    wl_bottom, eqe_bottom = read_scaps_two_columns(cdte_qe_path)
    wl_cdte_unfiltered, eqe_cdte_unfiltered = read_scaps_two_columns(cdte_unfiltered_qe_path)

    metrics_top = extract_pv_metrics(v_top, j_top, args.pin_mw_cm2)
    metrics_cdte_filtered = extract_pv_metrics(v_cdte_filtered, j_cdte_filtered, args.pin_mw_cm2)
    metrics_cdte_unfiltered = extract_pv_metrics(
        v_cdte_unfiltered, j_cdte_unfiltered, args.pin_mw_cm2
    )

    tandem = construct_tandem_2t(
        v_top=v_top,
        j_top=j_top,
        v_bottom=v_cdte_filtered,
        j_bottom=j_cdte_filtered,
        pin_mw_cm2=args.pin_mw_cm2,
    )
    v_tandem = tandem["voltage"]
    j_tandem = tandem["current"]
    metrics_tandem_2t = tandem["metrics"]
    tandem_status = tandem["status"]
    metrics_tandem_4t = construct_tandem_4t_equivalent(
        metrics_top, metrics_cdte_filtered, args.pin_mw_cm2
    )

    (
        wl_common,
        eqe_top_interp,
        eqe_bottom_interp,
        eqe_tandem_optical,
    ) = construct_tandem_optical_eqe(
        wl_top, eqe_top, wl_bottom, eqe_bottom
    )
    eqe_mask = (wl_common >= args.eqe_min_nm) & (wl_common <= args.eqe_max_nm)
    wl_eqe_plot = wl_common[eqe_mask]
    eqe_top_plot = eqe_top_interp[eqe_mask]
    eqe_bottom_plot = eqe_bottom_interp[eqe_mask]
    eqe_tandem_plot = eqe_tandem_optical[eqe_mask]
    wl_cdte_unfiltered_sorted, eqe_cdte_unfiltered_sorted = sort_unique(
        wl_cdte_unfiltered, np.clip(eqe_cdte_unfiltered, 0.0, 100.0)
    )
    eqe_cdte_unfiltered_plot = np.interp(
        wl_eqe_plot, wl_cdte_unfiltered_sorted, eqe_cdte_unfiltered_sorted
    )
    if wl_eqe_plot.size < 2:
        raise ValueError(
            "Insufficient EQE points inside selected wavelength window. "
            "Check --eqe-min-nm/--eqe-max-nm."
        )

    # JV plot
    fig_jv, ax_jv = plt.subplots(figsize=(7.2, 5.0))
    ax_jv.plot(v_top, j_top, lw=1.8, color="#1f4e79", label="PVSK top")
    ax_jv.plot(
        v_cdte_unfiltered,
        j_cdte_unfiltered,
        lw=1.6,
        ls="--",
        color="#a56a00",
        label="CdTe unfiltered",
    )
    ax_jv.plot(
        v_cdte_filtered, j_cdte_filtered, lw=1.8, color="#7a3e00", label="CdTe filtered"
    )
    ax_jv.plot(v_tandem, j_tandem, lw=2.3, color="#0f7f54", label="Constructed 2T tandem")
    ax_jv.axhline(0.0, color="black", lw=0.9, alpha=0.7)
    ax_jv.set_xlabel("Voltage (V)")
    ax_jv.set_ylabel("Current Density (mA cm$^{-2}$)")
    ax_jv.set_title("J-V Characteristics of Subcells and Constructed 2T Tandem")
    x_right = 1.02 * max(
        np.max(v_top), np.max(v_cdte_unfiltered), np.max(v_cdte_filtered), np.max(v_tandem)
    )
    ax_jv.set_xlim(0.0, x_right)
    y_all = np.concatenate([j_top, j_cdte_unfiltered, j_cdte_filtered, j_tandem])
    y_min = np.min(y_all)
    y_pad = 0.04 * abs(y_min) if y_min < 0 else 1.0
    # Show only power-generating quadrant in current: J <= 0.
    ax_jv.set_ylim(y_min - y_pad, 0.0)
    style_axis(ax_jv)
    ax_jv.legend(loc="best")
    fig_jv.tight_layout()
    jv_plot_path = outdir / "tandem_jv_curve.png"
    fig_jv.savefig(jv_plot_path, dpi=600)
    plt.close(fig_jv)

    # EQE plot
    fig_eqe, ax_eqe = plt.subplots(figsize=(7.2, 5.0))
    ax_eqe.plot(
        wl_eqe_plot,
        eqe_top_plot,
        lw=1.9,
        color="#1f4e79",
        label="PVSK top EQE",
    )
    ax_eqe.plot(
        wl_eqe_plot,
        eqe_bottom_plot,
        lw=1.9,
        color="#7a3e00",
        label="CdTe bottom EQE (filtered)",
    )
    ax_eqe.plot(
        wl_eqe_plot,
        eqe_cdte_unfiltered_plot,
        lw=1.7,
        ls=":",
        color="#a56a00",
        label="CdTe bottom EQE (unfiltered)",
    )
    ax_eqe.plot(
        wl_eqe_plot,
        eqe_tandem_plot,
        lw=2.5,
        color="#7d1f6f",
        label="Constructed tandem optical EQE",
    )
    ax_eqe.set_xlabel("Wavelength (nm)")
    ax_eqe.set_ylabel("EQE (%)")
    ax_eqe.set_title("Constructed Tandem Optical EQE")
    ax_eqe.set_xlim(args.eqe_min_nm, args.eqe_max_nm)
    ax_eqe.set_ylim(0.0, 100.0)
    style_axis(ax_eqe)
    ax_eqe.legend(loc="best")
    fig_eqe.tight_layout()
    eqe_plot_path = outdir / "tandem_eqe_curve.png"
    fig_eqe.savefig(eqe_plot_path, dpi=600)
    plt.close(fig_eqe)

    # 4T plot 1: subcell power curves and additive 4T MPP
    v_top_s, j_top_s = sort_unique(v_top, j_top)
    v_cdte_f_s, j_cdte_f_s = sort_unique(v_cdte_filtered, j_cdte_filtered)
    p_top_s = -v_top_s * j_top_s
    p_cdte_f_s = -v_cdte_f_s * j_cdte_f_s

    fig_4t_curves, ax_4t_curves = plt.subplots(figsize=(7.2, 5.0))
    ax_4t_curves.plot(v_top_s, p_top_s, lw=2.0, color="#1f4e79", label="PVSK top power")
    ax_4t_curves.plot(
        v_cdte_f_s, p_cdte_f_s, lw=2.0, color="#7a3e00", label="CdTe filtered power"
    )
    ax_4t_curves.axhline(
        metrics_tandem_4t["pmax"],
        lw=1.6,
        ls="--",
        color="#0f7f54",
        label=f"4T total Pmax = {metrics_tandem_4t['pmax']:.2f} mW cm$^{{-2}}$",
    )
    ax_4t_curves.scatter(
        [metrics_top["vmpp"]],
        [metrics_top["pmax"]],
        color="#1f4e79",
        s=28,
        zorder=3,
    )
    ax_4t_curves.scatter(
        [metrics_cdte_filtered["vmpp"]],
        [metrics_cdte_filtered["pmax"]],
        color="#7a3e00",
        s=28,
        zorder=3,
    )
    ax_4t_curves.set_xlabel("Subcell Voltage (V)")
    ax_4t_curves.set_ylabel("Power Density (mW cm$^{-2}$)")
    ax_4t_curves.set_title("4T Tandem: Independent Subcell Power and Total Power")
    ax_4t_curves.set_xlim(0.0, 1.02 * max(np.max(v_top_s), np.max(v_cdte_f_s)))
    p_all = np.concatenate([p_top_s, p_cdte_f_s, np.array([metrics_tandem_4t["pmax"]])])
    p_pad = 0.05 * (np.max(p_all) - np.min(p_all))
    ax_4t_curves.set_ylim(np.min(p_all) - p_pad, np.max(p_all) + p_pad)
    style_axis(ax_4t_curves)
    ax_4t_curves.legend(loc="best")
    fig_4t_curves.tight_layout()
    plot_4t_curves_path = outdir / "tandem_4t_power_curves.png"
    fig_4t_curves.savefig(plot_4t_curves_path, dpi=600)
    plt.close(fig_4t_curves)

    # 4T plot 2: total power map P_total(V_top, V_bottom_filtered)
    p_total_map = p_cdte_f_s[:, np.newaxis] + p_top_s[np.newaxis, :]
    fig_4t_map, ax_4t_map = plt.subplots(figsize=(7.4, 5.6))
    contour = ax_4t_map.contourf(
        v_top_s,
        v_cdte_f_s,
        p_total_map,
        levels=35,
        cmap="viridis",
    )
    cbar = fig_4t_map.colorbar(contour, ax=ax_4t_map, pad=0.02)
    cbar.set_label("Total 4T Power Density (mW cm$^{-2}$)")
    ax_4t_map.scatter(
        [metrics_top["vmpp"]],
        [metrics_cdte_filtered["vmpp"]],
        color="white",
        edgecolor="black",
        s=44,
        zorder=4,
        label=(
            "4T MPP operating point\n"
            f"(Vtop={metrics_top['vmpp']:.3f} V, Vcdte={metrics_cdte_filtered['vmpp']:.3f} V)"
        ),
    )
    ax_4t_map.set_xlabel("PVSK Top Voltage (V)")
    ax_4t_map.set_ylabel("CdTe Filtered Voltage (V)")
    ax_4t_map.set_title("4T Tandem Power Map")
    ax_4t_map.set_xlim(0.0, np.max(v_top_s))
    ax_4t_map.set_ylim(0.0, np.max(v_cdte_f_s))
    style_axis(ax_4t_map)
    ax_4t_map.legend(loc="upper right")
    fig_4t_map.tight_layout()
    plot_4t_map_path = outdir / "tandem_4t_power_map.png"
    fig_4t_map.savefig(plot_4t_map_path, dpi=600)
    plt.close(fig_4t_map)

    # Save curve CSVs
    np.savetxt(
        outdir / "tandem_jv_curve.csv",
        np.column_stack([v_tandem, j_tandem]),
        delimiter=",",
        header="Voltage_V,CurrentDensity_mAcm2",
        comments="",
        fmt="%.8f",
    )
    np.savetxt(
        outdir / "tandem_eqe_curve.csv",
        np.column_stack([wl_eqe_plot, eqe_top_plot, eqe_bottom_plot, eqe_tandem_plot]),
        delimiter=",",
        header=(
            "Wavelength_nm,"
            "PVSK_EQE_percent,"
            "CdTe_Filtered_EQE_percent,"
            "Tandem_Optical_EQE_percent"
        ),
        comments="",
        fmt="%.8f",
    )
    save_metrics_csv(
        outdir / "tandem_metrics.csv",
        metrics_cdte_unfiltered,
        metrics_cdte_filtered,
        metrics_top,
        metrics_tandem_2t,
        metrics_tandem_4t,
    )

    print("PV parameters (SCAPS sign convention converted to |Jsc|):")
    print("Cell              Voc (V)    Jsc (mA/cm2)   FF (%)    eta (%)")
    print(
        f"CdTe unfiltered   {metrics_cdte_unfiltered['voc']:.4f}      "
        f"{metrics_cdte_unfiltered['jsc']:.4f}       "
        f"{metrics_cdte_unfiltered['ff']:.2f}      {metrics_cdte_unfiltered['eta']:.2f}"
    )
    print(
        f"CdTe filtered     {metrics_cdte_filtered['voc']:.4f}      "
        f"{metrics_cdte_filtered['jsc']:.4f}       "
        f"{metrics_cdte_filtered['ff']:.2f}      {metrics_cdte_filtered['eta']:.2f}"
    )
    print(
        f"PVSK top          {metrics_top['voc']:.4f}      {metrics_top['jsc']:.4f}       "
        f"{metrics_top['ff']:.2f}      {metrics_top['eta']:.2f}"
    )
    print(
        f"Tandem 2T         {metrics_tandem_2t['voc']:.4f}      {metrics_tandem_2t['jsc']:.4f}       "
        f"{metrics_tandem_2t['ff']:.2f}      {metrics_tandem_2t['eta']:.2f}"
    )
    print(
        f"Tandem 4T*        {metrics_tandem_4t['voc']:.4f}      {metrics_tandem_4t['jsc']:.4f}       "
        f"{metrics_tandem_4t['ff']:.2f}      {metrics_tandem_4t['eta']:.2f}"
    )
    print("")
    print("Files saved:")
    print(f"  {jv_plot_path}")
    print(f"  {eqe_plot_path}")
    print(f"  {plot_4t_curves_path}")
    print(f"  {plot_4t_map_path}")
    print(f"  {outdir / 'tandem_metrics.csv'}")
    print(f"  {outdir / 'tandem_jv_curve.csv'}")
    print(f"  {outdir / 'tandem_eqe_curve.csv'}")
    print("")
    if tandem_status["jsc_source"] != "tandem_curve" or tandem_status["voc_source"] != "tandem_curve":
        print(
            "JV note: tandem metric fallback used "
            f"(Jsc source: {tandem_status['jsc_source']}, "
            f"Voc source: {tandem_status['voc_source']})."
        )
    else:
        print("JV note: tandem Voc/Jsc extracted directly from the series-combined tandem curve.")
    print("")
    print(
        "Tandem JV span note: minimum constructed tandem voltage is "
        f"{float(tandem_status['v_tandem_min']):.3f} V at "
        f"J = {float(tandem_status['j_at_v_tandem_min']):.3f} mA/cm2 "
        f"(Vtop = {float(tandem_status['v_top_at_v_tandem_min']):.3f} V, "
        f"Vbottom = {float(tandem_status['v_bottom_at_v_tandem_min']):.3f} V)."
    )
    if float(tandem_status["v_tandem_min"]) > 0.0:
        print(
            "The constructed tandem curve stays above 0 V because the available reverse-bias "
            "range does not yet drive Vtop + Vbottom below zero over the shared current range."
        )
    else:
        print(
            "The constructed tandem curve now crosses 0 V, so tandem Jsc can be extracted "
            "directly from the constructed 2T JV."
        )
    print("")
    print(
        "EQE note: the EQE figure shows PVSK EQE, filtered CdTe EQE, unfiltered CdTe "
        "EQE as a dotted comparison, and the constructed tandem optical response "
        "(top + filtered bottom, clipped to 100%) in the requested wavelength window."
    )
    print(
        "4T note: starred 4T metrics are equivalent aggregated values from independent "
        "power addition (PVSK + filtered CdTe), not a unique 4T terminal JV."
    )


if __name__ == "__main__":
    main()
