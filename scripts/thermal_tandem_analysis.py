#!/usr/bin/env python3
"""Analyze temperature-dependent PVSK/CdTe sub-cell and tandem behavior."""

from __future__ import annotations

import argparse
import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence

import matplotlib.pyplot as plt
import numpy as np

from tandem_analysis import (
    apply_publication_style,
    construct_tandem_2t,
    construct_tandem_4t_equivalent,
    extract_pv_metrics,
    interpolate_y_at_x,
    resolve_project_path,
    style_axis,
)


METRIC_META = {
    "voc": {"label": "Voc (V)", "unit": "V"},
    "jsc": {"label": "Jsc (mA cm$^{-2}$)", "unit": "mA/cm2"},
    "ff": {"label": "FF (%)", "unit": "%"},
    "eta": {"label": "Efficiency (%)", "unit": "%"},
}

SUMMARY_METRIC_META = {
    **METRIC_META,
    "pmax": {"label": "Pmax (mW cm$^{-2}$)", "unit": "mW/cm2"},
}

COEFFICIENT_META = {
    "voc": {"symbol": "beta", "absolute_scale": 1000.0, "absolute_unit": "mV/K"},
    "pmax": {"symbol": "gamma", "absolute_scale": 1.0, "absolute_unit": "mW/cm2/K"},
}

RECOMB_COLUMNS = [
    "j_total_rec(mA/cm2)",
    "j_total_gen(mA/cm2)",
    "jbulk(mA/cm2)",
    "jifr(mA/cm2)",
    "jminor_left(mA/cm2)",
    "jminor_right(mA/cm2)",
    "j_SRH(mA/cm2)",
    "j_Radiative(mA/cm2)",
    "j_Auger(mA/cm2)",
]


@dataclass
class SweepBlock:
    temperature_k: float
    headers: List[str]
    data: np.ndarray

    @property
    def column_index(self) -> Dict[str, int]:
        return {header: idx for idx, header in enumerate(self.headers)}

    def column(self, name: str) -> np.ndarray:
        return self.data[:, self.column_index[name]]


def parse_scaps_batch_iv(path: Path) -> List[SweepBlock]:
    """Parse one SCAPS batch IV export into one block per temperature."""
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    blocks: List[SweepBlock] = []
    current_temp: float | None = None
    headers: List[str] | None = None
    rows: List[List[float]] = []
    reading_table = False

    def flush_block() -> None:
        if current_temp is None or headers is None or not rows:
            return
        blocks.append(
            SweepBlock(
                temperature_k=current_temp,
                headers=headers.copy(),
                data=np.asarray(rows, dtype=float),
            )
        )

    for line in lines:
        temp_match = re.search(r"temperature \[K\]:\s*([0-9.eE+-]+)", line)
        if temp_match:
            flush_block()
            current_temp = float(temp_match.group(1))
            headers = None
            rows = []
            reading_table = False
            continue

        stripped = line.strip()
        if current_temp is not None and stripped.startswith("v(V)"):
            headers = re.split(r"\s+", stripped)
            reading_table = True
            continue

        if not reading_table or headers is None:
            continue

        parts = re.split(r"\s+", stripped)
        if len(parts) < len(headers):
            continue
        try:
            row = [float(value) for value in parts[: len(headers)]]
        except ValueError:
            continue
        rows.append(row)

    flush_block()

    if not blocks:
        raise ValueError(f"No batch JV blocks found in {path}")

    return blocks


def temperature_to_celsius(temperature_k: Sequence[float] | np.ndarray) -> np.ndarray:
    return np.asarray(temperature_k, dtype=float) - 273.15


def interpolate_series(x: Sequence[float], y: Sequence[float], x_target: float) -> float:
    x_arr = np.asarray(x, dtype=float)
    y_arr = np.asarray(y, dtype=float)
    if x_target < np.min(x_arr) or x_target > np.max(x_arr):
        raise ValueError(
            f"Target {x_target:.2f} K is outside the sweep ({np.min(x_arr):.2f}-{np.max(x_arr):.2f} K)."
        )
    return float(np.interp(x_target, x_arr, y_arr))


def linear_slope_per_k(x: Sequence[float], y: Sequence[float]) -> float:
    x_arr = np.asarray(x, dtype=float)
    y_arr = np.asarray(y, dtype=float)
    if x_arr.size < 2:
        return float("nan")
    return float(np.polyfit(x_arr, y_arr, 1)[0])


def format_value(value: float, digits: int = 6) -> str:
    if np.isnan(value):
        return ""
    return f"{value:.{digits}f}"


def sanitize_column_name(column_name: str) -> str:
    name = column_name.replace("(mA/cm2)", "_mAcm2")
    name = name.replace("(V)", "_V")
    name = name.replace("[", "_").replace("]", "")
    name = name.replace("(", "_").replace(")", "")
    name = name.replace("/", "")
    name = name.replace("%", "percent")
    name = name.replace(" ", "_")
    name = name.replace("-", "_")
    name = re.sub(r"_+", "_", name)
    return name.strip("_").lower()


def parse_scaps_layers(path: Path) -> List[Dict[str, float | str]]:
    """Read layer name, thickness and nominal Eg from a SCAPS settings file."""
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    layers: List[Dict[str, float | str]] = []
    current_name: str | None = None
    current_thickness_m: float | None = None

    for line in lines:
        stripped = line.strip()
        if stripped == "layer":
            current_name = None
            current_thickness_m = None
            continue
        if stripped.startswith("name :"):
            current_name = stripped.split(":", 1)[1].strip()
            continue
        if stripped.startswith("d :"):
            match = re.search(r"d :\s*([0-9.eE+-]+)", stripped)
            if match:
                current_thickness_m = float(match.group(1))
            continue
        if stripped.startswith("Eg :") and current_name is not None:
            match = re.search(r"Eg :\s*([0-9.eE+-]+)", stripped)
            if match:
                layers.append(
                    {
                        "name": current_name,
                        "thickness_m": current_thickness_m if current_thickness_m else float("nan"),
                        "eg_ev": float(match.group(1)),
                    }
                )
                current_name = None
                current_thickness_m = None

    return layers


def select_absorber_layer(path: Path) -> Dict[str, float | str] | None:
    """Pick the thickest layer with an absorber-scale bandgap."""
    candidates = []
    for layer in parse_scaps_layers(path):
        eg_ev = float(layer["eg_ev"])
        thickness_m = float(layer["thickness_m"])
        if 1.0 <= eg_ev <= 2.0 and np.isfinite(thickness_m):
            candidates.append(layer)
    if not candidates:
        return None
    return max(candidates, key=lambda layer: float(layer["thickness_m"]))


def analyze_single_cell_block(block: SweepBlock, cell_name: str) -> Dict[str, float | str]:
    voltage = block.column("v(V)")
    current = block.column("jtot(mA/cm2)")
    metrics = extract_pv_metrics(voltage, current)
    voc = metrics["voc"]

    row: Dict[str, float | str] = {
        "cell": cell_name,
        "temperature_k": block.temperature_k,
        "temperature_c": block.temperature_k - 273.15,
        **metrics,
    }

    for column_name in RECOMB_COLUMNS:
        if column_name not in block.column_index:
            continue
        series = block.column(column_name)
        key = sanitize_column_name(column_name)
        row[f"{key}_at_voc"] = interpolate_y_at_x(voltage, series, voc)
        row[f"{key}_at_0v"] = interpolate_y_at_x(voltage, series, 0.0)

    total_rec_voc = float(row.get("j_total_rec_macm2_at_voc", float("nan")))
    bulk_rec_voc = float(row.get("jbulk_macm2_at_voc", float("nan")))
    left_minor_voc = float(row.get("jminor_left_macm2_at_voc", float("nan")))
    right_minor_voc = float(row.get("jminor_right_macm2_at_voc", float("nan")))
    interface_voc = float(row.get("jifr_macm2_at_voc", 0.0))

    total_rec_sc = float(row.get("j_total_rec_macm2_at_0v", float("nan")))
    total_gen_sc = float(row.get("j_total_gen_macm2_at_0v", float("nan")))

    if np.isfinite(total_rec_voc) and total_rec_voc > 0:
        nonbulk_voc = max(total_rec_voc - bulk_rec_voc, 0.0)
        row["bulk_fraction_at_voc"] = bulk_rec_voc / total_rec_voc
        row["nonbulk_fraction_at_voc"] = nonbulk_voc / total_rec_voc
        row["left_minor_fraction_at_voc"] = left_minor_voc / total_rec_voc
        row["right_minor_fraction_at_voc"] = right_minor_voc / total_rec_voc
        row["interface_fraction_at_voc"] = interface_voc / total_rec_voc
    else:
        row["bulk_fraction_at_voc"] = float("nan")
        row["nonbulk_fraction_at_voc"] = float("nan")
        row["left_minor_fraction_at_voc"] = float("nan")
        row["right_minor_fraction_at_voc"] = float("nan")
        row["interface_fraction_at_voc"] = float("nan")

    if np.isfinite(total_rec_sc) and np.isfinite(total_gen_sc) and total_gen_sc > 0:
        row["collection_loss_fraction_at_0v"] = total_rec_sc / total_gen_sc
    else:
        row["collection_loss_fraction_at_0v"] = float("nan")

    return row


def align_temperature_blocks(
    top_blocks: Sequence[SweepBlock], bottom_blocks: Sequence[SweepBlock]
) -> None:
    if len(top_blocks) != len(bottom_blocks):
        raise ValueError("Top and bottom sweeps do not have the same number of temperature steps.")
    for top_block, bottom_block in zip(top_blocks, bottom_blocks):
        if abs(top_block.temperature_k - bottom_block.temperature_k) > 1e-6:
            raise ValueError("Top and bottom sweeps are not aligned at the same temperatures.")


def analyze_tandems(
    top_blocks: Sequence[SweepBlock],
    bottom_blocks: Sequence[SweepBlock],
    pin_mw_cm2: float,
) -> tuple[List[Dict[str, float | str]], List[Dict[str, float | str]]]:
    tandem_rows: List[Dict[str, float | str]] = []
    drift_rows: List[Dict[str, float | str]] = []

    for top_block, bottom_block in zip(top_blocks, bottom_blocks):
        temperature_k = top_block.temperature_k
        top_metrics = extract_pv_metrics(
            top_block.column("v(V)"), top_block.column("jtot(mA/cm2)"), pin_mw_cm2
        )
        bottom_metrics = extract_pv_metrics(
            bottom_block.column("v(V)"), bottom_block.column("jtot(mA/cm2)"), pin_mw_cm2
        )
        tandem_2t = construct_tandem_2t(
            v_top=top_block.column("v(V)"),
            j_top=top_block.column("jtot(mA/cm2)"),
            v_bottom=bottom_block.column("v(V)"),
            j_bottom=bottom_block.column("jtot(mA/cm2)"),
            pin_mw_cm2=pin_mw_cm2,
        )
        tandem_4t_metrics = construct_tandem_4t_equivalent(top_metrics, bottom_metrics, pin_mw_cm2)

        tandem_rows.append(
            {
                "cell": "Tandem_2T",
                "temperature_k": temperature_k,
                "temperature_c": temperature_k - 273.15,
                **tandem_2t["metrics"],
                **tandem_2t["status"],
            }
        )
        tandem_rows.append(
            {
                "cell": "Tandem_4T",
                "temperature_k": temperature_k,
                "temperature_c": temperature_k - 273.15,
                **tandem_4t_metrics,
                "jsc_source": "equivalent_power_addition",
                "voc_source": "equivalent_power_addition",
            }
        )

        drift_rows.append(
            {
                "temperature_k": temperature_k,
                "temperature_c": temperature_k - 273.15,
                "jtop_jsc": top_metrics["jsc"],
                "jbottom_jsc": bottom_metrics["jsc"],
                "j_gap_top_minus_bottom": top_metrics["jsc"] - bottom_metrics["jsc"],
                "j_ratio_top_over_bottom": (
                    top_metrics["jsc"] / bottom_metrics["jsc"]
                    if bottom_metrics["jsc"] > 0
                    else float("nan")
                ),
            }
        )

    return tandem_rows, drift_rows


def save_csv(path: Path, rows: Sequence[Dict[str, float | str]], fieldnames: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            serialized_row = {}
            for key in fieldnames:
                value = row.get(key, "")
                if isinstance(value, float):
                    serialized_row[key] = format_value(value)
                else:
                    serialized_row[key] = value
            writer.writerow(serialized_row)


def extract_series(rows: Sequence[Dict[str, float | str]], cell: str, field: str) -> tuple[np.ndarray, np.ndarray]:
    filtered_rows = [row for row in rows if row["cell"] == cell]
    temperatures_k = np.asarray([float(row["temperature_k"]) for row in filtered_rows], dtype=float)
    values = np.asarray([float(row[field]) for row in filtered_rows], dtype=float)
    return temperatures_k, values


def build_metric_summary_rows(
    rows: Sequence[Dict[str, float | str]],
    baseline_k: float,
    target_k: float,
) -> List[Dict[str, float | str]]:
    summary_rows: List[Dict[str, float | str]] = []

    for cell in sorted({str(row["cell"]) for row in rows}):
        temperatures_k, _ = extract_series(rows, cell, "voc")
        for metric_name, metric_meta in SUMMARY_METRIC_META.items():
            _, values = extract_series(rows, cell, metric_name)
            baseline_value = interpolate_series(temperatures_k, values, baseline_k)
            target_value = interpolate_series(temperatures_k, values, target_k)
            slope = linear_slope_per_k(temperatures_k, values)
            summary_rows.append(
                {
                    "cell": cell,
                    "metric": metric_name,
                    "unit": metric_meta["unit"],
                    "baseline_temperature_k": baseline_k,
                    "baseline_temperature_c": baseline_k - 273.15,
                    "baseline_value": baseline_value,
                    "target_temperature_k": target_k,
                    "target_temperature_c": target_k - 273.15,
                    "target_value": target_value,
                    "linear_slope_per_k": slope,
                    "normalized_coeff_percent_per_k": (
                        100.0 * slope / baseline_value if baseline_value != 0 else float("nan")
                    ),
                }
            )

    return summary_rows


def build_temperature_coefficient_rows(
    rows: Sequence[Dict[str, float | str]],
    baseline_k: float,
    target_k: float,
) -> List[Dict[str, float | str]]:
    coefficient_rows: List[Dict[str, float | str]] = []

    for cell in sorted({str(row["cell"]) for row in rows}):
        temperatures_k, _ = extract_series(rows, cell, "voc")
        for metric_name, coeff_meta in COEFFICIENT_META.items():
            _, values = extract_series(rows, cell, metric_name)
            baseline_value = interpolate_series(temperatures_k, values, baseline_k)
            target_value = interpolate_series(temperatures_k, values, target_k)
            slope = linear_slope_per_k(temperatures_k, values)
            coefficient_rows.append(
                {
                    "cell": cell,
                    "coefficient_symbol": coeff_meta["symbol"],
                    "metric": metric_name,
                    "metric_unit": SUMMARY_METRIC_META[metric_name]["unit"],
                    "baseline_temperature_k": baseline_k,
                    "baseline_temperature_c": baseline_k - 273.15,
                    "baseline_value": baseline_value,
                    "target_temperature_k": target_k,
                    "target_temperature_c": target_k - 273.15,
                    "target_value": target_value,
                    "linear_slope_per_k": slope,
                    "absolute_coefficient": coeff_meta["absolute_scale"] * slope,
                    "absolute_unit": coeff_meta["absolute_unit"],
                    "normalized_coeff_percent_per_k": (
                        100.0 * slope / baseline_value if baseline_value != 0 else float("nan")
                    ),
                }
            )

    return coefficient_rows


def build_mismatch_summary_rows(
    drift_rows: Sequence[Dict[str, float | str]],
    baseline_k: float,
    target_k: float,
) -> List[Dict[str, float | str]]:
    temperatures_k = np.asarray([float(row["temperature_k"]) for row in drift_rows], dtype=float)
    gap = np.asarray([float(row["j_gap_top_minus_bottom"]) for row in drift_rows], dtype=float)
    ratio = np.asarray([float(row["j_ratio_top_over_bottom"]) for row in drift_rows], dtype=float)

    gap_baseline = interpolate_series(temperatures_k, gap, baseline_k)
    gap_target = interpolate_series(temperatures_k, gap, target_k)
    ratio_baseline = interpolate_series(temperatures_k, ratio, baseline_k)
    ratio_target = interpolate_series(temperatures_k, ratio, target_k)

    return [
        {
            "metric": "Jtop_minus_Jbottom",
            "unit": "mA/cm2",
            "baseline_temperature_k": baseline_k,
            "baseline_temperature_c": baseline_k - 273.15,
            "baseline_value": gap_baseline,
            "target_temperature_k": target_k,
            "target_temperature_c": target_k - 273.15,
            "target_value": gap_target,
            "absolute_change": gap_target - gap_baseline,
            "relative_change_percent": (
                100.0 * (gap_target / gap_baseline - 1.0) if gap_baseline != 0 else float("nan")
            ),
            "trend": "widens" if gap_target > gap_baseline else "closes",
        },
        {
            "metric": "Jtop_over_Jbottom",
            "unit": "ratio",
            "baseline_temperature_k": baseline_k,
            "baseline_temperature_c": baseline_k - 273.15,
            "baseline_value": ratio_baseline,
            "target_temperature_k": target_k,
            "target_temperature_c": target_k - 273.15,
            "target_value": ratio_target,
            "absolute_change": ratio_target - ratio_baseline,
            "relative_change_percent": (
                100.0 * (ratio_target / ratio_baseline - 1.0) if ratio_baseline != 0 else float("nan")
            ),
            "trend": "widens" if ratio_target > ratio_baseline else "closes",
        },
    ]


def add_reference_vertical(ax: plt.Axes, target_c: float) -> None:
    ax.axvline(target_c, color="black", lw=1.0, ls="--", alpha=0.5)


def plot_single_cell_metrics(
    out_path: Path,
    rows: Sequence[Dict[str, float | str]],
    title: str,
    line_color: str,
    label: str,
    target_temp_c: float,
) -> None:
    temperatures_k = np.asarray([float(row["temperature_k"]) for row in rows], dtype=float)
    x_celsius = temperature_to_celsius(temperatures_k)

    fig, axes = plt.subplots(2, 2, figsize=(9.2, 6.4), sharex=True)
    for ax, metric_name in zip(axes.flat, METRIC_META.keys()):
        values = np.asarray([float(row[metric_name]) for row in rows], dtype=float)
        ax.plot(x_celsius, values, color=line_color, lw=2.1, label=label)
        ax.scatter(
            [target_temp_c],
            [np.interp(target_temp_c, x_celsius, values)],
            s=34,
            color=line_color,
            edgecolor="white",
            linewidth=0.7,
            zorder=3,
        )
        add_reference_vertical(ax, target_temp_c)
        ax.set_ylabel(METRIC_META[metric_name]["label"])
        style_axis(ax)

    axes[1, 0].set_xlabel("Temperature (°C)")
    axes[1, 1].set_xlabel("Temperature (°C)")
    axes[0, 0].legend(loc="best")
    fig.suptitle(title, y=0.99)
    fig.tight_layout()
    fig.savefig(out_path, dpi=600)
    plt.close(fig)


def plot_subcell_metrics(
    out_path: Path,
    top_rows: Sequence[Dict[str, float | str]],
    bottom_rows: Sequence[Dict[str, float | str]],
    target_temp_c: float,
) -> None:
    temperatures_k = np.asarray([float(row["temperature_k"]) for row in top_rows], dtype=float)
    x_celsius = temperature_to_celsius(temperatures_k)
    baseline_k = float(np.min(temperatures_k))

    fig, axes = plt.subplots(2, 2, figsize=(9.2, 6.4), sharex=True)
    for ax, metric_name in zip(axes.flat, METRIC_META.keys()):
        top_values = np.asarray([float(row[metric_name]) for row in top_rows], dtype=float)
        bottom_values = np.asarray([float(row[metric_name]) for row in bottom_rows], dtype=float)
        top_baseline = interpolate_series(temperatures_k, top_values, baseline_k)
        bottom_baseline = interpolate_series(temperatures_k, bottom_values, baseline_k)
        top_rel = 100.0 * top_values / top_baseline
        bottom_rel = 100.0 * bottom_values / bottom_baseline

        ax.plot(x_celsius, top_rel, color="#1f4e79", lw=2.0, label="PVSK top")
        ax.plot(x_celsius, bottom_rel, color="#7a3e00", lw=2.0, label="CdTe filtered bottom")
        ax.scatter(
            [target_temp_c],
            [np.interp(target_temp_c, x_celsius, top_rel)],
            s=30,
            color="#1f4e79",
            edgecolor="white",
            linewidth=0.7,
            zorder=3,
        )
        ax.scatter(
            [target_temp_c],
            [np.interp(target_temp_c, x_celsius, bottom_rel)],
            s=30,
            color="#7a3e00",
            edgecolor="white",
            linewidth=0.7,
            zorder=3,
        )
        add_reference_vertical(ax, target_temp_c)
        ax.set_ylabel(f"{METRIC_META[metric_name]['label']} rel. to 300 K (%)")
        style_axis(ax)

    axes[1, 0].set_xlabel("Temperature (°C)")
    axes[1, 1].set_xlabel("Temperature (°C)")
    axes[0, 0].legend(loc="best")
    fig.suptitle("Sub-cell Normalized Thermal Drift", y=0.99)
    fig.tight_layout()
    fig.savefig(out_path, dpi=600)
    plt.close(fig)


def plot_subcell_metric_separation(
    out_path: Path,
    top_rows: Sequence[Dict[str, float | str]],
    bottom_rows: Sequence[Dict[str, float | str]],
    target_temp_c: float,
) -> None:
    temperatures_k = np.asarray([float(row["temperature_k"]) for row in top_rows], dtype=float)
    x_celsius = temperature_to_celsius(temperatures_k)
    separation_meta = {
        "voc": {"label": "Voc top - bottom (V)"},
        "jsc": {"label": "Jsc top - bottom (mA cm$^{-2}$)"},
        "ff": {"label": "FF top - bottom (%)"},
        "eta": {"label": "Efficiency top - bottom (%)"},
    }

    fig, axes = plt.subplots(2, 2, figsize=(9.2, 6.4), sharex=True)
    for ax, metric_name in zip(axes.flat, separation_meta.keys()):
        top_values = np.asarray([float(row[metric_name]) for row in top_rows], dtype=float)
        bottom_values = np.asarray([float(row[metric_name]) for row in bottom_rows], dtype=float)
        gap = top_values - bottom_values
        ax.plot(x_celsius, gap, color="#0f7f54", lw=2.1, label="Top - bottom")
        ax.scatter(
            [target_temp_c],
            [np.interp(target_temp_c, x_celsius, gap)],
            s=34,
            color="#0f7f54",
            edgecolor="white",
            linewidth=0.7,
            zorder=3,
        )
        add_reference_vertical(ax, target_temp_c)
        ax.set_ylabel(separation_meta[metric_name]["label"])
        style_axis(ax)

    axes[1, 0].set_xlabel("Temperature (°C)")
    axes[1, 1].set_xlabel("Temperature (°C)")
    axes[0, 0].legend(loc="best")
    fig.suptitle("Sub-cell Metric Separation", y=0.99)
    fig.tight_layout()
    fig.savefig(out_path, dpi=600)
    plt.close(fig)


def plot_tandem_metrics(
    out_path: Path,
    tandem_rows: Sequence[Dict[str, float | str]],
    target_temp_c: float,
) -> None:
    temperatures_k, _ = extract_series(tandem_rows, "Tandem_2T", "voc")
    x_celsius = temperature_to_celsius(temperatures_k)
    baseline_k = float(np.min(temperatures_k))
    series_map = {
        "Tandem 2T": {"cell": "Tandem_2T", "color": "#0f7f54", "ls": "-"},
        "Tandem 4T": {"cell": "Tandem_4T", "color": "#7d1f6f", "ls": "--"},
    }

    fig, axes = plt.subplots(2, 2, figsize=(9.2, 6.4), sharex=True)
    for ax, metric_name in zip(axes.flat, METRIC_META.keys()):
        for label, meta in series_map.items():
            _, values = extract_series(tandem_rows, meta["cell"], metric_name)
            baseline_value = interpolate_series(temperatures_k, values, baseline_k)
            rel_values = 100.0 * values / baseline_value
            ax.plot(x_celsius, rel_values, color=meta["color"], lw=2.0, ls=meta["ls"], label=label)
            ax.scatter(
                [target_temp_c],
                [np.interp(target_temp_c, x_celsius, rel_values)],
                s=30,
                color=meta["color"],
                edgecolor="white",
                linewidth=0.7,
                zorder=3,
            )
        add_reference_vertical(ax, target_temp_c)
        ax.set_ylabel(f"{METRIC_META[metric_name]['label']} rel. to 300 K (%)")
        style_axis(ax)

    axes[1, 0].set_xlabel("Temperature (°C)")
    axes[1, 1].set_xlabel("Temperature (°C)")
    axes[0, 0].legend(loc="best")
    fig.suptitle("Tandem Normalized Thermal Drift", y=0.99)
    fig.tight_layout()
    fig.savefig(out_path, dpi=600)
    plt.close(fig)


def plot_voc_physics(
    out_path: Path,
    top_rows: Sequence[Dict[str, float | str]],
    bottom_rows: Sequence[Dict[str, float | str]],
    target_temp_c: float,
) -> None:
    temperatures_k = np.asarray([float(row["temperature_k"]) for row in top_rows], dtype=float)
    x_celsius = temperature_to_celsius(temperatures_k)
    top_voc = np.asarray([float(row["voc"]) for row in top_rows], dtype=float)
    bottom_voc = np.asarray([float(row["voc"]) for row in bottom_rows], dtype=float)
    top_nonbulk = 100.0 * np.asarray([float(row["nonbulk_fraction_at_voc"]) for row in top_rows], dtype=float)
    bottom_nonbulk = (
        100.0 * np.asarray([float(row["nonbulk_fraction_at_voc"]) for row in bottom_rows], dtype=float)
    )

    top_fit = np.poly1d(np.polyfit(x_celsius, top_voc, 1))(x_celsius)
    bottom_fit = np.poly1d(np.polyfit(x_celsius, bottom_voc, 1))(x_celsius)
    top_slope_mv = 1000.0 * linear_slope_per_k(temperatures_k, top_voc)
    bottom_slope_mv = 1000.0 * linear_slope_per_k(temperatures_k, bottom_voc)

    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.2))

    ax_voc = axes[0]
    ax_voc.plot(x_celsius, top_voc, color="#1f4e79", lw=2.0, label="PVSK top")
    ax_voc.plot(x_celsius, bottom_voc, color="#7a3e00", lw=2.0, label="CdTe filtered bottom")
    ax_voc.plot(x_celsius, top_fit, color="#1f4e79", lw=1.3, ls="--", alpha=0.75)
    ax_voc.plot(x_celsius, bottom_fit, color="#7a3e00", lw=1.3, ls="--", alpha=0.75)
    add_reference_vertical(ax_voc, target_temp_c)
    ax_voc.set_xlabel("Temperature (°C)")
    ax_voc.set_ylabel("Voc (V)")
    ax_voc.set_title("Voc Thermal Slopes")
    ax_voc.text(
        0.03,
        0.05,
        f"PVSK: {top_slope_mv:.2f} mV/K\nCdTe: {bottom_slope_mv:.2f} mV/K",
        transform=ax_voc.transAxes,
        va="bottom",
        ha="left",
        bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "0.7"},
    )
    style_axis(ax_voc)
    ax_voc.legend(loc="best")

    ax_loss = axes[1]
    ax_loss.plot(x_celsius, top_nonbulk, color="#1f4e79", lw=2.0, label="PVSK top")
    ax_loss.plot(x_celsius, bottom_nonbulk, color="#7a3e00", lw=2.0, label="CdTe filtered bottom")
    add_reference_vertical(ax_loss, target_temp_c)
    ax_loss.set_xlabel("Temperature (°C)")
    ax_loss.set_ylabel("Non-bulk share of total recombination at Voc (%)")
    ax_loss.set_title("Voc Loss Channel Severity")
    style_axis(ax_loss)
    ax_loss.legend(loc="best")

    fig.tight_layout()
    fig.savefig(out_path, dpi=600)
    plt.close(fig)


def plot_mismatch_drift(
    out_path: Path,
    drift_rows: Sequence[Dict[str, float | str]],
    top_rows: Sequence[Dict[str, float | str]],
    bottom_rows: Sequence[Dict[str, float | str]],
    target_temp_c: float,
) -> None:
    temperatures_k = np.asarray([float(row["temperature_k"]) for row in drift_rows], dtype=float)
    x_celsius = temperature_to_celsius(temperatures_k)
    jtop = np.asarray([float(row["jtop_jsc"]) for row in drift_rows], dtype=float)
    jbottom = np.asarray([float(row["jbottom_jsc"]) for row in drift_rows], dtype=float)
    gap = np.asarray([float(row["j_gap_top_minus_bottom"]) for row in drift_rows], dtype=float)
    ratio = np.asarray([float(row["j_ratio_top_over_bottom"]) for row in drift_rows], dtype=float)
    top_jrec_sc = np.asarray([float(row["j_total_rec_macm2_at_0v"]) for row in top_rows], dtype=float)
    bottom_jrec_sc = np.asarray([float(row["j_total_rec_macm2_at_0v"]) for row in bottom_rows], dtype=float)

    fig, axes = plt.subplots(1, 3, figsize=(14.0, 4.2))

    ax_curr = axes[0]
    ax_curr.plot(x_celsius, jtop, color="#1f4e79", lw=2.0, label="Jtop")
    ax_curr.plot(x_celsius, jbottom, color="#7a3e00", lw=2.0, label="Jbottom")
    ax_curr.fill_between(x_celsius, jbottom, jtop, color="#c4c4c4", alpha=0.25)
    add_reference_vertical(ax_curr, target_temp_c)
    ax_curr.set_xlabel("Temperature (°C)")
    ax_curr.set_ylabel("Sub-cell Jsc (mA cm$^{-2}$)")
    ax_curr.set_title("Current Matching Drift")
    style_axis(ax_curr)
    ax_curr.legend(loc="best")

    ax_gap = axes[1]
    ax_gap.plot(x_celsius, gap, color="#0f7f54", lw=2.0, label="Jtop - Jbottom")
    ax_gap.scatter(
        [target_temp_c],
        [np.interp(target_temp_c, x_celsius, gap)],
        color="#0f7f54",
        s=34,
        edgecolor="white",
        linewidth=0.7,
        zorder=3,
    )
    add_reference_vertical(ax_gap, target_temp_c)
    ax_gap.set_xlabel("Temperature (°C)")
    ax_gap.set_ylabel("Current gap (mA cm$^{-2}$)")
    ax_gap.set_title("Mismatch Gap and Ratio")
    style_axis(ax_gap)
    ax_ratio = ax_gap.twinx()
    ax_ratio.plot(x_celsius, ratio, color="#7d1f6f", lw=1.7, ls="--", label="Jtop / Jbottom")
    ax_ratio.set_ylabel("Mismatch ratio")
    ax_gap.legend(loc="upper left")
    ax_ratio.legend(loc="lower right")

    ax_rec = axes[2]
    ax_rec.plot(x_celsius, top_jrec_sc, color="#1f4e79", lw=2.0, label="PVSK top")
    ax_rec.plot(x_celsius, bottom_jrec_sc, color="#7a3e00", lw=2.0, label="CdTe filtered bottom")
    add_reference_vertical(ax_rec, target_temp_c)
    ax_rec.set_xlabel("Temperature (°C)")
    ax_rec.set_ylabel("Recombination current at 0 V (mA cm$^{-2}$)")
    ax_rec.set_title("Short-circuit Recombination Loss")
    style_axis(ax_rec)
    ax_rec.legend(loc="best")

    fig.tight_layout()
    fig.savefig(out_path, dpi=600)
    plt.close(fig)


def write_summary(
    out_path: Path,
    top_rows: Sequence[Dict[str, float | str]],
    bottom_rows: Sequence[Dict[str, float | str]],
    tandem_rows: Sequence[Dict[str, float | str]],
    mismatch_summary_rows: Sequence[Dict[str, float | str]],
    baseline_k: float,
    target_k: float,
    pvsk_absorber: Dict[str, float | str] | None,
    cdte_absorber: Dict[str, float | str] | None,
) -> None:
    target_c = target_k - 273.15
    top_temperatures_k = np.asarray([float(row["temperature_k"]) for row in top_rows], dtype=float)
    bottom_temperatures_k = np.asarray([float(row["temperature_k"]) for row in bottom_rows], dtype=float)

    top_voc = np.asarray([float(row["voc"]) for row in top_rows], dtype=float)
    bottom_voc = np.asarray([float(row["voc"]) for row in bottom_rows], dtype=float)
    top_nonbulk = np.asarray([float(row["nonbulk_fraction_at_voc"]) for row in top_rows], dtype=float)
    bottom_nonbulk = np.asarray([float(row["nonbulk_fraction_at_voc"]) for row in bottom_rows], dtype=float)
    top_jrec_sc = np.asarray([float(row["j_total_rec_macm2_at_0v"]) for row in top_rows], dtype=float)
    bottom_jrec_sc = np.asarray([float(row["j_total_rec_macm2_at_0v"]) for row in bottom_rows], dtype=float)

    tandem_2t_temps, tandem_2t_voc = extract_series(tandem_rows, "Tandem_2T", "voc")
    _, tandem_2t_eta = extract_series(tandem_rows, "Tandem_2T", "eta")
    _, tandem_4t_eta = extract_series(tandem_rows, "Tandem_4T", "eta")
    _, tandem_2t_pmax = extract_series(tandem_rows, "Tandem_2T", "pmax")
    _, tandem_4t_pmax = extract_series(tandem_rows, "Tandem_4T", "pmax")
    top_pmax = np.asarray([float(row["pmax"]) for row in top_rows], dtype=float)
    bottom_pmax = np.asarray([float(row["pmax"]) for row in bottom_rows], dtype=float)

    gap_row = next(row for row in mismatch_summary_rows if row["metric"] == "Jtop_minus_Jbottom")
    ratio_row = next(row for row in mismatch_summary_rows if row["metric"] == "Jtop_over_Jbottom")

    lines = [
        "Thermal PVSK/CdTe analysis summary",
        "=================================",
        "",
        f"Available simulation range: {np.min(top_temperatures_k):.2f}-{np.max(top_temperatures_k):.2f} K",
        f"Baseline reported here: {baseline_k:.2f} K ({baseline_k - 273.15:.2f} C)",
        f"Target reported here: {target_k:.2f} K ({target_c:.2f} C)",
        "The sweep starts at 300 K, so exact STC 25.00 C = 298.15 K is not present in the dataset.",
        "",
        "Sub-cell thermal coefficients:",
        (
            f"- PVSK top Voc slope = {1000.0 * linear_slope_per_k(top_temperatures_k, top_voc):.3f} mV/K "
            f"({100.0 * linear_slope_per_k(top_temperatures_k, top_voc) / interpolate_series(top_temperatures_k, top_voc, baseline_k):.3f} %/K)"
        ),
        (
            f"- CdTe filtered Voc slope = {1000.0 * linear_slope_per_k(bottom_temperatures_k, bottom_voc):.3f} mV/K "
            f"({100.0 * linear_slope_per_k(bottom_temperatures_k, bottom_voc) / interpolate_series(bottom_temperatures_k, bottom_voc, baseline_k):.3f} %/K)"
        ),
        f"- 2T tandem Voc slope = {1000.0 * linear_slope_per_k(tandem_2t_temps, tandem_2t_voc):.3f} mV/K",
        "",
        "Power thermal coefficients (gamma for Pmax):",
        (
            f"- PVSK top gamma_Pmax = {linear_slope_per_k(top_temperatures_k, top_pmax):.5f} mW/cm2/K "
            f"({100.0 * linear_slope_per_k(top_temperatures_k, top_pmax) / interpolate_series(top_temperatures_k, top_pmax, baseline_k):.3f} %/K)"
        ),
        (
            f"- CdTe filtered gamma_Pmax = {linear_slope_per_k(bottom_temperatures_k, bottom_pmax):.5f} mW/cm2/K "
            f"({100.0 * linear_slope_per_k(bottom_temperatures_k, bottom_pmax) / interpolate_series(bottom_temperatures_k, bottom_pmax, baseline_k):.3f} %/K)"
        ),
        (
            f"- 2T tandem gamma_Pmax = {linear_slope_per_k(tandem_2t_temps, tandem_2t_pmax):.5f} mW/cm2/K "
            f"({100.0 * linear_slope_per_k(tandem_2t_temps, tandem_2t_pmax) / interpolate_series(tandem_2t_temps, tandem_2t_pmax, baseline_k):.3f} %/K)"
        ),
        (
            f"- 4T equivalent gamma_Pmax = {linear_slope_per_k(tandem_2t_temps, tandem_4t_pmax):.5f} mW/cm2/K "
            f"({100.0 * linear_slope_per_k(tandem_2t_temps, tandem_4t_pmax) / interpolate_series(tandem_2t_temps, tandem_4t_pmax, baseline_k):.3f} %/K)"
        ),
        "- Because Pin is fixed at 100 mW/cm2 in this workflow, eta and Pmax are numerically identical, but gamma_Pmax is now reported explicitly to avoid ambiguity.",
        "",
        "Why Voc drops faster in the CdTe bottom cell:",
        (
            f"- At {baseline_k:.0f} K, Voc is {interpolate_series(top_temperatures_k, top_voc, baseline_k):.4f} V "
            f"for PVSK but only {interpolate_series(bottom_temperatures_k, bottom_voc, baseline_k):.4f} V for filtered CdTe."
        ),
        (
            f"- The non-bulk share of total recombination at Voc is "
            f"{100.0 * interpolate_series(top_temperatures_k, top_nonbulk, baseline_k):.2f}% for PVSK at {baseline_k:.0f} K "
            f"and {100.0 * interpolate_series(top_temperatures_k, top_nonbulk, target_k):.2f}% at {target_c:.2f} C."
        ),
        (
            f"- The same quantity is {100.0 * interpolate_series(bottom_temperatures_k, bottom_nonbulk, baseline_k):.2f}% "
            f"for CdTe at {baseline_k:.0f} K and {100.0 * interpolate_series(bottom_temperatures_k, bottom_nonbulk, target_k):.2f}% "
            f"at {target_c:.2f} C."
        ),
        "- Both cells are SRH-dominated in SCAPS, but the CdTe bottom cell carries a much larger thermally growing minority-contact leakage term.",
    ]

    if pvsk_absorber is not None:
        lines.append(
            (
                f"- The guessed PVSK absorber from the local device file is "
                f"'{pvsk_absorber['name']}' with nominal Eg ~ {float(pvsk_absorber['eg_ev']):.3f} eV."
            )
        )
    if cdte_absorber is not None:
        lines.append(
            (
                f"- The guessed CdTe absorber from the local device file is "
                f"'{cdte_absorber['name']}' with nominal Eg ~ {float(cdte_absorber['eg_ev']):.3f} eV."
            )
        )
    if pvsk_absorber is not None and cdte_absorber is not None:
        lines.append(
            "- The lower-gap CdTe absorber also has a smaller Voc headroom, so the same thermally activated dark-current rise costs it more voltage."
        )

    lines.extend(
        [
            "",
            "Current mismatch drift:",
            (
                f"- Jtop - Jbottom = {float(gap_row['baseline_value']):.4f} mA/cm2 at {baseline_k:.0f} K and "
                f"{float(gap_row['target_value']):.4f} mA/cm2 at {target_c:.2f} C. The gap {gap_row['trend']} "
                f"by {float(gap_row['absolute_change']):.4f} mA/cm2 "
                f"({float(gap_row['relative_change_percent']):.3f}%)."
            ),
            (
                f"- Jtop / Jbottom = {float(ratio_row['baseline_value']):.4f} at {baseline_k:.0f} K and "
                f"{float(ratio_row['target_value']):.4f} at {target_c:.2f} C. The ratio {ratio_row['trend']} "
                f"by {float(ratio_row['absolute_change']):.4f} "
                f"({float(ratio_row['relative_change_percent']):.3f}%)."
            ),
            (
                f"- Jrec at 0 V changes from {interpolate_series(top_temperatures_k, top_jrec_sc, baseline_k):.4f} to "
                f"{interpolate_series(top_temperatures_k, top_jrec_sc, target_k):.4f} mA/cm2 in PVSK, but from "
                f"{interpolate_series(bottom_temperatures_k, bottom_jrec_sc, baseline_k):.4f} to "
                f"{interpolate_series(bottom_temperatures_k, bottom_jrec_sc, target_k):.4f} mA/cm2 in CdTe."
            ),
            "- That is why the drift widens rather than closes: the top-cell current stays nearly flat, while the filtered bottom-cell current slowly erodes with temperature.",
            "",
            "Tandem performance at the target temperature:",
            (
                f"- 2T efficiency changes from {interpolate_series(tandem_2t_temps, tandem_2t_eta, baseline_k):.3f}% "
                f"to {interpolate_series(tandem_2t_temps, tandem_2t_eta, target_k):.3f}%."
            ),
            (
                f"- 4T equivalent efficiency changes from {interpolate_series(tandem_2t_temps, tandem_4t_eta, baseline_k):.3f}% "
                f"to {interpolate_series(tandem_2t_temps, tandem_4t_eta, target_k):.3f}%."
            ),
        ]
    )

    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Parse SCAPS temperature sweeps for PVSK top and filtered CdTe bottom cells, "
            "build 2T/4T tandems at each temperature, and export thermal figures."
        )
    )
    parser.add_argument("--pvsk-iv", default="data/scaps_exports/iv/PVSK_thermal_sweep.iv")
    parser.add_argument("--cdte-iv", default="data/scaps_exports/iv/CdTe_filtered_thermal_sweep.iv")
    parser.add_argument("--pvsk-device", default="devices/pvsk/scaps/Louis_PVSK.scaps")
    parser.add_argument("--cdte-device", default="devices/cdte/scaps/Group1_CdTe_with_PVSK_filter.scaps")
    parser.add_argument("--pin-mw-cm2", type=float, default=100.0)
    parser.add_argument("--baseline-temp-k", type=float, default=300.0)
    parser.add_argument("--target-temp-c", type=float, default=65.0)
    parser.add_argument("--outdir", default="results/thermal")
    args = parser.parse_args()

    target_temp_k = args.target_temp_c + 273.15
    apply_publication_style()

    outdir = resolve_project_path(args.outdir).resolve()
    figures_dir = outdir / "figures"
    csv_dir = outdir / "csv"
    figures_dir.mkdir(parents=True, exist_ok=True)
    csv_dir.mkdir(parents=True, exist_ok=True)

    top_blocks = parse_scaps_batch_iv(resolve_project_path(args.pvsk_iv))
    bottom_blocks = parse_scaps_batch_iv(resolve_project_path(args.cdte_iv))
    align_temperature_blocks(top_blocks, bottom_blocks)

    top_rows = [analyze_single_cell_block(block, "PVSK_top") for block in top_blocks]
    bottom_rows = [analyze_single_cell_block(block, "CdTe_filtered") for block in bottom_blocks]
    tandem_rows, drift_rows = analyze_tandems(top_blocks, bottom_blocks, args.pin_mw_cm2)

    summary_rows = build_metric_summary_rows(
        rows=top_rows + bottom_rows + tandem_rows,
        baseline_k=args.baseline_temp_k,
        target_k=target_temp_k,
    )
    coefficient_rows = build_temperature_coefficient_rows(
        rows=top_rows + bottom_rows + tandem_rows,
        baseline_k=args.baseline_temp_k,
        target_k=target_temp_k,
    )
    mismatch_summary_rows = build_mismatch_summary_rows(
        drift_rows=drift_rows,
        baseline_k=args.baseline_temp_k,
        target_k=target_temp_k,
    )

    top_fieldnames = [
        "cell",
        "temperature_k",
        "temperature_c",
        "voc",
        "jsc",
        "ff",
        "eta",
        "vmpp",
        "jmpp",
        "pmax",
        "j_total_rec_macm2_at_voc",
        "j_total_gen_macm2_at_voc",
        "jbulk_macm2_at_voc",
        "j_srh_macm2_at_voc",
        "jminor_left_macm2_at_voc",
        "jminor_right_macm2_at_voc",
        "bulk_fraction_at_voc",
        "nonbulk_fraction_at_voc",
        "left_minor_fraction_at_voc",
        "right_minor_fraction_at_voc",
        "j_total_rec_macm2_at_0v",
        "j_total_gen_macm2_at_0v",
        "collection_loss_fraction_at_0v",
    ]
    tandem_fieldnames = [
        "cell",
        "temperature_k",
        "temperature_c",
        "voc",
        "jsc",
        "ff",
        "eta",
        "vmpp",
        "jmpp",
        "pmax",
        "jsc_source",
        "voc_source",
        "v_tandem_min",
        "j_at_v_tandem_min",
        "v_top_at_v_tandem_min",
        "v_bottom_at_v_tandem_min",
    ]
    summary_fieldnames = [
        "cell",
        "metric",
        "unit",
        "baseline_temperature_k",
        "baseline_temperature_c",
        "baseline_value",
        "target_temperature_k",
        "target_temperature_c",
        "target_value",
        "linear_slope_per_k",
        "normalized_coeff_percent_per_k",
    ]
    coefficient_fieldnames = [
        "cell",
        "coefficient_symbol",
        "metric",
        "metric_unit",
        "baseline_temperature_k",
        "baseline_temperature_c",
        "baseline_value",
        "target_temperature_k",
        "target_temperature_c",
        "target_value",
        "linear_slope_per_k",
        "absolute_coefficient",
        "absolute_unit",
        "normalized_coeff_percent_per_k",
    ]
    mismatch_fieldnames = [
        "metric",
        "unit",
        "baseline_temperature_k",
        "baseline_temperature_c",
        "baseline_value",
        "target_temperature_k",
        "target_temperature_c",
        "target_value",
        "absolute_change",
        "relative_change_percent",
        "trend",
    ]
    drift_fieldnames = [
        "temperature_k",
        "temperature_c",
        "jtop_jsc",
        "jbottom_jsc",
        "j_gap_top_minus_bottom",
        "j_ratio_top_over_bottom",
    ]

    save_csv(csv_dir / "subcell_temperature_metrics.csv", top_rows + bottom_rows, top_fieldnames)
    save_csv(csv_dir / "tandem_temperature_metrics.csv", tandem_rows, tandem_fieldnames)
    save_csv(csv_dir / "temperature_metric_summary.csv", summary_rows, summary_fieldnames)
    save_csv(csv_dir / "temperature_coefficients.csv", coefficient_rows, coefficient_fieldnames)
    save_csv(csv_dir / "mismatch_drift_summary.csv", mismatch_summary_rows, mismatch_fieldnames)
    save_csv(csv_dir / "mismatch_drift_vs_temperature.csv", drift_rows, drift_fieldnames)

    plot_single_cell_metrics(
        figures_dir / "pvsk_metrics_vs_temperature.png",
        top_rows,
        "PVSK Top-Cell Metrics vs Temperature",
        "#1f4e79",
        "PVSK top",
        args.target_temp_c,
    )
    plot_single_cell_metrics(
        figures_dir / "cdte_filtered_metrics_vs_temperature.png",
        bottom_rows,
        "Filtered CdTe Bottom-Cell Metrics vs Temperature",
        "#7a3e00",
        "CdTe filtered bottom",
        args.target_temp_c,
    )
    plot_subcell_metrics(
        figures_dir / "subcell_metrics_vs_temperature.png",
        top_rows,
        bottom_rows,
        args.target_temp_c,
    )
    plot_subcell_metric_separation(
        figures_dir / "subcell_metric_separation_vs_temperature.png",
        top_rows,
        bottom_rows,
        args.target_temp_c,
    )
    tandem_2t_rows = [row for row in tandem_rows if row["cell"] == "Tandem_2T"]
    tandem_4t_rows = [row for row in tandem_rows if row["cell"] == "Tandem_4T"]

    plot_single_cell_metrics(
        figures_dir / "tandem_2t_metrics_vs_temperature.png",
        tandem_2t_rows,
        "2T Tandem Metrics vs Temperature",
        "#0f7f54",
        "Tandem 2T",
        args.target_temp_c,
    )
    plot_single_cell_metrics(
        figures_dir / "tandem_4t_metrics_vs_temperature.png",
        tandem_4t_rows,
        "4T Equivalent Tandem Metrics vs Temperature",
        "#7d1f6f",
        "Tandem 4T",
        args.target_temp_c,
    )
    plot_tandem_metrics(figures_dir / "tandem_metrics_vs_temperature.png", tandem_rows, args.target_temp_c)
    plot_voc_physics(figures_dir / "voc_physics_vs_temperature.png", top_rows, bottom_rows, args.target_temp_c)
    plot_mismatch_drift(
        figures_dir / "mismatch_drift_vs_temperature.png",
        drift_rows,
        top_rows,
        bottom_rows,
        args.target_temp_c,
    )

    pvsk_device_path = resolve_project_path(args.pvsk_device)
    cdte_device_path = resolve_project_path(args.cdte_device)
    pvsk_absorber = select_absorber_layer(pvsk_device_path) if pvsk_device_path.exists() else None
    cdte_absorber = select_absorber_layer(cdte_device_path) if cdte_device_path.exists() else None

    write_summary(
        outdir / "thermal_summary.txt",
        top_rows,
        bottom_rows,
        tandem_rows,
        mismatch_summary_rows,
        args.baseline_temp_k,
        target_temp_k,
        pvsk_absorber,
        cdte_absorber,
    )

    gap_row = next(row for row in mismatch_summary_rows if row["metric"] == "Jtop_minus_Jbottom")
    ratio_row = next(row for row in mismatch_summary_rows if row["metric"] == "Jtop_over_Jbottom")
    top_voc_slope = 1000.0 * linear_slope_per_k(
        [float(row["temperature_k"]) for row in top_rows],
        [float(row["voc"]) for row in top_rows],
    )
    bottom_voc_slope = 1000.0 * linear_slope_per_k(
        [float(row["temperature_k"]) for row in bottom_rows],
        [float(row["voc"]) for row in bottom_rows],
    )

    print("Thermal analysis complete.")
    print(f"Reference: {args.baseline_temp_k:.2f} K ({args.baseline_temp_k - 273.15:.2f} C)")
    print(f"Target:    {target_temp_k:.2f} K ({args.target_temp_c:.2f} C)")
    print(f"PVSK Voc slope: {top_voc_slope:.3f} mV/K")
    print(f"CdTe Voc slope: {bottom_voc_slope:.3f} mV/K")
    print(f"Jtop - Jbottom at target: {float(gap_row['target_value']):.4f} mA/cm2 ({str(gap_row['trend'])})")
    print(f"Jtop / Jbottom at target: {float(ratio_row['target_value']):.4f} ({str(ratio_row['trend'])})")
    print(f"Outputs written to: {outdir}")


if __name__ == "__main__":
    main()
