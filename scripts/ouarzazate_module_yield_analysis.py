#!/usr/bin/env python3
"""Estimate hourly tandem module performance from PVGIS weather data."""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List

PROJECT_ROOT = Path(__file__).resolve().parent.parent
STC_TEMPERATURE_C = 25.0
STC_TEMPERATURE_K = STC_TEMPERATURE_C + 273.15
STC_IRRADIANCE_W_M2 = 1000.0
THERMAL_VOLTAGE_STC_V = 8.617333262145e-5 * STC_TEMPERATURE_K


@dataclass(frozen=True)
class MetricPoint:
    cell: str
    metric: str
    unit: str
    baseline_temperature_k: float
    baseline_value: float
    linear_slope_per_k: float
    normalized_coeff_percent_per_k: float

    def value_at_temperature_k(self, temperature_k: float) -> float:
        return self.baseline_value + self.linear_slope_per_k * (
            temperature_k - self.baseline_temperature_k
        )


@dataclass(frozen=True)
class PvgisHour:
    timestamp: datetime
    irradiance_w_m2: float
    ambient_c: float
    wind_10m_ms: float
    reconstructed: bool


def resolve_project_path(path_str: str) -> Path:
    path = Path(path_str)
    if path.is_absolute():
        return path
    return PROJECT_ROOT / path


def read_temperature_metric_summary(path: Path) -> Dict[tuple[str, str], MetricPoint]:
    metrics: Dict[tuple[str, str], MetricPoint] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            key = (row["cell"], row["metric"])
            metrics[key] = MetricPoint(
                cell=row["cell"],
                metric=row["metric"],
                unit=row["unit"],
                baseline_temperature_k=float(row["baseline_temperature_k"]),
                baseline_value=float(row["baseline_value"]),
                linear_slope_per_k=float(row["linear_slope_per_k"]),
                normalized_coeff_percent_per_k=float(row["normalized_coeff_percent_per_k"]),
            )
    return metrics


def read_pvgis_hourly(path: Path) -> List[PvgisHour]:
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    try:
        header_index = next(index for index, line in enumerate(lines) if line.startswith("time,"))
    except StopIteration as exc:
        raise ValueError(f"Could not find PVGIS hourly header in {path}") from exc

    rows: List[PvgisHour] = []
    reader = csv.DictReader(lines[header_index:])
    for row in reader:
        time_token = row.get("time", "")
        if not time_token or len(time_token) < 13 or not time_token[:8].isdigit():
            break
        try:
            rows.append(
                PvgisHour(
                    timestamp=datetime.strptime(time_token, "%Y%m%d:%H%M"),
                    irradiance_w_m2=float(row["G(i)"]),
                    ambient_c=float(row["T2m"]),
                    wind_10m_ms=float(row["WS10m"]),
                    reconstructed=bool(float(row["Int"])),
                )
            )
        except (KeyError, ValueError) as exc:
            raise ValueError(f"Failed to parse PVGIS row: {row}") from exc

    if not rows:
        raise ValueError(f"No hourly PVGIS rows found in {path}")

    return rows


def ff_from_dimensionless_voc(voc_star: float) -> float:
    return (voc_star - math.log(voc_star + 0.72)) / (voc_star + 1.0)


def solve_dimensionless_voc(ff_percent: float) -> float:
    ff_target = ff_percent / 100.0
    if not 0.0 < ff_target < 1.0:
        raise ValueError(f"FF must be between 0 and 100%, got {ff_percent}")

    lower = 1.0
    upper = 200.0
    for _ in range(200):
        midpoint = 0.5 * (lower + upper)
        if ff_from_dimensionless_voc(midpoint) < ff_target:
            lower = midpoint
        else:
            upper = midpoint
    return 0.5 * (lower + upper)


def infer_n_ns(voc_stc_v: float, ff_stc_percent: float) -> float:
    voc_star = solve_dimensionless_voc(ff_stc_percent)
    return voc_stc_v / (voc_star * THERMAL_VOLTAGE_STC_V)


def thermal_voltage_v(temperature_c: float) -> float:
    temperature_k = temperature_c + 273.15
    return 8.617333262145e-5 * temperature_k


def wind_speed_at_module_height(
    wind_speed_10m_ms: float,
    module_height_m: float,
    wind_shear_exponent: float,
) -> float:
    if wind_speed_10m_ms <= 0.0:
        return 0.0
    return wind_speed_10m_ms * (module_height_m / 10.0) ** wind_shear_exponent


def faiman_module_temperature_c(
    irradiance_w_m2: float,
    ambient_c: float,
    wind_module_ms: float,
    u0: float,
    u1: float,
) -> float:
    if irradiance_w_m2 <= 0.0:
        return ambient_c
    return ambient_c + irradiance_w_m2 / (u0 + u1 * wind_module_ms)


def irradiance_adjusted_voc_v(
    voc_stc_v: float,
    beta_v_per_k: float,
    n_ns: float,
    module_temperature_c: float,
    irradiance_w_m2: float,
) -> float:
    if irradiance_w_m2 <= 0.0:
        return 0.0

    irradiance_ratio = irradiance_w_m2 / STC_IRRADIANCE_W_M2
    voc_v = voc_stc_v + beta_v_per_k * (module_temperature_c - STC_TEMPERATURE_C)
    voc_v += n_ns * thermal_voltage_v(module_temperature_c) * math.log(irradiance_ratio)
    return max(voc_v, 0.0)


def irradiance_adjusted_pmax_density_mw_cm2(
    pmax_stc_mw_cm2: float,
    gamma_per_k: float,
    module_temperature_c: float,
    irradiance_w_m2: float,
) -> float:
    if irradiance_w_m2 <= 0.0:
        return 0.0
    irradiance_ratio = irradiance_w_m2 / STC_IRRADIANCE_W_M2
    thermal_factor = 1.0 + gamma_per_k * (module_temperature_c - STC_TEMPERATURE_C)
    return max(pmax_stc_mw_cm2 * irradiance_ratio * thermal_factor, 0.0)


def density_to_module_power_w(pmax_density_mw_cm2: float, module_area_m2: float) -> float:
    return pmax_density_mw_cm2 * module_area_m2 * 10.0


def float_mean(values: Iterable[float]) -> float:
    values_list = list(values)
    if not values_list:
        return float("nan")
    return sum(values_list) / len(values_list)


def write_csv(path: Path, fieldnames: List[str], rows: Iterable[Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compute hourly 2T and 4T tandem module performance from PVGIS weather data."
    )
    parser.add_argument(
        "--pvgis-file",
        default="data/PVGIS/Timeseries_30.937_-6.906_SA3_32deg_-4deg_2023_2023.csv",
        help="PVGIS hourly CSV exported for the array plane of interest.",
    )
    parser.add_argument(
        "--thermal-summary",
        default="results/thermal/csv/temperature_metric_summary.csv",
        help="Thermal summary CSV with baseline values and extracted coefficients.",
    )
    parser.add_argument(
        "--output-dir",
        default="results/pvgis/csv",
        help="Directory for the hourly and annual summary CSV files.",
    )
    parser.add_argument("--cell-count", type=int, default=60, help="Series cell count per module.")
    parser.add_argument(
        "--module-area-m2",
        type=float,
        default=1.65,
        help="Total module area used to scale power density to module power.",
    )
    parser.add_argument(
        "--faiman-u0",
        type=float,
        default=25.0,
        help="Faiman U0 coefficient in W/m^2/K.",
    )
    parser.add_argument(
        "--faiman-u1",
        type=float,
        default=6.84,
        help="Faiman U1 coefficient in W s/m^3/K.",
    )
    parser.add_argument(
        "--module-height-m",
        type=float,
        default=2.0,
        help="Representative module height used to downscale PVGIS 10 m wind speed.",
    )
    parser.add_argument(
        "--wind-shear-exponent",
        type=float,
        default=0.14,
        help="Power-law exponent for adjusting 10 m wind speed to module height.",
    )
    args = parser.parse_args()

    pvgis_file = resolve_project_path(args.pvgis_file)
    thermal_summary = resolve_project_path(args.thermal_summary)
    output_dir = resolve_project_path(args.output_dir)

    metrics = read_temperature_metric_summary(thermal_summary)
    pvgis_rows = read_pvgis_hourly(pvgis_file)

    top_voc = metrics[("PVSK_top", "voc")]
    top_ff = metrics[("PVSK_top", "ff")]
    bottom_voc = metrics[("CdTe_filtered", "voc")]
    bottom_ff = metrics[("CdTe_filtered", "ff")]
    tandem_2t_pmax = metrics[("Tandem_2T", "pmax")]
    tandem_4t_pmax = metrics[("Tandem_4T", "pmax")]

    top_voc_stc_v = top_voc.value_at_temperature_k(STC_TEMPERATURE_K)
    bottom_voc_stc_v = bottom_voc.value_at_temperature_k(STC_TEMPERATURE_K)
    top_ff_stc_percent = top_ff.value_at_temperature_k(STC_TEMPERATURE_K)
    bottom_ff_stc_percent = bottom_ff.value_at_temperature_k(STC_TEMPERATURE_K)
    module_voc_stc_v = (top_voc_stc_v + bottom_voc_stc_v) * args.cell_count

    top_n_ns = infer_n_ns(top_voc_stc_v, top_ff_stc_percent)
    bottom_n_ns = infer_n_ns(bottom_voc_stc_v, bottom_ff_stc_percent)

    tandem_2t_pmax_stc_mw_cm2 = tandem_2t_pmax.value_at_temperature_k(STC_TEMPERATURE_K)
    tandem_4t_pmax_stc_mw_cm2 = tandem_4t_pmax.value_at_temperature_k(STC_TEMPERATURE_K)

    gamma_2t_per_k = tandem_2t_pmax.linear_slope_per_k / tandem_2t_pmax_stc_mw_cm2
    gamma_4t_per_k = tandem_4t_pmax.linear_slope_per_k / tandem_4t_pmax_stc_mw_cm2

    module_pmax_2t_stc_w = density_to_module_power_w(
        tandem_2t_pmax_stc_mw_cm2, args.module_area_m2
    )
    module_pmax_4t_stc_w = density_to_module_power_w(
        tandem_4t_pmax_stc_mw_cm2, args.module_area_m2
    )

    hourly_rows: List[Dict[str, object]] = []
    voc_module_series: List[float] = []
    voc_module_sunlit_series: List[float] = []
    pmax_2t_series_w: List[float] = []
    pmax_4t_series_w: List[float] = []

    for hour in pvgis_rows:
        wind_module_ms = wind_speed_at_module_height(
            hour.wind_10m_ms,
            args.module_height_m,
            args.wind_shear_exponent,
        )
        module_temperature_c = faiman_module_temperature_c(
            hour.irradiance_w_m2,
            hour.ambient_c,
            wind_module_ms,
            args.faiman_u0,
            args.faiman_u1,
        )

        voc_top_v = irradiance_adjusted_voc_v(
            top_voc_stc_v,
            top_voc.linear_slope_per_k,
            top_n_ns,
            module_temperature_c,
            hour.irradiance_w_m2,
        )
        voc_bottom_v = irradiance_adjusted_voc_v(
            bottom_voc_stc_v,
            bottom_voc.linear_slope_per_k,
            bottom_n_ns,
            module_temperature_c,
            hour.irradiance_w_m2,
        )
        voc_module_v = (voc_top_v + voc_bottom_v) * args.cell_count

        pmax_2t_w = density_to_module_power_w(
            irradiance_adjusted_pmax_density_mw_cm2(
                tandem_2t_pmax_stc_mw_cm2,
                gamma_2t_per_k,
                module_temperature_c,
                hour.irradiance_w_m2,
            ),
            args.module_area_m2,
        )
        pmax_4t_w = density_to_module_power_w(
            irradiance_adjusted_pmax_density_mw_cm2(
                tandem_4t_pmax_stc_mw_cm2,
                gamma_4t_per_k,
                module_temperature_c,
                hour.irradiance_w_m2,
            ),
            args.module_area_m2,
        )

        hourly_rows.append(
            {
                "timestamp": hour.timestamp.isoformat(timespec="minutes"),
                "irradiance_poa_w_m2": f"{hour.irradiance_w_m2:.2f}",
                "ambient_temperature_c": f"{hour.ambient_c:.2f}",
                "wind_10m_ms": f"{hour.wind_10m_ms:.3f}",
                "wind_module_ms": f"{wind_module_ms:.3f}",
                "module_temperature_c": f"{module_temperature_c:.3f}",
                "reconstructed_flag": int(hour.reconstructed),
                "voc_module_2t_v": f"{voc_module_v:.3f}",
                "voc_module_4t_eq_v": f"{voc_module_v:.3f}",
                "pmax_module_2t_w": f"{pmax_2t_w:.3f}",
                "pmax_module_4t_eq_w": f"{pmax_4t_w:.3f}",
                "energy_2t_wh": f"{pmax_2t_w:.3f}",
                "energy_4t_eq_wh": f"{pmax_4t_w:.3f}",
            }
        )

        voc_module_series.append(voc_module_v)
        pmax_2t_series_w.append(pmax_2t_w)
        pmax_4t_series_w.append(pmax_4t_w)
        if hour.irradiance_w_m2 > 0.0:
            voc_module_sunlit_series.append(voc_module_v)

    max_voc_index = max(range(len(voc_module_series)), key=voc_module_series.__getitem__)
    min_voc_index = min(range(len(voc_module_series)), key=voc_module_series.__getitem__)

    sunlit_indices = [idx for idx, hour in enumerate(pvgis_rows) if hour.irradiance_w_m2 > 0.0]
    min_sunlit_voc_index = min(sunlit_indices, key=voc_module_series.__getitem__)

    summary_rows = [
        {
            "architecture": "Tandem_2T",
            "cell_count": args.cell_count,
            "cell_area_m2": f"{args.module_area_m2 / args.cell_count:.6f}",
            "module_area_m2": f"{args.module_area_m2:.6f}",
            "stc_temperature_c": f"{STC_TEMPERATURE_C:.2f}",
            "stc_irradiance_w_m2": f"{STC_IRRADIANCE_W_M2:.1f}",
            "module_pmax_stc_w": f"{module_pmax_2t_stc_w:.3f}",
            "module_voc_stc_v": f"{module_voc_stc_v:.3f}",
            "beta_tandem_mV_per_k": f"{1000.0 * (top_voc.linear_slope_per_k + bottom_voc.linear_slope_per_k):.3f}",
            "gamma_pmax_percent_per_k": f"{100.0 * gamma_2t_per_k:.3f}",
            "annual_yield_kwh": f"{sum(pmax_2t_series_w) / 1000.0:.3f}",
            "average_module_voc_all_hours_v": f"{float_mean(voc_module_series):.3f}",
            "average_module_voc_sunlit_hours_v": f"{float_mean(voc_module_sunlit_series):.3f}",
            "max_module_voc_v": f"{max(voc_module_series):.3f}",
            "timestamp_of_max_voc": pvgis_rows[max_voc_index].timestamp.isoformat(timespec="minutes"),
            "min_module_voc_all_hours_v": f"{min(voc_module_series):.3f}",
            "timestamp_of_min_voc_all_hours": pvgis_rows[min_voc_index].timestamp.isoformat(
                timespec="minutes"
            ),
            "min_module_voc_sunlit_hours_v": f"{voc_module_series[min_sunlit_voc_index]:.3f}",
            "timestamp_of_min_voc_sunlit_hours": pvgis_rows[
                min_sunlit_voc_index
            ].timestamp.isoformat(timespec="minutes"),
            "faiman_u0_w_m2k": f"{args.faiman_u0:.3f}",
            "faiman_u1_wsm3k": f"{args.faiman_u1:.3f}",
            "module_height_m": f"{args.module_height_m:.2f}",
            "wind_shear_exponent": f"{args.wind_shear_exponent:.3f}",
        },
        {
            "architecture": "Tandem_4T_eq",
            "cell_count": args.cell_count,
            "cell_area_m2": f"{args.module_area_m2 / args.cell_count:.6f}",
            "module_area_m2": f"{args.module_area_m2:.6f}",
            "stc_temperature_c": f"{STC_TEMPERATURE_C:.2f}",
            "stc_irradiance_w_m2": f"{STC_IRRADIANCE_W_M2:.1f}",
            "module_pmax_stc_w": f"{module_pmax_4t_stc_w:.3f}",
            "module_voc_stc_v": f"{module_voc_stc_v:.3f}",
            "beta_tandem_mV_per_k": f"{1000.0 * (top_voc.linear_slope_per_k + bottom_voc.linear_slope_per_k):.3f}",
            "gamma_pmax_percent_per_k": f"{100.0 * gamma_4t_per_k:.3f}",
            "annual_yield_kwh": f"{sum(pmax_4t_series_w) / 1000.0:.3f}",
            "average_module_voc_all_hours_v": f"{float_mean(voc_module_series):.3f}",
            "average_module_voc_sunlit_hours_v": f"{float_mean(voc_module_sunlit_series):.3f}",
            "max_module_voc_v": f"{max(voc_module_series):.3f}",
            "timestamp_of_max_voc": pvgis_rows[max_voc_index].timestamp.isoformat(timespec="minutes"),
            "min_module_voc_all_hours_v": f"{min(voc_module_series):.3f}",
            "timestamp_of_min_voc_all_hours": pvgis_rows[min_voc_index].timestamp.isoformat(
                timespec="minutes"
            ),
            "min_module_voc_sunlit_hours_v": f"{voc_module_series[min_sunlit_voc_index]:.3f}",
            "timestamp_of_min_voc_sunlit_hours": pvgis_rows[
                min_sunlit_voc_index
            ].timestamp.isoformat(timespec="minutes"),
            "faiman_u0_w_m2k": f"{args.faiman_u0:.3f}",
            "faiman_u1_wsm3k": f"{args.faiman_u1:.3f}",
            "module_height_m": f"{args.module_height_m:.2f}",
            "wind_shear_exponent": f"{args.wind_shear_exponent:.3f}",
        },
    ]

    write_csv(
        output_dir / "ouarzazate_2023_tandem_module_hourly.csv",
        [
            "timestamp",
            "irradiance_poa_w_m2",
            "ambient_temperature_c",
            "wind_10m_ms",
            "wind_module_ms",
            "module_temperature_c",
            "reconstructed_flag",
            "voc_module_2t_v",
            "voc_module_4t_eq_v",
            "pmax_module_2t_w",
            "pmax_module_4t_eq_w",
            "energy_2t_wh",
            "energy_4t_eq_wh",
        ],
        hourly_rows,
    )
    write_csv(
        output_dir / "ouarzazate_2023_tandem_module_summary.csv",
        list(summary_rows[0].keys()),
        summary_rows,
    )


if __name__ == "__main__":
    main()
