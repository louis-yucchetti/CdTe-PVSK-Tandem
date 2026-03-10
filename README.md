# PVSK/CdTe Tandem Analysis from SCAPS Exports

This repository collects the SCAPS input files, exported JV and EQE data, and Python post-processing scripts used for a perovskite/CdTe tandem solar-cell project. The code does not run SCAPS itself. It takes SCAPS exports that already exist in the repository, reconstructs tandem metrics, analyzes temperature trends, and turns the extracted thermal coefficients into simple module-level estimates using PVGIS weather data.

The workflow is built around one specific tandem pairing:

- Top cell: `HTM (SpiroOMeTAD) / IDL2 / (Cs-FA-MA)(4-81-15)Pb(I-Br)_3(83-17) / IDL1 / BL (SnO2) / ITO`
- Bottom cell: `CdTe / CdS / SnOx`

In the tandem case, the CdTe device is evaluated under the spectrum transmitted through the top perovskite cell, so the filtered CdTe current is much lower than the standalone CdTe current. That current bottleneck drives most of the 2T tandem behavior in this project.

## What This Repository Contains

```text
.
|-- data/
|   |-- PVGIS/                  # Hourly weather and irradiance data
|   `-- scaps_exports/          # SCAPS JV and QE exports used by the scripts
|       |-- iv/
|       `-- qe/
|-- devices/
|   |-- cdte/
|   |   |-- absorption/
|   |   `-- scaps/
|   |-- pvsk/
|   |   |-- band_diagrams/
|   |   `-- scaps/
|   `-- tandem/
|       `-- scaps/
|-- results/
|   |-- csv/                    # Single-temperature tandem outputs
|   |-- figures/
|   |-- pvgis/
|   |   `-- csv/                # Module-yield outputs
|   `-- thermal/
|       |-- csv/
|       `-- figures/
|-- scripts/
|   |-- tandem_analysis.py
|   |-- thermal_tandem_analysis.py
|   `-- ouarzazate_module_yield_analysis.py
|-- OUTDATED/                   # Archived material kept for reference
`-- Solar Cells Modelling kits/ # SCAPS tutorials and course material, left untouched
```

## Project Scope

The repository is organized around three linked questions:

1. What do the isolated PVSK top cell, filtered CdTe bottom cell, and reconstructed tandems look like at one temperature?
2. How do `Voc`, `Jsc`, `FF`, `Pmax`, and efficiency change between `300 K` and `350 K`?
3. What do those thermal trends imply for a simple module-level projection in Ouarzazate using hourly PVGIS weather data?

## Main Scripts

### `scripts/tandem_analysis.py`

Purpose:

- Reads single-temperature SCAPS `.iv` and `.qe` exports.
- Extracts standard PV metrics for the PVSK top cell, the filtered CdTe bottom cell, and the unfiltered CdTe comparison case.
- Reconstructs a `2T` tandem by combining the two sub-cells in the current domain.
- Computes an equivalent `4T` summary by adding the independent sub-cell powers.
- Saves publication-style JV, EQE, and 4T power plots plus CSV tables.

Default inputs:

- `data/scaps_exports/iv/PVSK_louis.iv`
- `data/scaps_exports/iv/CdTe_filtered_louis.iv`
- `data/scaps_exports/iv/CdTe_louis.iv`
- `data/scaps_exports/qe/PVSK_louis.qe`
- `data/scaps_exports/qe/CdTe_filtered_louis.qe`
- `data/scaps_exports/qe/CdTe_louis.qe`

Key outputs:

- `results/figures/tandem_jv_curve.png`
- `results/figures/tandem_eqe_curve.png`
- `results/figures/tandem_4t_power_curves.png`
- `results/figures/tandem_4t_power_map.png`
- `results/csv/tandem_metrics.csv`
- `results/csv/tandem_jv_curve.csv`
- `results/csv/tandem_eqe_curve.csv`

### `scripts/thermal_tandem_analysis.py`

Purpose:

- Parses SCAPS batch JV exports for the perovskite top cell and filtered CdTe bottom cell across temperature.
- Extracts `Voc`, `Jsc`, `FF`, `eta`, `Vmpp`, `Jmpp`, and `Pmax` at every temperature step.
- Builds a reconstructed `2T` tandem and an equivalent `4T` tandem at every temperature.
- Uses the SCAPS recombination columns, when present, to track bulk and non-bulk contributions at `Voc` and collection loss at `0 V`.
- Writes figures, summary tables, a text summary, and a compact standardized reporting matrix.

Default inputs:

- `data/scaps_exports/iv/PVSK_thermal_sweep.iv`
- `data/scaps_exports/iv/CdTe_filtered_thermal_sweep.iv`
- `devices/pvsk/scaps/Louis_PVSK.scaps`
- `devices/cdte/scaps/Group1_CdTe_with_PVSK_filter.scaps`

The two `.scaps` device files are only read to identify absorber-layer names and nominal bandgaps for the plain-text summary. The thermal calculations themselves come from the exported SCAPS JV sweep.

Key outputs:

- `results/thermal/figures/pvsk_metrics_vs_temperature.png`
- `results/thermal/figures/cdte_filtered_metrics_vs_temperature.png`
- `results/thermal/figures/subcell_metrics_vs_temperature.png`
- `results/thermal/figures/subcell_metric_separation_vs_temperature.png`
- `results/thermal/figures/tandem_2t_metrics_vs_temperature.png`
- `results/thermal/figures/tandem_4t_metrics_vs_temperature.png`
- `results/thermal/figures/tandem_metrics_vs_temperature.png`
- `results/thermal/figures/voc_physics_vs_temperature.png`
- `results/thermal/figures/mismatch_drift_vs_temperature.png`
- `results/thermal/figures/jtop_over_jbottom_vs_temperature.png`
- `results/thermal/csv/subcell_temperature_metrics.csv`
- `results/thermal/csv/tandem_temperature_metrics.csv`
- `results/thermal/csv/temperature_metric_summary.csv`
- `results/thermal/csv/temperature_coefficients.csv`
- `results/thermal/csv/mismatch_drift_vs_temperature.csv`
- `results/thermal/csv/mismatch_drift_summary.csv`
- `results/thermal/csv/pillar5_standardized_data_matrix.csv`
- `results/thermal/pillar5_standardized_data_matrix.md`
- `results/thermal/thermal_summary.txt`

### `scripts/ouarzazate_module_yield_analysis.py`

Purpose:

- Reads hourly PVGIS weather data and the thermal summary generated by `thermal_tandem_analysis.py`.
- Extrapolates the thermal sweep from the `300 K` baseline to exact STC `25 C = 298.15 K`.
- Estimates module temperature with the Faiman model.
- Applies the extracted `Voc` and `Pmax` temperature coefficients to a simple `60`-cell tandem module.
- Writes hourly and annual summary CSV files for Ouarzazate `2023`.

Default inputs:

- `data/PVGIS/Timeseries_30.937_-6.906_SA3_32deg_-4deg_2023_2023.csv`
- `results/thermal/csv/temperature_metric_summary.csv`

Default assumptions:

- `60` series-connected tandem cells
- module area `1.65 m^2`
- `U0 = 25 W/m^2/K`
- `U1 = 6.84 W s/m^3/K`
- module-height wind speed estimated at `2 m`
- wind-shear exponent `0.14`

Key outputs:

- `results/pvgis/csv/ouarzazate_2023_tandem_module_hourly.csv`
- `results/pvgis/csv/ouarzazate_2023_tandem_module_summary.csv`

## Equations Used in the Analysis

### 1. PV metrics extracted from JV data

The scripts read SCAPS voltage-current exports, interpolate where needed, and compute:

$$
P(V) = -VJ
$$

The minus sign is required because the SCAPS export uses the generator sign convention with negative current in the power-producing quadrant.

From that curve the code extracts:

$$
FF = 100 \times \frac{P_{max}}{V_{oc}J_{sc}}
$$

$$
\eta = 100 \times \frac{P_{max}}{P_{in}}
$$

where `P_in = 100 mW/cm^2` by default.

### 2. Current from EQE

The README figures compare measured or simulated EQE curves, but the physical link to current is still:

$$
J_{sc} = q \int \Phi_{ph}(\lambda)\,EQE(\lambda)\,d\lambda
$$

For the filtered CdTe cell, the relevant photon flux is the transmitted spectrum after the perovskite top cell, not the original AM1.5G spectrum.

### 3. Constructed optical tandem EQE

The optical tandem EQE shown by `tandem_analysis.py` is a constructed quantity:

$$
EQE_{tandem,opt}(\lambda) = \min\left(EQE_{top}(\lambda) + EQE_{bottom,filtered}(\lambda), 100\%\right)
$$

This is useful for visualization, but it is not a direct electrical 2T measurement.

### 4. Reconstructed `2T` tandem

A series-connected tandem must carry the same current through both sub-cells, so the script works in the current domain and builds:

$$
V_{2T}(J) = V_{top}(J) + V_{bottom}(J)
$$

This is why `tandem_analysis.py` interpolates each sub-cell as `V(J)` on a shared current grid rather than simply adding voltages point-by-point at equal applied voltage.

### 5. Equivalent `4T` reporting

The repository also reports an equivalent `4T` case in which each sub-cell operates independently at its own maximum-power point:

$$
P_{max,4T} = P_{max,top} + P_{max,bottom}
$$

$$
\eta_{4T} = \eta_{top} + \eta_{bottom}
$$

For table convenience it also reports:

$$
V_{oc,4T(eq)} = V_{oc,top} + V_{oc,bottom}
$$

$$
J_{sc,4T(eq)} = J_{sc,top} + J_{sc,bottom}
$$

Those `4T` values are bookkeeping quantities, not a unique four-terminal JV curve.

### 6. Thermal coefficients

The temperature-sweep script fits linear slopes to the metrics exported by SCAPS:

$$
x(T) \approx x(T_{ref}) + \frac{dx}{dT}(T - T_{ref})
$$

It then reports normalized coefficients in the usual form:

$$
\beta_{Voc} = \frac{1}{V_{oc,ref}}\frac{dV_{oc}}{dT}\times 100 \;\; [\%/K]
$$

$$
\gamma_{Pmax} = \frac{1}{P_{max,ref}}\frac{dP_{max}}{dT}\times 100 \;\; [\%/K]
$$

The scripts also keep the absolute slopes:

- `beta_Voc` in `mV/K`
- `gamma_Pmax` in `mW/cm^2/K`

Important clarification: for the temperature sweep, temperature was varied in SCAPS and the Python code post-processed those exports. The repository does not inject separate external `E_g(T)` or mobility-versus-temperature datasets on top of the SCAPS sweep.

### 7. Recombination and collection-loss indicators

When the SCAPS batch export includes recombination channels, `thermal_tandem_analysis.py` derives a few simple indicators:

$$
f_{nonbulk,Voc} = \frac{J_{total,rec}(V_{oc}) - J_{bulk}(V_{oc})}{J_{total,rec}(V_{oc})}
$$

$$
f_{collection,0V} = \frac{J_{total,rec}(0)}{J_{total,gen}(0)}
$$

These are not new physical models. They are compact ways to summarize the SCAPS recombination output at `Voc` and short-circuit conditions.

### 8. Current mismatch in the tandem

To track whether the current match improves or worsens with temperature, the script uses:

$$
\Delta J(T) = J_{top}(T) - J_{bottom}(T)
$$

$$
M(T) = \frac{J_{top}(T)}{J_{bottom}(T)}
$$

If either quantity increases, the current mismatch is widening.

### 9. Module temperature and yield model

The Ouarzazate workflow first converts PVGIS `10 m` wind speed to an approximate module-height wind speed:

$$
u_{module} = u_{10m}\left(\frac{h_{module}}{10}\right)^{0.14}
$$

Module temperature is then estimated with the Faiman relation:

$$
T_{mod} = T_{air} + \frac{G_{POA}}{U_0 + U_1 u_{module}}
$$

Power is scaled from the STC baseline with irradiance and the extracted `Pmax` coefficient:

$$
P_{max}(G,T) = P_{max,STC}\left(\frac{G}{1000}\right)\left[1 + \gamma_{Pmax}(T_{mod} - 25)\right]
$$

For `Voc`, the script combines the linear temperature slope with a logarithmic irradiance term:

$$
V_{oc}(G,T) = V_{oc,STC} + \beta_{Voc}(T_{mod} - 25) + nN_sV_T(T)\ln\left(\frac{G}{1000}\right)
$$

where:

- `V_T(T) = kT/q`
- `nN_s` is inferred from the STC `Voc` and `FF`

The inference uses the standard dimensionless-`Voc` approximation:

$$
FF \approx \frac{v_{oc}^* - \ln(v_{oc}^* + 0.72)}{v_{oc}^* + 1}
$$

followed by:

$$
nN_s = \frac{V_{oc,STC}}{v_{oc}^*V_{T,STC}}
$$

At module level, the script assumes `60` series-connected tandem cells and scales power density to power with the module area.

## Running the Scripts

The repository does not currently include a `requirements.txt` or `pyproject.toml`. The active scripts use:

- Python `3.10+`
- `numpy`
- `matplotlib`

Typical usage:

```bash
python scripts/tandem_analysis.py
python scripts/thermal_tandem_analysis.py --target-temp-c 65
python scripts/ouarzazate_module_yield_analysis.py
```

Useful CLI options:

- `tandem_analysis.py`: `--outdir`, `--pvsk-iv`, `--cdte-iv`, `--pvsk-qe`, `--cdte-qe`
- `thermal_tandem_analysis.py`: `--baseline-temp-k`, `--target-temp-c`, `--outdir`
- `ouarzazate_module_yield_analysis.py`: `--pvgis-file`, `--thermal-summary`, `--output-dir`, `--cell-count`, `--module-area-m2`

## Current Headline Results in the Repository

The tables below summarize the outputs already checked into `results/`.

### Single-temperature metrics (`results/csv/tandem_metrics.csv`)

| Case | `Voc` (V) | `Jsc` (mA/cm²) | `FF` (%) | `eta` (%) |
| --- | ---: | ---: | ---: | ---: |
| `CdTe_unfiltered` | 0.8529 | 26.0557 | 79.22 | 17.61 |
| `CdTe_filtered` | 0.7915 | 5.4030 | 74.33 | 3.18 |
| `PVSK_top` | 1.2854 | 20.5430 | 81.79 | 21.60 |
| `Tandem_2T` | 2.0769 | 5.5714 | 81.04 | 9.38 |
| `Tandem_4T` | 2.0769 | 25.9460 | 45.98 | 24.78 |

The main point is immediate: once the CdTe cell is filtered by the top perovskite, it becomes strongly current-limited, so the reconstructed `2T` tandem is bottom-cell-limited.

### Thermal coefficients from `300 K` to `65 C` (`results/thermal/csv/temperature_coefficients.csv`)

| Case | `beta_Voc` (mV/K) | `beta_Voc` (%/K) | `gamma_Pmax` (mW/cm²/K) | `gamma_Pmax` (%/K) |
| --- | ---: | ---: | ---: | ---: |
| `PVSK_top` | -1.373 | -0.1068 | -0.02990 | -0.1385 |
| `CdTe_filtered` | -2.761 | -0.3488 | -0.01377 | -0.4330 |
| `Tandem_2T` | -4.134 | -0.1990 | -0.02079 | -0.2217 |
| `Tandem_4T` | -4.134 | -0.1990 | -0.04367 | -0.1763 |

The filtered CdTe bottom cell is much more temperature-sensitive in voltage than the PVSK top cell.

### Current mismatch drift (`results/thermal/csv/mismatch_drift_summary.csv`)

From the generated thermal outputs:

- `Jtop - Jbottom` increases from `15.1400 mA/cm²` at `300 K` to `15.1692 mA/cm²` at `65 C`.
- `Jtop / Jbottom` increases from `3.8021` to `3.8195`.

So the mismatch widens slightly with temperature rather than closing.

### Ouarzazate `2023` module projection (`results/pvgis/csv/ouarzazate_2023_tandem_module_summary.csv`)

| Architecture | `Pmax @ STC` (W) | `Voc @ STC` (V) | Annual yield (kWh) |
| --- | ---: | ---: | ---: |
| `Tandem_2T` | 155.366 | 125.072 | 376.823 |
| `Tandem_4T_eq` | 410.133 | 125.072 | 1004.461 |

These are simple projections based on the extracted thermal slopes and the Faiman temperature model, not a full system simulation.

## Reading the Output Folders

- `results/csv/` and `results/figures/` contain the one-temperature tandem reconstruction.
- `results/thermal/csv/` and `results/thermal/figures/` contain the temperature sweep, tandem comparison, mismatch tracking, and reporting tables.
- `results/pvgis/csv/` contains the hourly Ouarzazate module estimates and the annual summary.

If you only want the main deliverables, start with:

- `results/csv/tandem_metrics.csv`
- `results/thermal/csv/temperature_coefficients.csv`
- `results/thermal/csv/mismatch_drift_summary.csv`
- `results/thermal/thermal_summary.txt`
- `results/pvgis/csv/ouarzazate_2023_tandem_module_summary.csv`

## Notes and Limits

- The `4T` numbers are equivalent reporting quantities, not a measured four-terminal JV trace.
- The thermal sweep starts at `300 K = 26.85 C`, not exact STC `25 C = 298.15 K`.
- The module-yield script therefore extrapolates from the `300 K` baseline to STC before applying weather-driven irradiance and temperature changes.
- The Ouarzazate workflow uses default Faiman coefficients rather than a module-specific thermal fit.
- `Solar Cells Modelling kits/` and `OUTDATED/` are kept for reference and archive value; they are not part of the active Python workflow.
