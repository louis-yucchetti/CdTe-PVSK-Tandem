# PVSK/CdTe Tandem Analysis From SCAPS Exports

## Overview

This repository turns SCAPS exports and PVGIS weather data into reproducible figures, tables, and physics interpretations for a perovskite/CdTe tandem solar-cell project.

It supports three complementary workflows:

- `scripts/tandem_analysis.py` for single-temperature JV, EQE, 2T tandem, and equivalent 4T analysis.
- `scripts/thermal_tandem_analysis.py` for the `300-350 K` thermal sweep, including sub-cell trends, 2T/4T tandem trends, and current-mismatch drift at `65 C`.
- `scripts/ouarzazate_module_yield_analysis.py` for PVGIS-driven hourly module temperature, `Voc`, `Pmax`, and annual yield estimates for a `60`-cell tandem module in Ouarzazate using the extracted thermal coefficients and the Faiman model.

The analysis is built around your specific simulated devices:

- Top cell: `HTM (SpiroOMeTAD) / IDL2 / (Cs-FA-MA)(4-81-15)Pb(I-Br)_3(83-17) / IDL1 / BL (SnO2) / ITO`
- Bottom cell: `CdTe / CdS / SnOx`, operated under the spectral filter transmitted by the top cell

This means the repository is not only a plotting utility. It is also a compact scientific notebook for explaining:

- why the filtered CdTe bottom cell carries much less current than the PVSK top cell,
- how 2T and 4T tandem definitions differ physically,
- why `Voc` drops faster in one sub-cell than the other,
- and whether current mismatch closes or widens at elevated temperature.

## Repository Layout

```text
.
|-- data/
|   |-- PVGIS/
|   `-- scaps_exports/
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
|   |-- csv/
|   |-- figures/
|   |-- pvgis/
|   |   `-- csv/
|   `-- thermal/
|       |-- csv/
|       `-- figures/
`-- scripts/
    |-- ouarzazate_module_yield_analysis.py
    |-- tandem_analysis.py
    `-- thermal_tandem_analysis.py
```

`Solar Cells Modelling kits` is intentionally left untouched. Any folder named `OUTDATED` is also intentionally kept as archive material rather than active analysis.

## Project Context

For the course, the graded component is the oral defense on **March 11, 2026**. The written manuscript is not graded, but it matters as scientific project documentation and possible paper material.

Because of that, the priority is not curve fitting for its own sake. The priority is a physically defensible story:

- what each layer does,
- how the tandem shares the spectrum,
- where current limitation comes from,
- how recombination changes with temperature,
- and what those trends imply for tandem operation.

## Main Scientific Questions

This repository is organized around six questions:

1. What are the electrical metrics of the isolated PVSK top cell and the filtered CdTe bottom cell?
2. What happens when those two sub-cells are combined into a `2T` series tandem?
3. What is the corresponding `4T` equivalent when each sub-cell operates independently?
4. How do `Voc`, `Jsc`, `FF`, and efficiency vary between `300 K` and `350 K`?
5. At `65 C = 338.15 K`, does the current mismatch between top and bottom sub-cells shrink or widen?
6. What module-level `Voc`, `Pmax`, and annual energy yield do the extracted tandem coefficients predict for `2023` hourly weather in Ouarzazate?

## Device Physics

## Top Perovskite Cell

The top cell is a wide-gap perovskite device with a nominal absorber bandgap close to `1.60 eV`. In the SCAPS stack:

- `SpiroOMeTAD` acts as the hole-transport material.
- `SnO2` acts as the electron-selective transport or blocking layer.
- The perovskite absorber captures the high-energy part of the AM1.5G spectrum.
- The thin `IDL` layers model non-ideal interfaces and their recombination activity.

A wide-gap top cell is useful in tandem because it converts short-wavelength photons efficiently while transmitting part of the red and near-infrared spectrum to the bottom cell.

## Filtered CdTe Bottom Cell

The bottom device is a CdTe heterojunction operated under the transmitted spectrum of the top cell rather than under full AM1.5G.

- `CdTe` is the main absorber, with nominal `Eg ~ 1.50 eV`.
- `CdS` is the buffer layer.
- `SnOx` is the contact/selective layer.

In a real tandem context, the bottom cell does not see the full solar spectrum. It sees only the photons not already absorbed by the top cell. That is why the filtered CdTe `Jsc` is much lower than the unfiltered CdTe `Jsc`, and also why current matching becomes the central 2T design constraint.

## Why the Bottom Cell Current Is Small

The filtered CdTe current is low for a simple physical reason: the top perovskite cell harvests most of the photons above its bandgap first.

The short-circuit current density is related to EQE by

$$
J_{sc} = q \int \Phi_{ph}(\lambda)\, EQE(\lambda)\, d\lambda
$$

where:

- `q` is the elementary charge,
- `\Phi_{ph}(\lambda)` is the incident photon flux,
- `EQE(\lambda)` is the external quantum efficiency.

For the filtered bottom cell, the relevant photon flux is not the original AM1.5G spectrum but the transmitted spectrum after top-cell absorption and parasitic optical losses. That is why a bottom cell that looks good in single-junction mode can become current-starved in tandem mode.

## Core Solar-Cell Equations

## Single-Junction JV Model

A practical illuminated diode model is

$$
J(V) = J_{ph} - J_0 \left[\exp\left(\frac{q(V + J R_s)}{n k T}\right) - 1\right] - \frac{V + J R_s}{R_{sh}}
$$

where:

- `J_ph` is the photocurrent density,
- `J_0` is the dark saturation current density,
- `n` is the ideality factor,
- `R_s` is the series resistance,
- `R_sh` is the shunt resistance.

The SCAPS outputs are more detailed than this compact equation because they also separate recombination channels, but this expression is still the best starting point for interpretation.

## Key Performance Metrics

The repository extracts the standard photovoltaic quantities:

- `Voc`: voltage where `J = 0`
- `Jsc`: current magnitude at `V = 0`
- `P(V) = -VJ` using the SCAPS sign convention
- `Pmax`: maximum generated power density
- `FF = 100 * Pmax / (Voc * Jsc)`
- `eta = 100 * Pmax / P_in`

Equivalently,

$$
\eta = \frac{V_{oc} J_{sc} FF}{P_{in}}
$$

with `FF` treated as a unitless fraction in that expression.

## Why Voc Is So Sensitive to Recombination

Ignoring series and shunt losses, a useful approximation is

$$
V_{oc} \approx \frac{n k T}{q}\ln\left(\frac{J_{ph}}{J_0} + 1\right)
$$

This equation explains several important facts at once:

- `Voc` increases when photocurrent increases.
- `Voc` decreases when dark current `J0` increases.
- `Voc` is very sensitive to recombination because stronger recombination usually means larger `J0`.
- A lower-gap absorber usually suffers a stronger `Voc` penalty because it has less thermodynamic voltage headroom.

That is the core reason temperature analysis is so important: temperature usually increases recombination-related dark current faster than it improves carrier collection.

## 2T and 4T Tandem Physics

## 2T Tandem

A `2T` tandem is a series-connected device. The same current must pass through both sub-cells, so the correct electrical construction is performed in the current domain:

$$
V_{2T}(J) = V_{top}(J) + V_{bottom}(J)
$$

The tandem current is therefore constrained by the weaker current-producing sub-cell over the operating range. In practice, for this repository, that means the filtered CdTe bottom cell strongly controls the current side of the `2T` tandem.

## 4T Tandem

A true `4T` tandem has independent terminals for each sub-cell. There is no unique single JV curve for the full device because each sub-cell can operate at its own maximum-power point.

The repository therefore reports an equivalent aggregated `4T` definition:

$$
P_{max,4T} = P_{max,top} + P_{max,bottom}
$$

$$
\eta_{4T} = \eta_{top} + \eta_{bottom}
$$

and for table convenience:

$$
V_{oc,4T(eq)} = V_{oc,top} + V_{oc,bottom}
$$

$$
J_{sc,4T(eq)} = J_{sc,top} + J_{sc,bottom}
$$

These are useful reporting quantities, but they should not be confused with a unique measured four-terminal JV trace.

## Temperature Physics

The thermal analysis assumes the temperature-dependent SCAPS models are active for at least:

- `Eg(T)`
- `mu(T)`
- `Nc(T)`
- `Nv(T)`

Several linked effects matter:

- The bandgap usually decreases with temperature.
- The effective density of states changes with temperature.
- Carrier mobility can change with phonon scattering and transport barriers.
- Recombination pathways become more thermally activated.
- Dark current therefore rises, which pushes `Voc` downward.

A common bandgap model is the Varshni form

$$
E_g(T) = E_g(0) - \frac{\alpha T^2}{T + \beta}
$$

You do not explicitly fit `\alpha` and `\beta` in this repository, but the equation is helpful for understanding why a lower-gap device often becomes more vulnerable to thermal voltage loss.

Thermal coefficients are reported in the usual normalized form:

$$
\beta_{Voc} = \frac{1}{V_{oc,ref}}\frac{dV_{oc}}{dT}\times 100 \;\; [\%/K]
$$

$$
\gamma_{Pmax} = \frac{1}{P_{max,ref}}\frac{dP_{max}}{dT}\times 100 \;\; [\%/K]
$$

The repository also stores the corresponding absolute slopes:

- `beta_Voc` in `mV/K`
- `gamma_Pmax` in `mW/cm^2/K`

## Module Temperature and Yield Model

For the Ouarzazate field-yield workflow, the repository converts hourly PVGIS weather data into module operating conditions.

Module temperature is estimated with the Faiman relation

$$
T_{mod} = T_{air} + \frac{G_{POA}}{U_0 + U_1 u}
$$

using the default coefficients `U0 = 25 W/m^2/K` and `U1 = 6.84 W s/m^3/K`.

The PVGIS file provides wind speed at `10 m`, so the script first downscales that to an approximate module-height wind speed with

$$
u_{2m} = u_{10m}\left(\frac{2}{10}\right)^{0.14}
$$

to avoid overestimating convective cooling.

The module scaling assumes:

- `60` series-connected tandem cells,
- total module area `1.65 m^2`,
- cell area `1.65 / 60 = 0.0275 m^2 = 275 cm^2`.

The script extrapolates the SCAPS-derived metrics from the `300 K` thermal baseline to exact STC `25 C = 298.15 K`, then applies:

$$
P_{max}(G,T) = P_{max,STC}\left(\frac{G}{1000}\right)\left[1 + \gamma_{Pmax}(T_{mod} - 25)\right]
$$

For `Voc`, the script combines the extracted linear temperature slope with a logarithmic irradiance term,

$$
V_{oc}(G,T) = V_{oc,STC} + \beta_{Voc}(T_{mod} - 25) + nN_sV_T(T)\ln\left(\frac{G}{1000}\right)
$$

where the effective `nN_s` is inferred from the STC `Voc` and `FF` so that low-irradiance `Voc` is not artificially flat.

At module level:

- `Voc,module = 60 * (Voc,top + Voc,bottom)`
- the `4T` module `Voc` reported here is the same equivalent sum as in the SCAPS-based `4T` workflow,
- `Pmax,module = Pmax,density * module_area`.

## Why Voc Can Drop Faster in One Sub-Cell Than the Other

For your specific data, the filtered CdTe bottom cell loses `Voc` faster than the PVSK top cell.

The physical interpretation is:

- CdTe has the lower absorber gap, so it starts with less voltage headroom.
- The filtered CdTe cell also operates with much lower photocurrent than the top PVSK cell.
- The SCAPS channel decomposition shows a much larger non-bulk contribution to total recombination at `Voc` in the filtered CdTe cell than in the PVSK cell.
- As temperature rises, that non-bulk or leakage-related contribution grows more strongly in CdTe.

So even though both devices are SRH-dominated in the present simulations, the CdTe bottom cell pays a larger `Voc` price per kelvin.

## Current Mismatch Drift

For tandem design, the most relevant thermal comparison is not only absolute `Jsc`, but the mismatch between sub-cells:

$$
\Delta J(T) = J_{top}(T) - J_{bottom}(T)
$$

$$
M(T) = \frac{J_{top}(T)}{J_{bottom}(T)}
$$

If `\Delta J` or `M` increases with temperature, current matching is getting worse for a `2T` tandem. If they decrease, current matching is improving.

In your present results, the mismatch widens slightly with temperature because the PVSK current stays almost constant while the filtered CdTe current decreases slowly.

## Scripts

## `scripts/tandem_analysis.py`

This script handles the single-temperature workflow:

- reads single-shot SCAPS `.iv` and `.qe` exports,
- extracts sub-cell metrics,
- constructs a full `2T` tandem JV curve in the current domain,
- computes equivalent `4T` metrics,
- generates JV, EQE, and 4T power plots,
- writes CSV summaries.

Default inputs:

- `data/scaps_exports/iv/PVSK_louis.iv`
- `data/scaps_exports/iv/CdTe_filtered_louis.iv`
- `data/scaps_exports/iv/CdTe_louis.iv`
- `data/scaps_exports/qe/PVSK_louis.qe`
- `data/scaps_exports/qe/CdTe_filtered_louis.qe`
- `data/scaps_exports/qe/CdTe_louis.qe`

## `scripts/thermal_tandem_analysis.py`

This script handles the temperature sweep:

- parses batch SCAPS JV files for each temperature,
- extracts `Voc`, `Jsc`, `FF`, `eta`, `Vmpp`, `Jmpp`, and `Pmax`,
- evaluates SCAPS recombination channels at `Voc` and `0 V`,
- constructs `2T` and equivalent `4T` tandems at every temperature,
- generates separate absolute plots for PVSK, CdTe, `2T`, and `4T`,
- generates normalized comparison plots so thermal sensitivity is easy to see,
- generates mismatch drift plots,
- writes CSV tables and a text summary.

Default inputs:

- `data/scaps_exports/iv/PVSK_thermal_sweep.iv`
- `data/scaps_exports/iv/CdTe_filtered_thermal_sweep.iv`
- `devices/pvsk/scaps/Louis_PVSK.scaps`
- `devices/cdte/scaps/Group1_CdTe_with_PVSK_filter.scaps`

## `scripts/ouarzazate_module_yield_analysis.py`

This script handles the PVGIS-driven module-yield workflow:

- reads hourly PVGIS data for `G(i)`, `T2m`, and `WS10m`,
- reads `results/thermal/csv/temperature_metric_summary.csv`,
- extrapolates thermal-sweep metrics from `300 K` to exact STC,
- computes hourly module temperature with the Faiman model,
- applies the extracted `beta` and `gamma` coefficients to `Voc` and `Pmax`,
- scales the tandem to a `60`-cell, `1.65 m^2` module,
- writes hourly and annual summary CSV files.

Default inputs and assumptions:

- `data/PVGIS/Timeseries_30.937_-6.906_SA3_32deg_-4deg_2023_2023.csv`
- `results/thermal/csv/temperature_metric_summary.csv`
- `60` cells
- `1.65 m^2` module area
- `U0 = 25`, `U1 = 6.84`
- module height `2 m`
- wind-shear exponent `0.14`

## Usage

Run the single-temperature analysis:

```bash
python scripts/tandem_analysis.py
```

Choose another output folder:

```bash
python scripts/tandem_analysis.py --outdir results
```

Run the thermal analysis:

```bash
python scripts/thermal_tandem_analysis.py
```

Set another target temperature for the drift discussion:

```bash
python scripts/thermal_tandem_analysis.py --target-temp-c 65
```

Run the Ouarzazate module-yield analysis:

```bash
python scripts/ouarzazate_module_yield_analysis.py
```

Point it to another PVGIS file or output folder:

```bash
python scripts/ouarzazate_module_yield_analysis.py --pvgis-file data/PVGIS/your_file.csv --output-dir results/pvgis/csv
```

## Outputs

## Single-Temperature Outputs

The single-temperature script writes to `results/figures/` and `results/csv/`:

- `tandem_jv_curve.png`
- `tandem_eqe_curve.png`
- `tandem_4t_power_curves.png`
- `tandem_4t_power_map.png`
- `tandem_metrics.csv`
- `tandem_jv_curve.csv`
- `tandem_eqe_curve.csv`

## Thermal Outputs

The thermal script writes to `results/thermal/figures/` and `results/thermal/csv/`.

Figures:

- `pvsk_metrics_vs_temperature.png`
- `cdte_filtered_metrics_vs_temperature.png`
- `subcell_metrics_vs_temperature.png`
- `subcell_metric_separation_vs_temperature.png`
- `tandem_2t_metrics_vs_temperature.png`
- `tandem_4t_metrics_vs_temperature.png`
- `tandem_metrics_vs_temperature.png`
- `voc_physics_vs_temperature.png`
- `mismatch_drift_vs_temperature.png`

CSV files:

- `subcell_temperature_metrics.csv`
- `tandem_temperature_metrics.csv`
- `temperature_metric_summary.csv`
- `temperature_coefficients.csv`
- `mismatch_drift_vs_temperature.csv`
- `mismatch_drift_summary.csv`

Text summary:

- `thermal_summary.txt`

## PVGIS Module-Yield Outputs

The PVGIS module-yield script writes to `results/pvgis/csv/`:

- `ouarzazate_2023_tandem_module_hourly.csv`
- `ouarzazate_2023_tandem_module_summary.csv`

## How To Read the Plots

The plotting strategy is intentionally split into absolute and normalized views.

- Use `pvsk_metrics_vs_temperature.png` and `cdte_filtered_metrics_vs_temperature.png` to see the absolute behavior of each sub-cell without wasting axis range.
- Use `subcell_metrics_vs_temperature.png` to compare thermal sensitivity directly because each curve is normalized to its own `300 K` value.
- Use `subcell_metric_separation_vs_temperature.png` to see how far apart the two sub-cells are in `Voc`, `Jsc`, `FF`, and efficiency.
- Use `tandem_2t_metrics_vs_temperature.png` and `tandem_4t_metrics_vs_temperature.png` for absolute tandem behavior.
- Use `tandem_metrics_vs_temperature.png` to compare normalized thermal stability of `2T` and `4T`.
- Use `voc_physics_vs_temperature.png` to connect the `Voc(T)` slope to recombination-channel severity.
- Use `mismatch_drift_vs_temperature.png` for the main paper result on current drift at elevated temperature.

## Current Results From the Present Dataset

Important caveat: the thermal sweep starts at `300 K = 26.85 C`, not exact STC `25 C = 298.15 K`. Therefore, the repository reports `300 K` directly and interpolates `65 C = 338.15 K` inside the available sweep.

### Sub-Cells

| Cell | `Voc @ 300 K` (V) | `Voc @ 65 C` (V) | `beta_Voc` (%/K) | `Jsc @ 300 K` (mA/cm^2) | `Jsc @ 65 C` (mA/cm^2) | `eta @ 300 K` (%) | `eta @ 65 C` (%) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `PVSK_top` | 1.2854 | 1.2336 | -0.1068 | 20.5430 | 20.5494 | 21.5970 | 20.4580 |
| `CdTe_filtered` | 0.7915 | 0.6864 | -0.3488 | 5.4030 | 5.3801 | 3.1788 | 2.6552 |

### Tandems

| Tandem | `Voc @ 300 K` (V) | `Voc @ 65 C` (V) | `beta_Voc` (%/K) | `Jsc @ 300 K` (mA/cm^2) | `Jsc @ 65 C` (mA/cm^2) | `eta @ 300 K` (%) | `eta @ 65 C` (%) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `Tandem_2T` | 2.0769 | 1.9200 | -0.1990 | 5.5714 | 5.5632 | 9.3777 | 8.5922 |
| `Tandem_4T(eq)` | 2.0769 | 1.9200 | -0.1990 | 25.9460 | 25.9295 | 24.7758 | 23.1132 |

### Ouarzazate `2023` Module Projection

Assumptions here are the default ones used by `scripts/ouarzazate_module_yield_analysis.py`: `60` series cells, total area `1.65 m^2`, and Faiman thermal modeling.

| Architecture | `Pmax @ STC` (W) | `Voc @ STC` (V) | Annual yield (kWh) | Average `Voc` over all `8760` h (V) | Average `Voc` over sunlit hours (V) | Max `Voc` (V) | Min `Voc` over all hours (V) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `Tandem_2T` | 155.366 | 125.072 | 376.823 | 56.715 | 114.871 | 124.965 | 0.000 |
| `Tandem_4T(eq)` | 410.133 | 125.072 | 1004.461 | 56.715 | 114.871 | 124.965 | 0.000 |

Notes for reading those numbers:

- The `Voc` values are identical for `2T` and `4T(eq)` because both are reported as the sum of top and bottom sub-cell `Voc`.
- The minimum `Voc` over all hours is `0 V` because night hours are included.
- The minimum sunlit `Voc` is `69.196 V`, occurring at low-irradiance hot conditions near sunset.

### Thermal Conclusions

- The filtered CdTe bottom cell loses `Voc` much faster than the PVSK top cell: about `-2.761 mV/K` versus `-1.373 mV/K`.
- The current mismatch widens slightly with temperature rather than closing.
- At `300 K`, `Jtop - Jbottom = 15.1400 mA/cm^2`.
- At `65 C`, `Jtop - Jbottom = 15.1692 mA/cm^2`.
- The mismatch ratio `Jtop / Jbottom` increases from `3.8021` to `3.8195`.
- The 2T tandem therefore remains bottom-cell-limited and becomes very slightly more mismatched at high temperature.

## Practical Interpretation for the Paper

The most defensible physics message from the present simulations is:

1. The perovskite top cell is thermally more robust in voltage than the filtered CdTe bottom cell.
2. The bottom cell is already current-starved by spectral filtering, so any additional thermal current loss directly worsens 2T matching.
3. The steeper CdTe `Voc` drop is consistent with lower gap, lower voltage headroom, and larger non-bulk recombination contribution at `Voc`.
4. The thermal drift does not rescue current matching at `65 C`; it slightly worsens it.

That is exactly the kind of result that matters for tandem engineering: not just whether each isolated cell is good, but whether their mismatch becomes easier or harder to manage under realistic operating temperatures.

## Important Caveats

- The `4T` numbers are equivalent aggregated metrics, not a unique measured four-terminal JV curve.
- The thermal baseline is `300 K`, not exact STC.
- The module-yield script linearly extrapolates the `300 K` thermal baseline to exact STC `25 C`.
- The module-yield script uses generic Faiman defaults (`U0 = 25`, `U1 = 6.84`) rather than a module-specific thermal fit.
- PVGIS supplies wind at `10 m`; the module-yield script converts that to `2 m` with a `0.14` power-law exponent before applying the Faiman model.
- Annual average `Voc` over all `8760` hours includes night-time zeros; the generated summary CSV also reports sunlit-only `Voc` averages and minima.
- The optical tandem EQE is a constructed optical sum using the filtered bottom response you requested.
- If the constructed `2T` JV does not cross `V = 0` over the overlapping current range, the code transparently falls back to `Jsc_2T = min(Jsc_top, Jsc_bottom)`.
- SCAPS sign convention is handled internally with `P = -VJ` for generated power.
