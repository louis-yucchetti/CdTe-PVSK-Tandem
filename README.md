# PVSK/CdTe Tandem Analysis (SCAPS Exports)

## What This Project Is For

This repository supports the final tandem-solar-cell project by turning SCAPS exports into plots, metrics, and comparisons that can be used in the oral defense and in the non-graded scientific manuscript.

The main script is `scripts/tandem_analysis.py`.

## Project Layout

```text
.
|-- data/
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
|   `-- figures/
`-- scripts/
    `-- tandem_analysis.py
```

`Solar Cells Modelling kits` is intentionally left untouched. Any folder named `outdated` is also intentionally left in place.

## What We Actually Have To Deliver

### 1. Graded Component: Oral Defense Only

- The project grade is based **only** on the oral defense on **March 11, 2026**.
- The oral defense counts for **100% of the project grade**, and the project itself is **40% of the course grade**.
- Format: **15-minute presentation + 10-minute Q&A**.
- There is **no graded written report**.

### 2. Non-Graded Component: Scientific Manuscript

- The report is **not** a course deliverable for grading.
- It is meant to serve as material for a possible joint benchmarking paper.
- Because of that, the report still needs strong bibliographic grounding and scientific interpretation.

### 3. Main Priority Shift

- The focus is now **bibliographic depth** and **physical interpretation**.
- The goal is to explain **why** the device behaves as it does, not to spend time on complex curve fitting.
- We are being evaluated as **R&D engineers**, not as fitting specialists.
- A weak result with a strong physical explanation is more valuable than a perfect-looking unexplained simulation.

## What The Defense Must Cover

### Pillar 1. Bibliographic Analysis and Device Physics

- Compare the chosen tandem technology with current efficiency records from the literature.
- Explain the role of each layer: absorber, buffer, HTL, and ETL.
- Justify why those materials are used and what their limitations are.
- Present and justify the band diagrams at `0 V` and `Voc`.
- Explain carrier transport and interface barriers using bibliographic support.

### Pillar 2. Tandem Dynamics and EQE Analysis

- Plot the J-V curves of each sub-cell and compare them with the full tandem response.
- Use the **filtered bottom-cell** response when discussing tandem operation.
- Use EQE curves to show spectral sharing between the top and bottom cells.
- Explain any mismatch physically: bandgap misalignment, parasitic absorption, transport limits, or interface losses.
- Compare extracted `Jsc` values with literature.

### Pillar 3. Thermal Modeling and Mismatch Drift

- Ensure `Eg(T)`, `mu(T)`, and `Nc, Nv(T)` are active in SCAPS.
- Analyze how temperature affects the specific materials in each sub-cell.
- Explain why `Voc` may fall faster in one sub-cell than in the other.
- Study the drift at `65°C`: does the gap between `Jtop` and `Jbottom` shrink or widen?
- This temperature-dependent mismatch drift is one of the main scientific results expected for the paper.

### Pillar 4. Energy Yield Estimation (PVGIS)

- Use **PVGIS hourly data** for a high-irradiance site such as **Ouarzazate, Morocco**.
- Scale the result to a **60-cell module** using Python or Excel.
- Include thermal losses using extracted coefficients such as `beta` for `Voc` and `gamma` for `Pmax`.

### Pillar 5. Standardized Data Matrix

The following table is mandatory in both the report and the oral defense:

| Metric (Tandem) | Value @ 25°C (STC) | Value @ 65°C | Temp. Coeff. (% / K) |
| --- | --- | --- | --- |
| `Voc` (V) | TBD | TBD | TBD |
| `Jsc` (mA/cm^2) | TBD | TBD | TBD |
| `FF` (%) | TBD | TBD | TBD |
| `PCE` (%) | TBD | TBD | TBD |
| Mismatch Ratio (`Jtop / Jbot`) | TBD | TBD | TBD |

Mismatch drift should be discussed separately as the change in `Jtop / Jbot` between `25°C` and `65°C`; it is not a literal entry in the temperature-coefficient column.

## Workload Relief From The Professor

- **No curve fitting:** extract parameters directly from SCAPS and explain their physical meaning.
- **NGSPICE can be skipped:** simple Python or Excel scaling is enough for the yield study.
- **Phase 4 is theoretical only:** about `0.5-1 page`, based on bibliography and Phase 2 results.

## What This Repository Produces

This project contains a Python workflow to:

1. Parse SCAPS `.iv` and `.qe` text exports.
2. Plot PVSK top-cell JV, CdTe unfiltered JV, CdTe filtered JV, and a constructed 2T tandem JV.
3. Extract `Voc`, `Jsc`, `FF`, `Vmpp`, `Jmpp`, and `efficiency` for all requested devices.
4. Plot PVSK EQE, filtered CdTe EQE, and constructed tandem optical EQE.
5. Compute equivalent 4T tandem metrics and generate dedicated 4T power plots.

## How To Use This Repo For The Project

- Use the generated J-V and EQE plots to support the oral defense.
- Use the extracted metrics to populate the standardized data matrix at `25°C` and `65°C`.
- Tie every figure back to literature, layer selection, interface physics, spectral sharing, and thermal behavior.
- Keep the discussion centered on **physical interpretation**, not on numerical fitting.

## Physics Model Used

### 1) JV construction for a 2T series tandem

The tandem is constructed in the **current domain**:

- Interpolate each subcell voltage as a function of current: `V_top(J)`, `V_bottom(J)`.
- Sum voltages at the same current:
  - `V_tandem(J) = V_top(J) + V_bottom(J)`

This is the correct series-connection rule for a 2-terminal tandem.

### 2) PV metrics extraction

From each JV curve:

- `Voc`: voltage at `J = 0`
- `Jsc`: absolute current magnitude at `V = 0`
- `P(V) = -V * J` (SCAPS sign convention)
- `Pmax`: maximum generated power density
- `FF = 100 * Pmax / (Voc * Jsc)`
- `eta = 100 * Pmax / Pin`, with default `Pin = 100 mW/cm^2`

### 3) Tandem optical EQE definition

This script constructs the optical tandem EQE as:

- `EQE_tandem_optical(lambda) = clip(EQE_top(lambda) + EQE_bottom_filtered(lambda), 0, 100)`

where bottom EQE is already the filtered response from your CdTe file.
The EQE figure includes all three curves: top, filtered bottom, and tandem optical.

### 4) 4T tandem definition used in this script

In a true 4T tandem, the two subcells have independent terminals, so there is no unique single JV curve for the whole tandem.

This workflow therefore reports an **equivalent 4T aggregate**:

- `Pmax_4T = Pmax_top + Pmax_bottom_filtered`
- `eta_4T = eta_top + eta_bottom_filtered`

and also reports equivalent `Voc`, `Jsc`, `FF`, `Vmpp`, `Jmpp` values for table convenience:

- `Voc_4T(eq) = Voc_top + Voc_bottom_filtered`
- `Jsc_4T(eq) = Jsc_top + Jsc_bottom_filtered`
- `Vmpp_4T(eq) = Vmpp_top + Vmpp_bottom_filtered`
- `Jmpp_4T(eq) = Pmax_4T / Vmpp_4T(eq)`
- `FF_4T(eq) = 100 * Pmax_4T / (Voc_4T(eq) * Jsc_4T(eq))`

These equivalent values are useful for comparative reporting, but they do not represent a unique measured 4T terminal JV.

## Why the Constructed Tandem JV Can Start Above 0 V

You may see the tandem JV start at a positive voltage (for example ~`0.26 V`) even when the x-axis starts at `0 V`.

Reason:

1. The series tandem is constructed only over the **shared current range** of both subcell JV datasets.
2. At the most negative common current, the bottom cell may be limited by the reverse-bias endpoint in your data.
3. If `V_top(J)` is still strongly positive there, `V_top(J) + V_bottom(J)` remains positive.
4. Then the constructed tandem curve never reaches `V = 0` inside the available overlap.

In that case, the script reports a transparent fallback for tandem `Jsc`:

- `Jsc_tandem = min(Jsc_top, Jsc_bottom)`

To force direct extraction of tandem `Jsc` from the constructed curve, extend reverse-bias sweep (usually bottom cell, often both cells) until the summed curve crosses `V = 0`.

## Plot Formatting (Publication-Oriented)

The script applies a journal-style theme:

- Serif font stack (`Times New Roman` fallback to `Times`, `DejaVu Serif`)
- Inward ticks on all sides
- Minor ticks enabled
- Clean light dotted major grid
- High-resolution export (`600 dpi`)
- Tight margins

## Inputs

Default expected files:

- `data/scaps_exports/iv/PVSK_louis.iv`
- `data/scaps_exports/iv/CdTe_filtered_louis.iv`
- `data/scaps_exports/iv/CdTe_louis.iv`
- `data/scaps_exports/qe/PVSK_louis.qe`
- `data/scaps_exports/qe/CdTe_filtered_louis.qe`
- `data/scaps_exports/qe/CdTe_louis.qe`

SCAPS headers are automatically ignored; first two numeric columns are used.

## Usage

Run with defaults:

```bash
python scripts/tandem_analysis.py
```

Set EQE wavelength window explicitly:

```bash
python scripts/tandem_analysis.py --eqe-min-nm 250 --eqe-max-nm 900
```

Override unfiltered CdTe file path:

```bash
python scripts/tandem_analysis.py --cdte-unfiltered-iv data/scaps_exports/iv/CdTe_louis.iv
```

Choose another output folder:

```bash
python scripts/tandem_analysis.py --outdir results
```

## Outputs

The script writes into `results/figures/` and `results/csv/` by default:

- `results/figures/tandem_jv_curve.png`
- `results/figures/tandem_eqe_curve.png`
- `results/figures/tandem_4t_power_curves.png`
- `results/figures/tandem_4t_power_map.png`
- `results/csv/tandem_metrics.csv`
- `results/csv/tandem_jv_curve.csv`
- `results/csv/tandem_eqe_curve.csv`

### Output CSV details

- `tandem_metrics.csv`: extracted metrics for `CdTe_unfiltered`, `CdTe_filtered`, `PVSK_top`, `Tandem_2T`, and `Tandem_4T` with a `model` column that states the metric definition.
- `tandem_jv_curve.csv`: constructed tandem JV points (`Voltage_V`, `CurrentDensity_mAcm2`).
- `tandem_eqe_curve.csv`: top EQE, filtered bottom EQE, and tandem optical EQE in selected wavelength window.

## Practical Notes

1. SCAPS current sign is handled internally (power generated uses `-V*J`).
2. If `Jsc` appears as fallback, check printed JV span diagnostics.
3. For strict 2T electrical spectral response studies, monochromatic EQE treatment depends on measurement bias-light protocol; this script intentionally reports the optical-sum construction you requested.
4. 4T plots are power-based and physically appropriate for independent terminals:
   - `tandem_4t_power_curves.png`: top/bottom power curves + total 4T `Pmax`.
   - `tandem_4t_power_map.png`: `P_total(V_top, V_bottom_filtered)` map with 4T MPP point.
