# PVSK/CdTe Tandem Analysis (SCAPS Exports)

This project contains a Python workflow to:

1. Parse SCAPS `.iv` and `.qe` text exports.
2. Plot PVSK top-cell JV, CdTe unfiltered JV, CdTe filtered JV, and a constructed 2T tandem JV.
3. Extract `Voc`, `Jsc`, `FF`, `Vmpp`, `Jmpp`, and `efficiency` for all requested devices.
4. Plot PVSK EQE, filtered CdTe EQE, and constructed tandem optical EQE.
5. Compute equivalent 4T tandem metrics and generate dedicated 4T power plots.

The main script is: `tandem_analysis.py`

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

Default expected files in project root:

- `PVSK_louis.iv`
- `CdTe_filtered_louis.iv`
- `CdTe_louis.iv`
- `PVSK_louis.qe`
- `CdTe_filtered_louis.qe`

SCAPS headers are automatically ignored; first two numeric columns are used.

## Usage

Run with defaults:

```bash
python tandem_analysis.py
```

Set EQE wavelength window explicitly:

```bash
python tandem_analysis.py --eqe-min-nm 250 --eqe-max-nm 900
```

Override unfiltered CdTe file path:

```bash
python tandem_analysis.py --cdte-unfiltered-iv CdTe_louis.iv
```

Choose another output folder:

```bash
python tandem_analysis.py --outdir results
```

## Outputs

The script writes:

- `tandem_jv_curve.png`
- `tandem_eqe_curve.png`
- `tandem_4t_power_curves.png`
- `tandem_4t_power_map.png`
- `tandem_metrics.csv`
- `tandem_jv_curve.csv`
- `tandem_eqe_curve.csv`

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
