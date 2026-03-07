# Pillar 5 Standardized Data Matrix

Note: the thermal sweep baseline is `300.00 K = 26.85 C`, not exact STC `25.00 C = 298.15 K`.

## PVSK Top Cell

| Metric | Value @ 300 K (26.85 C) | Value @ 65.00 C (338.15 K) | Temp. Coeff. (%/K) | Delta Drift |
| --- | ---: | ---: | ---: | ---: |
| Voc (V) | 1.2854 | 1.2336 | -0.1068 | — |
| Jsc (mA/cm2) | 20.5430 | 20.5494 | 0.0007 | — |
| FF (%) | 81.7884 | 80.7050 | -0.0342 | — |
| PCE (%) | 21.5970 | 20.4580 | -0.1385 | — |

## 2T Tandem

| Metric | Value @ 300 K (26.85 C) | Value @ 65.00 C (338.15 K) | Temp. Coeff. (%/K) | Delta Drift |
| --- | ---: | ---: | ---: | ---: |
| Voc (V) | 2.0769 | 1.9200 | -0.1990 | — |
| Jsc (mA/cm2) | 5.5714 | 5.5632 | -0.0039 | — |
| FF (%) | 81.0438 | 80.4406 | -0.0212 | — |
| PCE (%) | 9.3777 | 8.5922 | -0.2217 | — |
| Mismatch Ratio (Jtop/Jbottom) (ratio) | 3.8021 | 3.8195 | — | 0.0174 (0.457%, widens) |

Mismatch ratio is derived from the filtered top-cell and bottom-cell `Jsc` values, so it is the same design constraint for both `2T` and equivalent `4T` reporting.

## 4T Tandem (equivalent)

| Metric | Value @ 300 K (26.85 C) | Value @ 65.00 C (338.15 K) | Temp. Coeff. (%/K) | Delta Drift |
| --- | ---: | ---: | ---: | ---: |
| Voc (V) | 2.0769 | 1.9200 | -0.1990 | — |
| Jsc (mA/cm2) | 25.9460 | 25.9295 | -0.0018 | — |
| FF (%) | 45.9773 | 46.4264 | 0.0270 | — |
| PCE (%) | 24.7758 | 23.1132 | -0.1763 | — |
| Mismatch Ratio (Jtop/Jbottom) (ratio) | 3.8021 | 3.8195 | — | 0.0174 (0.457%, widens) |

Mismatch ratio is derived from the filtered top-cell and bottom-cell `Jsc` values, so it is the same design constraint for both `2T` and equivalent `4T` reporting.

