# Wafer-scale Mach-Zehnder Modulator Analysis

Automated extraction and statistical analysis of Mach-Zehnder modulator (MZM)
parameters from wafer-level optical and electrical measurements.  Targets the
HY202103 device set but the methodology generalizes to any depletion-mode
silicon MZM data delivered in the same XML schema.

---

## 1. Introduction

### 1.1 Purpose

Wafer-level silicon photonics fabrication produces hundreds of nominally
identical devices, but process variation, measurement instrumentation drift,
and probe-station contact quality all introduce per-die scatter.  Engineers
need to:

1. Extract operationally meaningful parameters (extinction ratio,
   insertion loss, half-wave voltage) for every die.
2. Distinguish broken or unreliable measurements from real device variation.
3. Decompose remaining variation into systematic spatial trends (which a
   process engineer can act on) and random scatter (which requires
   redundancy strategies).

This codebase performs all three tasks in a single reproducible pipeline.

### 1.2 Scope

Input is a folder of XML files (HY202103 schema), each containing one die's
optical spectra (six bias points) and current-voltage characteristic.

Output is a per-die table (CSV and XLSX), per-metric wafer maps, 1D
distributions, a Robust Z-score map, a per-date breakdown CSV with a
red-highlighted Excel version, and a systematic/random variation
decomposition figure.

### 1.3 Pipeline overview

For each die the pipeline computes three primary metrics (ER, IL, V_pi)
and several derived parameters (FSR, dlambda/dV, linearity R^2, coupler
imbalance, MZM section loss).  Per-die rows are then aggregated, outliers
are flagged using physical bounds and a per-wafer-per-band Robust Z-score,
and the resulting tables and figures are written to disk and pushed to
the GitHub remote.

```
XML files
  -> per-die extraction (ER, IL, V_pi, passive params)
  -> outlier flagging (physical bounds + Robust Z)
  -> aggregation
  -> wafer maps / 1D plots / Z-score maps / decomposition / by-date
  -> CSV + XLSX + PNG written to res/
  -> auto-commit and push to GitHub
```

---

## 2. Theory

### 2.1 Mach-Zehnder modulator operation

A MZM splits light into two arms, applies a voltage-controlled phase shift
in one arm, and recombines the two paths.  The output transmission is

```
T(lambda, V) = |a * exp(j*phi1) + b * exp(j*phi2(V))|^2
             = a^2 + b^2 + 2*a*b*cos(Delta_phi(V))
```

where `a` and `b` are the amplitude coefficients after the splitter and the
phase difference `Delta_phi` depends on the bias `V`.  Maximum transmission
occurs when the two arms are in phase, minimum when they are pi out of phase.

The wavelength-dependent transmission therefore shows a periodic null
pattern.  The period (`FSR`, free spectral range) and the rate at which
each null moves with bias (`dlambda/dV`) together fix the half-wave voltage:

```
V_pi = FSR / (2 * |dlambda/dV|)
```

### 2.2 Plasma dispersion

In a silicon depletion-mode MZM the phase modulation arises from the
plasma dispersion effect (Soref and Bennett, 1987):

```
Delta_n = -8.8e-22 * Delta_N_e - 8.5e-18 * (Delta_N_h)^0.8
```

where `Delta_N_e` and `Delta_N_h` are the changes in free electron and hole
density.  Reverse bias on a PN junction widens the depletion region, which
*removes* carriers from the optical mode and *increases* the effective
index.  Forward bias *injects* carriers, *decreases* the effective index.

### 2.3 Why reverse bias is the operating regime

Commercial high-speed silicon MZMs operate exclusively in reverse bias.
The comparison is summarized in Table 1.

**Table 1.** Reverse (depletion) versus forward (injection) bias regimes.

| Property | Reverse (depletion) | Forward (injection) |
|---|---|---|
| Mechanism | depletion-region width modulation | direct carrier injection |
| Response time | RC time constant (~5 ps) | minority-carrier lifetime (~1-100 ns) |
| Bandwidth | 30 GHz and above | 10-100 MHz |
| Steady-state power | sub-pJ/bit (capacitive only) | mW (continuous current) |
| Linearity in V | nearly linear (`W proportional to sqrt(V_bi-V)`) | exponential (Shockley equation) |
| Suitability for telecom | standard | unusable above ~100 MHz |

Reverse bias modulates the *spatial distribution* of carriers without
significant steady-state current.  This redistribution responds on a
capacitive timescale, giving the bandwidth required for modern
communication links.  Forward bias requires actual carrier transport,
which is slow and dissipates substantial power.

The extraction code therefore uses only reverse biases (V <= 0) for the
`dlambda/dV` fit; forward bias points are retained in the raw XML but are
not included in the slope estimation.

### 2.4 Linearity as a quality metric

A perfect depletion-mode MZM produces a strictly linear `Delta_lambda` versus
`V_bias` relationship.  Deviation from linearity indicates one of:

- Onset of forward conduction below the expected threshold (high series
  resistance or unintended forward operation).
- Significant carrier injection at moderate reverse bias (process defect).
- Measurement noise or null-tracking failure.

The coefficient of determination R^2 of the linear fit therefore acts as a
per-die signal-quality metric.  The pipeline reports the median R^2 over
all nulls successfully tracked for each die in the column `R2_dlam_vs_V`.

### 2.5 Coupler split ratio from extinction ratio

For an imbalanced two-arm interferometer with amplitude split ratio
`k = b/a` (k <= 1) the linear extinction ratio is

```
ER_linear = ((1 + k) / (1 - k))^2
```

so

```
k = (sqrt(ER_linear) - 1) / (sqrt(ER_linear) + 1)
```

This closes the loop between the directly measurable ER and the underlying
splitter imbalance.  The pipeline reports `amplitude_ratio_k`,
`power_split_ratio` (`k^2`), and `imbalance_dB` (`-20*log10(k)`) per die,
plus the imbalance evaluated at each individual bias point.

### 2.6 Outlier detection strategy

Two independent layers are combined:

1. **Physical bounds.**  Each metric is checked against a hard range
   derived from device physics and validated against the observed
   distribution (Section 5.1).  Out-of-bounds values represent extraction
   failures rather than real device variation and are flagged
   `is_problematic = True` in the by-date CSV.

2. **Robust Z-score.**  Within each (Wafer, Band) group the modified
   Z-score `z' = (x - median) / (1.4826 * MAD)` is computed.  Dies with
   `|z'| > 3` are flagged.  The Robust Z is independent of the physical
   bound; either trigger marks the die as an outlier.

The Z-score is grouped by (Wafer, Band) rather than by spatial
neighborhood because per-wafer sample size (14 dies) is too small for
neighborhood comparison to be stable.

### 2.7 Variation decomposition

Following Xing et al. (ACS Photonics 2023), the per-die metric is
decomposed into a smooth spatial component and a residual:

```
observed(R, C) = systematic(R, C) + random(R, C)
```

where `systematic` is a degree-2 polynomial in (Row, Col) fit by least
squares.  The ratio

```
R^2 = Var(systematic) / Var(observed)
```

indicates the fraction of die-to-die variation that is location-dependent.
A high R^2 motivates process-uniformity work; a low R^2 indicates that the
remaining variation is essentially random and requires statistical
mitigation (more dies, redundancy).

---

## 3. Methodology

### 3.1 Input schema

Each XML file describes a single (Wafer, DieRow, DieColumn, Band) entry.
The relevant fields are:

- `TestSiteInfo` attributes Wafer, DieRow, DieColumn.
- `Modulator` named `ALIGN...`: passive reference spectrum (wavelength L
  and transmitted power IL).
- `Modulator` named `MZM...`: device under test.  Provides up to six
  `WavelengthSweep` entries (one per bias point) plus one `IVMeasurement`.

Files with suffix `_DCM_LMZC` are C-band (lambda_c = 1550 nm); `_DCM_LMZO`
are O-band (lambda_c = 1310 nm).  Width is parsed from the modulator name
pattern `MZMCTE_LULAB_<width>_<length>`.

### 3.2 Per-die parameter extraction

#### 3.2.1 Extinction ratio (`ER_dB`)

For each bias point the ALIGN-referenced transmission is

```
T_dev(lambda) = IL_mzm(lambda) - IL_ref(lambda)
```

The extinction ratio is computed within a fixed wavelength window:

```
ER_dB = max over all biases of max_lambda T_dev   -
        min over all biases of min_lambda T_dev
```

Window boundaries are chosen to comfortably contain at least one full
FSR per band (Section 5.2).  Using the best peak and the deepest null
over the entire bias sweep gives a fixed-bias-independent ER value that
correlates with the device's worst-case modulation depth.

#### 3.2.2 Insertion loss (`IL_dB`)

`IL_dB` is the peak `IL_mzm` value at V = -1 V within a +- 5 nm window
around the band center.  No ALIGN subtraction is performed because the
XML's IL field already represents device-level transmission.

#### 3.2.3 Half-wave voltage (`Vpi_V`)

The `extract_vpi` module:

1. Identifies the V = -2 V spectrum nulls with `scipy.signal.find_peaks`.
   Nulls deeper than -25 dB are retained.
2. Computes FSR as the median spacing between deep nulls.
3. For each deep null and each reverse bias point, refines the null
   wavelength by a parabolic fit on the five nearest samples.
4. Linear fits the null wavelength versus bias to obtain `dlambda/dV`.
5. Linear fits with total bias-induced shift exceeding `1.5 * window`
   are rejected as null tracking failures.
6. Surviving slopes are averaged after a 3 * MAD outlier trim.
7. Applies the slope filter: if `|dlambda/dV| < MIN_SLOPE_PM_PER_V`
   (default 10 pm/V) the measurement is declared broken and `Vpi_V`
   is set to NaN.  The slope value itself is preserved for diagnostic
   purposes.

The output dictionary also carries an explicit status field
(`vpi_status`) so downstream code can distinguish among
`ok | slope_filter | no_nulls | no_slopes | few_biases | no_sweeps`.

#### 3.2.4 Linearity R^2 (`R2_dlam_vs_V`)

For each tracked null the R^2 of the linear fit is computed.  After the
3 * MAD outlier trim used for slope averaging, the median R^2 of the
surviving fits is reported.  This is the primary signal-quality metric
for the V_pi extraction.

#### 3.2.5 Quality grade (`quality_grade`)

A letter grade is assigned from `R2_dlam_vs_V` and `vpi_status`:

| Grade | Condition |
|---|---|
| A | R^2 >= 0.99 |
| B | 0.95 <= R^2 < 0.99 |
| C | 0.90 <= R^2 < 0.95 |
| D | 0.50 <= R^2 < 0.90 |
| F | R^2 < 0.50 or vpi_status != 'ok' |

Grades A-B correspond to dies suitable for high-fidelity signal
extraction; C-D indicate increasing measurement uncertainty; F dies
should not be used in downstream statistics.

#### 3.2.6 Passive parameters

From the measured ER the splitter amplitude ratio is computed via the
closed-form expression of Section 2.5.  The MZM-section propagation loss
is the negative of the maximum T_dev over all biases within the ER
window (the best transmission achieved at any operating point relative
to the ALIGN reference).

Per-bias splitter imbalance is reported in columns
`imbalance_V<V>_dB`, allowing the user to observe whether the splitter
imbalance varies with operating point (it should be approximately
constant for an ideal MMI).

### 3.3 Physical bounds

Bounds are derived from published reviews and validated against the
HY202103 distribution (Section 5.1):

| Metric | Lower | Upper | Source |
|---|---|---|---|
| ER_dB | 10 | 45 | Witzens 2018 plus empirical margin |
| IL_dB | -15 | -1 | Reed et al. 2010 |
| Vpi_V | 2 | 60 | V_pi*L = 1 to 3 V*cm, L = 0.5 to 5 mm |

A value outside its bound is flagged in `reason_<metric>` with the
detected violation (for example `over: 5478.34 > 80.0`).  The
`is_problematic` boolean is the OR of the three per-metric flags.

### 3.4 Robust Z-score

Per (Wafer, Band):

```
sigma_robust = 1.4826 * MAD(values within group)
z'           = (value - median) / sigma_robust
outlier      = |z'| > 3
```

For dies with NaN metric value the outlier flag is forced True.  The
constant 1.4826 makes `sigma_robust` equal to the standard deviation of
a Gaussian sample with the same MAD, so the threshold 3 corresponds to
a standard 3-sigma rule but is insensitive to the extreme values it is
meant to detect (Iglewicz and Hoaglin, 1993).

### 3.5 Variation decomposition

For each (Wafer, Band) group, the function
`decompose_variation.fit_systematic` fits a degree-2 bivariate
polynomial in (Row, Col).  Six free parameters, evaluated against
fourteen data points, leaves eight degrees of freedom: adequate for a
stable fit while permitting curvature without overfitting.

The decomposition figure shows three panels per group: the observed
field, the polynomial fit, and the residual, with R^2 annotated on the
fit panel.

---

## 4. Implementation

### 4.1 Folder layout

```
.
|-- run.py                   entry point
|-- requirements.txt
|-- src/
|   |-- xml_loader.py        XML parsing
|   |-- extract_er.py        ER per die
|   |-- extract_il.py        IL per die
|   |-- extract_vpi.py       V_pi, FSR, dlambda/dV, linearity R^2, grade
|   |-- extract_passive_params.py   coupler imbalance, MZM section loss
|   |-- outlier_detect.py    physical bounds + Robust Z
|   |-- csv_export.py        run-folder creation
|   |-- plot_common.py       shared color and ordering helpers
|   |-- wafer_map.py         continuous-surface wafer maps
|   |-- plot_1d.py           IQR box plots
|   |-- plot_1d_mad.py       MAD box plots
|   |-- zscore_map.py        Robust Z-score grid
|   |-- decompose_variation.py    systematic/random spatial decomposition
|   |-- analyze_by_date.py   per-measurement-date analysis with XLSX
|   |-- investigate.py       9-point diagnostic report
|   `-- sensitivity_test.py  ER-window and slope-filter sensitivity sweep
|-- data/                    (placeholder; HY202103 currently referenced
|                              externally via DATA_ROOT in run.py)
|-- doc/                     methodology figures and investigation outputs
`-- res/
    |-- csv/                 tracked CSV and XLSX outputs (latest run)
    |-- figures/             tracked PNG outputs (latest run)
    `-- <timestamp>/         every run's raw outputs (gitignored)
```

### 4.2 Module responsibilities

`run.py` orchestrates the pipeline.  It enumerates XML files via
`xml_loader.find_all_xmls` (which retains only the most recent
measurement for each `(Wafer, Row, Col, Band)`), distributes per-die
extraction across CPU cores using `multiprocessing.Pool`, applies
`outlier_detect.mark_outliers` to the aggregated DataFrame, writes the
per-run results, mirrors the latest run into `res/csv` and `res/figures`,
invokes `analyze_by_date.export_and_plot` for the per-date breakdown,
and finally commits and pushes the tracked outputs to the GitHub remote.

Plot generation is fanned out via `concurrent.futures.ProcessPoolExecutor`
across the four plot types and three metrics for a total of twelve
parallel tasks.

### 4.3 Execution

```
python3 run.py
```

Two standalone scripts are also available:

```
python3 src/investigate.py        # 9-point diagnostic
python3 src/sensitivity_test.py   # sweep ER window and slope filter
```

### 4.4 Output

Two CSVs are produced (along with their PNG counterparts) on every run:

`res/csv/data.csv` contains one row per deduplicated die.  Columns are
grouped as identifier, primary metrics, V_pi diagnostics, splitter
parameters, and outlier flags.

`res/csv/data_by_date.csv` (and `data_by_date.xlsx`) preserves all
measurement dates without deduplication.  The XLSX version applies a red
background fill to cells whose `reason_<metric>` is non-empty and to the
`is_problematic` column.

Per-run snapshots are stored under `res/<YYYY-MM-DD_HH-MM-SS>/`, which is
gitignored.  Only the latest snapshot is mirrored into `res/csv` and
`res/figures` and tracked by git.

### 4.5 CSV column dictionary

`data.csv`:

| Column | Description |
|---|---|
| Wafer, Band, Row, Col, Width_nm | identifier |
| ER_dB | extinction ratio (Section 3.2.1) |
| IL_dB | insertion loss (Section 3.2.2) |
| Vpi_V | half-wave voltage (Section 3.2.3) |
| FSR_nm | free spectral range from V=-2 V spectrum |
| dlam_dV_pm_per_V | mean reverse-bias null shift rate |
| R2_dlam_vs_V | median R^2 of the linear `Delta_lambda` versus V fit |
| quality_grade | A through F derived from R2_dlam_vs_V and vpi_status |
| vpi_status | one of ok / slope_filter / no_nulls / no_slopes / few_biases / no_sweeps |
| amplitude_ratio_k, power_split_ratio, imbalance_dB | derived splitter parameters |
| mzm_loss_dB | MZM-section propagation loss (positive dB) |
| imbalance_V<bias>_dB | per-bias splitter imbalance |
| is_outlier_ER_dB, is_outlier_IL_dB, is_outlier_Vpi_V | physical-bound or Robust-Z outlier |
| robust_z_ER_dB, robust_z_IL_dB, robust_z_Vpi_V | Robust Z-score values |
| is_trusted | none of the three outlier flags is set |

`data_by_date.csv` adds the `Date` column and replaces the trust columns
with `reason_<metric>`, `out_of_bound_<metric>`, and `is_problematic`.

---

## 5. Results (HY202103)

### 5.1 Empirical validation of physical bounds

The initial ER upper bound was set to 35 dB based on Witzens (2018)
Table II.  This rejected all 28 C-band dies because the observed C-band
ER median is 37 dB.  Two factors drive the difference:

1. The HY202103 splitter (MMI) is more balanced than the population
   average, yielding higher achievable ER.
2. The ER definition used here is the peak-to-null span over all biases
   (Section 3.2.1), which is naturally larger than a single-bias ER.

The upper bound was therefore relaxed to 45 dB.  This value retains all
real device data while still rejecting unphysical artifacts in the
2019-05-31 broken measurement set.

### 5.2 ER window selection

`src/sensitivity_test.py` sweeps the window width.  Results for C-band:

| Window | Median ER | sigma | Out-of-bound (n>45) |
|---|---|---|---|
| 14 nm | 36.87 | 1.47 | 0 |
| 16 nm | 37.11 | 1.24 | 0 |
| 18 nm | 37.31 | 1.40 | 0 |
| 22 nm | 38.03 | 1.36 | 0 |
| 36 nm | 39.53 | 2.48 | 2 |

A 16-nm window minimizes the standard deviation across dies (better
peak/null statistics with one full FSR comfortably inside the window)
without introducing artifacts from band-edge regions where the ALIGN
reference is less reliable.  The C-band window is therefore set to
[1545, 1561] nm.  The O-band FSR is ~9.8 nm so the original 14 nm
window already contains 1.4 FSR and was kept unchanged.

### 5.3 Diagnostic findings

The 2019-05-31 measurement set (28 O-band dies, D23 and D24 wafers)
exhibits `dlambda/dV` values within +-5 pm/V.  Normal measurements of
the same dies on 2019-06-03 yield `|dlambda/dV| > 100 pm/V`.  The slope
filter at 10 pm/V cleanly separates the two populations: all 28 broken
measurements are marked broken (`vpi_status = slope_filter`) while no
normal measurement is incorrectly rejected.

The same dies measured later on 2019-06-03 produced clean spectra and
the 06-03 measurement is what the deduplicated `data.csv` retains.

### 5.4 Linearity distribution

`R2_dlam_vs_V` per group (mean +- standard deviation):

| Group | Mean R^2 | Min R^2 |
|---|---|---|
| D08, C-band | 0.987 | 0.974 |
| D08, O-band | 0.989 | 0.978 |
| D07, C-band | 0.981 | 0.965 |
| D23, O-band | 0.977 | 0.932 |
| D24, O-band | 0.958 | 0.862 |

The C-band measurements show the cleanest linear depletion behavior.
The D24 O-band has a few dies with R^2 below 0.90 (grade C-D) that
extract successfully but should be down-weighted in pooled statistics.

### 5.5 Variation decomposition

Fraction of metric variance explained by location (degree-2 polynomial
in Row, Col):

| Group | R^2 (ER) | R^2 (IL) | R^2 (V_pi) |
|---|---|---|---|
| D08, O-band | 0.48 | 0.90 | 0.72 |
| D23, O-band | 0.64 | 0.83 | 0.85 |
| D24, O-band | 0.31 | 0.81 | 0.40 |
| D07, C-band | 0.38 | 0.55 | 0.10 |
| D08, C-band | 0.51 | 0.67 | 0.29 |

IL is the metric most strongly governed by spatial process variation
across all groups, consistent with thickness and width gradients
affecting waveguide loss.  V_pi is dominated by random variation in the
C-band wafers but is strongly position-dependent in the O-band wafers,
suggesting that the C-band fabrication run had better local control
of the phase-shifter geometry while the O-band run shows a measurable
wafer-scale gradient.

---

## 6. References

1. R. A. Soref and B. R. Bennett, "Electrooptical effects in silicon",
   *IEEE Journal of Quantum Electronics* **23**, 123-129 (1987).
2. G. T. Reed et al., "Silicon optical modulators", *Nature Photonics*
   **4**, 518-526 (2010).
3. A. H. Patel et al., "Design, analysis, and transmission system
   performance of a 41 GHz silicon photonic modulator", *Optics
   Express* **23**, 14263 (2015).
4. J. Witzens, "High-speed silicon photonics modulators",
   *Proceedings of the IEEE* **106**, 2158-2182 (2018).
5. S. K. Selvaraja et al., "Process variation in silicon photonic
   devices", *Applied Optics* **52**, 7638 (2013).
6. Y. Xing, J. Dong, U. Khan, and W. Bogaerts, "Capturing the effects of
   spatial process variations in silicon photonic circuits",
   *ACS Photonics* **10**, 928 (2023).
7. P. Xu et al., "Optical and geometric parameter extraction across
   300-mm photonic integrated circuit wafers", *APL Photonics* **9**,
   016104 (2024).
8. B. Iglewicz and D. Hoaglin, *How to detect and handle outliers*,
   ASQC Quality Press, 1993.
