# Glacier Area Change in the Ala-Archa Catchment, Kyrgyzstan

### A Landsat-based assessment with explicit observation-quality screening

**Author:** Aydina Rysbaeva  
**Program:** Lumiere Research Scholar Program

---

## Overview

This repository contains the reproducible Google Earth Engine workflow developed for a study of glacier-area change in the Ala-Archa catchment, northern Tien Shan, Kyrgyzstan.

The study uses multi-decadal Landsat observations, the Randolph Glacier Inventory (RGI) Version 7, published glacier-area benchmarks, and ERA5 temperature data to investigate both long-term glacier change and the reliability of annual glacier-area estimates derived from the optical satellite archive.

A central objective of the project is to distinguish between **satellite-derived mapped area** and **physically interpretable glacier area**. Annual Landsat estimates are therefore evaluated using observation-quality criteria defined independently of the resulting mapped glacier area.

---

## Research Question

**To what extent can glacier-area change in the Ala-Archa catchment, northern Tien Shan, be quantified between 1994 and 2025 using satellite observations, and to what extent is annual glacier-area change associated with summer air temperature between 2002 and 2025?**

---

## Study Area

The analysis uses a reproducible HydroSHEDS Level 10 catchment:

- **HYBAS_ID:** 4101277960
- **Catchment area:** 183.18 km²
- **Elevation range:** 2,183–4,799 m a.s.l.
- **RGI Version 7 glacier outlines:** 95
- **RGI glacier footprint within the catchment:** approximately 37.60 km²

The catchment represents the glacierized headwaters of the Ala-Archa river system in the northern Tien Shan of Kyrgyzstan.

---

## Data

The workflow integrates:

- Landsat Collection 2 Level-2 Surface Reflectance imagery
- Landsat 5 TM
- Landsat 7 ETM+
- Landsat 8 OLI
- Landsat 9 OLI-2
- Randolph Glacier Inventory (RGI) Version 7
- HydroSHEDS Level 10 catchment boundaries
- ERA5 2 m air temperature
- Independent published glacier-area benchmarks

---

## Glacier Mapping Workflow

Glacier-like snow and clean-ice surfaces are identified using the Normalized Difference Snow Index (NDSI):

**NDSI = (Green - SWIR1) / (Green + SWIR1)**

The principal processing parameters are:

| Parameter | Specification |
|---|---|
| Seasonal window | 15 August–15 September |
| Scene cloud threshold | <30% |
| Spatial resolution | 30 m |
| Primary NDSI threshold | 0.40 |
| Surface reflectance scaling | DN × 0.0000275 − 0.2 |
| Minimum diagnostic observations | ≥2 valid observations |
| Climate period | June–August |
| Climate variable | ERA5 2 m air temperature |

Pixel-level quality masking and sensor-specific processing are applied before annual NDSI composites are generated.

---

## Observation-Quality Framework

A major component of the study is an observation-quality screen that is independent of mapped glacier area.

For each year, the workflow evaluates:

1. the number of suitable Landsat scenes;
2. the number of valid observations received by each pixel;
3. the percentage of the independent RGI glacier footprint receiving sufficient valid observations.

The primary year-level quality rule requires:

**≥3 suitable scenes AND ≥90% of the RGI footprint receiving ≥2 valid observations.**

This rule prevents apparently plausible glacier-area values from being used to determine whether the same observations are considered reliable.

Sensitivity analyses additionally evaluate:

- ≥2, ≥3, and ≥5 valid observations per pixel;
- RGI coverage thresholds between 50% and 90%;
- alternative NDSI thresholds;
- alternative temporal observation windows.

---

## Key Findings

Independent published observations indicate substantial long-term glacier recession in the Ala-Archa region.

Published glacier-area benchmarks include:

| Year | Glacier area |
|---|---:|
| 1994 | 35.8 ± 3.2 km² |
| 2003 | 34.6 ± 1.7 km² |
| 2010 | 33.4 ± 0.8 km² |
| 2021 | 31.45 km² |

The difference between the 1994 and 2021 published estimates is **4.35 km²**, equivalent to approximately **12.15%** of the 1994 estimate.

Because these estimates originate from different published mapping products and spatial definitions, this change is treated as an independent benchmark comparison rather than a fitted homogeneous annual trend.

---

## Observation Completeness

The Landsat analysis demonstrates that nominal satellite availability does not necessarily imply complete observation of glacierized terrain.

Under the primary observation-quality rule, only **2022** satisfies the required criteria.

Using an 80% RGI-coverage sensitivity threshold retains:

- 2022
- 2024
- 2025

These years are not consecutive and therefore cannot provide a defensible annual glacier-area change series.

A 2020 observation-window experiment further demonstrates the importance of observation completeness. Using the full 15 August–15 September window produced a diagnostic mapped area of approximately **25.26 km²**, whereas restricting the window to 1–15 September reduced the mapped area to approximately **7.83 km²**.

This difference is far too large to represent physical glacier loss over such a short interval and demonstrates the sensitivity of the mapped area to observation availability.

---

## Climate Analysis

The study evaluates annual glacier-area change against June–August ERA5 2 m air temperature.

The pre-screening diagnostic produced:

- **Pearson r = −0.0518**
- **p = 0.8143**
- **Spearman ρ = +0.0455**

These statistics are retained for transparency but are **not interpreted as physical climate relationships**.

After application of the independent observation-quality screen, there are no consecutive annual observations suitable for calculating screened annual glacier-area change.

The final screened annual temperature association is therefore **not estimable**.

---

## Main Interpretation

The study supports two distinct conclusions:

1. Independent observations provide evidence of substantial multi-decadal glacier recession in the Ala-Archa region.

2. The available Landsat optical archive does not support a defensible continuous annual glacier-area trajectory through 2025 under the declared observation-quality framework.

The principal methodological result is that a long satellite archive should not automatically be interpreted as a long physical annual record. Observation completeness over the glacierized footprint must first be evaluated independently of the mapped glacier area.

---

## Repository Structure

```text
ala-archa-glacier-research/
│
├── code/
│   └── Ala_Archa_Glacier_Analysis_FINAL.js
│
└── README.md
