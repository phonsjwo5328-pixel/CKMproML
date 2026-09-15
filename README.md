# CKMproML — Reproducibility Package

Reproducibility materials for:

> **Development of an Interpretable Machine Learning Framework for Cardiovascular–Kidney–Metabolic Syndrome Progression: SHAP-Guided Prediction with Phenotype Stratification**
> (submitted to *Diabetes, Metabolic Syndrome and Obesity: Targets and Therapy*, Dove Medical Press)

The pipeline derives a 12-predictor logistic regression model for 4-year progression
from CKM stages 0–2 to 3–4 in CHARLS 2011–2015 (n = 3,874), then applies the **locked
model without refitting or recalibration** to NHANES 2017–March 2020 (n = 2,983 US
adults) as external validation.

---

## Repository structure

```
R/
  01_nhanes_external_validation.R   # Full external-validation analysis (metrics layer)
python/
  03_build_validation.py           # NHANES 2017-2020.3 -> analysis-ready dataset (from CDC XPT)
data_public/
  nhanes_predictions.csv            # SEQN + outcome + predicted probability (n = 2,983)
  nhanes_calibration_deciles.csv    # Decile calibration table
  platform_metrics/                 # Archived platform metric outputs (incl. an exploratory all-algorithm external comparison not reported in the paper)
LICENSE (MIT for code)
```

## What can be reproduced here, and to what tolerance

This repository is designed around a **two-layer reproducibility contract**. Please read
this table before comparing your output against the paper.

| Layer | What you need | Expected agreement |
|---|---|---|
| **Metrics layer (fully reproducible here)** — external-validation AUROC/PR-AUC/Brier, calibration slope and calibration-in-the-large, decile table, 0.5-threshold classification metrics | Only the files in `data_public/` + R ≥ 4.x with `pROC` | **Exact to machine precision** (≤ 1e-12 in predicted probabilities; identical metrics to 4 dp). The final model is a deterministic Weka Logistic fit; predicted probabilities are provided, so every reported external-validation number in Table 4 / Supplementary Table S8 is recomputable from `nhanes_predictions.csv` alone. |
| **Modeling layer (data not redistributable)** — refitting the final model, internal test-set predictions, 16-algorithm benchmark | CHARLS microdata (application required) + the FreeStatistics/Weka environment | The final-model numbers (internal AUC 0.7508 etc.) are deterministic given the same pipeline inputs and seed (123), and were verified to machine precision against archived platform predictions. Cross-validated benchmark entries (Supplementary Table S6) may vary by **±0.003** across package versions/platforms; the model-ranking conclusions are unaffected. |

Note on calibration: the paper reports the external calibration slope (0.40)
and the calibration-in-the-large (−1.37, the intercept with the linear predictor fixed
as an offset). The script prints both that quantity and the 2-parameter logistic-recalibration
intercept (alpha = −1.0245), which is not reported in the paper, so that either convention
can be verified.

If your reproduction differs from the paper, check which layer you are in. Differences
inside the stated tolerances are expected environment behavior, not errors.

## Quick start (metrics layer, ~1 minute)

```r
install.packages("pROC")   # once
source("R/01_nhanes_external_validation.R")
```

Expected console output (compare with paper Table 4 / Supplementary Table S8):

```
n = 2983 | events = 777 (26.0%)
AUC   = 0.7357 (95% CI 0.7150 - 0.7564)   # paper: 0.736 (0.715-0.756)
PR-AUC = 0.5180 (stepwise average precision; platform estimator 0.5175; both round to 0.518 in Figure 11B)
Brier = 0.2271
Calibration: slope = 0.4002 | recalibration intercept (alpha) = -1.0245 | calibration-in-the-large = -1.3746   # paper: slope 0.40, intercept -1.37 (calibration-in-the-large)
At threshold 0.5: sens = 0.6873  spec = 0.6714  ppv = 0.4241  npv = 0.8590  acc = 0.6755
```

All reported external-validation point estimates (Table 4 / Supplementary Table S8)
are recomputed from `data_public/nhanes_predictions.csv` alone; the script ends with
an automated cross-check (9/9 reported metrics agree with the paper to the reported precision).

## Rebuilding the NHANES dataset from CDC source files

1. Download the P-cycle (2017–March 2020) XPT files listed in
   `python/03_build_validation.py` from
   <https://wwwn.cdc.gov/nchs/nhanes/ContinuousNhanes/Default.aspx?BeginYear=2017>
2. Place them in a `raw/` folder and point `RAW` in the script to it.
3. `python python/03_build_validation.py` regenerates the analysis dataset
   (n = 2,983; 777 events, 26.0%) — this is a deterministic transformation of public
   CDC files, so your output should match `data_public/` bit-for-bit apart from
   row-ordering of identical SEQNs (none expected).

## CHARLS microdata (not included, by design)

CHARLS user agreement **prohibits redistribution** of microdata. To reproduce the
modeling layer:

1. Register at <https://charls.pku.edu.cn/en/> and download Harmonized/Waves 1 & 3.
2. Follow Supplementary Methods S1–S3 of the paper (staging criteria, cohort assembly
   17,708 → 12,434 → 6,075 → 4,041 → 3,927 → 3,875 → 3,874; MICE imputation m = 20,
   seed 123; LASSO λ1se = 0.0201; 12 predictors).
3. Final classifier: Weka `Logistic` (RWeka), default hyperparameters (C = FALSE,
   R = 1e-8, M = −1), with kknn imputation (k = 7) + one-hot encoding in the pipeline.

The NHANES predictor operationalization (including the composite bodily-pain mapping)
is fully specified in `python/03_build_validation.py` and Supplementary Table S7.

## Software environment

| Component | Version used in the paper |
|---|---|
| R | 4.5.3 |
| RWeka | (package bundled Weka 3) |
| pROC | 1.18.5 |
| Python (dataset build) | 3.x with pandas, numpy |
| Development platform | FreeStatistics (FengRui) for model development/benchmarking |

## Data availability statements

- **NHANES**: public-use data, U.S. CDC/NCHS — no redistribution restrictions; derived
  analysis dataset and per-participant predictions are included in `data_public/`.
- **CHARLS**: public upon registration; **not** redistributable under its data-use
  agreement.

## License

Code: MIT (see `LICENSE`). No claim is made over CDC/NHANES data (public domain) or
CHARLS data (governed by its own agreement).
