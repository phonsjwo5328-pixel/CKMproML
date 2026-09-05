# ============================================================
# 01_nhanes_external_validation.R
# Metrics-layer reproduction of the NHANES external validation.
#
# Input : data_public/nhanes_predictions.csv  (SEQN, ckm_y, pred_prob; n = 2,983)
# Output: all external-validation metrics of manuscript Table 4 /
#         Supplementary Table S8 (AUC + CI, PR-AUC, Brier,
#         logistic-recalibration slope/intercept, decile calibration,
#         0.5-threshold classification metrics).
#
# Expected agreement (see README): identical to the paper to 4 dp,
# because the predicted probabilities themselves are provided and
# every downstream statistic is a deterministic function of them.
#
# Model provenance: final Weka Logistic fit (RWeka, C = FALSE,
# R = 1e-8, M = -1) locked after development on CHARLS training data
# (LASSO-selected 12 predictors, kknn k = 7 imputation + one-hot
# pipeline, seed 123); applied here WITHOUT refitting or recalibration.
# ============================================================
suppressPackageStartupMessages({
  library(pROC)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
here <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
pred_file <- file.path(here, "data_public", "nhanes_predictions.csv")
if (!file.exists(pred_file))
  pred_file <- file.path(here, "..", "data_public", "nhanes_predictions.csv")
stopifnot(file.exists(pred_file))

d <- read.csv(pred_file)
stopifnot(all(c("SEQN", "ckm_y", "pred_prob") %in% names(d)),
          nrow(d) == 2983, sum(d$ckm_y) == 777)

pe <- d$pred_prob
y  <- d$ckm_y

cat(sprintf("n = %d | events = %d (%.1f%%)\n", length(y), sum(y), 100 * mean(y)))

## ---- discrimination ----
r   <- roc(y, pe, quiet = TRUE, direction = "<")
ci  <- as.numeric(ci.auc(r))
auc <- as.numeric(auc(r))
cat(sprintf("AUC   = %.4f (95%% CI %.4f - %.4f)\n", auc, ci[1], ci[3]))

# PR-AUC (average precision, exact stepwise sum; no package needed)
ord <- order(pe, decreasing = TRUE)
yy  <- y[ord]
tp   <- cumsum(yy); fp <- cumsum(1 - yy)
rec  <- tp / sum(yy)
prec <- tp / (tp + fp)
prauc <- sum(diff(c(0, rec)) * prec)
cat(sprintf("PR-AUC = %.4f\n", prauc))

## ---- calibration ----
brier <- mean((pe - y)^2)
cat(sprintf("Brier = %.4f\n", brier))

lp <- log(pe / (1 - pe))
calmod <- glm(y ~ lp, family = binomial())
cat(sprintf("Calibration: intercept = %.4f, slope = %.4f\n",
            unname(coef(calmod)[1]), unname(coef(calmod)[2])))

q    <- quantile(pe, seq(0, 1, 0.1))
bins <- cut(pe, breaks = unique(q), include.lowest = TRUE)
cal  <- aggregate(list(pred = pe, obs = y), list(decile = bins), mean)
cat("\nDecile calibration:\n"); print(cal, digits = 3)

## ---- classification at prespecified 0.5 threshold ----
cls <- ifelse(pe >= 0.5, 1, 0)
tab <- table(y = y, pred = cls)
sens <- tab["1", "1"] / sum(tab["1", ])
spec <- tab["0", "0"] / sum(tab["0", ])
ppv  <- tab["1", "1"] / sum(tab[, "1"])
npv  <- tab["0", "0"] / sum(tab[, "0"])
acc  <- sum(diag(tab)) / sum(tab)
dor  <- (tab["1", "1"] * tab["0", "0"]) / (tab["0", "1"] * tab["1", "0"])
cat(sprintf("\nAt threshold 0.5: sens = %.4f  spec = %.4f  ppv = %.4f  npv = %.4f  acc = %.4f  dor = %.3f\n",
            sens, spec, ppv, npv, acc, dor))

cat("\n--- Cross-check against paper (Table 4 / Supp Table S8) ---\n")
expected <- c(AUC = 0.7357, Brier = 0.2271, slope = 0.4002, intercept = -1.0245,
              sens = 0.6873, spec = 0.6714, ppv = 0.4241, npv = 0.8590, acc = 0.6755)
got <- c(AUC = auc, Brier = brier, slope = unname(coef(calmod)[2]),
         intercept = unname(coef(calmod)[1]),
         sens = sens, spec = spec, ppv = ppv, npv = npv, acc = acc)
chk <- data.frame(metric = names(expected), paper = expected,
                  reproduced = round(got[names(expected)], 4))
chk$match <- abs(chk$paper - chk$reproduced) < 5e-4
print(chk, row.names = FALSE)
cat(sprintf("\n%d/%d metrics agree with the paper to 4 dp.\n", sum(chk$match), nrow(chk)))
