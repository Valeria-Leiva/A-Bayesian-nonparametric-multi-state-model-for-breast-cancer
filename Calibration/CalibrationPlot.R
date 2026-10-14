# ==============================================================================
# CALIBRATION PLOTS ACROSS RISK STRATA AT 5 AND 10 YEARS (TRANSITION 5)
# METABRIC Cohort Analysis: BNP (H = 25) vs. BP (H = 1)
# ==============================================================================
# Note:
# 1. The data splits (e.g., Split1) are identical to those used across all other 
#    evaluations (such as C-index and Brier Score scripts). To evaluate another fold,
#    simply update the folder path/directory.
# 2. This script computes predicted individual survival probabilities and compares
#    them against observed Kaplan-Meier estimates across risk strata (quintiles).
# ==============================================================================

# Clear workspace
rm(list = ls())

# Required Libraries -----------------------------------------------------------
library(survival)
library(tidyverse)
library(splines)
library(ggplot2)

# 1. Load Data -----------------------------------------------------------------
# Note: Update paths according to local directory or split fold
load("Split1/id_res_up.RData")
load("Split1/resp_up_II.RData")

col_er <- grep("^ER", colnames(res), ignore.case = TRUE, value = TRUE)[1]

res2 <- res %>%
  mutate(
    time = Tstop - Tstart,
    trans2 = case_when(
      trans == 10 ~ 1, trans == 11 ~ 2, trans == 12 ~ 3,
      trans == 13 ~ 4, trans == 14 ~ 5, trans == 15 ~ 6,
      trans == 16 ~ 7, trans == 17 ~ 8, trans == 18 ~ 9,
      TRUE ~ as.numeric(as.character(trans))
    ),
    ER.status = ifelse(.data[[col_er]] == "ER+" | .data[[col_er]] == 1, 1, 0)
  )

cols_req <- c("id", "AGE", "LN_O", "Size", "Grade", "ER.status", 
              "hospital", "trans2", "Tstart", "Tstop", "time", "status")

res2.na  <- na.omit(res2[, cols_req])
res.test <- res2.na %>% filter(id %in% id_res_up$id.test)
test_tr2 <- res.test %>% filter(trans2 == 5)

# 2. Individual Survival Prediction Function at Target Time t ------------------
get_indiv_pred_surv <- function(file_rdata, test_data, tr_select = 5, target_time = 5) {
  load(file_rdata) # Loads out.DDP2
  out_obj <- out.DDP2$BUGSoutput$sims.list
  
  betah2 <- out_obj$betah    # [Nmc, p, tra, Hdp]
  bk2    <- out_obj$bk       # [Nmc, n_hosp]
  wh2    <- out_obj$BetaProb # [Nmc, Hdp, tra]
  sig2   <- out_obj$sigma2   # [Nmc, tra]
  
  Nmc <- dim(betah2)[1]
  Hdp <- dim(betah2)[4]
  
  data_j <- test_data %>% filter(trans2 == tr_select)
  X_j    <- model.matrix(~ AGE + LN_O + bs(Size, df = 3) + as.factor(Grade) + ER.status, data = data_j)
  
  beta_tr <- betah2[, , tr_select, , drop = FALSE]
  w_tr    <- if (Hdp == 1) matrix(1, nrow = Nmc, ncol = 1) else wh2[, , tr_select]
  sig_tr  <- sqrt(sig2[, tr_select])
  
  pred_surv <- numeric(nrow(data_j))
  
  for (i in 1:nrow(data_j)) {
    hp_id <- data_j$hospital[i]
    x_i   <- X_j[i, ]
    bk_i  <- bk2[, hp_id]
    
    mc_cdf <- numeric(Nmc)
    for (jdp in 1:Hdp) {
      beta_jdp   <- beta_tr[, , 1, jdp]
      mean_log_j <- as.numeric(beta_jdp %*% x_i) + bk_i
      weight_j   <- w_tr[, jdp]
      mc_cdf     <- mc_cdf + weight_j * plnorm(target_time, meanlog = mean_log_j, sdlog = sig_tr)
    }
    pred_surv[i] <- 1 - mean(mc_cdf)
  }
  
  data_j$pred_surv <- pred_surv
  return(data_j)
}

# 3. Calibration Data Extractor Across Risk Groups ------------------------------
get_calibration_data <- function(df_preds, target_time, model_name, n_g = 5) {
  # Assign risk groups (quintiles by default)
  df_preds <- df_preds %>%
    mutate(risk_group = ntile(pred_surv, n_g))
  
  calib_list <- list()
  
  for (g in 1:n_g) {
    sub_data <- df_preds %>% filter(risk_group == g)
    mean_pred <- mean(sub_data$pred_surv)
    
    # Estimate observed survival using Kaplan-Meier at target time
    km_fit <- survfit(Surv(Tstart, Tstop, status) ~ 1, data = sub_data)
    
    # Extract S(target_time) from KM curve
    km_sum <- summary(km_fit, times = target_time, extend = TRUE)
    
    calib_list[[g]] <- tibble(
      horizon   = paste(target_time, "Years"),
      Model     = model_name,
      group     = g,
      pred_mean = mean_pred,
      obs_km    = km_sum$surv[1],
      lower     = km_sum$lower[1],
      upper     = km_sum$upper[1]
    )
  }
  return(bind_rows(calib_list))
}

# 4. Generate Calibration Data at 5 and 10 Years -------------------------------
file_bnp <- "Split1/Res_bsSizeLNO_BNP_25_zellnerTrain.RData"
file_bp  <- "Split1/Res_bsSizeLNO_BNP_1_zellnerTrain.RData"

# Compute predictions
pred_bnp_5  <- get_indiv_pred_surv(file_bnp, res.test, tr_select = 5, target_time = 5)
pred_bnp_10 <- get_indiv_pred_surv(file_bnp, res.test, tr_select = 5, target_time = 10)

pred_bp_5   <- get_indiv_pred_surv(file_bp, res.test, tr_select = 5, target_time = 5)
pred_bp_10  <- get_indiv_pred_surv(file_bp, res.test, tr_select = 5, target_time = 10)

# Grouped calibration (using 5 risk strata)
calib_bnp_5  <- get_calibration_data(pred_bnp_5, 5, "BNP (H=25)", n_g = 5)
calib_bnp_10 <- get_calibration_data(pred_bnp_10, 10, "BNP (H=25)", n_g = 5)

calib_bp_5   <- get_calibration_data(pred_bp_5, 5, "BP (H=1)", n_g = 5)
calib_bp_10  <- get_calibration_data(pred_bp_10, 10, "BP (H=1)", n_g = 5)

# Concatenate and enforce factor ordering for time horizons
df_calib_all <- bind_rows(calib_bnp_5, calib_bnp_10, calib_bp_5, calib_bp_10) %>%
  mutate(horizon = factor(horizon, levels = c("5 Years", "10 Years")))

# 5. Render Calibration Plots --------------------------------------------------
p_calib <- ggplot(df_calib_all, aes(x = pred_mean, y = obs_km, color = Model, shape = Model)) +
  # Reference line for perfect calibration (45-degree diagonal)
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  # Vertical confidence intervals for observed survival
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.02, alpha = 0.6) +
  # Calibration curves and point estimates per model
  geom_line(linewidth = 0.8) +
  geom_point(size = 3) +
  # Facet panel by evaluation horizon
  facet_wrap(~ horizon) +
  scale_color_manual(values = c("BNP (H=25)" = "#2b5c8f", "BP (H=1)" = "#d95f02")) +
  scale_x_continuous(limits = c(0, 1), name = "Predicted Survival Probability") +
  scale_y_continuous(limits = c(0, 1), name = "Observed Survival Probability (Kaplan-Meier)") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "gray92"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Calibration Plots across Risk Strata (Transition 5)",
    subtitle = "Assessment of agreement between predicted and observed survival at 5 and 10 years"
  )

print(p_calib)
