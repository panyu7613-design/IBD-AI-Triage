# ============================================================
# 01_Shapley_AUC_Selection.R
# Shapley-AUC contribution + multi-method consistency feature selection
# ============================================================
# Purpose:
#   1) Use TRAINING data only to enumerate candidate feature subsets.
#   2) Estimate each subset's cross-validated AUC using Logistic regression.
#   3) Calculate exact Shapley-style marginal AUC contribution for each feature.
#   4) Combine Shapley-AUC contribution with multi-method consistency evidence
#      (LASSO, RF, XGBoost, SVM).
#   5) Export selected feature sets for subsequent models.
#
# GitHub Demo Version: 
#   Includes an auto-fallback to a generated dummy dataset to ensure 
#   100% reproducibility for reviewers without compromising patient privacy.
# ============================================================

message("==== Running 01_Shapley_AUC_Selection.R ====")

if (basename(getwd()) == "scripts") {
  setwd("..")
} else if (!dir.exists("scripts") && dir.exists("IBD_data&result/scripts")) {
  setwd("IBD_data&result")
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(caret)
  library(pROC)
  if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")
  library(ggrepel)
})

# ----------------------------
# Fallback Helpers for GitHub Standalone Demo
# ----------------------------
if (!exists("log_message")) log_message <- function(...) cat(paste0(..., "\n"))
if (!exists("log_warning")) log_warning <- function(...) warning(..., call. = FALSE)
if (!exists("pretty_var")) pretty_var <- function(x) x
if (!exists("save_pdf")) save_pdf <- function(p, file, ...) { pdf(file, ...); print(p); dev.off() }
if (!exists("theme_paper")) theme_paper <- function(base_size=11) ggplot2::theme_minimal(base_size=base_size)
if (!exists("safe_rank_desc")) safe_rank_desc <- function(x) rank(-x, ties.method = "min")
if (!exists("ensure_dirs")) ensure_dirs <- function() {}

if (!exists("resolve_predictors")) {
  resolve_predictors <- function(df, preds) {
    tibble::tibble(
      standard_name = preds,
      resolved_name = preds,
      found = preds %in% names(df)
    )
  }
}

helper_path <- if (file.exists("scripts/helpers_pipeline.R")) "scripts/helpers_pipeline.R" else "IBD_data&result/scripts/helpers_pipeline.R"
if (file.exists(helper_path)) source(helper_path, local = .GlobalEnv)

options(IBD_MAKE_FIGURES_ON_SOURCE = FALSE)
helper_fig_path <- "scripts/helpers_figures.R"
if (file.exists(helper_fig_path)) source(helper_fig_path, local = .GlobalEnv)

if (!exists("apply_unified_split")) {
  apply_unified_split <- function(df, split_obj, patient_id_col = NA_character_) {
    if (split_obj$method == "grouped_patient_id" && !is.na(patient_id_col) && (patient_id_col %in% names(df))) {
      train <- df %>% filter(as.character(.data[[patient_id_col]]) %in% as.character(split_obj$train_ids[[1]]))
      test  <- df %>% filter(as.character(.data[[patient_id_col]]) %in% as.character(split_obj$test_ids[[1]]))
    } else {
      train <- df[1:length(split_obj$train_ids[[1]]), ]
      test <- df[(length(split_obj$train_ids[[1]])+1):nrow(df), ]
    }
    list(train = train, test = test)
  }
}

# ----------------------------
# User-adjustable parameters
# ----------------------------
SEED_FEATURE_SELECTION <- 2026
CV_K_MAX <- 10
TOP_K_FINAL <- 4
TOP_K_EVIDENCE <- 4
MIN_CONSISTENCY_FOR_PRIMARY <- 0.50

set.seed(SEED_FEATURE_SELECTION)

# ----------------------------
# Required inputs (Adapted for GitHub Demo Mode)
# ----------------------------
real_data_exists <- file.exists("data/clean/data_model_features.rds") && file.exists("results/tables/screening/screening_evidence_combined.csv")

if (real_data_exists) {
  # --- REAL ENVIRONMENT MODE ---
  data_features <- readRDS("data/clean/data_model_features.rds")
  UNIFIED_SPLITS <- readRDS("data/clean/unified_splits.rds")
  task_cfg <- get_task_config()
  screening_evidence <- readr::read_csv("results/tables/screening/screening_evidence_combined.csv", show_col_types = FALSE)
  log_message("Real internal data loaded successfully.")
  
} else {
  # --- GITHUB DEMO MODE ---
  log_message("Real data not found. Initializing GitHub DEMO mode to protect patient privacy...")
  
  if (!file.exists("dummy_data.csv")) {
    log_message("Auto-generating synthetic dummy dataset (dummy_data.csv)...")
    set.seed(SEED_FEATURE_SELECTION)
    n_samples <- 200
    dummy <- tibble::tibble(
      patient_id = paste0("PT", 1:n_samples),
      age = rnorm(n_samples, 40, 15),
      FCP = runif(n_samples, 10, 1000),
      ASCA_IgG = runif(n_samples, 0, 100),
      PR3 = runif(n_samples, 0, 50),
      pANCA = runif(n_samples, 0, 50),
      IBD_status = sample(c("IBD", "non-IBD"), n_samples, replace = TRUE),
      Subtype = ifelse(IBD_status == "IBD", sample(c("CD", "UC"), sum(IBD_status == "IBD"), replace = TRUE), NA_character_),
      group2 = Subtype
    )
    readr::write_csv(dummy, "dummy_data.csv")
  }
  
  data_features <- readr::read_csv("dummy_data.csv", show_col_types = FALSE)
  
  train_ids <- data_features$patient_id[1:floor(nrow(data_features)*0.7)]
  test_ids <- data_features$patient_id[(floor(nrow(data_features)*0.7)+1):nrow(data_features)]
  UNIFIED_SPLITS <- list(
    patient_id_col = "patient_id",
    model1 = list(method = "grouped_patient_id", train_ids = list(train_ids), test_ids = list(test_ids)),
    model2 = list(method = "grouped_patient_id", train_ids = list(train_ids), test_ids = list(test_ids))
  )
  
  screening_evidence <- tibble::tibble(
    task = rep(c("Model1_IBD_vs_nonIBD", "Model2_CD_vs_UC"), each = 5),
    variable = rep(c("FCP", "ASCA_IgG", "age", "PR3", "pANCA"), 2),
    lasso_lambda_1se_selected = TRUE, 
    rf_shap_rank = 1L, 
    xgb_shap_rank = 1L, 
    svm_permutation_rank = 1L
  )
  
  task_cfg <- list(
    Model1_IBD_vs_nonIBD = list(outcome_col = "IBD_status", positive_class = "IBD", predictors_std = c("FCP", "ASCA_IgG", "age", "PR3", "pANCA")),
    Model2_CD_vs_UC = list(outcome_col = "Subtype", positive_class = "CD", predictors_std = c("ASCA_IgG", "age", "PR3", "pANCA"))
  )
}

# ----------------------------
# Output directories
# ----------------------------
dir.create("results/tables/feature_selection", recursive = TRUE, showWarnings = FALSE)
dir.create("models/feature_selection", recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Helper functions
# ============================================================

as_logical_safe <- function(x) {
  if (is.logical(x)) return(replace_na(x, FALSE))
  if (is.numeric(x)) return(replace_na(x != 0, FALSE))
  x_chr <- tolower(trimws(as.character(x)))
  replace_na(x_chr %in% c("true", "t", "1", "yes", "y"), FALSE)
}

make_subset_key <- function(vars) {
  if (length(vars) == 0) return("__EMPTY__")
  paste(sort(vars), collapse = "|")
}

all_feature_subsets <- function(vars) {
  p <- length(vars)
  out <- vector("list", 2^p)
  idx <- 1L
  out[[idx]] <- character(0)
  idx <- idx + 1L
  for (m in seq_len(p)) {
    cmb <- utils::combn(vars, m, simplify = FALSE)
    for (s in cmb) {
      out[[idx]] <- s
      idx <- idx + 1L
    }
  }
  out
}

make_stratified_folds <- function(y_factor, k_max = 10, seed = 2026) {
  set.seed(seed)
  class_counts <- table(y_factor)
  k <- min(k_max, as.integer(min(class_counts)))
  k <- max(2L, k)
  caret::createFolds(y_factor, k = k, list = TRUE, returnTrain = FALSE)
}

safe_auc <- function(y_factor, prob_positive, positive_class) {
  y_chr <- as.character(y_factor)
  ok <- !is.na(y_chr) & !is.na(prob_positive)
  y_chr <- y_chr[ok]
  prob_positive <- as.numeric(prob_positive[ok])
  if (length(unique(y_chr)) < 2 || length(unique(prob_positive)) < 2) return(0.5)
  neg_class <- setdiff(levels(factor(y_chr)), positive_class)
  if (length(neg_class) == 0) neg_class <- setdiff(unique(y_chr), positive_class)
  neg_class <- neg_class[[1]]
  auc_val <- tryCatch({
    as.numeric(pROC::auc(pROC::roc(
      response = factor(y_chr, levels = c(neg_class, positive_class)),
      predictor = prob_positive,
      levels = c(neg_class, positive_class),
      direction = "<",
      quiet = TRUE
    )))
  }, error = function(e) NA_real_)
  ifelse(is.na(auc_val), 0.5, auc_val)
}

cv_auc_logistic_subset <- function(df, outcome_col, positive_class, subset_vars, folds) {
  if (length(subset_vars) == 0) return(0.5)
  
  y_all <- df[[outcome_col]]
  pred_all <- rep(NA_real_, nrow(df))
  
  for (fold_idx in seq_along(folds)) {
    val_idx <- folds[[fold_idx]]
    tr_idx <- setdiff(seq_len(nrow(df)), val_idx)
    
    tr <- df[tr_idx, , drop = FALSE]
    va <- df[val_idx, , drop = FALSE]
    
    if (length(unique(as.character(tr[[outcome_col]]))) < 2) {
      pred_all[val_idx] <- mean(as.character(tr[[outcome_col]]) == positive_class, na.rm = TRUE)
      next
    }
    
    form <- stats::reformulate(subset_vars, response = outcome_col)
    fit <- tryCatch(
      suppressWarnings(stats::glm(form, data = tr, family = binomial())),
      error = function(e) NULL
    )
    
    if (is.null(fit)) {
      pred_all[val_idx] <- mean(as.character(tr[[outcome_col]]) == positive_class, na.rm = TRUE)
    } else {
      pp <- tryCatch(
        suppressWarnings(as.numeric(stats::predict(fit, newdata = va, type = "response"))),
        error = function(e) rep(mean(as.character(tr[[outcome_col]]) == positive_class, na.rm = TRUE), nrow(va))
      )
      pred_all[val_idx] <- pp
    }
  }
  safe_auc(y_all, pred_all, positive_class)
}

calculate_subset_auc_table <- function(df, outcome_col, positive_class, predictors, folds, task_name) {
  subsets <- all_feature_subsets(predictors)
  n_subsets <- length(subsets)
  log_message("[Feature selection][", task_name, "] Calculating CV AUC for ", n_subsets, " subsets...")
  
  purrr::map_dfr(seq_along(subsets), function(i) {
    vars_i <- subsets[[i]]
    if (i %% 50 == 0 || i == n_subsets) {
      log_message("[Feature selection][", task_name, "] subset ", i, "/", n_subsets)
    }
    auc_i <- cv_auc_logistic_subset(df, outcome_col, positive_class, vars_i, folds)
    tibble::tibble(
      task = task_name,
      subset_key = make_subset_key(vars_i),
      subset_size = length(vars_i),
      variables = paste(vars_i, collapse = ";"),
      cv_auc = auc_i
    )
  })
}

calculate_shapley_auc <- function(subset_auc_tbl, predictors, task_name) {
  p <- length(predictors)
  auc_lookup <- setNames(subset_auc_tbl$cv_auc, subset_auc_tbl$subset_key)
  
  purrr::map_dfr(predictors, function(j) {
    other_vars <- setdiff(predictors, j)
    subsets_without_j <- all_feature_subsets(other_vars)
    
    contrib_tbl <- purrr::map_dfr(subsets_without_j, function(S) {
      S_with_j <- c(S, j)
      key0 <- make_subset_key(S)
      key1 <- make_subset_key(S_with_j)
      s <- length(S)
      weight <- factorial(s) * factorial(p - s - 1) / factorial(p)
      delta <- unname(auc_lookup[[key1]]) - unname(auc_lookup[[key0]])
      tibble::tibble(
        task = task_name,
        variable_resolved = j,
        background_subset = make_subset_key(S),
        background_size = s,
        auc_without = unname(auc_lookup[[key0]]),
        auc_with = unname(auc_lookup[[key1]]),
        marginal_delta_auc = delta,
        shapley_weight = weight,
        weighted_delta_auc = weight * delta
      )
    })
    
    contrib_tbl %>%
      summarise(
        task = task_name,
        variable_resolved = j,
        phi_auc = sum(weighted_delta_auc, na.rm = TRUE),
        mean_marginal_delta_auc = mean(marginal_delta_auc, na.rm = TRUE),
        median_marginal_delta_auc = median(marginal_delta_auc, na.rm = TRUE),
        min_marginal_delta_auc = min(marginal_delta_auc, na.rm = TRUE),
        max_marginal_delta_auc = max(marginal_delta_auc, na.rm = TRUE),
        n_marginal_comparisons = dplyr::n(),
        .groups = "drop"
      )
  })
}

# [HOTFIX]: Safely handle missing columns to prevent dplyr .data$ missing column errors
extract_consistency_evidence <- function(task_name, pred_map, screening_evidence, top_k = 4) {
  ev <- screening_evidence %>%
    filter(task == task_name) %>%
    mutate(variable = as.character(variable))
  
  out <- pred_map %>%
    filter(found) %>%
    transmute(variable = standard_name, variable_resolved = resolved_name) %>%
    left_join(ev, by = "variable")
  
  # Ensure target columns exist (fallback for dummy or missing data)
  if (!"lasso_lambda_1se_selected" %in% names(out)) out$lasso_lambda_1se_selected <- FALSE
  if (!"rf_shap_rank" %in% names(out)) out$rf_shap_rank <- NA_integer_
  if (!"xgb_shap_rank" %in% names(out)) out$xgb_shap_rank <- NA_integer_
  if (!"svm_permutation_rank" %in% names(out)) out$svm_permutation_rank <- NA_integer_
  
  if (!"rf_shap_top4" %in% names(out)) out$rf_shap_top4 <- out$rf_shap_rank <= top_k
  if (!"xgb_shap_top4" %in% names(out)) out$xgb_shap_top4 <- out$xgb_shap_rank <= top_k
  if (!"svm_permutation_top4" %in% names(out)) out$svm_permutation_top4 <- out$svm_permutation_rank <= top_k
  
  out %>%
    mutate(
      lasso_selected = as_logical_safe(lasso_lambda_1se_selected),
      rf_shap_rank = suppressWarnings(as.integer(rf_shap_rank)),
      xgb_shap_rank = suppressWarnings(as.integer(xgb_shap_rank)),
      svm_permutation_rank = suppressWarnings(as.integer(svm_permutation_rank)),
      rf_shap_top4 = as_logical_safe(rf_shap_top4),
      xgb_shap_top4 = as_logical_safe(xgb_shap_top4),
      svm_permutation_top4 = as_logical_safe(svm_permutation_top4),
      method_votes = as.integer(lasso_selected) + as.integer(rf_shap_top4) + as.integer(xgb_shap_top4) + as.integer(svm_permutation_top4),
      consistency = method_votes / 4
    ) %>%
    select(variable, variable_resolved, lasso_selected, rf_shap_rank, rf_shap_top4, xgb_shap_rank, xgb_shap_top4, svm_permutation_rank, svm_permutation_top4, method_votes, consistency, everything())
}

select_final_features <- function(final_tbl, top_k = 4, min_consistency = 0.50) {
  ranked <- final_tbl %>%
    arrange(desc(F_score), desc(phi_auc_positive), desc(consistency), variable) %>%
    mutate(F_rank = row_number())
  
  primary <- ranked %>%
    filter(consistency >= min_consistency, phi_auc_positive > 0) %>%
    arrange(desc(F_score), desc(phi_auc_positive), desc(consistency)) %>%
    slice_head(n = top_k)
  
  if (nrow(primary) < top_k) {
    fill <- ranked %>%
      filter(!(variable %in% primary$variable), phi_auc_positive > 0) %>%
      arrange(desc(F_score), desc(phi_auc_positive), desc(consistency)) %>%
      slice_head(n = top_k - nrow(primary))
    primary <- bind_rows(primary, fill)
  }
  
  selected_vars <- primary$variable
  ranked %>%
    mutate(
      selected_final = variable %in% selected_vars,
      selection_note = case_when(
        selected_final & consistency >= min_consistency ~ "selected: high F_score with acceptable consistency",
        selected_final ~ "selected: fallback by positive Shapley-AUC contribution",
        TRUE ~ "not selected"
      )
    )
}

plot_feature_selection_plane <- function(final_tbl, task_name) {
  p <- final_tbl %>%
    mutate(variable = pretty_var(variable)) %>%
    ggplot(aes(x = phi_auc, y = consistency)) +
    geom_hline(yintercept = MIN_CONSISTENCY_FOR_PRIMARY, linetype = "dashed", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.3) +
    geom_point(aes(size = F_score, shape = selected_final, color = selected_final), alpha = 0.8) +
    ggrepel::geom_text_repel(
      aes(label = variable), size = 3.5, max.overlaps = Inf,
      box.padding = 0.8, point.padding = 0.5, min.segment.length = 0,
      segment.size = 0.3, segment.alpha = 0.65, segment.color = "grey40",
      force = 3, force_pull = 0.5, seed = 123, show.legend = FALSE
    ) +
    scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16)) +
    scale_color_manual(values = c(`FALSE` = "grey40", `TRUE` = "#D55E00")) +
    scale_size_continuous(name = "F score", range = c(1.8, 5.0)) +
    scale_x_continuous(name = "Shapley-AUC contribution", expand = expansion(mult = c(0.08, 0.18))) +
    scale_y_continuous(name = "Method consistency", limits = c(0, 1), expand = expansion(mult = c(0.08, 0.18))) +
    coord_cartesian(clip = "off") +
    labs(title = paste0(ifelse(task_name == "Model1_IBD_vs_nonIBD", "Model1", "Model2"), ": Shapley-AUC contribution vs multi-method consistency"), shape = "Selected", color = "Selected") +
    theme_paper(base_size = 11) +
    theme(legend.position = "right", plot.margin = ggplot2::margin(10, 25, 10, 10))
  
  save_pdf(p, paste0("results/tables/feature_selection/", task_name, "_feature_selection_plane.pdf"), width = 8, height = 5.5)
  invisible(p)
}

run_shapley_auc_feature_selection_for_task <- function(task_name, train_data, outcome_col, positive_class, predictors_std) {
  log_message("[Feature selection] Task: ", task_name)
  
  pred_map <- resolve_predictors(train_data, predictors_std)
  missing_vars <- pred_map %>% filter(!found) %>% pull(standard_name)
  if (length(missing_vars) > 0) log_warning("[Feature selection][", task_name, "] Missing predictors: ", paste(missing_vars, collapse = ", "))
  
  usable <- pred_map %>% filter(found)
  if (nrow(usable) == 0) stop("No usable predictors for ", task_name)
  
  df <- train_data %>%
    select(all_of(c(outcome_col, usable$resolved_name))) %>%
    tidyr::drop_na()
  
  df[[outcome_col]] <- factor(df[[outcome_col]])
  n_class <- table(df[[outcome_col]])
  if (length(n_class) < 2 || min(n_class) < 2) stop("Insufficient class counts for CV AUC in ", task_name)
  
  predictors_resolved <- usable$resolved_name
  names(predictors_resolved) <- usable$standard_name
  
  folds <- make_stratified_folds(df[[outcome_col]], k_max = CV_K_MAX, seed = SEED_FEATURE_SELECTION)
  
  subset_auc_tbl <- calculate_subset_auc_table(df, outcome_col, positive_class, unname(predictors_resolved), folds, task_name)
  shapley_tbl_resolved <- calculate_shapley_auc(subset_auc_tbl, unname(predictors_resolved), task_name)
  
  name_map <- usable %>% transmute(variable = standard_name, variable_resolved = resolved_name)
  
  shapley_tbl <- shapley_tbl_resolved %>%
    left_join(name_map, by = "variable_resolved") %>%
    relocate(variable, .before = variable_resolved) %>%
    mutate(phi_auc_positive = pmax(phi_auc, 0), phi_auc_rank = safe_rank_desc(phi_auc))
  
  consistency_tbl <- extract_consistency_evidence(task_name, pred_map, screening_evidence, TOP_K_EVIDENCE) %>%
    select(variable, variable_resolved, lasso_selected, rf_shap_rank, rf_shap_top4, xgb_shap_rank, xgb_shap_top4, svm_permutation_rank, svm_permutation_top4, method_votes, consistency)
  
  final_tbl <- shapley_tbl %>%
    left_join(consistency_tbl, by = c("variable", "variable_resolved")) %>%
    mutate(
      across(c(lasso_selected, rf_shap_top4, xgb_shap_top4, svm_permutation_top4), ~ replace_na(.x, FALSE)),
      method_votes = replace_na(method_votes, 0L),
      consistency = replace_na(consistency, 0),
      F_score = phi_auc_positive * consistency,
      task = task_name
    ) %>%
    select_final_features(top_k = TOP_K_FINAL, min_consistency = MIN_CONSISTENCY_FOR_PRIMARY)
  
  readr::write_csv(subset_auc_tbl, file.path("results/tables/feature_selection", paste0(task_name, "_subset_cv_auc.csv")))
  readr::write_csv(shapley_tbl, file.path("results/tables/feature_selection", paste0(task_name, "_shapley_auc_contribution.csv")))
  readr::write_csv(consistency_tbl, file.path("results/tables/feature_selection", paste0(task_name, "_method_consistency.csv")))
  readr::write_csv(final_tbl, file.path("results/tables/feature_selection", paste0(task_name, "_feature_selection_final.csv")))
  
  plot_feature_selection_plane(final_tbl, task_name)
  
  selected <- final_tbl %>% filter(selected_final) %>% arrange(F_rank)
  log_message("[Feature selection][", task_name, "] Selected variables: ", paste(selected$variable, collapse = "; "))
  
  final_tbl
}

# ============================================================
# Run both tasks
# ============================================================

split_model1 <- apply_unified_split(data_features, UNIFIED_SPLITS$model1, UNIFIED_SPLITS$patient_id_col %||% NA_character_)
split_model2 <- apply_unified_split(filter(data_features, !is.na(group2)), UNIFIED_SPLITS$model2, UNIFIED_SPLITS$patient_id_col %||% NA_character_)

res_model1 <- run_shapley_auc_feature_selection_for_task(
  task_name = "Model1_IBD_vs_nonIBD", train_data = split_model1$train,
  outcome_col = task_cfg$Model1_IBD_vs_nonIBD$outcome_col, positive_class = task_cfg$Model1_IBD_vs_nonIBD$positive_class,
  predictors_std = task_cfg$Model1_IBD_vs_nonIBD$predictors_std
)

res_model2 <- run_shapley_auc_feature_selection_for_task(
  task_name = "Model2_CD_vs_UC", train_data = split_model2$train,
  outcome_col = task_cfg$Model2_CD_vs_UC$outcome_col, positive_class = task_cfg$Model2_CD_vs_UC$positive_class,
  predictors_std = task_cfg$Model2_CD_vs_UC$predictors_std
)

combined_final <- bind_rows(res_model1, res_model2)
readr::write_csv(combined_final, "results/tables/feature_selection/feature_selection_final_combined.csv")

dir.create(file.path("results", "figures", "feature_selection"), recursive = TRUE, showWarnings = FALSE)

plot_fscore_bar <- function(df, task_name, out_file, title_text) {
  p <- df %>%
    mutate(
      variable_label = pretty_var(variable),
      variable_label = forcats::fct_reorder(variable_label, F_score),
      selected_final = as.logical(selected_final)
    ) %>%
    ggplot(aes(x = variable_label, y = F_score, fill = selected_final)) +
    geom_col(width = 0.65, alpha = 0.85) +
    coord_flip() +
    scale_fill_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "#BDBDBD"), guide = "none") +
    labs(x = NULL, y = "Composite screening index (F-score)", title = title_text) +
    theme_paper(12)
  save_pdf(p, out_file, width = 8, height = 5.5)
}

plot_fscore_bar(res_model1, "Model1_IBD_vs_nonIBD", "results/figures/feature_selection/Model1_IBD_vs_nonIBD_F_score_bar.pdf", "Model1")
plot_fscore_bar(res_model2, "Model2_CD_vs_UC", "results/figures/feature_selection/Model2_CD_vs_UC_F_score_bar.pdf", "Model2")

selected_feature_sets <- combined_final %>%
  filter(selected_final) %>%
  arrange(task, F_rank) %>%
  group_by(task) %>%
  summarise(
    feature_set = "selected",
    selected_variables = paste(variable, collapse = ";"),
    n_selected = n(),
    selection_rule = paste0("Top ", TOP_K_FINAL, " by F_score"),
    .groups = "drop"
  )
readr::write_csv(selected_feature_sets, "results/tables/feature_selection/selected_feature_sets.csv")

log_message("[Feature selection] Successfully completed and results saved!")