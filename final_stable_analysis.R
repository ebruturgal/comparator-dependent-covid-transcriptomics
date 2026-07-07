############################################################
# FINAL STABLE SCRIPT
# GSE282464
# Revised tasks + GPBoost LightGBM part kept (gpb.train)
############################################################

rm(list = ls())
gc()

############################################################
# 0) PACKAGES
############################################################
cran_pkgs <- c(
  "data.table", "dplyr", "stringr", "ggplot2", "tidyr", "pROC",
  "glmnet", "xgboost", "randomForest", "pheatmap", "scales",
  "gpboost", "GEOquery"
)
bioc_pkgs <- c("edgeR")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
for (p in bioc_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    BiocManager::install(p, ask = FALSE, update = FALSE)
  }
}

library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)
library(pROC)
library(glmnet)
library(xgboost)
library(randomForest)
library(pheatmap)
library(scales)
library(gpboost)
library(GEOquery)
library(edgeR)

set.seed(123)

############################################################
# 1) BUILD count_mat2 and pheno2
############################################################
gse_id <- "GSE282464"
untar_dir <- file.path(gse_id, "untarred_raw")

if (!dir.exists(untar_dir)) {
  stop("GSE282464/untarred_raw klas??r?? bulunamad??.")
}

count_files <- list.files(untar_dir, pattern = "\\.counts\\.gz$", full.names = TRUE)
if (length(count_files) == 0) stop(".counts.gz dosyas?? bulunamad??.")

read_one_count <- function(f) {
  dt <- fread(f, header = FALSE)
  dt <- dt[, 1:2]
  sample_id <- sub("_.*$", "", basename(f))
  colnames(dt) <- c("gene_id", sample_id)
  dt
}

count_list <- lapply(count_files, read_one_count)
count_merged <- Reduce(function(x, y) merge(x, y, by = "gene_id", all = TRUE), count_list)
count_merged[is.na(count_merged)] <- 0

count_mat <- as.data.frame(count_merged)
rownames(count_mat) <- count_mat$gene_id
count_mat <- count_mat[, -1, drop = FALSE]
count_mat <- as.matrix(sapply(count_mat, as.numeric))
rownames(count_mat) <- count_merged$gene_id

rm(count_list)
gc()

gset_list <- getGEO(gse_id, GSEMatrix = TRUE)
gset <- gset_list[[1]]
pheno <- pData(gset)

pheno$sample_id <- as.character(pheno$geo_accession)
common_samples <- intersect(colnames(count_mat), pheno$sample_id)

if (length(common_samples) == 0) stop("Metadata ile count matrix e??le??medi.")

pheno2 <- pheno %>% filter(sample_id %in% common_samples)
count_mat2 <- count_mat[, pheno2$sample_id, drop = FALSE]
pheno2 <- pheno2[match(colnames(count_mat2), pheno2$sample_id), , drop = FALSE]

stopifnot(all(colnames(count_mat2) == pheno2$sample_id))

cat("pheno2 dim:", dim(pheno2), "\n")
cat("count_mat2 dim:", dim(count_mat2), "\n")
print(table(pheno2$`type of_infection:ch1`, useNA = "ifany"))

############################################################
# 2) METRICS
############################################################
calc_auc <- function(y_true, prob) {
  if (length(unique(y_true)) < 2) return(NA_real_)
  as.numeric(pROC::roc(y_true, prob, quiet = TRUE)$auc)
}

calc_pr_auc <- function(y_true, prob) {
  if (length(unique(y_true)) < 2) return(NA_real_)
  if (sum(y_true == 1) == 0) return(NA_real_)
  
  ord <- order(prob, decreasing = TRUE)
  y_sorted <- y_true[ord]
  
  tp_cum <- cumsum(y_sorted == 1)
  fp_cum <- cumsum(y_sorted == 0)
  
  precision_curve <- tp_cum / (tp_cum + fp_cum)
  recall_curve <- tp_cum / sum(y_true == 1)
  
  precision_curve <- c(1, precision_curve)
  recall_curve <- c(0, recall_curve)
  
  sum(
    (recall_curve[-1] - recall_curve[-length(recall_curve)]) *
      (precision_curve[-1] + precision_curve[-length(precision_curve)]) / 2
  )
}

calc_metrics <- function(y_true, prob, threshold = 0.5) {
  pred <- ifelse(prob >= threshold, 1, 0)
  
  tp <- sum(pred == 1 & y_true == 1)
  tn <- sum(pred == 0 & y_true == 0)
  fp <- sum(pred == 1 & y_true == 0)
  fn <- sum(pred == 0 & y_true == 1)
  
  accuracy <- (tp + tn) / length(y_true)
  sensitivity <- ifelse((tp + fn) == 0, NA_real_, tp / (tp + fn))
  specificity <- ifelse((tn + fp) == 0, NA_real_, tn / (tn + fp))
  precision <- ifelse((tp + fp) == 0, NA_real_, tp / (tp + fp))
  recall <- sensitivity
  
  f1 <- ifelse(
    is.na(precision) || is.na(recall) || (precision + recall) == 0,
    NA_real_,
    2 * precision * recall / (precision + recall)
  )
  
  balanced_accuracy <- mean(c(sensitivity, specificity), na.rm = TRUE)
  auc <- calc_auc(y_true, prob)
  pr_auc <- calc_pr_auc(y_true, prob)
  
  mcc_denom <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  mcc <- ifelse(mcc_denom == 0, NA_real_,
                ((tp * tn) - (fp * fn)) / mcc_denom)
  
  brier <- mean((prob - y_true)^2)
  
  data.frame(
    auc = auc,
    pr_auc = pr_auc,
    accuracy = accuracy,
    balanced_accuracy = balanced_accuracy,
    sensitivity = sensitivity,
    specificity = specificity,
    precision = precision,
    f1 = f1,
    mcc = mcc,
    brier = brier
  )
}

############################################################
# 3) SAFE FOLDS
############################################################
make_subject_folds <- function(subject_table,
                               label_col = "label",
                               subject_col = "subject_id2",
                               k = 5,
                               seed = 123) {
  set.seed(seed)
  
  classes <- unique(subject_table[[label_col]])
  if (length(classes) != 2) stop("Binary classification expected.")
  
  out <- vector("list", k)
  for (i in seq_len(k)) out[[i]] <- character(0)
  
  for (cl in classes) {
    subj <- subject_table %>%
      filter(.data[[label_col]] == cl) %>%
      pull(.data[[subject_col]]) %>%
      unique()
    
    subj <- sample(subj)
    
    for (j in seq_along(subj)) {
      fold_id <- ((j - 1) %% k) + 1
      out[[fold_id]] <- c(out[[fold_id]], subj[j])
    }
  }
  
  out
}

############################################################
# 4) REVISED TASKS
############################################################
prepare_gse282464_task <- function(pheno2, count_mat2,
                                   task_name = c("covid_vs_influenza",
                                                 "covid_vs_noninfluenza_viral",
                                                 "covid_vs_sepsis")) {
  task_name <- match.arg(task_name)
  
  ph <- pheno2
  ph$subject_id2 <- ph$sample_id
  inf <- as.character(ph$`type of_infection:ch1`)
  
  if (task_name == "covid_vs_influenza") {
    ph$label <- case_when(
      inf == "COVID-19" ~ "COVID",
      inf %in% c("Influenza_A", "Influenza_B") ~ "Influenza",
      TRUE ~ NA_character_
    )
  }
  
  if (task_name == "covid_vs_noninfluenza_viral") {
    ph$label <- case_when(
      inf == "COVID-19" ~ "COVID",
      inf %in% c("Coronavirus",
                 "Co-infection (Viral/Viral)",
                 "Co-infection (viral/viral)",
                 "Co-infection (Viral/Fungal)") ~ "NonInfluenza_Viral",
      TRUE ~ NA_character_
    )
  }
  
  if (task_name == "covid_vs_sepsis") {
    ph$label <- case_when(
      inf == "COVID-19" ~ "COVID",
      inf %in% c("Sepsis", "Septic_shock") ~ "Sepsis_Group",
      TRUE ~ NA_character_
    )
  }
  
  cat("\nRaw label distribution for", task_name, ":\n")
  print(table(ph$label, useNA = "ifany"))
  
  ph_task <- ph %>% filter(!is.na(label))
  count_task <- count_mat2[, ph_task$sample_id, drop = FALSE]
  ph_task <- ph_task[match(colnames(count_task), ph_task$sample_id), , drop = FALSE]
  
  stopifnot(all(colnames(count_task) == ph_task$sample_id))
  
  y <- DGEList(counts = count_task)
  keep <- filterByExpr(y, group = ph_task$label)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")
  xmat <- cpm(y, log = TRUE, prior.count = 1)
  
  list(
    pheno_task = ph_task,
    xmat = xmat
  )
}

############################################################
# 5) MODEL FITTER
# GPBoost package kept, but only LightGBM part via gpb.train
############################################################
fit_one_model <- function(model_name, x_train_fs, x_test_fs, y_train, y_test) {
  x_train_fs <- as.matrix(x_train_fs)
  x_test_fs  <- as.matrix(x_test_fs)
  
  if (nrow(x_train_fs) == 0 || nrow(x_test_fs) == 0) return(NULL)
  if (length(unique(y_train)) < 2 || length(unique(y_test)) < 2) return(NULL)
  
  if (model_name == "GPBoost_LightGBM") {
    prob <- tryCatch({
      dtrain_plain <- gpb.Dataset(data = x_train_fs, label = y_train)
      
      params_plain <- list(
        objective = "binary",
        learning_rate = 0.03,
        max_depth = 4,
        num_leaves = 31,
        min_data_in_leaf = 5,
        feature_fraction = 0.8,
        bagging_fraction = 1.0,
        bagging_freq = 0,
        verbose = 0
      )
      
      bst_plain <- gpb.train(
        params = params_plain,
        data = dtrain_plain,
        nrounds = 120,
        verbose = 0
      )
      
      out_prob <- as.numeric(predict(bst_plain, x_test_fs))
      rm(dtrain_plain, bst_plain)
      gc()
      out_prob
    }, error = function(e) {
      cat("GPBoost_LightGBM error:", e$message, "\n")
      NULL
    })
    return(prob)
  }
  
  if (model_name == "ElasticNet") {
    prob <- tryCatch({
      nfolds_use <- max(2, min(5, min(table(y_train))))
      
      cv_fit <- cv.glmnet(
        x = as.matrix(x_train_fs),
        y = y_train,
        family = "binomial",
        alpha = 0.5,
        nfolds = nfolds_use,
        type.measure = "auc"
      )
      
      out_prob <- as.numeric(
        predict(cv_fit, newx = as.matrix(x_test_fs), s = "lambda.min", type = "response")
      )
      
      rm(cv_fit)
      gc()
      out_prob
    }, error = function(e) {
      cat("ElasticNet error:", e$message, "\n")
      NULL
    })
    return(prob)
  }
  
  if (model_name == "XGBoost") {
    prob <- tryCatch({
      dtrain_xgb <- xgb.DMatrix(data = x_train_fs, label = y_train)
      dtest_xgb  <- xgb.DMatrix(data = x_test_fs, label = y_test)
      
      params_xgb <- list(
        objective = "binary:logistic",
        eval_metric = "auc",
        eta = 0.03,
        max_depth = 4,
        subsample = 1.0,
        colsample_bytree = 0.8,
        min_child_weight = 1
      )
      
      bst_xgb <- xgb.train(
        params = params_xgb,
        data = dtrain_xgb,
        nrounds = 120,
        verbose = 0
      )
      
      out_prob <- as.numeric(predict(bst_xgb, dtest_xgb))
      rm(dtrain_xgb, dtest_xgb, bst_xgb)
      gc()
      out_prob
    }, error = function(e) {
      cat("XGBoost error:", e$message, "\n")
      NULL
    })
    return(prob)
  }
  
  if (model_name == "RandomForest") {
    prob <- tryCatch({
      df_train_rf <- as.data.frame(x_train_fs)
      df_test_rf  <- as.data.frame(x_test_fs)
      
      y_train_fac <- factor(ifelse(y_train == 1, "Case", "Control"),
                            levels = c("Control", "Case"))
      
      rf_model <- randomForest(
        x = df_train_rf,
        y = y_train_fac,
        ntree = 300,
        importance = TRUE
      )
      
      out_prob <- as.numeric(predict(rf_model, df_test_rf, type = "prob")[, "Case"])
      rm(df_train_rf, df_test_rf, y_train_fac, rf_model)
      gc()
      out_prob
    }, error = function(e) {
      cat("RandomForest error:", e$message, "\n")
      NULL
    })
    return(prob)
  }
  
  stop("Unknown model name.")
}

############################################################
# 6) RUN BENCHMARK
############################################################
run_task_full_benchmark <- function(task_obj,
                                    dataset_name,
                                    task_name,
                                    seed = 123,
                                    top_k = 100,
                                    k = 5) {
  pheno_task <- task_obj$pheno_task
  xmat <- task_obj$xmat
  
  subject_table <- pheno_task %>%
    distinct(subject_id2, label)
  
  fold_subjects <- make_subject_folds(
    subject_table,
    label_col = "label",
    subject_col = "subject_id2",
    k = k,
    seed = seed
  )
  
  positive_class <- "COVID"
  model_list <- c("GPBoost_LightGBM", "ElasticNet", "XGBoost", "RandomForest")
  
  cv_results <- data.frame(
    dataset = character(),
    task = character(),
    fold = integer(),
    model = character(),
    n_train = integer(),
    n_test = integer(),
    n_test_case = integer(),
    n_test_control = integer(),
    auc = numeric(),
    pr_auc = numeric(),
    accuracy = numeric(),
    balanced_accuracy = numeric(),
    sensitivity = numeric(),
    specificity = numeric(),
    precision = numeric(),
    f1 = numeric(),
    mcc = numeric(),
    brier = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (fold_i in seq_len(k)) {
    cat("\n====================\n")
    cat("Task:", task_name, "| Fold:", fold_i, "\n")
    cat("====================\n")
    
    test_subjects <- fold_subjects[[fold_i]]
    if (length(test_subjects) == 0) next
    
    test_idx <- pheno_task$subject_id2 %in% test_subjects
    train_idx <- !test_idx
    
    if (sum(test_idx) == 0 || sum(train_idx) == 0) next
    
    x_train <- t(xmat[, train_idx, drop = FALSE])
    x_test  <- t(xmat[, test_idx, drop = FALSE])
    
    y_train <- ifelse(pheno_task$label[train_idx] == positive_class, 1, 0)
    y_test  <- ifelse(pheno_task$label[test_idx] == positive_class, 1, 0)
    
    cat("Train:", nrow(x_train), " Test:", nrow(x_test), "\n")
    print(table(y_train))
    print(table(y_test))
    
    if (length(unique(y_train)) < 2) next
    if (length(unique(y_test)) < 2) next
    
    gene_var <- apply(x_train, 2, var, na.rm = TRUE)
    top_genes <- names(sort(gene_var, decreasing = TRUE))[1:min(top_k, length(gene_var))]
    
    x_train_fs <- x_train[, top_genes, drop = FALSE]
    x_test_fs  <- x_test[, top_genes, drop = FALSE]
    
    for (m in model_list) {
      cat("Running:", m, "\n")
      
      prob <- fit_one_model(
        model_name = m,
        x_train_fs = x_train_fs,
        x_test_fs = x_test_fs,
        y_train = y_train,
        y_test = y_test
      )
      
      if (is.null(prob)) next
      if (length(prob) != nrow(x_test_fs)) next
      if (!all(is.finite(prob))) next
      
      met <- calc_metrics(y_test, prob)
      
      cv_results <- rbind(
        cv_results,
        data.frame(
          dataset = dataset_name,
          task = task_name,
          fold = fold_i,
          model = m,
          n_train = nrow(x_train_fs),
          n_test = nrow(x_test_fs),
          n_test_case = sum(y_test == 1),
          n_test_control = sum(y_test == 0),
          met,
          stringsAsFactors = FALSE
        )
      )
    }
    
    rm(x_train, x_test, x_train_fs, x_test_fs, y_train, y_test, gene_var, top_genes)
    gc()
  }
  
  summary_results <- cv_results %>%
    group_by(dataset, task, model) %>%
    summarise(
      mean_auc = mean(auc, na.rm = TRUE),
      sd_auc = sd(auc, na.rm = TRUE),
      mean_pr_auc = mean(pr_auc, na.rm = TRUE),
      sd_pr_auc = sd(pr_auc, na.rm = TRUE),
      mean_acc = mean(accuracy, na.rm = TRUE),
      sd_acc = sd(accuracy, na.rm = TRUE),
      mean_bal_acc = mean(balanced_accuracy, na.rm = TRUE),
      sd_bal_acc = sd(balanced_accuracy, na.rm = TRUE),
      mean_f1 = mean(f1, na.rm = TRUE),
      sd_f1 = sd(f1, na.rm = TRUE),
      mean_mcc = mean(mcc, na.rm = TRUE),
      sd_mcc = sd(mcc, na.rm = TRUE),
      mean_brier = mean(brier, na.rm = TRUE),
      sd_brier = sd(brier, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_auc))
  
  list(
    cv_results = cv_results,
    summary_results = summary_results
  )
}

############################################################
# 7) PLOTS
############################################################
plot_summary_bar <- function(summary_df, title_text, outfile) {
  plot_df <- summary_df %>%
    select(model, mean_auc, mean_pr_auc, mean_acc, mean_bal_acc, mean_f1) %>%
    pivot_longer(
      cols = -model,
      names_to = "metric",
      values_to = "value"
    ) %>%
    mutate(
      metric = recode(metric,
                      mean_auc = "AUC",
                      mean_pr_auc = "PR AUC",
                      mean_acc = "Accuracy",
                      mean_bal_acc = "Balanced Accuracy",
                      mean_f1 = "F1 Score"),
      model = factor(model, levels = summary_df$model[order(summary_df$mean_auc, decreasing = TRUE)])
    )
  
  p <- ggplot(plot_df, aes(x = model, y = value, fill = metric)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.72) +
    scale_fill_manual(values = c(
      "AUC" = "#1b9e77",
      "PR AUC" = "#d95f02",
      "Accuracy" = "#7570b3",
      "Balanced Accuracy" = "#e7298a",
      "F1 Score" = "#66a61e"
    )) +
    scale_y_continuous(limits = c(0, 1.05), labels = number_format(accuracy = 0.01)) +
    labs(
      title = title_text,
      x = NULL,
      y = "Performance metric value",
      fill = NULL
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 20, hjust = 1, face = "bold"),
      legend.position = "top",
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  print(p)
  ggsave(outfile, p, width = 11, height = 6, dpi = 600)
}

plot_foldwise <- function(cv_df, title_text, outfile) {
  plot_df <- cv_df %>%
    select(fold, model, auc, pr_auc, accuracy, balanced_accuracy, f1) %>%
    pivot_longer(
      cols = c(auc, pr_auc, accuracy, balanced_accuracy, f1),
      names_to = "metric",
      values_to = "value"
    ) %>%
    mutate(
      metric = recode(metric,
                      auc = "AUC",
                      pr_auc = "PR AUC",
                      accuracy = "Accuracy",
                      balanced_accuracy = "Balanced Accuracy",
                      f1 = "F1 Score")
    )
  
  p <- ggplot(plot_df, aes(x = factor(fold), y = value, fill = model)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.72) +
    facet_wrap(~ metric, scales = "free_y") +
    scale_y_continuous(labels = number_format(accuracy = 0.01)) +
    labs(
      title = title_text,
      x = "Fold",
      y = "Value",
      fill = "Model"
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
  
  print(p)
  ggsave(outfile, p, width = 12, height = 8, dpi = 600)
}

plot_heatmap_metrics <- function(summary_df, title_text, outfile) {
  heat_df <- summary_df %>%
    select(model, mean_auc, mean_pr_auc, mean_acc, mean_bal_acc, mean_f1, mean_brier) %>%
    as.data.frame()
  
  rownames(heat_df) <- heat_df$model
  heat_df <- heat_df[, -1, drop = FALSE]
  
  png(outfile, width = 2200, height = 1600, res = 300)
  pheatmap(
    as.matrix(heat_df),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    scale = "none",
    color = colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(100),
    display_numbers = TRUE,
    number_format = "%.3f",
    fontsize = 12,
    fontsize_number = 10,
    border_color = "grey85",
    main = title_text
  )
  dev.off()
}

############################################################
# 8) FEATURE IMPORTANCE
############################################################
extract_model_importance <- function(xmat, pheno_task,
                                     model_name = c("ElasticNet", "XGBoost", "RandomForest"),
                                     top_k = 100,
                                     positive_class = "COVID") {
  model_name <- match.arg(model_name)
  
  x <- t(xmat)
  y <- ifelse(pheno_task$label == positive_class, 1, 0)
  
  gene_var <- apply(x, 2, var, na.rm = TRUE)
  top_genes <- names(sort(gene_var, decreasing = TRUE))[1:min(top_k, length(gene_var))]
  x_fs <- x[, top_genes, drop = FALSE]
  
  if (model_name == "ElasticNet") {
    nfolds_use <- max(2, min(5, min(table(y))))
    cv_fit <- cv.glmnet(
      x = as.matrix(x_fs),
      y = y,
      family = "binomial",
      alpha = 0.5,
      nfolds = nfolds_use,
      type.measure = "auc"
    )
    
    coef_mat <- as.matrix(coef(cv_fit, s = "lambda.min"))
    coef_df <- data.frame(
      gene = rownames(coef_mat),
      coefficient = as.numeric(coef_mat[, 1]),
      stringsAsFactors = FALSE
    ) %>%
      filter(gene != "(Intercept)") %>%
      mutate(abs_coef = abs(coefficient)) %>%
      arrange(desc(abs_coef))
    
    return(coef_df)
  }
  
  if (model_name == "XGBoost") {
    dtrain <- xgb.DMatrix(data = as.matrix(x_fs), label = y)
    
    params_xgb <- list(
      objective = "binary:logistic",
      eval_metric = "auc",
      eta = 0.03,
      max_depth = 4,
      subsample = 1.0,
      colsample_bytree = 0.8,
      min_child_weight = 1
    )
    
    bst_xgb <- xgb.train(
      params = params_xgb,
      data = dtrain,
      nrounds = 120,
      verbose = 0
    )
    
    imp <- xgb.importance(model = bst_xgb)
    imp <- as.data.frame(imp)
    colnames(imp)[colnames(imp) == "Feature"] <- "gene"
    return(imp)
  }
  
  if (model_name == "RandomForest") {
    df_train <- as.data.frame(x_fs, check.names = FALSE)
    y_fac <- factor(ifelse(y == 1, "Case", "Control"), levels = c("Control", "Case"))
    
    rf_model <- randomForest(
      x = df_train,
      y = y_fac,
      ntree = 300,
      importance = TRUE
    )
    
    imp <- importance(rf_model)
    imp_df <- data.frame(
      gene = rownames(imp),
      MeanDecreaseGini = imp[, "MeanDecreaseGini"],
      stringsAsFactors = FALSE
    ) %>%
      arrange(desc(MeanDecreaseGini))
    
    return(imp_df)
  }
}

############################################################
# 9) ENRICHMENT
############################################################
run_enrichment <- function(gene_symbols, out_prefix = "enrichment") {
  gene_symbols <- unique(gene_symbols)
  gene_symbols <- gene_symbols[!is.na(gene_symbols) & gene_symbols != ""]
  
  eg <- tryCatch(
    clusterProfiler::bitr(
      gene_symbols,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    ),
    error = function(e) NULL
  )
  
  if (is.null(eg) || nrow(eg) == 0) {
    cat("No genes mapped for enrichment.\n")
    return(NULL)
  }
  
  go_res <- clusterProfiler::enrichGO(
    gene = unique(eg$ENTREZID),
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    readable = TRUE
  )
  
  write.csv(as.data.frame(go_res), paste0(out_prefix, "_GO.csv"), row.names = FALSE)
  go_res
}

############################################################
# 10) RUN TASKS
############################################################
task_names_new <- c(
  "covid_vs_influenza",
  "covid_vs_noninfluenza_viral",
  "covid_vs_sepsis"
)

task_objects <- lapply(task_names_new, function(tt) {
  prepare_gse282464_task(pheno2, count_mat2, task_name = tt)
})
names(task_objects) <- task_names_new

res_influenza <- run_task_full_benchmark(
  task_obj = task_objects[["covid_vs_influenza"]],
  dataset_name = "GSE282464",
  task_name = "covid_vs_influenza",
  seed = 123, top_k = 100, k = 5
)

res_noninfluenza <- run_task_full_benchmark(
  task_obj = task_objects[["covid_vs_noninfluenza_viral"]],
  dataset_name = "GSE282464",
  task_name = "covid_vs_noninfluenza_viral",
  seed = 123, top_k = 100, k = 5
)

res_sepsis <- run_task_full_benchmark(
  task_obj = task_objects[["covid_vs_sepsis"]],
  dataset_name = "GSE282464",
  task_name = "covid_vs_sepsis",
  seed = 123, top_k = 100, k = 5
)

############################################################
# 11) SAVE TABLES
############################################################
write.csv(res_influenza$cv_results,
          "GSE282464_COVID_vs_Influenza_Revised_All_CV_Results.csv",
          row.names = FALSE)
write.csv(res_influenza$summary_results,
          "GSE282464_COVID_vs_Influenza_Revised_Summary.csv",
          row.names = FALSE)

write.csv(res_noninfluenza$cv_results,
          "GSE282464_COVID_vs_nonInfluenzaViral_Revised_All_CV_Results.csv",
          row.names = FALSE)
write.csv(res_noninfluenza$summary_results,
          "GSE282464_COVID_vs_nonInfluenzaViral_Revised_Summary.csv",
          row.names = FALSE)

write.csv(res_sepsis$cv_results,
          "GSE282464_COVID_vs_Sepsis_Revised_All_CV_Results.csv",
          row.names = FALSE)
write.csv(res_sepsis$summary_results,
          "GSE282464_COVID_vs_Sepsis_Revised_Summary.csv",
          row.names = FALSE)

task_char_table <- data.frame(
  task = c("COVID vs Influenza", "COVID vs Non-Influenza Viral", "COVID vs Sepsis/Septic shock"),
  comparator_definition = c(
    "Influenza_A + Influenza_B only",
    "Coronavirus + Viral/Viral + viral/viral + Viral/Fungal; influenza excluded",
    "Sepsis + Septic_shock"
  ),
  note = c(
    "Pure influenza comparator",
    "Non-influenza viral comparator; Bac/Viral mixed groups excluded",
    "Clinically defined sepsis comparator; not interpreted as purely viral"
  ),
  stringsAsFactors = FALSE
)
write.csv(task_char_table, "Table1_Task_Characteristics_Revised.csv", row.names = FALSE)

############################################################
# 12) PLOTS
############################################################
plot_summary_bar(
  res_influenza$summary_results,
  "GSE282464: COVID-19 vs Influenza",
  "Figure_COVID_vs_Influenza_Barplot_Revised.png"
)

plot_summary_bar(
  res_noninfluenza$summary_results,
  "GSE282464: COVID-19 vs Non-Influenza Viral Infections",
  "Figure_COVID_vs_nonInfluenzaViral_Barplot_Revised.png"
)

plot_summary_bar(
  res_sepsis$summary_results,
  "GSE282464: COVID-19 vs Sepsis / Septic Shock",
  "Figure_COVID_vs_Sepsis_Barplot_Revised.png"
)

plot_foldwise(
  res_influenza$cv_results,
  "GSE282464: Fold-wise Performance (COVID-19 vs Influenza)",
  "Figure_COVID_vs_Influenza_Foldwise_Revised.png"
)

plot_foldwise(
  res_noninfluenza$cv_results,
  "GSE282464: Fold-wise Performance (COVID-19 vs Non-Influenza Viral)",
  "Figure_COVID_vs_nonInfluenzaViral_Foldwise_Revised.png"
)

plot_foldwise(
  res_sepsis$cv_results,
  "GSE282464: Fold-wise Performance (COVID-19 vs Sepsis)",
  "Figure_COVID_vs_Sepsis_Foldwise_Revised.png"
)

plot_heatmap_metrics(
  res_influenza$summary_results,
  "GSE282464: COVID-19 vs Influenza",
  "Figure_COVID_vs_Influenza_Heatmap_Revised.png"
)

plot_heatmap_metrics(
  res_noninfluenza$summary_results,
  "GSE282464: COVID-19 vs Non-Influenza Viral Infections",
  "Figure_COVID_vs_nonInfluenzaViral_Heatmap_Revised.png"
)

plot_heatmap_metrics(
  res_sepsis$summary_results,
  "GSE282464: COVID-19 vs Sepsis / Septic Shock",
  "Figure_COVID_vs_Sepsis_Heatmap_Revised.png"
)

############################################################
# 13) FEATURE IMPORTANCE + TOP GENE PLOTS
############################################################
enet_importance_list <- lapply(names(task_objects), function(tt) {
  obj <- task_objects[[tt]]
  imp <- extract_model_importance(
    xmat = obj$xmat,
    pheno_task = obj$pheno_task,
    model_name = "ElasticNet",
    top_k = 100,
    positive_class = "COVID"
  )
  imp$task <- tt
  imp$model <- "ElasticNet"
  imp
})
names(enet_importance_list) <- names(task_objects)
enet_importance_all <- bind_rows(enet_importance_list)
write.csv(enet_importance_all, "ElasticNet_Importance_Revised_AllTasks.csv", row.names = FALSE)

rf_importance_list <- lapply(names(task_objects), function(tt) {
  obj <- task_objects[[tt]]
  imp <- extract_model_importance(
    xmat = obj$xmat,
    pheno_task = obj$pheno_task,
    model_name = "RandomForest",
    top_k = 100,
    positive_class = "COVID"
  )
  imp$task <- tt
  imp$model <- "RandomForest"
  imp
})
names(rf_importance_list) <- names(task_objects)
rf_importance_all <- bind_rows(rf_importance_list)
write.csv(rf_importance_all, "RandomForest_Importance_Revised_AllTasks.csv", row.names = FALSE)

plot_top_enet <- enet_importance_all %>%
  group_by(task) %>%
  slice_max(order_by = abs_coef, n = 15) %>%
  ungroup() %>%
  mutate(gene = reorder(gene, abs_coef))

p_enet <- ggplot(plot_top_enet, aes(x = gene, y = abs_coef, fill = task)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ task, scales = "free_y") +
  labs(
    title = "Top Elastic Net Features Across Revised Tasks",
    x = NULL,
    y = "Absolute coefficient"
  ) +
  theme_bw(base_size = 13)

print(p_enet)
ggsave("TopGenes_ElasticNet_Revised_AllTasks.png", p_enet, width = 12, height = 8, dpi = 600)

plot_top_rf <- rf_importance_all %>%
  group_by(task) %>%
  slice_max(order_by = MeanDecreaseGini, n = 15) %>%
  ungroup() %>%
  mutate(gene = reorder(gene, MeanDecreaseGini))

p_rf <- ggplot(plot_top_rf, aes(x = gene, y = MeanDecreaseGini, fill = task)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ task, scales = "free_y") +
  labs(
    title = "Top Random Forest Features Across Revised Tasks",
    x = NULL,
    y = "Mean Decrease Gini"
  ) +
  theme_bw(base_size = 13)

print(p_rf)
ggsave("TopGenes_RandomForest_Revised_AllTasks.png", p_rf, width = 12, height = 8, dpi = 600)

############################################################
# 14) SIMPLE ENRICHMENT EXAMPLE
############################################################
common_enet <- Reduce(intersect, lapply(task_names_new, function(tt) {
  enet_importance_all %>% filter(task == tt) %>% slice_head(n = 20) %>% pull(gene)
}))

enrich_common_enet <- run_enrichment(common_enet, out_prefix = "CommonGenes_ElasticNet_Revised")

############################################################
# 15) CROSS-TASK COMPARISON
############################################################
cross_task_df <- bind_rows(
  res_influenza$summary_results %>% mutate(task_label = "COVID vs Influenza"),
  res_noninfluenza$summary_results %>% mutate(task_label = "COVID vs Non-Influenza Viral"),
  res_sepsis$summary_results %>% mutate(task_label = "COVID vs Sepsis")
) %>%
  select(task_label, model, mean_auc, mean_pr_auc, mean_acc, mean_bal_acc, mean_f1)

p_cross_auc <- ggplot(cross_task_df, aes(x = task_label, y = mean_auc, fill = model)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72) +
  scale_y_continuous(limits = c(0, 1.05)) +
  labs(
    title = "Cross-task model comparison",
    x = NULL,
    y = "Mean AUC"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 20, hjust = 1)
  )

print(p_cross_auc)
ggsave("CrossTask_ModelComparison_AUC_Revised.png", p_cross_auc, width = 10, height = 6, dpi = 600)

cross_heat <- cross_task_df %>%
  select(task_label, model, mean_auc) %>%
  pivot_wider(names_from = model, values_from = mean_auc) %>%
  as.data.frame()

rownames(cross_heat) <- cross_heat$task_label
cross_heat <- cross_heat[, -1, drop = FALSE]

png("CrossTask_AUC_Heatmap_Revised.png", width = 2200, height = 1600, res = 300)
pheatmap(
  as.matrix(cross_heat),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(100),
  display_numbers = TRUE,
  number_format = "%.3f",
  main = "Cross-task AUC comparison"
)
dev.off()

############################################################
# 16) MASTER FILES
############################################################
all_cv_results_long <- bind_rows(
  res_influenza$cv_results,
  res_noninfluenza$cv_results,
  res_sepsis$cv_results
)
write.csv(all_cv_results_long, "All_CV_Results_Long_Revised.csv", row.names = FALSE)

paper_ready_summary <- bind_rows(
  res_influenza$summary_results %>% mutate(task_label = "COVID vs Influenza"),
  res_noninfluenza$summary_results %>% mutate(task_label = "COVID vs Non-Influenza Viral"),
  res_sepsis$summary_results %>% mutate(task_label = "COVID vs Sepsis")
) %>%
  select(dataset, task_label, model,
         mean_auc, sd_auc,
         mean_pr_auc, sd_pr_auc,
         mean_acc, sd_acc,
         mean_bal_acc, sd_bal_acc,
         mean_f1, sd_f1,
         mean_mcc, sd_mcc,
         mean_brier, sd_brier)

write.csv(paper_ready_summary, "Paper_Ready_Summary_Revised.csv", row.names = FALSE)

############################################################
# 17) DONE
############################################################
cat("\nSaved revised files.\n")





class_col_ext <- "source_name_ch1"

library(dplyr)
library(tibble)
library(stringr)
library(ggplot2)


############################################################
# 18) EXTERNAL VALIDATION: GSE199816
# Locked Elastic Net model from GSE282464 COVID vs Sepsis
############################################################

safe_col_lower <- function(df) {
  names(df) <- tolower(names(df))
  df
}

collapse_duplicate_genes <- function(expr_mat) {
  expr_df <- as.data.frame(expr_mat, check.names = FALSE)
  expr_df$gene <- rownames(expr_df)
  
  expr_df <- expr_df %>%
    group_by(gene) %>%
    summarise(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
  
  out <- as.matrix(expr_df[, -1, drop = FALSE])
  rownames(out) <- expr_df$gene
  out
}

fit_train_scaler <- function(x_train) {
  mu <- colMeans(x_train, na.rm = TRUE)
  sdv <- apply(x_train, 2, sd, na.rm = TRUE)
  sdv[is.na(sdv) | sdv == 0] <- 1
  list(center = mu, scale = sdv)
}

apply_scaler <- function(x, scaler) {
  sweep(sweep(x, 2, scaler$center, "-"), 2, scaler$scale, "/")
}

calc_external_metrics <- function(y_true_fac, prob, threshold = 0.5) {
  y_true <- ifelse(y_true_fac == "COVID", 1, 0)
  pred <- ifelse(prob >= threshold, 1, 0)
  
  tp <- sum(pred == 1 & y_true == 1)
  tn <- sum(pred == 0 & y_true == 0)
  fp <- sum(pred == 1 & y_true == 0)
  fn <- sum(pred == 0 & y_true == 1)
  
  accuracy <- (tp + tn) / length(y_true)
  sensitivity <- ifelse((tp + fn) == 0, NA_real_, tp / (tp + fn))
  specificity <- ifelse((tn + fp) == 0, NA_real_, tn / (tn + fp))
  precision <- ifelse((tp + fp) == 0, NA_real_, tp / (tp + fp))
  recall <- sensitivity
  
  f1 <- ifelse(
    is.na(precision) || is.na(recall) || (precision + recall) == 0,
    NA_real_,
    2 * precision * recall / (precision + recall)
  )
  
  balanced_accuracy <- mean(c(sensitivity, specificity), na.rm = TRUE)
  
  auc_obj <- tryCatch(
    pROC::roc(y_true, prob, quiet = TRUE),
    error = function(e) NULL
  )
  
  auc <- if (is.null(auc_obj)) NA_real_ else as.numeric(pROC::auc(auc_obj))
  auc_ci <- if (is.null(auc_obj)) c(NA_real_, NA_real_, NA_real_) else as.numeric(pROC::ci.auc(auc_obj))
  
  mcc_denom <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  mcc <- ifelse(mcc_denom == 0, NA_real_,
                ((tp * tn) - (fp * fn)) / mcc_denom)
  
  brier <- mean((prob - y_true)^2)
  
  data.frame(
    auc = auc,
    auc_ci_low = auc_ci[1],
    auc_ci_mid = auc_ci[2],
    auc_ci_high = auc_ci[3],
    accuracy = accuracy,
    balanced_accuracy = balanced_accuracy,
    sensitivity = sensitivity,
    specificity = specificity,
    precision = precision,
    f1 = f1,
    mcc = mcc,
    brier = brier,
    stringsAsFactors = FALSE
  )
}

label_covid_sepsis_external <- function(pheno_df, class_col) {
  ph <- pheno_df
  x <- tolower(as.character(ph[[class_col]]))
  
  ph$label_binary <- dplyr::case_when(
    str_detect(x, "covid|covid-19|sars-cov-2|sars cov 2|coronavirus disease") ~ "COVID",
    str_detect(x, "sepsis|septic shock|septic_shock|septic") ~ "Sepsis_Group",
    TRUE ~ NA_character_
  )
  ph
}

############################################################
# 19) BUILD FINAL LOCKED ELASTIC NET MODEL
# using GSE282464 COVID vs Sepsis full development cohort
############################################################

cat("\n====================\n")
cat("FINAL LOCKED MODEL: GSE282464 COVID vs Sepsis\n")
cat("====================\n")

dev_obj <- task_objects[["covid_vs_sepsis"]]
dev_pheno <- dev_obj$pheno_task
dev_xmat <- dev_obj$xmat

positive_class_ext <- "COVID"

# full development design
x_dev_full <- t(dev_xmat)
y_dev_full <- ifelse(dev_pheno$label == positive_class_ext, 1, 0)

cat("Development samples:", nrow(x_dev_full), "\n")
print(table(dev_pheno$label))

# same feature strategy style as main script: variance-based top genes
top_k_external <- 100
gene_var_dev <- apply(x_dev_full, 2, var, na.rm = TRUE)
selected_genes_external <- names(sort(gene_var_dev, decreasing = TRUE))[1:min(top_k_external, length(gene_var_dev))]

x_dev_fs <- x_dev_full[, selected_genes_external, drop = FALSE]

# scale using development only
scaler_external <- fit_train_scaler(x_dev_fs)
x_dev_fs_sc <- apply_scaler(x_dev_fs, scaler_external)

# final Elastic Net on full development cohort
nfolds_use_final <- max(2, min(5, min(table(y_dev_full))))

cv_fit_external <- cv.glmnet(
  x = as.matrix(x_dev_fs_sc),
  y = y_dev_full,
  family = "binomial",
  alpha = 0.5,
  nfolds = nfolds_use_final,
  type.measure = "auc"
)

final_enet_external <- glmnet(
  x = as.matrix(x_dev_fs_sc),
  y = y_dev_full,
  family = "binomial",
  alpha = 0.5,
  lambda = cv_fit_external$lambda.min
)

final_coef_external <- as.matrix(coef(final_enet_external))
final_coef_external_df <- data.frame(
  gene = rownames(final_coef_external),
  coefficient = as.numeric(final_coef_external[, 1]),
  stringsAsFactors = FALSE
) %>%
  filter(gene != "(Intercept)") %>%
  filter(coefficient != 0) %>%
  mutate(abs_coef = abs(coefficient)) %>%
  arrange(desc(abs_coef))

write.csv(
  final_coef_external_df,
  "GSE282464_COVID_vs_Sepsis_Final_ElasticNet_Coefficients_For_ExternalValidation.csv",
  row.names = FALSE
)




























if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("AnnotationDbi", quietly = TRUE)) BiocManager::install("AnnotationDbi", ask = FALSE, update = FALSE)
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db", ask = FALSE, update = FALSE)

library(AnnotationDbi)
library(org.Hs.eg.db)

############################################################
# 20) DOWNLOAD AND PREPARE GSE199816
############################################################

safe_col_lower <- function(df) {
  names(df) <- tolower(names(df))
  df
}

cat("\n====================\n")
cat("DOWNLOADING / PREPARING GSE199816\n")
cat("====================\n")

gse_ext_id <- "GSE199816"
dir.create(gse_ext_id, showWarnings = FALSE)

# Download GEO metadata
gset_ext_list <- getGEO(gse_ext_id, GSEMatrix = TRUE, getGPL = FALSE)
gset_ext <- gset_ext_list[[1]]

pheno_ext <- pData(gset_ext) %>%
  as.data.frame()

pheno_ext$sample_id <- rownames(pheno_ext)
pheno_ext <- pheno_ext[, c("sample_id", setdiff(names(pheno_ext), "sample_id"))]
pheno_ext <- safe_col_lower(pheno_ext)

cat("External pheno dim:", dim(pheno_ext), "\n")
cat("External metadata columns:\n")
print(names(pheno_ext))

############################################################
# 21) DOWNLOAD RAW SUPPLEMENTARY FILES
############################################################

getGEOSuppFiles(gse_ext_id, baseDir = ".", fetch_files = TRUE)

raw_tar_candidates <- c(
  file.path(gse_ext_id, "GSE199816_RAW.tar"),
  file.path(gse_ext_id, "suppl", "GSE199816_RAW.tar")
)

raw_tar <- raw_tar_candidates[file.exists(raw_tar_candidates)][1]

if (is.na(raw_tar) || !file.exists(raw_tar)) {
  stop("GSE199816_RAW.tar bulunamad??.")
}

cat("Using RAW tar:", raw_tar, "\n")

untar_dir_ext <- file.path(gse_ext_id, "untarred_raw")
dir.create(untar_dir_ext, showWarnings = FALSE)

untar(raw_tar, exdir = untar_dir_ext)

raw_files_ext <- list.files(untar_dir_ext, full.names = TRUE)
cat("Number of extracted raw files:", length(raw_files_ext), "\n")
print(head(basename(raw_files_ext), 10))

if (length(raw_files_ext) == 0) {
  stop("RAW tar a????ld?? ama i??inde dosya bulunamad??.")
}

############################################################
# 22) BUILD EXTERNAL MATRIX FROM RAW FILES
############################################################

read_one_external_file <- function(f) {
  dt <- tryCatch(
    fread(f, header = TRUE, data.table = FALSE),
    error = function(e) NULL
  )
  
  if (is.null(dt) || ncol(dt) < 2) {
    dt <- tryCatch(
      fread(f, header = FALSE, data.table = FALSE),
      error = function(e) NULL
    )
  }
  
  if (is.null(dt) || ncol(dt) < 2) {
    stop(paste("Dosya okunamad?? veya yetersiz kolon var:", basename(f)))
  }
  
  cn <- tolower(colnames(dt))
  
  gene_col_idx <- which(cn %in% c(
    "gene_id", "geneid", "gene", "ensembl_gene_id", "id", "symbol"
  ))
  if (length(gene_col_idx) == 0) {
    gene_col_idx <- 1
  } else {
    gene_col_idx <- gene_col_idx[1]
  }
  
  value_col_idx <- which(cn %in% c(
    "count", "counts", "raw_count", "raw_counts", "expected_count"
  ))
  if (length(value_col_idx) == 0) {
    other_cols <- setdiff(seq_len(ncol(dt)), gene_col_idx)
    if (length(other_cols) == 0) {
      stop(paste("Value column bulunamad??:", basename(f)))
    }
    value_col_idx <- other_cols[1]
  } else {
    value_col_idx <- value_col_idx[1]
  }
  
  sample_id <- sub("_.*$", "", basename(f))
  
  out <- dt[, c(gene_col_idx, value_col_idx), drop = FALSE]
  colnames(out) <- c("gene_id", sample_id)
  out
}

ext_list <- lapply(raw_files_ext, read_one_external_file)
ext_merged <- Reduce(function(x, y) merge(x, y, by = "gene_id", all = TRUE), ext_list)
ext_merged[is.na(ext_merged)] <- 0

expr_ext <- as.data.frame(ext_merged, check.names = FALSE)
rownames(expr_ext) <- expr_ext$gene_id
expr_ext <- expr_ext[, -1, drop = FALSE]
expr_ext <- as.matrix(sapply(expr_ext, as.numeric))
rownames(expr_ext) <- ext_merged$gene_id

cat("\nExternal raw matrix dim:", dim(expr_ext), "\n")
cat("Example external gene IDs:\n")
print(head(rownames(expr_ext), 10))

############################################################
# 23) DEFINE EXTERNAL COVID VS SEPSIS LABELS
############################################################

label_covid_sepsis_external <- function(pheno_df, class_col) {
  ph <- pheno_df
  x <- tolower(as.character(ph[[class_col]]))
  
  ph$label_binary <- dplyr::case_when(
    str_detect(x, "covid|covid-19|sars-cov-2|sars cov 2|coronavirus disease") ~ "COVID",
    str_detect(x, "sepsis|septic shock|septic_shock|septic") ~ "Sepsis_Group",
    TRUE ~ NA_character_
  )
  ph
}

# Prefer diagnosis:ch1 if available
if ("diagnosis:ch1" %in% names(pheno_ext)) {
  class_col_ext <- "diagnosis:ch1"
} else if ("source_name_ch1" %in% names(pheno_ext)) {
  class_col_ext <- "source_name_ch1"
} else {
  stop("Ne diagnosis:ch1 ne de source_name_ch1 bulundu.")
}

cat("Using phenotype column for labels:", class_col_ext, "\n")

pheno_ext <- label_covid_sepsis_external(pheno_ext, class_col = class_col_ext)

cat("\nRaw external label distribution:\n")
print(table(pheno_ext$label_binary, useNA = "ifany"))

pheno_ext_bin <- pheno_ext %>%
  filter(!is.na(label_binary)) %>%
  distinct(sample_id, .keep_all = TRUE)

common_ext_samples <- intersect(colnames(expr_ext), pheno_ext_bin$sample_id)

if (length(common_ext_samples) == 0) {
  stop("External expression matrix ile phenotype e??le??medi.")
}

expr_ext2 <- expr_ext[, common_ext_samples, drop = FALSE]
pheno_ext_bin <- pheno_ext_bin[match(colnames(expr_ext2), pheno_ext_bin$sample_id), , drop = FALSE]

stopifnot(all(colnames(expr_ext2) == pheno_ext_bin$sample_id))

cat("\nFiltered external matrix dim:", dim(expr_ext2), "\n")
cat("Filtered external label distribution:\n")
print(table(pheno_ext_bin$label_binary))











############################################################
# 24) HARMONIZE IDS AND MATCH FEATURES
############################################################

mode(expr_ext2) <- "numeric"

# External rownames already SYMBOL-like, just collapse duplicates
expr_ext2 <- collapse_duplicate_genes(expr_ext2)

# Development selected genes are ENSG -> map to SYMBOL
selected_genes_external_clean <- sub("\\..*$", "", selected_genes_external)

map_dev <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(selected_genes_external_clean),
  keytype = "ENSEMBL",
  columns = c("ENSEMBL", "SYMBOL")
)

map_dev <- map_dev[!is.na(map_dev$SYMBOL) & map_dev$SYMBOL != "", ]
map_dev <- map_dev[!duplicated(map_dev$ENSEMBL), ]

selected_symbols_external <- map_dev$SYMBOL[match(selected_genes_external_clean, map_dev$ENSEMBL)]

map_table_external <- data.frame(
  ensembl = selected_genes_external,
  symbol = selected_symbols_external,
  stringsAsFactors = FALSE
)

write.csv(
  map_table_external,
  "Development_SelectedGenes_ENSEMBL_to_SYMBOL.csv",
  row.names = FALSE
)

selected_symbols_external <- selected_symbols_external[!is.na(selected_symbols_external) & selected_symbols_external != ""]
selected_symbols_external <- unique(selected_symbols_external)

cat("Mapped development selected genes to symbols:", length(selected_symbols_external), "\n")

common_genes_external <- intersect(selected_symbols_external, rownames(expr_ext2))
cat("Overlapping genes between locked model and GSE199816 after ENSG->SYMBOL mapping:", length(common_genes_external), "\n")

if (length(common_genes_external) < 20) {
  cat("Example mapped development symbols:\n")
  print(head(selected_symbols_external, 20))
  cat("Example external rownames:\n")
  print(head(rownames(expr_ext2), 20))
  stop("Too few overlapping genes (<20) after ENSG-to-SYMBOL mapping.")
}

selected_genes_external_ordered <- selected_symbols_external[selected_symbols_external %in% common_genes_external]


























############################################################
# 25) PREPARE EXTERNAL MATRIX AND PREDICT
############################################################

# Original model feature order (must stay exactly the same)
model_feature_order <- colnames(x_dev_fs)

# Clean ENSG ids
model_feature_order_clean <- sub("\\..*$", "", model_feature_order)

# Map development ENSG -> SYMBOL
map_dev_full <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(model_feature_order_clean),
  keytype = "ENSEMBL",
  columns = c("ENSEMBL", "SYMBOL")
)

map_dev_full <- map_dev_full[!is.na(map_dev_full$SYMBOL) & map_dev_full$SYMBOL != "", ]
map_dev_full <- map_dev_full[!duplicated(map_dev_full$ENSEMBL), ]

model_feature_symbols <- map_dev_full$SYMBOL[match(model_feature_order_clean, map_dev_full$ENSEMBL)]

feature_map_df <- data.frame(
  model_ensembl = model_feature_order,
  model_ensembl_clean = model_feature_order_clean,
  symbol = model_feature_symbols,
  stringsAsFactors = FALSE
)

write.csv(
  feature_map_df,
  "ExternalValidation_ModelFeatureMap_ENSEMBL_to_SYMBOL.csv",
  row.names = FALSE
)

# Start with full external matrix in model feature space
n_ext <- ncol(expr_ext2)
x_ext_full <- matrix(
  NA_real_,
  nrow = n_ext,
  ncol = length(model_feature_order)
)

colnames(x_ext_full) <- model_feature_order
rownames(x_ext_full) <- colnames(expr_ext2)

# Fill columns that can be matched through SYMBOL
matched_count <- 0

for (i in seq_len(nrow(feature_map_df))) {
  ens_col <- feature_map_df$model_ensembl[i]
  sym_row <- feature_map_df$symbol[i]
  
  if (!is.na(sym_row) && sym_row != "" && sym_row %in% rownames(expr_ext2)) {
    x_ext_full[, ens_col] <- as.numeric(expr_ext2[sym_row, rownames(x_ext_full)])
    matched_count <- matched_count + 1
  }
}

cat("Number of model features matched in external data:", matched_count, "of", ncol(x_ext_full), "\n")

# Impute missing features with training mean
# After scaling, these become 0
for (j in seq_len(ncol(x_ext_full))) {
  feat <- colnames(x_ext_full)[j]
  miss_idx <- is.na(x_ext_full[, j])
  
  if (any(miss_idx)) {
    x_ext_full[miss_idx, j] <- scaler_external$center[feat]
  }
}

# Safety checks
stopifnot(ncol(x_ext_full) == length(model_feature_order))
stopifnot(all(colnames(x_ext_full) == model_feature_order))

# Scale using locked development scaler
scaler_locked <- list(
  center = scaler_external$center[model_feature_order],
  scale  = scaler_external$scale[model_feature_order]
)

x_ext_sc <- apply_scaler(x_ext_full, scaler_locked)

# Final labels
y_ext_fac <- factor(pheno_ext_bin$label_binary, levels = c("Sepsis_Group", "COVID"))

# Predict
ext_prob <- as.numeric(
  predict(
    final_enet_external,
    newx = as.matrix(x_ext_sc),
    type = "response"
  )
)

external_predictions_df <- data.frame(
  sample_id = pheno_ext_bin$sample_id,
  truth = as.character(y_ext_fac),
  prob_covid = ext_prob,
  pred_label_05 = ifelse(ext_prob >= 0.5, "COVID", "Sepsis_Group"),
  stringsAsFactors = FALSE
)

write.csv(
  external_predictions_df,
  "GSE199816_ExternalValidation_Predictions.csv",
  row.names = FALSE
)

external_metrics_df <- calc_external_metrics(y_ext_fac, ext_prob, threshold = 0.5)

write.csv(
  external_metrics_df,
  "GSE199816_ExternalValidation_Metrics.csv",
  row.names = FALSE
)

cat("\nExternal validation metrics:\n")
print(external_metrics_df)






















############################################################
# 29) OPTIMIZE THRESHOLD BY YOUDEN INDEX
############################################################

roc_obj_ext <- pROC::roc(
  response = ifelse(y_ext_fac == "COVID", 1, 0),
  predictor = ext_prob,
  quiet = TRUE
)

youden_coords <- pROC::coords(
  roc_obj_ext,
  x = "best",
  best.method = "youden",
  ret = c("threshold", "sensitivity", "specificity"),
  transpose = FALSE
)

best_threshold <- as.numeric(youden_coords["threshold"])
cat("Best Youden threshold:", best_threshold, "\n")

external_metrics_youden_df <- calc_external_metrics(
  y_true_fac = y_ext_fac,
  prob = ext_prob,
  threshold = best_threshold
)

external_metrics_youden_df$threshold <- best_threshold

write.csv(
  external_metrics_youden_df,
  "GSE199816_ExternalValidation_Metrics_YoudenThreshold.csv",
  row.names = FALSE
)

cat("\nExternal validation metrics at Youden threshold:\n")
print(external_metrics_youden_df)

external_predictions_df$pred_label_youden <- ifelse(
  external_predictions_df$prob_covid >= best_threshold,
  "COVID",
  "Sepsis_Group"
)

write.csv(
  external_predictions_df,
  "GSE199816_ExternalValidation_Predictions_WithYouden.csv",
  row.names = FALSE
)











############################################################
# 30) ROC FIGURE WITH YOUDEN THRESHOLD POINT
############################################################

# Rebuild ROC-related objects safely
truth01_ext <- ifelse(y_ext_fac == "COVID", 1, 0)

roc_obj_ext <- pROC::roc(
  response = truth01_ext,
  predictor = ext_prob,
  quiet = TRUE
)

auc_ext <- as.numeric(pROC::auc(roc_obj_ext))
auc_ci_ext <- as.numeric(pROC::ci.auc(roc_obj_ext))

youden_coords <- pROC::coords(
  roc_obj_ext,
  x = "best",
  best.method = "youden",
  ret = c("threshold", "sensitivity", "specificity"),
  transpose = FALSE
)

best_threshold <- as.numeric(youden_coords["threshold"])
best_sens <- as.numeric(youden_coords["sensitivity"])
best_spec <- as.numeric(youden_coords["specificity"])
best_fpr <- 1 - best_spec

roc_df_ext <- data.frame(
  fpr = 1 - roc_obj_ext$specificities,
  tpr = roc_obj_ext$sensitivities
)

p_ext_roc_youden <- ggplot(roc_df_ext, aes(x = fpr, y = tpr)) +
  geom_line(linewidth = 1.1, color = "#1b9e77") +
  geom_abline(intercept = 0, slope = 1, linetype = 2, color = "grey40") +
  geom_point(
    aes(x = best_fpr, y = best_sens),
    size = 3,
    color = "#d95f02"
  ) +
  annotate(
    "text",
    x = best_fpr,
    y = best_sens,
    label = paste0("Youden threshold = ", sprintf("%.3f", best_threshold)),
    vjust = -1,
    size = 4
  ) +
  labs(
    title = "External validation ROC: GSE199816",
    subtitle = paste0(
      "Elastic Net | AUC = ",
      sprintf("%.3f", auc_ext),
      " (95% CI ",
      sprintf("%.3f", auc_ci_ext[1]), " - ",
      sprintf("%.3f", auc_ci_ext[3]), ")"
    ),
    x = "False Positive Rate",
    y = "True Positive Rate"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

print(p_ext_roc_youden)

ggsave(
  "Figure_ExternalValidation_ROC_GSE199816_Youden.png",
  p_ext_roc_youden,
  width = 7,
  height = 5.5,
  dpi = 600
)

############################################################
# 31) YOUDEN THRESHOLD METRICS
############################################################

calc_external_metrics <- function(y_true_fac, prob, threshold = 0.5) {
  y_true <- ifelse(y_true_fac == "COVID", 1, 0)
  pred <- ifelse(prob >= threshold, 1, 0)
  
  tp <- sum(pred == 1 & y_true == 1)
  tn <- sum(pred == 0 & y_true == 0)
  fp <- sum(pred == 1 & y_true == 0)
  fn <- sum(pred == 0 & y_true == 1)
  
  accuracy <- (tp + tn) / length(y_true)
  sensitivity <- ifelse((tp + fn) == 0, NA_real_, tp / (tp + fn))
  specificity <- ifelse((tn + fp) == 0, NA_real_, tn / (tn + fp))
  precision <- ifelse((tp + fp) == 0, NA_real_, tp / (tp + fp))
  recall <- sensitivity
  
  f1 <- ifelse(
    is.na(precision) || is.na(recall) || (precision + recall) == 0,
    NA_real_,
    2 * precision * recall / (precision + recall)
  )
  
  balanced_accuracy <- mean(c(sensitivity, specificity), na.rm = TRUE)
  
  auc_obj <- tryCatch(
    pROC::roc(y_true, prob, quiet = TRUE),
    error = function(e) NULL
  )
  
  auc <- if (is.null(auc_obj)) NA_real_ else as.numeric(pROC::auc(auc_obj))
  auc_ci <- if (is.null(auc_obj)) c(NA_real_, NA_real_, NA_real_) else as.numeric(pROC::ci.auc(auc_obj))
  
  mcc_denom <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  mcc <- ifelse(mcc_denom == 0, NA_real_,
                ((tp * tn) - (fp * fn)) / mcc_denom)
  
  brier <- mean((prob - y_true)^2)
  
  data.frame(
    auc = auc,
    auc_ci_low = auc_ci[1],
    auc_ci_mid = auc_ci[2],
    auc_ci_high = auc_ci[3],
    accuracy = accuracy,
    balanced_accuracy = balanced_accuracy,
    sensitivity = sensitivity,
    specificity = specificity,
    precision = precision,
    f1 = f1,
    mcc = mcc,
    brier = brier,
    stringsAsFactors = FALSE
  )
}

external_metrics_youden_df <- calc_external_metrics(
  y_true_fac = y_ext_fac,
  prob = ext_prob,
  threshold = best_threshold
)

external_metrics_youden_df$threshold <- best_threshold

write.csv(
  external_metrics_youden_df,
  "GSE199816_ExternalValidation_Metrics_YoudenThreshold.csv",
  row.names = FALSE
)

cat("\nExternal validation metrics at Youden threshold:\n")
print(external_metrics_youden_df)

external_predictions_df$pred_label_youden <- ifelse(
  external_predictions_df$prob_covid >= best_threshold,
  "COVID",
  "Sepsis_Group"
)

write.csv(
  external_predictions_df,
  "GSE199816_ExternalValidation_Predictions_WithYouden.csv",
  row.names = FALSE
)

############################################################
# 32) FINAL ELASTIC NET COEFFICIENT BIOLOGY
############################################################

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("AnnotationDbi", quietly = TRUE)) BiocManager::install("AnnotationDbi", ask = FALSE, update = FALSE)
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db", ask = FALSE, update = FALSE)

library(AnnotationDbi)
library(org.Hs.eg.db)

coef_tbl <- final_coef_external_df
coef_tbl$ensembl_clean <- sub("\\..*$", "", coef_tbl$gene)

coef_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(coef_tbl$ensembl_clean),
  keytype = "ENSEMBL",
  columns = c("ENSEMBL", "SYMBOL", "GENENAME")
)

coef_map <- coef_map[!duplicated(coef_map$ENSEMBL), ]

coef_tbl$symbol <- coef_map$SYMBOL[match(coef_tbl$ensembl_clean, coef_map$ENSEMBL)]
coef_tbl$genename <- coef_map$GENENAME[match(coef_tbl$ensembl_clean, coef_map$ENSEMBL)]

coef_tbl <- coef_tbl %>%
  mutate(direction = ifelse(coefficient > 0, "Higher in COVID", "Higher in Sepsis")) %>%
  arrange(desc(abs_coef))

write.csv(
  coef_tbl,
  "Final_ElasticNet_Coefficients_Annotated.csv",
  row.names = FALSE
)

top_coef_tbl <- coef_tbl %>%
  filter(!is.na(symbol), symbol != "") %>%
  slice_head(n = 20)

write.csv(
  top_coef_tbl,
  "Top20_Final_ElasticNet_Coefficients_Annotated.csv",
  row.names = FALSE
)

cat("\nTop annotated coefficients:\n")
print(top_coef_tbl[, c("gene", "symbol", "genename", "coefficient", "abs_coef", "direction")])

############################################################
# 33) KEY BIOLOGICAL GENES TABLE
############################################################

key_genes <- c("IFI27", "CXCL10", "SIGLEC14", "MMP8", "CTSG", "LAMP3", "MDK", "MOV10")

key_gene_tbl <- coef_tbl %>%
  filter(symbol %in% key_genes) %>%
  arrange(desc(abs_coef))

write.csv(
  key_gene_tbl,
  "Key_Biological_Genes_In_FinalModel.csv",
  row.names = FALSE
)

cat("\nKey biological genes found in final model:\n")
print(key_gene_tbl[, c("symbol", "genename", "coefficient", "abs_coef", "direction")])

############################################################
# 34) TOP COEFFICIENT FIGURE
############################################################

plot_coef_tbl <- coef_tbl %>%
  filter(!is.na(symbol), symbol != "") %>%
  slice_head(n = 20) %>%
  mutate(symbol = factor(symbol, levels = rev(symbol)))

p_coef <- ggplot(plot_coef_tbl, aes(x = symbol, y = coefficient, fill = direction)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Top Elastic Net coefficients in the locked development model",
    x = NULL,
    y = "Coefficient"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

print(p_coef)

ggsave(
  "Figure_TopElasticNetCoefficients_Annotated.png",
  p_coef,
  width = 8,
  height = 6,
  dpi = 600
)

############################################################
# 35) DONE
############################################################

cat("\nSaved files:\n")
cat("- Figure_ExternalValidation_ROC_GSE199816_Youden.png\n")
cat("- GSE199816_ExternalValidation_Metrics_YoudenThreshold.csv\n")
cat("- GSE199816_ExternalValidation_Predictions_WithYouden.csv\n")
cat("- Final_ElasticNet_Coefficients_Annotated.csv\n")
cat("- Top20_Final_ElasticNet_Coefficients_Annotated.csv\n")
cat("- Key_Biological_Genes_In_FinalModel.csv\n")
cat("- Figure_TopElasticNetCoefficients_Annotated.png\n")










############################################################
# DECISION CURVE ANALYSIS
############################################################

if (!requireNamespace("rmda", quietly = TRUE)) {
  install.packages("rmda")
}

library(rmda)
library(ggplot2)

dca_df <- data.frame(
  outcome = ifelse(y_ext_fac == "COVID", 1, 0),
  pred = ext_prob
)

dca_model <- decision_curve(
  outcome ~ pred,
  data = dca_df,
  family = binomial(link = "logit"),
  thresholds = seq(0.01, 0.75, by = 0.01),
  confidence.intervals = 0.95,
  study.design = "cohort"
)

pdf("Figure_DecisionCurve_ExternalValidation.pdf", width = 8, height = 6)
plot_decision_curve(
  dca_model,
  curve.names = "Elastic Net",
  xlab = "Threshold Probability",
  ylab = "Net Benefit",
  standardize = FALSE
)
dev.off()

















############################################################
# FULL WORKING SCRIPT
# External validation with GSE161731
# Development: GSE282464
# Task: COVID vs Influenza
############################################################

rm(list = ls())
gc()

############################################################
# 0) PACKAGES
############################################################
cran_pkgs <- c(
  "data.table", "dplyr", "stringr", "ggplot2", "pROC",
  "glmnet", "GEOquery", "tibble", "rmda"
)
bioc_pkgs <- c("edgeR", "AnnotationDbi", "org.Hs.eg.db")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
for (p in bioc_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    BiocManager::install(p, ask = FALSE, update = FALSE)
  }
}

library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(pROC)
library(glmnet)
library(GEOquery)
library(tibble)
library(rmda)
library(edgeR)
library(AnnotationDbi)
library(org.Hs.eg.db)

set.seed(123)

############################################################
# 1) HELPERS
############################################################

safe_col_lower <- function(df) {
  names(df) <- tolower(names(df))
  df
}

collapse_duplicate_genes <- function(expr_mat) {
  expr_df <- as.data.frame(expr_mat, check.names = FALSE)
  expr_df$gene <- rownames(expr_df)
  
  expr_df <- expr_df %>%
    group_by(gene) %>%
    summarise(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
  
  out <- as.matrix(expr_df[, -1, drop = FALSE])
  rownames(out) <- expr_df$gene
  out
}

fit_train_scaler <- function(x_train) {
  mu <- colMeans(x_train, na.rm = TRUE)
  sdv <- apply(x_train, 2, sd, na.rm = TRUE)
  sdv[is.na(sdv) | sdv == 0] <- 1
  list(center = mu, scale = sdv)
}

apply_scaler <- function(x, scaler) {
  sweep(sweep(x, 2, scaler$center, "-"), 2, scaler$scale, "/")
}

calc_external_metrics <- function(y_true_fac, prob, threshold = 0.5, positive_class = "COVID") {
  y_true <- ifelse(y_true_fac == positive_class, 1, 0)
  pred <- ifelse(prob >= threshold, 1, 0)
  
  tp <- sum(pred == 1 & y_true == 1)
  tn <- sum(pred == 0 & y_true == 0)
  fp <- sum(pred == 1 & y_true == 0)
  fn <- sum(pred == 0 & y_true == 1)
  
  accuracy <- (tp + tn) / length(y_true)
  sensitivity <- ifelse((tp + fn) == 0, NA_real_, tp / (tp + fn))
  specificity <- ifelse((tn + fp) == 0, NA_real_, tn / (tn + fp))
  precision <- ifelse((tp + fp) == 0, NA_real_, tp / (tp + fp))
  recall <- sensitivity
  
  f1 <- ifelse(
    is.na(precision) || is.na(recall) || (precision + recall) == 0,
    NA_real_,
    2 * precision * recall / (precision + recall)
  )
  
  balanced_accuracy <- mean(c(sensitivity, specificity), na.rm = TRUE)
  
  auc_obj <- tryCatch(
    pROC::roc(y_true, prob, quiet = TRUE),
    error = function(e) NULL
  )
  
  auc <- if (is.null(auc_obj)) NA_real_ else as.numeric(pROC::auc(auc_obj))
  auc_ci <- if (is.null(auc_obj)) c(NA_real_, NA_real_, NA_real_) else as.numeric(pROC::ci.auc(auc_obj))
  
  mcc_denom <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  mcc <- ifelse(mcc_denom == 0, NA_real_,
                ((tp * tn) - (fp * fn)) / mcc_denom)
  
  brier <- mean((prob - y_true)^2)
  
  data.frame(
    auc = auc,
    auc_ci_low = auc_ci[1],
    auc_ci_mid = auc_ci[2],
    auc_ci_high = auc_ci[3],
    accuracy = accuracy,
    balanced_accuracy = balanced_accuracy,
    sensitivity = sensitivity,
    specificity = specificity,
    precision = precision,
    f1 = f1,
    mcc = mcc,
    brier = brier,
    stringsAsFactors = FALSE
  )
}

############################################################
# 2) BUILD GSE282464 DEVELOPMENT DATA
############################################################

gse_id <- "GSE282464"
untar_dir <- file.path(gse_id, "untarred_raw")

if (!dir.exists(untar_dir)) {
  stop("GSE282464/untarred_raw klas??r?? bulunamad??.")
}

count_files <- list.files(untar_dir, pattern = "\\.counts\\.gz$", full.names = TRUE)
if (length(count_files) == 0) stop(".counts.gz dosyas?? bulunamad??.")

read_one_count <- function(f) {
  dt <- fread(f, header = FALSE)
  dt <- dt[, 1:2]
  sample_id <- sub("_.*$", "", basename(f))
  colnames(dt) <- c("gene_id", sample_id)
  dt
}

count_list <- lapply(count_files, read_one_count)
count_merged <- Reduce(function(x, y) merge(x, y, by = "gene_id", all = TRUE), count_list)
count_merged[is.na(count_merged)] <- 0

count_mat <- as.data.frame(count_merged)
rownames(count_mat) <- count_mat$gene_id
count_mat <- count_mat[, -1, drop = FALSE]
count_mat <- as.matrix(sapply(count_mat, as.numeric))
rownames(count_mat) <- count_merged$gene_id

gset_list <- getGEO(gse_id, GSEMatrix = TRUE)
gset <- gset_list[[1]]
pheno <- pData(gset)

pheno$sample_id <- as.character(pheno$geo_accession)
common_samples <- intersect(colnames(count_mat), pheno$sample_id)

if (length(common_samples) == 0) stop("Metadata ile count matrix e??le??medi.")

pheno2 <- pheno %>% filter(sample_id %in% common_samples)
count_mat2 <- count_mat[, pheno2$sample_id, drop = FALSE]
pheno2 <- pheno2[match(colnames(count_mat2), pheno2$sample_id), , drop = FALSE]

stopifnot(all(colnames(count_mat2) == pheno2$sample_id))

############################################################
# 3) PREPARE DEVELOPMENT TASK: COVID VS INFLUENZA
############################################################

prepare_gse282464_task <- function(pheno2, count_mat2) {
  ph <- pheno2
  inf <- as.character(ph$`type of_infection:ch1`)
  
  ph$label <- dplyr::case_when(
    inf == "COVID-19" ~ "COVID",
    inf %in% c("Influenza_A", "Influenza_B") ~ "Influenza",
    TRUE ~ NA_character_
  )
  
  ph_task <- ph %>% filter(!is.na(label))
  count_task <- count_mat2[, ph_task$sample_id, drop = FALSE]
  ph_task <- ph_task[match(colnames(count_task), ph_task$sample_id), , drop = FALSE]
  
  y <- DGEList(counts = count_task)
  keep <- filterByExpr(y, group = ph_task$label)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")
  xmat <- cpm(y, log = TRUE, prior.count = 1)
  
  list(pheno_task = ph_task, xmat = xmat)
}

dev_obj <- prepare_gse282464_task(pheno2, count_mat2)
dev_pheno <- dev_obj$pheno_task
dev_xmat <- dev_obj$xmat

cat("Development label distribution:\n")
print(table(dev_pheno$label))

############################################################
# 4) FINAL LOCKED ELASTIC NET MODEL
############################################################

x_dev_full <- t(dev_xmat)
y_dev_full <- ifelse(dev_pheno$label == "COVID", 1, 0)

top_k_external <- 100
gene_var_dev <- apply(x_dev_full, 2, var, na.rm = TRUE)
selected_genes_external <- names(sort(gene_var_dev, decreasing = TRUE))[1:min(top_k_external, length(gene_var_dev))]

x_dev_fs <- x_dev_full[, selected_genes_external, drop = FALSE]

scaler_external <- fit_train_scaler(x_dev_fs)
x_dev_fs_sc <- apply_scaler(x_dev_fs, scaler_external)

nfolds_use_final <- max(2, min(5, min(table(y_dev_full))))

cv_fit_external <- cv.glmnet(
  x = as.matrix(x_dev_fs_sc),
  y = y_dev_full,
  family = "binomial",
  alpha = 0.5,
  nfolds = nfolds_use_final,
  type.measure = "auc"
)

final_enet_external <- glmnet(
  x = as.matrix(x_dev_fs_sc),
  y = y_dev_full,
  family = "binomial",
  alpha = 0.5,
  lambda = cv_fit_external$lambda.min
)

final_coef_external <- as.matrix(coef(final_enet_external))
final_coef_external_df <- data.frame(
  gene = rownames(final_coef_external),
  coefficient = as.numeric(final_coef_external[, 1]),
  stringsAsFactors = FALSE
) %>%
  filter(gene != "(Intercept)") %>%
  filter(coefficient != 0) %>%
  mutate(abs_coef = abs(coefficient)) %>%
  arrange(desc(abs_coef))

write.csv(
  final_coef_external_df,
  "GSE282464_COVID_vs_Influenza_Final_ElasticNet_Coefficients.csv",
  row.names = FALSE
)

############################################################
# 5) DOWNLOAD GSE161731 SUPPLEMENTARY FILES
############################################################

gse_ext_id <- "GSE161731"
dir.create(gse_ext_id, showWarnings = FALSE)

getGEOSuppFiles(gse_ext_id, baseDir = ".", fetch_files = TRUE)

supp_dir_candidates <- c(
  file.path(gse_ext_id),
  file.path(gse_ext_id, "suppl")
)

supp_dir <- supp_dir_candidates[dir.exists(supp_dir_candidates)][1]
if (is.na(supp_dir)) stop("Supplementary directory bulunamad??.")

supp_files <- list.files(supp_dir, full.names = TRUE)
cat("Supplementary files:\n")
print(basename(supp_files))

xpr_file <- supp_files[grepl("xpr_nlcpm\\.csv\\.gz$", basename(supp_files), ignore.case = TRUE)][1]
if (is.na(xpr_file) || !file.exists(xpr_file)) {
  stop("GSE161731_xpr_nlcpm.csv.gz bulunamad??.")
}

cat("Using expression file:", basename(xpr_file), "\n")

############################################################
# 6) READ EXPRESSION MATRIX
############################################################

expr_df <- fread(xpr_file, data.table = FALSE)

cat("Imported expression dim:", dim(expr_df), "\n")
cat("First columns:\n")
print(colnames(expr_df)[1:min(10, ncol(expr_df))])

gene_col_idx <- which(tolower(colnames(expr_df)) %in% c(
  "gene", "gene_id", "geneid", "symbol", "ensembl_gene_id", "x"
))

if (length(gene_col_idx) == 0) gene_col_idx <- 1 else gene_col_idx <- gene_col_idx[1]

gene_ids <- as.character(expr_df[[gene_col_idx]])
expr_df2 <- expr_df[, -gene_col_idx, drop = FALSE]

expr_ext <- as.matrix(sapply(expr_df2, as.numeric))
rownames(expr_ext) <- gene_ids

keep_rows <- !is.na(rownames(expr_ext)) & rownames(expr_ext) != ""
expr_ext <- expr_ext[keep_rows, , drop = FALSE]

cat("External matrix dim:", dim(expr_ext), "\n")
cat("Example external rownames:\n")
print(head(rownames(expr_ext), 10))

############################################################
# 7) DEFINE EXTERNAL COVID VS INFLUENZA LABELS FROM COLUMN NAMES
############################################################

col_txt <- tolower(colnames(expr_ext))

external_label <- dplyr::case_when(
  str_detect(col_txt, "covid") ~ "COVID",
  str_detect(col_txt, "influenza") ~ "Influenza",
  TRUE ~ NA_character_
)

cat("Raw external label distribution from column names:\n")
print(table(external_label, useNA = "ifany"))

keep_cols <- !is.na(external_label)
expr_ext2 <- expr_ext[, keep_cols, drop = FALSE]
external_label2 <- external_label[keep_cols]

if (ncol(expr_ext2) == 0) {
  stop("Column names i??inde COVID / influenza etiketleri bulunamad??.")
}

cat("Filtered external matrix dim:", dim(expr_ext2), "\n")
print(table(external_label2))

############################################################
# 8) HARMONIZE IDS
############################################################

mode(expr_ext2) <- "numeric"
rownames(expr_ext2) <- sub("\\..*$", "", rownames(expr_ext2))
expr_ext2 <- collapse_duplicate_genes(expr_ext2)

external_is_ensg <- sum(grepl("^ENSG", rownames(expr_ext2))) > 0
cat("External rownames are ENSG:", external_is_ensg, "\n")

model_feature_order <- colnames(x_dev_fs)
model_feature_order_clean <- sub("\\..*$", "", model_feature_order)

if (external_is_ensg) {
  feature_map_df <- data.frame(
    model_ensembl = model_feature_order,
    model_key = model_feature_order_clean,
    stringsAsFactors = FALSE
  )
} else {
  map_dev_full <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = unique(model_feature_order_clean),
    keytype = "ENSEMBL",
    columns = c("ENSEMBL", "SYMBOL")
  )
  
  map_dev_full <- map_dev_full[!is.na(map_dev_full$SYMBOL) & map_dev_full$SYMBOL != "", ]
  map_dev_full <- map_dev_full[!duplicated(map_dev_full$ENSEMBL), ]
  
  model_feature_symbols <- map_dev_full$SYMBOL[match(model_feature_order_clean, map_dev_full$ENSEMBL)]
  
  feature_map_df <- data.frame(
    model_ensembl = model_feature_order,
    model_key = model_feature_symbols,
    stringsAsFactors = FALSE
  )
}

write.csv(
  feature_map_df,
  "GSE161731_ExternalValidation_ModelFeatureMap.csv",
  row.names = FALSE
)

############################################################
# 9) BUILD FULL EXTERNAL MODEL MATRIX
############################################################

n_ext <- ncol(expr_ext2)
x_ext_full <- matrix(
  NA_real_,
  nrow = n_ext,
  ncol = length(model_feature_order)
)

colnames(x_ext_full) <- model_feature_order
rownames(x_ext_full) <- colnames(expr_ext2)

matched_count <- 0

for (i in seq_len(nrow(feature_map_df))) {
  ens_col <- feature_map_df$model_ensembl[i]
  key_row <- feature_map_df$model_key[i]
  
  if (!is.na(key_row) && key_row != "" && key_row %in% rownames(expr_ext2)) {
    x_ext_full[, ens_col] <- as.numeric(expr_ext2[key_row, rownames(x_ext_full)])
    matched_count <- matched_count + 1
  }
}

cat("Matched model features in external data:", matched_count, "of", ncol(x_ext_full), "\n")

for (j in seq_len(ncol(x_ext_full))) {
  feat <- colnames(x_ext_full)[j]
  miss_idx <- is.na(x_ext_full[, j])
  
  if (any(miss_idx)) {
    x_ext_full[miss_idx, j] <- scaler_external$center[feat]
  }
}

scaler_locked <- list(
  center = scaler_external$center[model_feature_order],
  scale  = scaler_external$scale[model_feature_order]
)

x_ext_sc <- apply_scaler(x_ext_full, scaler_locked)

############################################################
# 10) PREDICT
############################################################

y_ext_fac <- factor(external_label2, levels = c("Influenza", "COVID"))

ext_prob <- as.numeric(
  predict(
    final_enet_external,
    newx = as.matrix(x_ext_sc),
    type = "response"
  )
)

external_predictions_df <- data.frame(
  sample_id = colnames(expr_ext2),
  truth = as.character(y_ext_fac),
  prob_covid = ext_prob,
  pred_label_05 = ifelse(ext_prob >= 0.5, "COVID", "Influenza"),
  stringsAsFactors = FALSE
)

write.csv(
  external_predictions_df,
  "GSE161731_ExternalValidation_Predictions.csv",
  row.names = FALSE
)

############################################################
# 11) METRICS
############################################################

external_metrics_df <- calc_external_metrics(
  y_true_fac = y_ext_fac,
  prob = ext_prob,
  threshold = 0.5,
  positive_class = "COVID"
)

write.csv(
  external_metrics_df,
  "GSE161731_ExternalValidation_Metrics.csv",
  row.names = FALSE
)

cat("\nExternal validation metrics at 0.5:\n")
print(external_metrics_df)

############################################################
# 12) YOUDEN THRESHOLD
############################################################

truth01_ext <- ifelse(y_ext_fac == "COVID", 1, 0)

roc_obj_ext <- pROC::roc(
  response = truth01_ext,
  predictor = ext_prob,
  quiet = TRUE
)

auc_ext <- as.numeric(pROC::auc(roc_obj_ext))
auc_ci_ext <- as.numeric(pROC::ci.auc(roc_obj_ext))

youden_coords <- pROC::coords(
  roc_obj_ext,
  x = "best",
  best.method = "youden",
  ret = c("threshold", "sensitivity", "specificity"),
  transpose = FALSE
)

best_threshold <- as.numeric(youden_coords["threshold"])

external_metrics_youden_df <- calc_external_metrics(
  y_true_fac = y_ext_fac,
  prob = ext_prob,
  threshold = best_threshold,
  positive_class = "COVID"
)

external_metrics_youden_df$threshold <- best_threshold

write.csv(
  external_metrics_youden_df,
  "GSE161731_ExternalValidation_Metrics_YoudenThreshold.csv",
  row.names = FALSE
)

############################################################
# 13) ROC FIGURE
############################################################

best_sens <- as.numeric(youden_coords["sensitivity"])
best_spec <- as.numeric(youden_coords["specificity"])
best_fpr <- 1 - best_spec

roc_df_ext <- data.frame(
  fpr = 1 - roc_obj_ext$specificities,
  tpr = roc_obj_ext$sensitivities
)

p_ext_roc_youden <- ggplot(roc_df_ext, aes(x = fpr, y = tpr)) +
  geom_line(linewidth = 1.1) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, color = "grey40") +
  geom_point(aes(x = best_fpr, y = best_sens), size = 3) +
  annotate(
    "text",
    x = best_fpr,
    y = best_sens,
    label = paste0("Youden threshold = ", sprintf("%.3f", best_threshold)),
    vjust = -1,
    size = 4
  ) +
  labs(
    title = "External validation ROC: GSE161731",
    subtitle = paste0(
      "Elastic Net | AUC = ",
      sprintf("%.3f", auc_ext),
      " (95% CI ",
      sprintf("%.3f", auc_ci_ext[1]), " - ",
      sprintf("%.3f", auc_ci_ext[3]), ")"
    ),
    x = "False Positive Rate",
    y = "True Positive Rate"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

ggsave(
  "Figure_ExternalValidation_ROC_GSE161731_Youden.png",
  p_ext_roc_youden,
  width = 7,
  height = 5.5,
  dpi = 600
)

############################################################
# 14) DCA
############################################################

dca_df <- data.frame(
  outcome = ifelse(y_ext_fac == "COVID", 1, 0),
  pred = ext_prob
)

dca_model <- decision_curve(
  outcome ~ pred,
  data = dca_df,
  family = binomial(link = "logit"),
  thresholds = seq(0.01, 0.75, by = 0.01),
  confidence.intervals = 0.95,
  study.design = "cohort"
)

pdf("Figure_DecisionCurve_GSE161731.pdf", width = 8, height = 6)
plot_decision_curve(
  dca_model,
  curve.names = "Elastic Net",
  xlab = "Threshold Probability",
  ylab = "Net Benefit",
  standardize = FALSE
)
dev.off()

############################################################
# 15) COEFFICIENT BIOLOGY
############################################################

coef_tbl <- final_coef_external_df
coef_tbl$ensembl_clean <- sub("\\..*$", "", coef_tbl$gene)

coef_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(coef_tbl$ensembl_clean),
  keytype = "ENSEMBL",
  columns = c("ENSEMBL", "SYMBOL", "GENENAME")
)

coef_map <- coef_map[!duplicated(coef_map$ENSEMBL), ]

coef_tbl$symbol <- coef_map$SYMBOL[match(coef_tbl$ensembl_clean, coef_map$ENSEMBL)]
coef_tbl$genename <- coef_map$GENENAME[match(coef_tbl$ensembl_clean, coef_map$ENSEMBL)]

coef_tbl <- coef_tbl %>%
  mutate(direction = ifelse(coefficient > 0, "Higher in COVID", "Higher in Influenza")) %>%
  arrange(desc(abs_coef))

write.csv(
  coef_tbl,
  "GSE161731_Final_ElasticNet_Coefficients_Annotated.csv",
  row.names = FALSE
)

############################################################
# 16) SUMMARY
############################################################

summary_df <- data.frame(
  development_dataset = "GSE282464",
  development_task = "COVID vs Influenza",
  external_dataset = "GSE161731",
  model = "ElasticNet",
  n_external = nrow(external_predictions_df),
  n_external_covid = sum(external_predictions_df$truth == "COVID"),
  n_external_influenza = sum(external_predictions_df$truth == "Influenza"),
  auc = external_metrics_df$auc,
  auc_ci_low = external_metrics_df$auc_ci_low,
  auc_ci_high = external_metrics_df$auc_ci_high,
  accuracy = external_metrics_df$accuracy,
  balanced_accuracy = external_metrics_df$balanced_accuracy,
  sensitivity = external_metrics_df$sensitivity,
  specificity = external_metrics_df$specificity,
  precision = external_metrics_df$precision,
  f1 = external_metrics_df$f1,
  mcc = external_metrics_df$mcc,
  brier = external_metrics_df$brier,
  youden_threshold = best_threshold,
  stringsAsFactors = FALSE
)

write.csv(
  summary_df,
  "Paper_Ready_ExternalValidation_Summary_GSE161731.csv",
  row.names = FALSE
)

cat("\nFinished external validation for GSE161731\n")
print(summary_df)





rm(prepare_gse282464_task)










prepare_gse282464_task <- function(pheno2, count_mat2,
                                   task_name = c("covid_vs_influenza",
                                                 "covid_vs_noninfluenza_viral",
                                                 "covid_vs_sepsis")) {
  task_name <- match.arg(task_name)
  
  ph <- pheno2
  ph$subject_id2 <- ph$sample_id
  inf <- as.character(ph$`type of_infection:ch1`)
  
  if (task_name == "covid_vs_influenza") {
    ph$label <- dplyr::case_when(
      inf == "COVID-19" ~ "COVID",
      inf %in% c("Influenza_A", "Influenza_B") ~ "Influenza",
      TRUE ~ NA_character_
    )
  }
  
  if (task_name == "covid_vs_noninfluenza_viral") {
    ph$label <- dplyr::case_when(
      inf == "COVID-19" ~ "COVID",
      inf %in% c("Coronavirus",
                 "Co-infection (Viral/Viral)",
                 "Co-infection (viral/viral)",
                 "Co-infection (Viral/Fungal)") ~ "NonInfluenza_Viral",
      TRUE ~ NA_character_
    )
  }
  
  if (task_name == "covid_vs_sepsis") {
    ph$label <- dplyr::case_when(
      inf == "COVID-19" ~ "COVID",
      inf %in% c("Sepsis", "Septic_shock") ~ "Sepsis_Group",
      TRUE ~ NA_character_
    )
  }
  
  ph_task <- ph %>% dplyr::filter(!is.na(label))
  count_task <- count_mat2[, ph_task$sample_id, drop = FALSE]
  ph_task <- ph_task[match(colnames(count_task), ph_task$sample_id), , drop = FALSE]
  
  stopifnot(all(colnames(count_task) == ph_task$sample_id))
  
  y <- edgeR::DGEList(counts = count_task)
  keep <- edgeR::filterByExpr(y, group = ph_task$label)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- edgeR::calcNormFactors(y, method = "TMM")
  xmat <- edgeR::cpm(y, log = TRUE, prior.count = 1)
  
  list(
    pheno_task = ph_task,
    xmat = xmat
  )
}

args(prepare_gse282464_task)









task_names_new <- c(
  "covid_vs_influenza",
  "covid_vs_noninfluenza_viral",
  "covid_vs_sepsis"
)

task_objects <- lapply(task_names_new, function(tt) {
  prepare_gse282464_task(pheno2, count_mat2, task_name = tt)
})
names(task_objects) <- task_names_new







############################################################
# MINIMUM FUNCTIONS TO REBUILD res_sepsis AND SAVE HEATMAP
############################################################

calc_auc <- function(y_true, prob) {
  if (length(unique(y_true)) < 2) return(NA_real_)
  as.numeric(pROC::roc(y_true, prob, quiet = TRUE)$auc)
}

calc_pr_auc <- function(y_true, prob) {
  if (length(unique(y_true)) < 2) return(NA_real_)
  if (sum(y_true == 1) == 0) return(NA_real_)
  
  ord <- order(prob, decreasing = TRUE)
  y_sorted <- y_true[ord]
  
  tp_cum <- cumsum(y_sorted == 1)
  fp_cum <- cumsum(y_sorted == 0)
  
  precision_curve <- tp_cum / (tp_cum + fp_cum)
  recall_curve <- tp_cum / sum(y_true == 1)
  
  precision_curve <- c(1, precision_curve)
  recall_curve <- c(0, recall_curve)
  
  sum(
    (recall_curve[-1] - recall_curve[-length(recall_curve)]) *
      (precision_curve[-1] + precision_curve[-length(precision_curve)]) / 2
  )
}

calc_metrics <- function(y_true, prob, threshold = 0.5) {
  pred <- ifelse(prob >= threshold, 1, 0)
  
  tp <- sum(pred == 1 & y_true == 1)
  tn <- sum(pred == 0 & y_true == 0)
  fp <- sum(pred == 1 & y_true == 0)
  fn <- sum(pred == 0 & y_true == 1)
  
  accuracy <- (tp + tn) / length(y_true)
  sensitivity <- ifelse((tp + fn) == 0, NA_real_, tp / (tp + fn))
  specificity <- ifelse((tn + fp) == 0, NA_real_, tn / (tn + fp))
  precision <- ifelse((tp + fp) == 0, NA_real_, tp / (tp + fp))
  recall <- sensitivity
  
  f1 <- ifelse(
    is.na(precision) || is.na(recall) || (precision + recall) == 0,
    NA_real_,
    2 * precision * recall / (precision + recall)
  )
  
  balanced_accuracy <- mean(c(sensitivity, specificity), na.rm = TRUE)
  auc <- calc_auc(y_true, prob)
  pr_auc <- calc_pr_auc(y_true, prob)
  
  mcc_denom <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  mcc <- ifelse(mcc_denom == 0, NA_real_,
                ((tp * tn) - (fp * fn)) / mcc_denom)
  
  brier <- mean((prob - y_true)^2)
  
  data.frame(
    auc = auc,
    pr_auc = pr_auc,
    accuracy = accuracy,
    balanced_accuracy = balanced_accuracy,
    sensitivity = sensitivity,
    specificity = specificity,
    precision = precision,
    f1 = f1,
    mcc = mcc,
    brier = brier
  )
}

make_subject_folds <- function(subject_table,
                               label_col = "label",
                               subject_col = "subject_id2",
                               k = 5,
                               seed = 123) {
  set.seed(seed)
  
  classes <- unique(subject_table[[label_col]])
  if (length(classes) != 2) stop("Binary classification expected.")
  
  out <- vector("list", k)
  for (i in seq_len(k)) out[[i]] <- character(0)
  
  for (cl in classes) {
    subj <- subject_table %>%
      dplyr::filter(.data[[label_col]] == cl) %>%
      dplyr::pull(.data[[subject_col]]) %>%
      unique()
    
    subj <- sample(subj)
    
    for (j in seq_along(subj)) {
      fold_id <- ((j - 1) %% k) + 1
      out[[fold_id]] <- c(out[[fold_id]], subj[j])
    }
  }
  
  out
}

fit_one_model <- function(model_name, x_train_fs, x_test_fs, y_train, y_test) {
  x_train_fs <- as.matrix(x_train_fs)
  x_test_fs  <- as.matrix(x_test_fs)
  
  if (nrow(x_train_fs) == 0 || nrow(x_test_fs) == 0) return(NULL)
  if (length(unique(y_train)) < 2 || length(unique(y_test)) < 2) return(NULL)
  
  if (model_name == "ElasticNet") {
    prob <- tryCatch({
      nfolds_use <- max(2, min(5, min(table(y_train))))
      
      cv_fit <- glmnet::cv.glmnet(
        x = as.matrix(x_train_fs),
        y = y_train,
        family = "binomial",
        alpha = 0.5,
        nfolds = nfolds_use,
        type.measure = "auc"
      )
      
      as.numeric(
        predict(cv_fit, newx = as.matrix(x_test_fs), s = "lambda.min", type = "response")
      )
    }, error = function(e) NULL)
    return(prob)
  }
  
  if (model_name == "XGBoost") {
    prob <- tryCatch({
      dtrain_xgb <- xgboost::xgb.DMatrix(data = x_train_fs, label = y_train)
      dtest_xgb  <- xgboost::xgb.DMatrix(data = x_test_fs, label = y_test)
      
      params_xgb <- list(
        objective = "binary:logistic",
        eval_metric = "auc",
        eta = 0.03,
        max_depth = 4,
        subsample = 1.0,
        colsample_bytree = 0.8,
        min_child_weight = 1
      )
      
      bst_xgb <- xgboost::xgb.train(
        params = params_xgb,
        data = dtrain_xgb,
        nrounds = 120,
        verbose = 0
      )
      
      as.numeric(predict(bst_xgb, dtest_xgb))
    }, error = function(e) NULL)
    return(prob)
  }
  
  if (model_name == "RandomForest") {
    prob <- tryCatch({
      df_train_rf <- as.data.frame(x_train_fs)
      df_test_rf  <- as.data.frame(x_test_fs)
      
      y_train_fac <- factor(ifelse(y_train == 1, "Case", "Control"),
                            levels = c("Control", "Case"))
      
      rf_model <- randomForest::randomForest(
        x = df_train_rf,
        y = y_train_fac,
        ntree = 300,
        importance = TRUE
      )
      
      as.numeric(predict(rf_model, df_test_rf, type = "prob")[, "Case"])
    }, error = function(e) NULL)
    return(prob)
  }
  
  if (model_name == "GPBoost_LightGBM") {
    prob <- tryCatch({
      dtrain_plain <- gpboost::gpb.Dataset(data = x_train_fs, label = y_train)
      
      params_plain <- list(
        objective = "binary",
        learning_rate = 0.03,
        max_depth = 4,
        num_leaves = 31,
        min_data_in_leaf = 5,
        feature_fraction = 0.8,
        bagging_fraction = 1.0,
        bagging_freq = 0,
        verbose = 0
      )
      
      bst_plain <- gpboost::gpb.train(
        params = params_plain,
        data = dtrain_plain,
        nrounds = 120,
        verbose = 0
      )
      
      as.numeric(predict(bst_plain, x_test_fs))
    }, error = function(e) NULL)
    return(prob)
  }
  
  stop("Unknown model name.")
}

run_task_full_benchmark <- function(task_obj,
                                    dataset_name,
                                    task_name,
                                    seed = 123,
                                    top_k = 100,
                                    k = 5) {
  pheno_task <- task_obj$pheno_task
  xmat <- task_obj$xmat
  
  subject_table <- pheno_task %>%
    dplyr::distinct(subject_id2, label)
  
  fold_subjects <- make_subject_folds(
    subject_table,
    label_col = "label",
    subject_col = "subject_id2",
    k = k,
    seed = seed
  )
  
  positive_class <- "COVID"
  model_list <- c("GPBoost_LightGBM", "ElasticNet", "XGBoost", "RandomForest")
  
  cv_results <- data.frame(
    dataset = character(),
    task = character(),
    fold = integer(),
    model = character(),
    n_train = integer(),
    n_test = integer(),
    n_test_case = integer(),
    n_test_control = integer(),
    auc = numeric(),
    pr_auc = numeric(),
    accuracy = numeric(),
    balanced_accuracy = numeric(),
    sensitivity = numeric(),
    specificity = numeric(),
    precision = numeric(),
    f1 = numeric(),
    mcc = numeric(),
    brier = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (fold_i in seq_len(k)) {
    test_subjects <- fold_subjects[[fold_i]]
    if (length(test_subjects) == 0) next
    
    test_idx <- pheno_task$subject_id2 %in% test_subjects
    train_idx <- !test_idx
    
    if (sum(test_idx) == 0 || sum(train_idx) == 0) next
    
    x_train <- t(xmat[, train_idx, drop = FALSE])
    x_test  <- t(xmat[, test_idx, drop = FALSE])
    
    y_train <- ifelse(pheno_task$label[train_idx] == positive_class, 1, 0)
    y_test  <- ifelse(pheno_task$label[test_idx] == positive_class, 1, 0)
    
    if (length(unique(y_train)) < 2) next
    if (length(unique(y_test)) < 2) next
    
    gene_var <- apply(x_train, 2, var, na.rm = TRUE)
    top_genes <- names(sort(gene_var, decreasing = TRUE))[1:min(top_k, length(gene_var))]
    
    x_train_fs <- x_train[, top_genes, drop = FALSE]
    x_test_fs  <- x_test[, top_genes, drop = FALSE]
    
    for (m in model_list) {
      prob <- fit_one_model(
        model_name = m,
        x_train_fs = x_train_fs,
        x_test_fs = x_test_fs,
        y_train = y_train,
        y_test = y_test
      )
      
      if (is.null(prob)) next
      if (length(prob) != nrow(x_test_fs)) next
      if (!all(is.finite(prob))) next
      
      met <- calc_metrics(y_test, prob)
      
      cv_results <- rbind(
        cv_results,
        data.frame(
          dataset = dataset_name,
          task = task_name,
          fold = fold_i,
          model = m,
          n_train = nrow(x_train_fs),
          n_test = nrow(x_test_fs),
          n_test_case = sum(y_test == 1),
          n_test_control = sum(y_test == 0),
          met,
          stringsAsFactors = FALSE
        )
      )
    }
  }
  
  summary_results <- cv_results %>%
    dplyr::group_by(dataset, task, model) %>%
    dplyr::summarise(
      mean_auc = mean(auc, na.rm = TRUE),
      mean_pr_auc = mean(pr_auc, na.rm = TRUE),
      mean_acc = mean(accuracy, na.rm = TRUE),
      mean_bal_acc = mean(balanced_accuracy, na.rm = TRUE),
      mean_f1 = mean(f1, na.rm = TRUE),
      mean_brier = mean(brier, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(desc(mean_auc))
  
  list(
    cv_results = cv_results,
    summary_results = summary_results
  )
}

plot_heatmap_metrics <- function(summary_df, title_text, outfile) {
  heat_df <- summary_df %>%
    dplyr::select(model, mean_auc, mean_pr_auc, mean_acc, mean_bal_acc, mean_f1, mean_brier) %>%
    as.data.frame()
  
  rownames(heat_df) <- heat_df$model
  heat_df <- heat_df[, -1, drop = FALSE]
  
  png(outfile, width = 2200, height = 1600, res = 300)
  pheatmap::pheatmap(
    as.matrix(heat_df),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    scale = "none",
    color = colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(100),
    display_numbers = TRUE,
    number_format = "%.3f",
    fontsize = 12,
    fontsize_number = 10,
    border_color = "grey85",
    main = title_text
  )
  dev.off()
}

############################################################
# REBUILD res_sepsis AND SAVE HEATMAP
############################################################

task_names_new <- c(
  "covid_vs_influenza",
  "covid_vs_noninfluenza_viral",
  "covid_vs_sepsis"
)

task_objects <- lapply(task_names_new, function(tt) {
  prepare_gse282464_task(pheno2, count_mat2, task_name = tt)
})
names(task_objects) <- task_names_new

res_sepsis <- run_task_full_benchmark(
  task_obj = task_objects[["covid_vs_sepsis"]],
  dataset_name = "GSE282464",
  task_name = "covid_vs_sepsis",
  seed = 123,
  top_k = 100,
  k = 5
)

print(res_sepsis$summary_results)

plot_heatmap_metrics(
  res_sepsis$summary_results,
  "GSE282464: COVID-19 vs Sepsis / Septic Shock",
  "Figure_COVID_vs_Sepsis_Heatmap_Revised.png"
)

file.exists("Figure_COVID_vs_Sepsis_Heatmap_Revised.png")











############################################################
# FINAL WORKING SCRIPT
# Development: GSE282464
# External validation: GSE161731
# Task: COVID vs Influenza
# COUNTS-BASED VERSION
############################################################

rm(list = ls())
gc()

############################################################
# 0) PACKAGES
############################################################
cran_pkgs <- c(
  "data.table", "dplyr", "stringr", "ggplot2", "pROC",
  "glmnet", "GEOquery", "tibble", "rmda"
)
bioc_pkgs <- c("edgeR", "AnnotationDbi", "org.Hs.eg.db")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
for (p in bioc_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    BiocManager::install(p, ask = FALSE, update = FALSE)
  }
}

library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(pROC)
library(glmnet)
library(GEOquery)
library(tibble)
library(rmda)
library(edgeR)
library(AnnotationDbi)
library(org.Hs.eg.db)

set.seed(123)

############################################################
# 1) HELPERS
############################################################

safe_col_lower <- function(df) {
  names(df) <- tolower(names(df))
  df
}

clean_id <- function(x) {
  x <- as.character(x)
  x <- gsub("^X", "", x)
  x <- gsub("\\.gz$", "", x)
  x <- gsub("\\.csv$", "", x)
  x <- gsub("\\.txt$", "", x)
  x <- gsub("[^A-Za-z0-9]", "", x)
  tolower(x)
}

collapse_duplicate_genes <- function(expr_mat) {
  expr_df <- as.data.frame(expr_mat, check.names = FALSE)
  expr_df$gene <- rownames(expr_df)
  
  expr_df <- expr_df %>%
    group_by(gene) %>%
    summarise(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
  
  out <- as.matrix(expr_df[, -1, drop = FALSE])
  rownames(out) <- expr_df$gene
  out
}

fit_train_scaler <- function(x_train) {
  mu <- colMeans(x_train, na.rm = TRUE)
  sdv <- apply(x_train, 2, sd, na.rm = TRUE)
  sdv[is.na(sdv) | sdv == 0] <- 1
  list(center = mu, scale = sdv)
}

apply_scaler <- function(x, scaler) {
  sweep(sweep(x, 2, scaler$center, "-"), 2, scaler$scale, "/")
}

calc_external_metrics <- function(y_true_fac, prob, threshold = 0.5, positive_class = "COVID") {
  y_true <- ifelse(y_true_fac == positive_class, 1, 0)
  pred <- ifelse(prob >= threshold, 1, 0)
  
  tp <- sum(pred == 1 & y_true == 1)
  tn <- sum(pred == 0 & y_true == 0)
  fp <- sum(pred == 1 & y_true == 0)
  fn <- sum(pred == 0 & y_true == 1)
  
  accuracy <- (tp + tn) / length(y_true)
  sensitivity <- ifelse((tp + fn) == 0, NA_real_, tp / (tp + fn))
  specificity <- ifelse((tn + fp) == 0, NA_real_, tn / (tn + fp))
  precision <- ifelse((tp + fp) == 0, NA_real_, tp / (tp + fp))
  recall <- sensitivity
  
  f1 <- ifelse(
    is.na(precision) || is.na(recall) || (precision + recall) == 0,
    NA_real_,
    2 * precision * recall / (precision + recall)
  )
  
  balanced_accuracy <- mean(c(sensitivity, specificity), na.rm = TRUE)
  
  auc_obj <- tryCatch(
    pROC::roc(y_true, prob, quiet = TRUE),
    error = function(e) NULL
  )
  
  auc <- if (is.null(auc_obj)) NA_real_ else as.numeric(pROC::auc(auc_obj))
  auc_ci <- if (is.null(auc_obj)) c(NA_real_, NA_real_, NA_real_) else as.numeric(pROC::ci.auc(auc_obj))
  
  mcc_denom <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  mcc <- ifelse(mcc_denom == 0, NA_real_,
                ((tp * tn) - (fp * fn)) / mcc_denom)
  
  brier <- mean((prob - y_true)^2)
  
  data.frame(
    auc = auc,
    auc_ci_low = auc_ci[1],
    auc_ci_mid = auc_ci[2],
    auc_ci_high = auc_ci[3],
    accuracy = accuracy,
    balanced_accuracy = balanced_accuracy,
    sensitivity = sensitivity,
    specificity = specificity,
    precision = precision,
    f1 = f1,
    mcc = mcc,
    brier = brier,
    stringsAsFactors = FALSE
  )
}

label_covid_influenza_from_text <- function(x) {
  x <- tolower(as.character(x))
  dplyr::case_when(
    str_detect(x, "covid") ~ "COVID",
    str_detect(x, "influenza") ~ "Influenza",
    TRUE ~ NA_character_
  )
}

############################################################
# 2) BUILD GSE282464 DEVELOPMENT DATA
############################################################

gse_id <- "GSE282464"
untar_dir <- file.path(gse_id, "untarred_raw")

if (!dir.exists(untar_dir)) {
  stop("GSE282464/untarred_raw klas??r?? bulunamad??.")
}

count_files <- list.files(untar_dir, pattern = "\\.counts\\.gz$", full.names = TRUE)
if (length(count_files) == 0) stop(".counts.gz dosyas?? bulunamad??.")

read_one_count <- function(f) {
  dt <- fread(f, header = FALSE)
  dt <- dt[, 1:2]
  sample_id <- sub("_.*$", "", basename(f))
  colnames(dt) <- c("gene_id", sample_id)
  dt
}

count_list <- lapply(count_files, read_one_count)
count_merged <- Reduce(function(x, y) merge(x, y, by = "gene_id", all = TRUE), count_list)
count_merged[is.na(count_merged)] <- 0

count_mat <- as.data.frame(count_merged)
rownames(count_mat) <- count_mat$gene_id
count_mat <- count_mat[, -1, drop = FALSE]
count_mat <- as.matrix(sapply(count_mat, as.numeric))
rownames(count_mat) <- count_merged$gene_id

gset_list <- getGEO(gse_id, GSEMatrix = TRUE)
gset <- gset_list[[1]]
pheno <- pData(gset)

pheno$sample_id <- as.character(pheno$geo_accession)
common_samples <- intersect(colnames(count_mat), pheno$sample_id)
if (length(common_samples) == 0) stop("Metadata ile count matrix e??le??medi.")

pheno2 <- pheno %>% filter(sample_id %in% common_samples)
count_mat2 <- count_mat[, pheno2$sample_id, drop = FALSE]
pheno2 <- pheno2[match(colnames(count_mat2), pheno2$sample_id), , drop = FALSE]

inf <- as.character(pheno2$`type of_infection:ch1`)
pheno2$label <- dplyr::case_when(
  inf == "COVID-19" ~ "COVID",
  inf %in% c("Influenza_A", "Influenza_B") ~ "Influenza",
  TRUE ~ NA_character_
)

dev_pheno <- pheno2 %>% filter(!is.na(label))
dev_counts <- count_mat2[, dev_pheno$sample_id, drop = FALSE]
dev_pheno <- dev_pheno[match(colnames(dev_counts), dev_pheno$sample_id), , drop = FALSE]

y <- DGEList(counts = dev_counts)
keep <- filterByExpr(y, group = dev_pheno$label)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y, method = "TMM")
dev_xmat <- cpm(y, log = TRUE, prior.count = 1)

############################################################
# 3) FINAL LOCKED ELASTIC NET MODEL
############################################################

x_dev_full <- t(dev_xmat)
y_dev_full <- ifelse(dev_pheno$label == "COVID", 1, 0)

top_k_external <- 100
gene_var_dev <- apply(x_dev_full, 2, var, na.rm = TRUE)
selected_genes_external <- names(sort(gene_var_dev, decreasing = TRUE))[1:min(top_k_external, length(gene_var_dev))]
x_dev_fs <- x_dev_full[, selected_genes_external, drop = FALSE]

scaler_external <- fit_train_scaler(x_dev_fs)
x_dev_fs_sc <- apply_scaler(x_dev_fs, scaler_external)

nfolds_use_final <- max(2, min(5, min(table(y_dev_full))))

cv_fit_external <- cv.glmnet(
  x = as.matrix(x_dev_fs_sc),
  y = y_dev_full,
  family = "binomial",
  alpha = 0.5,
  nfolds = nfolds_use_final,
  type.measure = "auc"
)

final_enet_external <- glmnet(
  x = as.matrix(x_dev_fs_sc),
  y = y_dev_full,
  family = "binomial",
  alpha = 0.5,
  lambda = cv_fit_external$lambda.min
)

############################################################
# 4) DOWNLOAD GSE161731 FILES
############################################################

gse_ext_id <- "GSE161731"
dir.create(gse_ext_id, showWarnings = FALSE)
getGEOSuppFiles(gse_ext_id, baseDir = ".", fetch_files = TRUE)

all_files <- list.files(gse_ext_id, recursive = TRUE, full.names = TRUE)

counts_file <- all_files[grepl("GSE161731_counts\\.csv\\.gz$", all_files)][1]
counts_key_file <- all_files[grepl("GSE161731_counts_key\\.csv\\.gz$", all_files)][1]
key_file <- all_files[grepl("GSE161731_key\\.csv\\.gz$", all_files)][1]

if (is.na(counts_file) || !file.exists(counts_file)) stop("GSE161731_counts.csv.gz bulunamad??.")
if (is.na(counts_key_file) || !file.exists(counts_key_file)) stop("GSE161731_counts_key.csv.gz bulunamad??.")
if (is.na(key_file) || !file.exists(key_file)) stop("GSE161731_key.csv.gz bulunamad??.")

cat("Using files:\n")
print(c(counts_file, counts_key_file, key_file))

############################################################
# 5) READ COUNTS MATRIX
############################################################

counts_ext <- fread(counts_file, data.table = FALSE)
cat("counts_ext dim:", dim(counts_ext), "\n")
print(colnames(counts_ext)[1:min(10, ncol(counts_ext))])

gene_col_idx <- 1
gene_ids <- as.character(counts_ext[[gene_col_idx]])
counts_ext2 <- counts_ext[, -gene_col_idx, drop = FALSE]

expr_ext_counts <- as.matrix(sapply(counts_ext2, as.numeric))
rownames(expr_ext_counts) <- gene_ids

keep_rows <- !is.na(rownames(expr_ext_counts)) & rownames(expr_ext_counts) != ""
expr_ext_counts <- expr_ext_counts[keep_rows, , drop = FALSE]

cat("External counts matrix dim:", dim(expr_ext_counts), "\n")
cat("Example external count colnames:\n")
print(head(colnames(expr_ext_counts), 20))

############################################################
# 6) READ KEY FILES AND BUILD LABEL MAP
############################################################

counts_key <- fread(counts_key_file, data.table = FALSE)
counts_key <- safe_col_lower(as.data.frame(counts_key))

meta_key <- fread(key_file, data.table = FALSE)
meta_key <- safe_col_lower(as.data.frame(meta_key))

cat("counts_key columns:\n")
print(names(counts_key))
cat("meta_key columns:\n")
print(names(meta_key))

# Determine sample ID column in counts_key
if ("rna_id" %in% names(counts_key)) {
  sample_col_ck <- "rna_id"
} else {
  sample_col_candidates <- names(counts_key)[grepl("rna|sample|id", names(counts_key))]
  sample_col_ck <- sample_col_candidates[1]
}

if (is.na(sample_col_ck) || !sample_col_ck %in% names(counts_key)) {
  stop("counts_key i??inde sample/rna_id kolonu bulunamad??.")
}

# Determine label column
label_col_ck <- names(counts_key)[grepl("cohort|group|diagnosis|disease|status|infection", names(counts_key))][1]

if (!is.na(label_col_ck) && label_col_ck %in% names(counts_key)) {
  counts_key$label_binary <- label_covid_influenza_from_text(counts_key[[label_col_ck]])
} else {
  counts_key$label_binary <- NA_character_
}

# If counts_key labels weak, try joining meta_key through subject_id
if (sum(!is.na(counts_key$label_binary)) == 0 && "subject_id" %in% names(counts_key) && "subject_id" %in% names(meta_key)) {
  meta_label_col <- names(meta_key)[grepl("cohort|group|diagnosis|disease|status|infection", names(meta_key))][1]
  if (!is.na(meta_label_col) && meta_label_col %in% names(meta_key)) {
    meta_key$label_binary <- label_covid_influenza_from_text(meta_key[[meta_label_col]])
    counts_key <- counts_key %>%
      left_join(meta_key %>% select(subject_id, label_binary), by = "subject_id", suffix = c("", ".meta")) %>%
      mutate(label_binary = ifelse(is.na(label_binary), label_binary.meta, label_binary)) %>%
      select(-label_binary.meta)
  }
}

cat("Label distribution in counts_key after mapping:\n")
print(table(counts_key$label_binary, useNA = "ifany"))

counts_key$sample_clean <- clean_id(counts_key[[sample_col_ck]])
counts_key_bin <- counts_key %>%
  filter(!is.na(label_binary)) %>%
  distinct(sample_clean, .keep_all = TRUE)

############################################################
# 7) MATCH COUNTS MATRIX COLUMNS TO COUNTS_KEY
############################################################

expr_col_raw <- colnames(expr_ext_counts)
expr_col_clean <- clean_id(expr_col_raw)

match_idx <- match(expr_col_clean, counts_key_bin$sample_clean)

cat("Matched external columns:", sum(!is.na(match_idx)), "of", length(expr_col_raw), "\n")

if (sum(!is.na(match_idx)) == 0) {
  cat("Example expr columns:\n")
  print(head(expr_col_raw, 20))
  cat("Example cleaned expr columns:\n")
  print(head(expr_col_clean, 20))
  cat("Example counts_key sample IDs:\n")
  print(head(counts_key_bin[[sample_col_ck]], 20))
  cat("Example cleaned counts_key sample IDs:\n")
  print(head(counts_key_bin$sample_clean, 20))
  stop("GSE161731 counts matrix ile counts_key e??le??medi.")
}

keep_cols <- !is.na(match_idx)
ext_counts_matched <- expr_ext_counts[, keep_cols, drop = FALSE]
counts_key_matched <- counts_key_bin[match_idx[keep_cols], , drop = FALSE]
external_label2 <- counts_key_matched$label_binary

cat("Matched external counts dim:", dim(ext_counts_matched), "\n")
cat("Matched external label distribution:\n")
print(table(external_label2))

############################################################
# 8) NORMALIZE EXTERNAL COUNTS WITH TMM + logCPM
############################################################

dge_ext <- DGEList(counts = ext_counts_matched)
keep_ext <- filterByExpr(dge_ext, group = external_label2)
dge_ext <- dge_ext[keep_ext, , keep.lib.sizes = FALSE]
dge_ext <- calcNormFactors(dge_ext, method = "TMM")
expr_ext2 <- cpm(dge_ext, log = TRUE, prior.count = 1)

mode(expr_ext2) <- "numeric"
rownames(expr_ext2) <- sub("\\..*$", "", rownames(expr_ext2))
expr_ext2 <- collapse_duplicate_genes(expr_ext2)

############################################################
# 9) HARMONIZE IDS
############################################################

external_is_ensg <- sum(grepl("^ENSG", rownames(expr_ext2))) > 0

model_feature_order <- colnames(x_dev_fs)
model_feature_order_clean <- sub("\\..*$", "", model_feature_order)

if (external_is_ensg) {
  feature_map_df <- data.frame(
    model_ensembl = model_feature_order,
    model_key = model_feature_order_clean,
    stringsAsFactors = FALSE
  )
} else {
  map_dev_full <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = unique(model_feature_order_clean),
    keytype = "ENSEMBL",
    columns = c("ENSEMBL", "SYMBOL")
  )
  map_dev_full <- map_dev_full[!is.na(map_dev_full$SYMBOL) & map_dev_full$SYMBOL != "", ]
  map_dev_full <- map_dev_full[!duplicated(map_dev_full$ENSEMBL), ]
  
  model_feature_symbols <- map_dev_full$SYMBOL[match(model_feature_order_clean, map_dev_full$ENSEMBL)]
  
  feature_map_df <- data.frame(
    model_ensembl = model_feature_order,
    model_key = model_feature_symbols,
    stringsAsFactors = FALSE
  )
}

############################################################
# 10) BUILD FULL EXTERNAL MODEL MATRIX
############################################################

n_ext <- ncol(expr_ext2)
x_ext_full <- matrix(NA_real_, nrow = n_ext, ncol = length(model_feature_order))
colnames(x_ext_full) <- model_feature_order
rownames(x_ext_full) <- colnames(expr_ext2)

matched_count <- 0
for (i in seq_len(nrow(feature_map_df))) {
  ens_col <- feature_map_df$model_ensembl[i]
  key_row <- feature_map_df$model_key[i]
  
  if (!is.na(key_row) && key_row != "" && key_row %in% rownames(expr_ext2)) {
    x_ext_full[, ens_col] <- as.numeric(expr_ext2[key_row, rownames(x_ext_full)])
    matched_count <- matched_count + 1
  }
}

cat("Matched model features in external data:", matched_count, "of", ncol(x_ext_full), "\n")

for (j in seq_len(ncol(x_ext_full))) {
  feat <- colnames(x_ext_full)[j]
  miss_idx <- is.na(x_ext_full[, j])
  if (any(miss_idx)) {
    x_ext_full[miss_idx, j] <- scaler_external$center[feat]
  }
}

scaler_locked <- list(
  center = scaler_external$center[model_feature_order],
  scale  = scaler_external$scale[model_feature_order]
)

x_ext_sc <- apply_scaler(x_ext_full, scaler_locked)

############################################################
# 11) PREDICT
############################################################

y_ext_fac <- factor(external_label2, levels = c("Influenza", "COVID"))

ext_prob <- as.numeric(
  predict(final_enet_external, newx = as.matrix(x_ext_sc), type = "response")
)

external_predictions_df <- data.frame(
  sample_id = colnames(ext_counts_matched),
  truth = as.character(y_ext_fac),
  prob_covid = ext_prob,
  pred_label_05 = ifelse(ext_prob >= 0.5, "COVID", "Influenza"),
  stringsAsFactors = FALSE
)

write.csv(external_predictions_df, "GSE161731_ExternalValidation_Predictions.csv", row.names = FALSE)

############################################################
# 12) METRICS
############################################################

external_metrics_df <- calc_external_metrics(
  y_true_fac = y_ext_fac,
  prob = ext_prob,
  threshold = 0.5,
  positive_class = "COVID"
)

write.csv(external_metrics_df, "GSE161731_ExternalValidation_Metrics.csv", row.names = FALSE)

cat("\nExternal validation metrics at 0.5:\n")
print(external_metrics_df)

############################################################
# 13) YOUDEN
############################################################

truth01_ext <- ifelse(y_ext_fac == "COVID", 1, 0)
roc_obj_ext <- pROC::roc(response = truth01_ext, predictor = ext_prob, quiet = TRUE)
auc_ext <- as.numeric(pROC::auc(roc_obj_ext))
auc_ci_ext <- as.numeric(pROC::ci.auc(roc_obj_ext))

youden_coords <- pROC::coords(
  roc_obj_ext,
  x = "best",
  best.method = "youden",
  ret = c("threshold", "sensitivity", "specificity"),
  transpose = FALSE
)

best_threshold <- as.numeric(youden_coords["threshold"])
best_sens <- as.numeric(youden_coords["sensitivity"])
best_spec <- as.numeric(youden_coords["specificity"])
best_fpr <- 1 - best_spec

external_metrics_youden_df <- calc_external_metrics(
  y_true_fac = y_ext_fac,
  prob = ext_prob,
  threshold = best_threshold,
  positive_class = "COVID"
)
external_metrics_youden_df$threshold <- best_threshold

write.csv(external_metrics_youden_df, "GSE161731_ExternalValidation_Metrics_YoudenThreshold.csv", row.names = FALSE)

############################################################
# 14) ROC FIGURE
############################################################

roc_df_ext <- data.frame(
  fpr = 1 - roc_obj_ext$specificities,
  tpr = roc_obj_ext$sensitivities
)

p_ext_roc_youden <- ggplot(roc_df_ext, aes(x = fpr, y = tpr)) +
  geom_line(linewidth = 1.1) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, color = "grey40") +
  geom_point(aes(x = best_fpr, y = best_sens), size = 3) +
  annotate(
    "text",
    x = best_fpr,
    y = best_sens,
    label = paste0("Youden threshold = ", sprintf("%.3f", best_threshold)),
    vjust = -1,
    size = 4
  ) +
  labs(
    title = "External validation ROC: GSE161731",
    subtitle = paste0(
      "Elastic Net | AUC = ",
      sprintf("%.3f", auc_ext),
      " (95% CI ",
      sprintf("%.3f", auc_ci_ext[1]), " - ",
      sprintf("%.3f", auc_ci_ext[3]), ")"
    ),
    x = "False Positive Rate",
    y = "True Positive Rate"
  ) +
  theme_bw(base_size = 14)

ggsave("Figure_ExternalValidation_ROC_GSE161731_Youden.png", p_ext_roc_youden, width = 7, height = 5.5, dpi = 600)

############################################################
# 15) DCA
############################################################

dca_df <- data.frame(
  outcome = ifelse(y_ext_fac == "COVID", 1, 0),
  pred = ext_prob
)

dca_model <- decision_curve(
  outcome ~ pred,
  data = dca_df,
  family = binomial(link = "logit"),
  thresholds = seq(0.01, 0.75, by = 0.01),
  confidence.intervals = 0.95,
  study.design = "cohort"
)

pdf("Figure_DecisionCurve_GSE161731.pdf", width = 8, height = 6)
plot_decision_curve(
  dca_model,
  curve.names = "Elastic Net",
  xlab = "Threshold Probability",
  ylab = "Net Benefit",
  standardize = FALSE
)
dev.off()

############################################################
# 16) SUMMARY
############################################################

summary_df <- data.frame(
  development_dataset = "GSE282464",
  development_task = "COVID vs Influenza",
  external_dataset = "GSE161731",
  model = "ElasticNet",
  n_external = nrow(external_predictions_df),
  n_external_covid = sum(external_predictions_df$truth == "COVID"),
  n_external_influenza = sum(external_predictions_df$truth == "Influenza"),
  auc = external_metrics_df$auc,
  auc_ci_low = external_metrics_df$auc_ci_low,
  auc_ci_high = external_metrics_df$auc_ci_high,
  accuracy = external_metrics_df$accuracy,
  balanced_accuracy = external_metrics_df$balanced_accuracy,
  sensitivity = external_metrics_df$sensitivity,
  specificity = external_metrics_df$specificity,
  precision = external_metrics_df$precision,
  f1 = external_metrics_df$f1,
  mcc = external_metrics_df$mcc,
  brier = external_metrics_df$brier,
  youden_threshold = best_threshold,
  stringsAsFactors = FALSE
)

write.csv(summary_df, "Paper_Ready_ExternalValidation_Summary_GSE161731.csv", row.names = FALSE)

cat("\nFinished external validation for GSE161731\n")
print(summary_df)







































































































































































##sil


############################################################
# MASTER SCRIPT
# Recreate external predictions + Calibration + DeLong
############################################################

rm(list = ls())
gc()

############################################################
# 0) PACKAGES
############################################################
cran_pkgs <- c(
  "data.table", "dplyr", "stringr", "ggplot2", "pROC",
  "glmnet", "GEOquery", "tibble"
)
bioc_pkgs <- c("edgeR", "AnnotationDbi", "org.Hs.eg.db")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
for (p in bioc_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    BiocManager::install(p, ask = FALSE, update = FALSE)
  }
}

library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(pROC)
library(glmnet)
library(GEOquery)
library(tibble)
library(edgeR)
library(AnnotationDbi)
library(org.Hs.eg.db)

set.seed(123)

############################################################
# 1) HELPERS
############################################################

safe_col_lower <- function(df) {
  names(df) <- tolower(names(df))
  df
}

clean_id <- function(x) {
  x <- as.character(x)
  x <- gsub("^X", "", x)
  x <- gsub("\\.gz$", "", x)
  x <- gsub("\\.csv$", "", x)
  x <- gsub("\\.txt$", "", x)
  x <- gsub("[^A-Za-z0-9]", "", x)
  tolower(x)
}

collapse_duplicate_genes <- function(expr_mat) {
  expr_df <- as.data.frame(expr_mat, check.names = FALSE)
  expr_df$gene <- rownames(expr_df)
  
  expr_df <- expr_df %>%
    group_by(gene) %>%
    summarise(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
  
  out <- as.matrix(expr_df[, -1, drop = FALSE])
  rownames(out) <- expr_df$gene
  out
}

fit_train_scaler <- function(x_train) {
  mu <- colMeans(x_train, na.rm = TRUE)
  sdv <- apply(x_train, 2, sd, na.rm = TRUE)
  sdv[is.na(sdv) | sdv == 0] <- 1
  list(center = mu, scale = sdv)
}

apply_scaler <- function(x, scaler) {
  sweep(sweep(x, 2, scaler$center, "-"), 2, scaler$scale, "/")
}

calc_external_metrics <- function(y_true_fac, prob, threshold = 0.5, positive_class = "COVID") {
  y_true <- ifelse(y_true_fac == positive_class, 1, 0)
  pred <- ifelse(prob >= threshold, 1, 0)
  
  tp <- sum(pred == 1 & y_true == 1)
  tn <- sum(pred == 0 & y_true == 0)
  fp <- sum(pred == 1 & y_true == 0)
  fn <- sum(pred == 0 & y_true == 1)
  
  accuracy <- (tp + tn) / length(y_true)
  sensitivity <- ifelse((tp + fn) == 0, NA_real_, tp / (tp + fn))
  specificity <- ifelse((tn + fp) == 0, NA_real_, tn / (tn + fp))
  precision <- ifelse((tp + fp) == 0, NA_real_, tp / (tp + fp))
  recall <- sensitivity
  
  f1 <- ifelse(
    is.na(precision) || is.na(recall) || (precision + recall) == 0,
    NA_real_,
    2 * precision * recall / (precision + recall)
  )
  
  balanced_accuracy <- mean(c(sensitivity, specificity), na.rm = TRUE)
  
  auc_obj <- tryCatch(
    pROC::roc(y_true, prob, quiet = TRUE),
    error = function(e) NULL
  )
  
  auc <- if (is.null(auc_obj)) NA_real_ else as.numeric(pROC::auc(auc_obj))
  auc_ci <- if (is.null(auc_obj)) c(NA_real_, NA_real_, NA_real_) else as.numeric(pROC::ci.auc(auc_obj))
  
  mcc_denom <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  mcc <- ifelse(mcc_denom == 0, NA_real_,
                ((tp * tn) - (fp * fn)) / mcc_denom)
  
  brier <- mean((prob - y_true)^2)
  
  data.frame(
    auc = auc,
    auc_ci_low = auc_ci[1],
    auc_ci_mid = auc_ci[2],
    auc_ci_high = auc_ci[3],
    accuracy = accuracy,
    balanced_accuracy = balanced_accuracy,
    sensitivity = sensitivity,
    specificity = specificity,
    precision = precision,
    f1 = f1,
    mcc = mcc,
    brier = brier,
    stringsAsFactors = FALSE
  )
}

prepare_dev_model <- function(task = c("covid_vs_sepsis", "covid_vs_influenza")) {
  task <- match.arg(task)
  
  gse_id <- "GSE282464"
  untar_dir <- file.path(gse_id, "untarred_raw")
  if (!dir.exists(untar_dir)) stop("GSE282464/untarred_raw klas??r?? bulunamad??.")
  
  count_files <- list.files(untar_dir, pattern = "\\.counts\\.gz$", full.names = TRUE)
  if (length(count_files) == 0) stop(".counts.gz dosyas?? bulunamad??.")
  
  read_one_count <- function(f) {
    dt <- fread(f, header = FALSE)
    dt <- dt[, 1:2]
    sample_id <- sub("_.*$", "", basename(f))
    colnames(dt) <- c("gene_id", sample_id)
    dt
  }
  
  count_list <- lapply(count_files, read_one_count)
  count_merged <- Reduce(function(x, y) merge(x, y, by = "gene_id", all = TRUE), count_list)
  count_merged[is.na(count_merged)] <- 0
  
  count_mat <- as.data.frame(count_merged)
  rownames(count_mat) <- count_mat$gene_id
  count_mat <- count_mat[, -1, drop = FALSE]
  count_mat <- as.matrix(sapply(count_mat, as.numeric))
  rownames(count_mat) <- count_merged$gene_id
  
  gset_list <- getGEO(gse_id, GSEMatrix = TRUE)
  gset <- gset_list[[1]]
  pheno <- pData(gset)
  
  pheno$sample_id <- as.character(pheno$geo_accession)
  common_samples <- intersect(colnames(count_mat), pheno$sample_id)
  
  pheno2 <- pheno %>% filter(sample_id %in% common_samples)
  count_mat2 <- count_mat[, pheno2$sample_id, drop = FALSE]
  pheno2 <- pheno2[match(colnames(count_mat2), pheno2$sample_id), , drop = FALSE]
  
  inf <- as.character(pheno2$`type of_infection:ch1`)
  
  if (task == "covid_vs_sepsis") {
    pheno2$label <- dplyr::case_when(
      inf == "COVID-19" ~ "COVID",
      inf %in% c("Sepsis", "Septic_shock") ~ "Sepsis_Group",
      TRUE ~ NA_character_
    )
    negative_label <- "Sepsis_Group"
  } else {
    pheno2$label <- dplyr::case_when(
      inf == "COVID-19" ~ "COVID",
      inf %in% c("Influenza_A", "Influenza_B") ~ "Influenza",
      TRUE ~ NA_character_
    )
    negative_label <- "Influenza"
  }
  
  dev_pheno <- pheno2 %>% filter(!is.na(label))
  dev_counts <- count_mat2[, dev_pheno$sample_id, drop = FALSE]
  dev_pheno <- dev_pheno[match(colnames(dev_counts), dev_pheno$sample_id), , drop = FALSE]
  
  y <- DGEList(counts = dev_counts)
  keep <- filterByExpr(y, group = dev_pheno$label)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")
  dev_xmat <- cpm(y, log = TRUE, prior.count = 1)
  
  x_dev_full <- t(dev_xmat)
  y_dev_full <- ifelse(dev_pheno$label == "COVID", 1, 0)
  
  top_k <- 100
  gene_var_dev <- apply(x_dev_full, 2, var, na.rm = TRUE)
  selected_genes <- names(sort(gene_var_dev, decreasing = TRUE))[1:min(top_k, length(gene_var_dev))]
  x_dev_fs <- x_dev_full[, selected_genes, drop = FALSE]
  
  scaler <- fit_train_scaler(x_dev_fs)
  x_dev_fs_sc <- apply_scaler(x_dev_fs, scaler)
  
  nfolds_use <- max(2, min(5, min(table(y_dev_full))))
  
  cv_fit <- cv.glmnet(
    x = as.matrix(x_dev_fs_sc),
    y = y_dev_full,
    family = "binomial",
    alpha = 0.5,
    nfolds = nfolds_use,
    type.measure = "auc"
  )
  
  final_enet <- glmnet(
    x = as.matrix(x_dev_fs_sc),
    y = y_dev_full,
    family = "binomial",
    alpha = 0.5,
    lambda = cv_fit$lambda.min
  )
  
  list(
    final_enet = final_enet,
    x_dev_fs = x_dev_fs,
    scaler = scaler,
    selected_genes = selected_genes,
    negative_label = negative_label
  )
}

############################################################
# 2) REBUILD GSE199816 PREDICTIONS
############################################################

dev_sepsis <- prepare_dev_model("covid_vs_sepsis")

gse_ext_id <- "GSE199816"
getGEOSuppFiles(gse_ext_id, baseDir = ".", fetch_files = TRUE)

raw_tar_candidates <- c(
  file.path(gse_ext_id, "GSE199816_RAW.tar"),
  file.path(gse_ext_id, "suppl", "GSE199816_RAW.tar")
)

raw_tar <- raw_tar_candidates[file.exists(raw_tar_candidates)][1]
if (is.na(raw_tar) || !file.exists(raw_tar)) stop("GSE199816 RAW tar bulunamad??.")

untar_dir_ext <- file.path(gse_ext_id, "untarred_raw")
dir.create(untar_dir_ext, showWarnings = FALSE)
untar(raw_tar, exdir = untar_dir_ext)

raw_files_ext <- list.files(untar_dir_ext, full.names = TRUE)

read_one_external_file <- function(f) {
  dt <- tryCatch(fread(f, header = TRUE, data.table = FALSE), error = function(e) NULL)
  if (is.null(dt) || ncol(dt) < 2) {
    dt <- tryCatch(fread(f, header = FALSE, data.table = FALSE), error = function(e) NULL)
  }
  if (is.null(dt) || ncol(dt) < 2) stop(paste("Dosya okunamad??:", basename(f)))
  
  cn <- tolower(colnames(dt))
  gene_col_idx <- which(cn %in% c("gene_id", "geneid", "gene", "ensembl_gene_id", "id", "symbol"))
  if (length(gene_col_idx) == 0) gene_col_idx <- 1 else gene_col_idx <- gene_col_idx[1]
  
  value_col_idx <- which(cn %in% c("count", "counts", "raw_count", "raw_counts", "expected_count"))
  if (length(value_col_idx) == 0) {
    other_cols <- setdiff(seq_len(ncol(dt)), gene_col_idx)
    value_col_idx <- other_cols[1]
  } else {
    value_col_idx <- value_col_idx[1]
  }
  
  sample_id <- sub("_.*$", "", basename(f))
  out <- dt[, c(gene_col_idx, value_col_idx), drop = FALSE]
  colnames(out) <- c("gene_id", sample_id)
  out
}

ext_list <- lapply(raw_files_ext, read_one_external_file)
ext_merged <- Reduce(function(x, y) merge(x, y, by = "gene_id", all = TRUE), ext_list)
ext_merged[is.na(ext_merged)] <- 0

expr_ext <- as.data.frame(ext_merged, check.names = FALSE)
rownames(expr_ext) <- expr_ext$gene_id
expr_ext <- expr_ext[, -1, drop = FALSE]
expr_ext <- as.matrix(sapply(expr_ext, as.numeric))
rownames(expr_ext) <- ext_merged$gene_id

gset_ext_list <- getGEO(gse_ext_id, GSEMatrix = TRUE, getGPL = FALSE)
gset_ext <- gset_ext_list[[1]]

pheno_ext <- pData(gset_ext) %>% as.data.frame()
pheno_ext$sample_id <- rownames(pheno_ext)
pheno_ext <- safe_col_lower(pheno_ext)

class_col_ext <- if ("diagnosis:ch1" %in% names(pheno_ext)) "diagnosis:ch1" else "source_name_ch1"
txt <- tolower(as.character(pheno_ext[[class_col_ext]]))

pheno_ext$label_binary <- dplyr::case_when(
  str_detect(txt, "covid|covid-19|sars-cov-2") ~ "COVID",
  str_detect(txt, "sepsis|septic") ~ "Sepsis_Group",
  TRUE ~ NA_character_
)

pheno_ext_bin <- pheno_ext %>% filter(!is.na(label_binary))
common_ext_samples <- intersect(colnames(expr_ext), pheno_ext_bin$sample_id)

expr_ext2 <- expr_ext[, common_ext_samples, drop = FALSE]
pheno_ext_bin <- pheno_ext_bin[match(colnames(expr_ext2), pheno_ext_bin$sample_id), , drop = FALSE]

mode(expr_ext2) <- "numeric"
rownames(expr_ext2) <- sub("\\..*$", "", rownames(expr_ext2))
expr_ext2 <- collapse_duplicate_genes(expr_ext2)

model_feature_order <- colnames(dev_sepsis$x_dev_fs)
model_feature_order_clean <- sub("\\..*$", "", model_feature_order)

external_is_ensg <- sum(grepl("^ENSG", rownames(expr_ext2))) > 0

if (external_is_ensg) {
  feature_map_df <- data.frame(model_ensembl = model_feature_order, model_key = model_feature_order_clean)
} else {
  map_dev_full <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = unique(model_feature_order_clean),
    keytype = "ENSEMBL",
    columns = c("ENSEMBL", "SYMBOL")
  )
  map_dev_full <- map_dev_full[!is.na(map_dev_full$SYMBOL) & map_dev_full$SYMBOL != "", ]
  map_dev_full <- map_dev_full[!duplicated(map_dev_full$ENSEMBL), ]
  model_feature_symbols <- map_dev_full$SYMBOL[match(model_feature_order_clean, map_dev_full$ENSEMBL)]
  feature_map_df <- data.frame(model_ensembl = model_feature_order, model_key = model_feature_symbols)
}

x_ext_full <- matrix(NA_real_, nrow = ncol(expr_ext2), ncol = length(model_feature_order))
colnames(x_ext_full) <- model_feature_order
rownames(x_ext_full) <- colnames(expr_ext2)

for (i in seq_len(nrow(feature_map_df))) {
  ens_col <- feature_map_df$model_ensembl[i]
  key_row <- feature_map_df$model_key[i]
  if (!is.na(key_row) && key_row != "" && key_row %in% rownames(expr_ext2)) {
    x_ext_full[, ens_col] <- as.numeric(expr_ext2[key_row, rownames(x_ext_full)])
  }
}

for (j in seq_len(ncol(x_ext_full))) {
  feat <- colnames(x_ext_full)[j]
  miss_idx <- is.na(x_ext_full[, j])
  if (any(miss_idx)) x_ext_full[miss_idx, j] <- dev_sepsis$scaler$center[feat]
}

x_ext_sc <- apply_scaler(x_ext_full, list(
  center = dev_sepsis$scaler$center[model_feature_order],
  scale = dev_sepsis$scaler$scale[model_feature_order]
))

y_ext_fac_199816 <- factor(pheno_ext_bin$label_binary, levels = c("Sepsis_Group", "COVID"))
prob_199816 <- as.numeric(predict(dev_sepsis$final_enet, newx = as.matrix(x_ext_sc), type = "response"))

pred_199816 <- data.frame(
  sample_id = pheno_ext_bin$sample_id,
  truth = as.character(y_ext_fac_199816),
  prob_covid = prob_199816,
  stringsAsFactors = FALSE
)

write.csv(pred_199816, "GSE199816_ExternalValidation_Predictions.csv", row.names = FALSE)

############################################################
# 3) REBUILD GSE161731 PREDICTIONS
############################################################

dev_influ <- prepare_dev_model("covid_vs_influenza")

gse_ext_id <- "GSE161731"
getGEOSuppFiles(gse_ext_id, baseDir = ".", fetch_files = TRUE)
all_files <- list.files(gse_ext_id, recursive = TRUE, full.names = TRUE)

counts_file <- all_files[grepl("GSE161731_counts\\.csv\\.gz$", all_files)][1]
counts_key_file <- all_files[grepl("GSE161731_counts_key\\.csv\\.gz$", all_files)][1]
key_file <- all_files[grepl("GSE161731_key\\.csv\\.gz$", all_files)][1]

counts_ext <- fread(counts_file, data.table = FALSE)
gene_ids <- as.character(counts_ext[[1]])
counts_ext2 <- counts_ext[, -1, drop = FALSE]

expr_ext_counts <- as.matrix(sapply(counts_ext2, as.numeric))
rownames(expr_ext_counts) <- gene_ids
expr_ext_counts <- expr_ext_counts[!is.na(rownames(expr_ext_counts)) & rownames(expr_ext_counts) != "", , drop = FALSE]

counts_key <- fread(counts_key_file, data.table = FALSE) %>% safe_col_lower()
meta_key <- fread(key_file, data.table = FALSE) %>% safe_col_lower()

sample_col_ck <- if ("rna_id" %in% names(counts_key)) "rna_id" else names(counts_key)[grepl("rna|sample|id", names(counts_key))][1]
label_col_ck <- names(counts_key)[grepl("cohort|group|diagnosis|disease|status|infection", names(counts_key))][1]

if (!is.na(label_col_ck) && label_col_ck %in% names(counts_key)) {
  counts_key$label_binary <- dplyr::case_when(
    str_detect(tolower(as.character(counts_key[[label_col_ck]])), "covid") ~ "COVID",
    str_detect(tolower(as.character(counts_key[[label_col_ck]])), "influenza") ~ "Influenza",
    TRUE ~ NA_character_
  )
} else {
  counts_key$label_binary <- NA_character_
}

if (sum(!is.na(counts_key$label_binary)) == 0 && "subject_id" %in% names(counts_key) && "subject_id" %in% names(meta_key)) {
  meta_label_col <- names(meta_key)[grepl("cohort|group|diagnosis|disease|status|infection", names(meta_key))][1]
  if (!is.na(meta_label_col) && meta_label_col %in% names(meta_key)) {
    meta_key$label_binary <- dplyr::case_when(
      str_detect(tolower(as.character(meta_key[[meta_label_col]])), "covid") ~ "COVID",
      str_detect(tolower(as.character(meta_key[[meta_label_col]])), "influenza") ~ "Influenza",
      TRUE ~ NA_character_
    )
    counts_key <- counts_key %>%
      left_join(meta_key %>% select(subject_id, label_binary), by = "subject_id", suffix = c("", ".meta")) %>%
      mutate(label_binary = ifelse(is.na(label_binary), label_binary.meta, label_binary)) %>%
      select(-label_binary.meta)
  }
}

counts_key$sample_clean <- clean_id(counts_key[[sample_col_ck]])
counts_key_bin <- counts_key %>% filter(!is.na(label_binary)) %>% distinct(sample_clean, .keep_all = TRUE)

expr_col_clean <- clean_id(colnames(expr_ext_counts))
match_idx <- match(expr_col_clean, counts_key_bin$sample_clean)

keep_cols <- !is.na(match_idx)
ext_counts_matched <- expr_ext_counts[, keep_cols, drop = FALSE]
counts_key_matched <- counts_key_bin[match_idx[keep_cols], , drop = FALSE]
external_label2 <- counts_key_matched$label_binary

dge_ext <- DGEList(counts = ext_counts_matched)
keep_ext <- filterByExpr(dge_ext, group = external_label2)
dge_ext <- dge_ext[keep_ext, , keep.lib.sizes = FALSE]
dge_ext <- calcNormFactors(dge_ext, method = "TMM")
expr_ext2 <- cpm(dge_ext, log = TRUE, prior.count = 1)

mode(expr_ext2) <- "numeric"
rownames(expr_ext2) <- sub("\\..*$", "", rownames(expr_ext2))
expr_ext2 <- collapse_duplicate_genes(expr_ext2)

model_feature_order <- colnames(dev_influ$x_dev_fs)
model_feature_order_clean <- sub("\\..*$", "", model_feature_order)
external_is_ensg <- sum(grepl("^ENSG", rownames(expr_ext2))) > 0

if (external_is_ensg) {
  feature_map_df <- data.frame(model_ensembl = model_feature_order, model_key = model_feature_order_clean)
} else {
  map_dev_full <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = unique(model_feature_order_clean),
    keytype = "ENSEMBL",
    columns = c("ENSEMBL", "SYMBOL")
  )
  map_dev_full <- map_dev_full[!is.na(map_dev_full$SYMBOL) & map_dev_full$SYMBOL != "", ]
  map_dev_full <- map_dev_full[!duplicated(map_dev_full$ENSEMBL), ]
  model_feature_symbols <- map_dev_full$SYMBOL[match(model_feature_order_clean, map_dev_full$ENSEMBL)]
  feature_map_df <- data.frame(model_ensembl = model_feature_order, model_key = model_feature_symbols)
}

x_ext_full <- matrix(NA_real_, nrow = ncol(expr_ext2), ncol = length(model_feature_order))
colnames(x_ext_full) <- model_feature_order
rownames(x_ext_full) <- colnames(expr_ext2)

for (i in seq_len(nrow(feature_map_df))) {
  ens_col <- feature_map_df$model_ensembl[i]
  key_row <- feature_map_df$model_key[i]
  if (!is.na(key_row) && key_row != "" && key_row %in% rownames(expr_ext2)) {
    x_ext_full[, ens_col] <- as.numeric(expr_ext2[key_row, rownames(x_ext_full)])
  }
}

for (j in seq_len(ncol(x_ext_full))) {
  feat <- colnames(x_ext_full)[j]
  miss_idx <- is.na(x_ext_full[, j])
  if (any(miss_idx)) x_ext_full[miss_idx, j] <- dev_influ$scaler$center[feat]
}

x_ext_sc <- apply_scaler(x_ext_full, list(
  center = dev_influ$scaler$center[model_feature_order],
  scale = dev_influ$scaler$scale[model_feature_order]
))

y_ext_fac_161731 <- factor(external_label2, levels = c("Influenza", "COVID"))
prob_161731 <- as.numeric(predict(dev_influ$final_enet, newx = as.matrix(x_ext_sc), type = "response"))

pred_161731 <- data.frame(
  sample_id = colnames(ext_counts_matched),
  truth = as.character(y_ext_fac_161731),
  prob_covid = prob_161731,
  stringsAsFactors = FALSE
)

write.csv(pred_161731, "GSE161731_ExternalValidation_Predictions.csv", row.names = FALSE)

############################################################
# 4) CALIBRATION + DELONG
############################################################

prepare_pred_df <- function(df) {
  df$truth01 <- ifelse(df$truth == "COVID", 1, 0)
  df
}

calc_calibration_stats <- function(df) {
  eps <- 1e-6
  p <- pmin(pmax(df$prob_covid, eps), 1 - eps)
  logit_p <- qlogis(p)
  fit <- glm(truth01 ~ logit_p, data = df, family = binomial())
  
  data.frame(
    calibration_intercept = coef(fit)[1],
    calibration_slope = coef(fit)[2],
    stringsAsFactors = FALSE
  )
}

make_calibration_plot <- function(df, dataset_name, outfile) {
  cal_df <- df %>%
    mutate(bin = ntile(prob_covid, 10)) %>%
    group_by(bin) %>%
    summarise(
      mean_pred = mean(prob_covid, na.rm = TRUE),
      obs_rate = mean(truth01, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )
  
  p <- ggplot(cal_df, aes(x = mean_pred, y = obs_rate)) +
    geom_point(size = 3) +
    geom_line() +
    geom_abline(intercept = 0, slope = 1, linetype = 2, color = "grey40") +
    xlim(0, 1) +
    ylim(0, 1) +
    labs(
      title = paste0("Calibration plot: ", dataset_name),
      x = "Mean predicted probability",
      y = "Observed event rate"
    ) +
    theme_bw(base_size = 14)
  
  ggsave(outfile, p, width = 6.5, height = 5.5, dpi = 600)
  cal_df
}

pred_199816 <- prepare_pred_df(pred_199816)
pred_161731 <- prepare_pred_df(pred_161731)

cal_199816 <- calc_calibration_stats(pred_199816)
cal_161731 <- calc_calibration_stats(pred_161731)

met_199816 <- data.frame(
  auc = as.numeric(pROC::auc(pROC::roc(pred_199816$truth01, pred_199816$prob_covid, quiet = TRUE))),
  brier = mean((pred_199816$prob_covid - pred_199816$truth01)^2)
)

met_161731 <- data.frame(
  auc = as.numeric(pROC::auc(pROC::roc(pred_161731$truth01, pred_161731$prob_covid, quiet = TRUE))),
  brier = mean((pred_161731$prob_covid - pred_161731$truth01)^2)
)

summary_calibration <- bind_rows(
  cbind(dataset = "GSE199816", met_199816, cal_199816),
  cbind(dataset = "GSE161731", met_161731, cal_161731)
)

write.csv(summary_calibration, "ExternalValidation_Calibration_Summary.csv", row.names = FALSE)

cal_plot_199816 <- make_calibration_plot(pred_199816, "GSE199816", "Figure_Calibration_GSE199816.png")
cal_plot_161731 <- make_calibration_plot(pred_161731, "GSE161731", "Figure_Calibration_GSE161731.png")

cal_plot_199816$dataset <- "GSE199816"
cal_plot_161731$dataset <- "GSE161731"
cal_combined <- bind_rows(cal_plot_199816, cal_plot_161731)

p_cal_combined <- ggplot(cal_combined, aes(x = mean_pred, y = obs_rate)) +
  geom_point(size = 3) +
  geom_line() +
  geom_abline(intercept = 0, slope = 1, linetype = 2, color = "grey40") +
  facet_wrap(~ dataset) +
  xlim(0, 1) +
  ylim(0, 1) +
  labs(
    title = "Calibration plots for external validation cohorts",
    x = "Mean predicted probability",
    y = "Observed event rate"
  ) +
  theme_bw(base_size = 14)

ggsave("Figure_Calibration_ExternalValidation_Combined.png", p_cal_combined, width = 10, height = 5.5, dpi = 600)

roc_199816 <- pROC::roc(pred_199816$truth01, pred_199816$prob_covid, quiet = TRUE)
roc_161731 <- pROC::roc(pred_161731$truth01, pred_161731$prob_covid, quiet = TRUE)

delong_external <- tryCatch(
  pROC::roc.test(roc_199816, roc_161731, method = "delong"),
  error = function(e) NULL
)

if (!is.null(delong_external)) {
  delong_external_df <- data.frame(
    comparison = "GSE199816_vs_GSE161731",
    auc_1 = as.numeric(pROC::auc(roc_199816)),
    auc_2 = as.numeric(pROC::auc(roc_161731)),
    p_value = delong_external$p.value,
    method = delong_external$method,
    stringsAsFactors = FALSE
  )
} else {
  delong_external_df <- data.frame(
    comparison = "GSE199816_vs_GSE161731",
    auc_1 = as.numeric(pROC::auc(roc_199816)),
    auc_2 = as.numeric(pROC::auc(roc_161731)),
    p_value = NA_real_,
    method = "DeLong failed",
    stringsAsFactors = FALSE
  )
}

write.csv(delong_external_df, "DeLong_ExternalValidation_Comparison.csv", row.names = FALSE)

cat("\nDone.\n")
print(summary_calibration)
print(delong_external_df)




# =========================================================
# COLORFUL EXTERNAL VALIDATION FIGURE
# FIXED VERSION: dplyr::select conflict solved
# =========================================================

if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr")
if (!requireNamespace("pROC", quietly = TRUE)) install.packages("pROC")
if (!requireNamespace("cowplot", quietly = TRUE)) install.packages("cowplot")

library(ggplot2)
library(dplyr)
library(readr)
library(pROC)
library(cowplot)

# Directory containing external prediction CSV files for the final combined figure.
# Keep this relative for GitHub/reproducible runs.
fig_dir <- "."

f161731 <- file.path(fig_dir, "GSE161731_ExternalValidation_Predictions.csv")
f199816 <- file.path(fig_dir, "GSE199816_ExternalValidation_Predictions_WithYouden.csv")

if (!file.exists(f161731)) stop("GSE161731 prediction csv bulunamad??.")
if (!file.exists(f199816)) stop("GSE199816 prediction csv bulunamad??.")

# ---------------------------------------------------------
# helper
# ---------------------------------------------------------
read_pred_file <- function(path, dataset_name) {
  df <- readr::read_csv(path, show_col_types = FALSE)
  
  names(df) <- tolower(names(df))
  
  truth_col <- names(df)[names(df) %in% c("truth", "label", "true_label", "observed")]
  if (length(truth_col) == 0) stop(paste("Truth kolonu bulunamad??:", path))
  truth_col <- truth_col[1]
  
  prob_col <- names(df)[names(df) %in% c("prob_covid", "prob", "prediction", "pred_prob", "probability")]
  if (length(prob_col) == 0) stop(paste("Probability kolonu bulunamad??:", path))
  prob_col <- prob_col[1]
  
  out <- df %>%
    dplyr::mutate(
      dataset = dataset_name,
      truth = .data[[truth_col]],
      prob = as.numeric(.data[[prob_col]])
    ) %>%
    dplyr::select(dataset, truth, prob)
  
  out <- out %>%
    dplyr::mutate(
      truth01 = ifelse(tolower(as.character(truth)) %in% c("covid", "covid-19", "1"), 1, 0)
    )
  
  return(out)
}

make_calibration_df <- function(df, n_bins = 10) {
  df <- df[is.finite(df$prob) & !is.na(df$truth01), , drop = FALSE]
  
  tmp <- df %>%
    dplyr::mutate(bin = dplyr::ntile(prob, n_bins)) %>%
    dplyr::group_by(bin) %>%
    dplyr::summarise(
      mean_pred = mean(prob, na.rm = TRUE),
      obs_rate = mean(truth01, na.rm = TRUE),
      n = dplyr::n(),
      .groups = "drop"
    )
  
  return(tmp)
}

plot_roc_panel <- function(df, title_text, curve_color) {
  roc_obj <- pROC::roc(df$truth01, df$prob, quiet = TRUE)
  auc_val <- as.numeric(pROC::auc(roc_obj))
  
  roc_df <- data.frame(
    fpr = 1 - roc_obj$specificities,
    tpr = roc_obj$sensitivities
  )
  
  ggplot(roc_df, aes(x = fpr, y = tpr)) +
    geom_line(linewidth = 1.4, color = curve_color) +
    geom_abline(intercept = 0, slope = 1, linetype = 2, color = "grey50") +
    annotate(
      "label",
      x = 0.68, y = 0.15,
      label = paste0("AUC = ", sprintf("%.3f", auc_val)),
      fill = "white",
      color = curve_color,
      fontface = "bold",
      size = 4
    ) +
    labs(
      title = title_text,
      x = "False Positive Rate",
      y = "True Positive Rate"
    ) +
    coord_equal() +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.grid.minor = element_blank()
    )
}

plot_cal_panel <- function(df, title_text, point_color, line_color) {
  cal_df <- make_calibration_df(df, n_bins = 10)
  
  ggplot(cal_df, aes(x = mean_pred, y = obs_rate)) +
    geom_abline(intercept = 0, slope = 1, linetype = 2, color = "grey50") +
    geom_line(linewidth = 1.2, color = line_color) +
    geom_point(size = 3.2, color = point_color) +
    labs(
      title = title_text,
      x = "Mean predicted probability",
      y = "Observed event rate"
    ) +
    xlim(0, 1) +
    ylim(0, 1) +
    coord_equal() +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.grid.minor = element_blank()
    )
}

# ---------------------------------------------------------
# data
# ---------------------------------------------------------
d161731 <- read_pred_file(f161731, "GSE161731")
d199816 <- read_pred_file(f199816, "GSE199816")

# ---------------------------------------------------------
# panels
# ---------------------------------------------------------
p1 <- plot_roc_panel(
  d161731,
  "A. ROC ??? GSE161731",
  curve_color = "#D81B60"
)

p2 <- plot_roc_panel(
  d199816,
  "B. ROC ??? GSE199816",
  curve_color = "#1E88E5"
)

p3 <- plot_cal_panel(
  d161731,
  "C. Calibration ??? GSE161731",
  point_color = "#43A047",
  line_color = "#43A047"
)

p4 <- plot_cal_panel(
  d199816,
  "D. Calibration ??? GSE199816",
  point_color = "#FB8C00",
  line_color = "#FB8C00"
)

final_plot <- cowplot::plot_grid(
  p1, p2, p3, p4,
  ncol = 2,
  align = "hv"
)

ggsave(
  file.path(fig_dir, "Figure_ExternalValidation_Colorful_Lines.png"),
  final_plot,
  width = 14,
  height = 10,
  dpi = 700,
  bg = "white"
)

if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
if (!requireNamespace("pROC", quietly = TRUE)) install.packages("pROC")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")

library(ggplot2)
library(patchwork)
library(pROC)
library(dplyr)

############################################################
# 1) READ EXTERNAL PREDICTION FILES
############################################################

pred_199816 <- read.csv("GSE199816_ExternalValidation_Predictions_WithYouden.csv")
pred_161731 <- read.csv("GSE161731_ExternalValidation_Predictions_WithYouden.csv")

############################################################
# 2) PREPARE DATA
############################################################

prepare_pred_df <- function(df) {
  df$truth01 <- ifelse(df$truth == "COVID", 1, 0)
  df
}

pred_199816 <- prepare_pred_df(pred_199816)
pred_161731 <- prepare_pred_df(pred_161731)

############################################################
# 3) ROC PLOT FUNCTION
############################################################

make_roc_plot <- function(df, panel_title, cohort_name) {
  roc_obj <- pROC::roc(df$truth01, df$prob_covid, quiet = TRUE)
  auc_val <- as.numeric(pROC::auc(roc_obj))
  auc_ci <- as.numeric(pROC::ci.auc(roc_obj))
  
  roc_df <- data.frame(
    fpr = 1 - roc_obj$specificities,
    tpr = roc_obj$sensitivities
  )
  
  ggplot(roc_df, aes(x = fpr, y = tpr)) +
    geom_line(linewidth = 1.1) +
    geom_abline(intercept = 0, slope = 1, linetype = 2, color = "grey50") +
    labs(
      title = paste0(panel_title, " ", cohort_name),
      subtitle = paste0(
        "AUC = ", sprintf("%.3f", auc_val),
        " (95% CI ", sprintf("%.3f", auc_ci[1]), "???", sprintf("%.3f", auc_ci[3]), ")"
      ),
      x = "False Positive Rate",
      y = "True Positive Rate"
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10)
    )
}

############################################################
# 4) CALIBRATION PLOT FUNCTION
############################################################

make_calibration_plot <- function(df, panel_title, cohort_name) {
  cal_df <- df %>%
    mutate(bin = ntile(prob_covid, 10)) %>%
    group_by(bin) %>%
    summarise(
      mean_pred = mean(prob_covid, na.rm = TRUE),
      obs_rate = mean(truth01, na.rm = TRUE),
      .groups = "drop"
    )
  
  ggplot(cal_df, aes(x = mean_pred, y = obs_rate)) +
    geom_point(size = 3) +
    geom_line() +
    geom_abline(intercept = 0, slope = 1, linetype = 2, color = "grey50") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(
      title = paste0(panel_title, " ", cohort_name),
      x = "Mean predicted probability",
      y = "Observed event rate"
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold")
    )
}

############################################################
# 5) BUILD PANELS
############################################################

pA <- make_roc_plot(pred_161731, "A.", "GSE161731")
pB <- make_roc_plot(pred_199816, "B.", "GSE199816")
pC <- make_calibration_plot(pred_161731, "C.", "GSE161731")
pD <- make_calibration_plot(pred_199816, "D.", "GSE199816")

############################################################
# 6) COMBINE FIGURE 4
############################################################

p_fig4 <- (pA + pB) / (pC + pD) +
  plot_annotation(
    title = "External validation performance of the locked Elastic Net model",
    theme = theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16)
    )
  )

print(p_fig4)

ggsave(
  "Figure_4_ExternalValidation_ROC_Calibration_Combined.png",
  p_fig4,
  width = 12,
  height = 10,
  dpi = 600
)
############################################################
# Build Table 2: External validation performance summary
#
# Run this script in the same working directory where the external
# validation CSV outputs were written by final_stable_analysis.R.
############################################################

required_files <- c(
  "GSE199816_ExternalValidation_Metrics.csv",
  "GSE199816_ExternalValidation_Metrics_YoudenThreshold.csv",
  "GSE161731_ExternalValidation_Metrics.csv",
  "GSE161731_ExternalValidation_Metrics_YoudenThreshold.csv",
  "ExternalValidation_Calibration_Summary.csv",
  "DeLong_ExternalValidation_Comparison.csv"
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop(
    "The following required files are missing:\n",
    paste0(" - ", missing_files, collapse = "\n"),
    "\n\nRun final_stable_analysis.R first, then rerun this script from the output directory."
  )
}

read_one <- function(path) {
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

add_missing_columns <- function(df, cols) {
  for (cc in cols) {
    if (!cc %in% names(df)) df[[cc]] <- NA
  }
  df
}

metric_cols <- c(
  "auc", "auc_ci_low", "auc_ci_mid", "auc_ci_high",
  "accuracy", "balanced_accuracy", "sensitivity", "specificity",
  "precision", "f1", "mcc", "brier", "threshold"
)

format_dataset_rows <- function(dataset, comparator, metrics_05_path, metrics_youden_path) {
  m05 <- read_one(metrics_05_path)
  my <- read_one(metrics_youden_path)
  
  m05 <- add_missing_columns(m05, metric_cols)
  my <- add_missing_columns(my, metric_cols)
  
  m05$threshold <- 0.5
  m05$threshold_type <- "0.5"
  my$threshold_type <- "Youden"
  
  out <- rbind(m05[, c(metric_cols, "threshold_type")], my[, c(metric_cols, "threshold_type")])
  out$external_dataset <- dataset
  out$comparator_setting <- comparator
  
  out[, c(
    "external_dataset", "comparator_setting", "threshold_type", "threshold",
    "auc", "auc_ci_low", "auc_ci_high",
    "accuracy", "balanced_accuracy", "sensitivity", "specificity",
    "precision", "f1", "mcc", "brier"
  )]
}

detail_table <- rbind(
  format_dataset_rows(
    dataset = "GSE199816",
    comparator = "COVID-19 vs sepsis",
    metrics_05_path = "GSE199816_ExternalValidation_Metrics.csv",
    metrics_youden_path = "GSE199816_ExternalValidation_Metrics_YoudenThreshold.csv"
  ),
  format_dataset_rows(
    dataset = "GSE161731",
    comparator = "COVID-19 vs influenza",
    metrics_05_path = "GSE161731_ExternalValidation_Metrics.csv",
    metrics_youden_path = "GSE161731_ExternalValidation_Metrics_YoudenThreshold.csv"
  )
)

calibration <- read_one("ExternalValidation_Calibration_Summary.csv")
calibration <- add_missing_columns(
  calibration,
  c("dataset", "calibration_intercept", "calibration_slope", "auc", "brier")
)

delong <- read_one("DeLong_ExternalValidation_Comparison.csv")
delong <- add_missing_columns(delong, c("comparison", "auc_1", "auc_2", "p_value", "method"))

detail_table <- merge(
  detail_table,
  calibration[, c("dataset", "calibration_intercept", "calibration_slope")],
  by.x = "external_dataset",
  by.y = "dataset",
  all.x = TRUE
)

detail_table$delong_comparison <- if (nrow(delong) >= 1) delong$comparison[1] else NA
detail_table$delong_p_value <- if (nrow(delong) >= 1) delong$p_value[1] else NA
detail_table$delong_method <- if (nrow(delong) >= 1) delong$method[1] else NA

detail_table <- detail_table[order(detail_table$external_dataset, detail_table$threshold_type), ]

# Manuscript-ready Table 2: one row per external cohort using Youden
# threshold-based metrics, plus discrimination, Brier score, and calibration.
table2 <- subset(detail_table, threshold_type == "Youden")
table2 <- table2[, c(
  "external_dataset", "comparator_setting", "threshold",
  "auc", "auc_ci_low", "auc_ci_high",
  "accuracy", "balanced_accuracy", "sensitivity", "specificity",
  "precision", "f1", "mcc", "brier",
  "calibration_intercept", "calibration_slope",
  "delong_comparison", "delong_p_value", "delong_method"
)]

write.csv(table2, "Table2_ExternalValidation_Performance.csv", row.names = FALSE)
write.csv(detail_table, "Table2_ExternalValidation_Performance_Detailed.csv", row.names = FALSE)

cat("\nWrote:\n")
cat("- Table2_ExternalValidation_Performance.csv\n")
cat("- Table2_ExternalValidation_Performance_Detailed.csv\n")
