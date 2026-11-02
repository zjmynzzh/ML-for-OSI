# =========================================================
# 宏基因组四分类：Stacking + RF（绘图美化 + 审稿人全指标版）
# 新增：Sensitivity / Specificity / PPV / NPV / LR+ / LR- /
#       Balanced Accuracy / F1（OvR、OvO、PerClass 三处均输出）
# =========================================================

rm(list = ls())
set.seed(42)

# ---- 运行配置 ----
RUN_ROBUSTNESS <- TRUE
N_REP          <- 5

pkgs <- c("tidyverse","caret","xgboost","pROC","ggplot2",
          "patchwork","RColorBrewer","reshape2","scales",
          "MLmetrics","ggrepel","ggbeeswarm","randomForest")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages({
  library(tidyverse); library(caret); library(xgboost); library(pROC)
  library(ggplot2); library(patchwork); library(RColorBrewer)
  library(reshape2); library(scales); library(MLmetrics)
  library(ggrepel); library(ggbeeswarm); library(randomForest)
})

# =========================================================
# ★ ggsave 全局包装：固定使用 pdf 输出
# =========================================================
PDF_DEVICE <- "pdf"
.gs_orig <- ggplot2::ggsave
ggsave <- function(filename, plot = ggplot2::last_plot(),
                   ..., width = NA, height = NA,
                   dpi = 600, device = NULL, useDingbats = FALSE) {
  if (grepl("\\.png$", filename, ignore.case = TRUE)) {
    filename <- sub("\\.png$", ".pdf", filename, ignore.case = TRUE)
  }
  if (is.null(device)) device <- PDF_DEVICE
  .gs_orig(filename = filename, plot = plot, ...,
           width = width, height = height, dpi = dpi,
           device = device, useDingbats = useDingbats)
}

PAL     <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")
PAL_SET <- c("Train" = "#377EB8", "Test" = "#E41A1C")
PAL_PN  <- c("Positive" = "#4DAF4A", "Negative" = "#984EA3")

safe_ci <- function(r) tryCatch(
  as.numeric(pROC::ci.auc(r, method = "delong")),
  error = function(e) c(NA_real_, NA_real_, NA_real_)
)
safe_se <- function(r) tryCatch(
  sqrt(pROC::var(r, method = "delong")),
  error = function(e) NA_real_
)

# 安全的 LR 计算（除零返回 Inf，符合审稿人惯例）
safe_lr_pos <- function(sens, spec) ifelse(spec < 1, sens / (1 - spec), Inf)
safe_lr_neg <- function(sens, spec) ifelse(spec > 0, (1 - sens) / spec, Inf)

theme_top <- function(base_size = 18) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid        = element_blank(),
      panel.border      = element_rect(color = "black", linewidth = 1.5, fill = NA),
      axis.line         = element_blank(),
      axis.ticks        = element_line(linewidth = 1.2, color = "black"),
      axis.ticks.length = unit(0.2, "cm"),
      axis.text         = element_text(face = "bold", size = base_size, color = "black"),
      axis.title        = element_text(face = "bold", size = base_size + 2, color = "black"),
      plot.title        = element_text(face = "bold", hjust = 0.5, size = base_size + 3,
                                       margin = unit(c(0, 0, 10, 0), "pt")),
      strip.background  = element_rect(fill = "gray90", color = "black", linewidth = 1.2),
      strip.text        = element_text(face = "bold", size = base_size + 1),
      legend.text       = element_text(face = "bold", size = base_size - 2),
      legend.title      = element_text(face = "bold", size = base_size),
      legend.background = element_blank(),
      legend.key        = element_blank()
    )
}

preferred_order <- c("bacterial", "fungal", "viral", "Others")
TOP_N_SPECIES   <- 30

# =========================================================
# 文件定位与读取
# =========================================================
desktop <- file.path(Sys.getenv("USERPROFILE"), "Desktop")
if (!dir.exists(desktop)) desktop <- file.path(Sys.getenv("USERPROFILE"), "OneDrive", "Desktop")
if (!dir.exists(desktop)) desktop <- file.path(Sys.getenv("USERPROFILE"), "桌面")
if (!dir.exists(desktop)) stop("找不到桌面路径")

all_files <- list.files(desktop, full.names = FALSE)
group_name <- all_files[grepl("分组|group", all_files, ignore.case = TRUE) &
                          grepl("\\.csv$", all_files, ignore.case = TRUE)]
if (length(group_name) == 0) stop("找不到分组 CSV")
group_file <- file.path(desktop, group_name[1])

all_csv <- all_files[grepl("\\.csv$", all_files, ignore.case = TRUE)]
abund_candidates <- all_csv[!grepl("分组|group", all_csv, ignore.case = TRUE)]
abund_match <- abund_candidates[grepl("abundance|final-ALL|ra|microbial|species|taxa",
                                      abund_candidates, ignore.case = TRUE)]
abundance_file <- if (length(abund_match) > 0) file.path(desktop, abund_match[1]) else
  file.path(desktop, abund_candidates[1])

outdir   <- file.path(desktop, "project_results_stacking_hybrid")
plot_dir <- file.path(outdir, "plots")
dir.create(outdir, showWarnings = FALSE)
dir.create(plot_dir, showWarnings = FALSE)

groups <- read.csv(group_file, fileEncoding = "UTF-8-BOM",
                   stringsAsFactors = FALSE, check.names = FALSE)
colnames(groups)[1:2] <- c("Sample", "Group")
groups$Sample <- trimws(groups$Sample)
groups$Group  <- trimws(groups$Group)
groups <- groups %>% distinct(Sample, .keep_all = TRUE) %>% column_to_rownames("Sample")

X_raw <- read.csv(abundance_file, row.names = 1, check.names = FALSE,
                  stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
X_raw <- X_raw[, !colnames(X_raw) %in% paste0("H", 1:20), drop = FALSE]
X_raw[] <- lapply(X_raw, function(x) suppressWarnings(as.numeric(as.character(x))))
X_raw[is.na(X_raw)] <- 0
X_raw <- X_raw[, colSums(X_raw) > 0, drop = FALSE]
X <- as.data.frame(t(as.matrix(X_raw)))
common <- intersect(rownames(X), rownames(groups))
X <- X[common, , drop = FALSE]

actual_classes <- sort(unique(as.character(groups[common, "Group"])))
present_levels <- intersect(preferred_order, actual_classes)
if (length(present_levels) < 2) present_levels <- actual_classes
y <- factor(groups[common, "Group"], levels = present_levels)

cnt <- table(y)
if (any(cnt < 5)) {
  keep_cls <- names(cnt)[cnt >= 5]
  keep_smp <- as.character(y) %in% keep_cls
  y <- droplevels(y[keep_smp])
  X <- X[keep_smp, , drop = FALSE]
}

prevalence <- colMeans(X > 0)
for (th in c(0.10, 0.05, 0.03, 0.02, 0.01, 0.005, 0.001, 0)) {
  keep <- prevalence >= th
  if (sum(keep) >= 5) { X <- X[, keep, drop = FALSE]; break }
}
if (ncol(X) > 200) X <- X[, order(colSums(X), decreasing = TRUE)[1:200], drop = FALSE]
X <- X + 1e-6
X <- X / rowSums(X)
X_log <- log(X)
X_clr <- X_log - rowMeans(X_log)
top_species <- names(sort(colSums(X), decreasing = TRUE))[1:min(TOP_N_SPECIES, ncol(X))]

fit_stacking <- function(X_clr, y, top_species, seed = 42, verbose = TRUE) {
  set.seed(seed)
  train_idx <- createDataPartition(y, p = 0.7, list = FALSE)
  X_train <- X_clr[train_idx, , drop = FALSE]
  X_test  <- X_clr[-train_idx, , drop = FALSE]
  y_train <- droplevels(y[train_idx])
  y_test  <- factor(y[-train_idx], levels = levels(y_train))
  
  classes    <- levels(y_train)
  pairs_list <- combn(classes, 2, simplify = FALSE)
  pair_names <- vapply(pairs_list, function(x) paste0(x[1], "_vs_", x[2]), character(1))
  
  gb_params <- list(
    objective = "binary:logistic", eval_metric = "logloss",
    eta = 0.05, max_depth = 2, min_child_weight = 8,
    subsample = 0.7, colsample_bytree = 0.5,
    gamma = 1, lambda = 10, alpha = 1
  )
  
  n_fold_all  <- min(5, min(table(y_train)))
  set.seed(seed)
  fold_id <- createFolds(y_train, k = n_fold_all, list = FALSE)
  
  oof_scores  <- matrix(NA, nrow(X_train), length(pairs_list))
  test_scores <- matrix(NA, nrow(X_test),  length(pairs_list))
  pair_models <- list()
  
  for (j in seq_along(pairs_list)) {
    cls1 <- pairs_list[[j]][1]; cls2 <- pairs_list[[j]][2]
    pair_name <- pair_names[j]
    in_pair_all <- as.character(y_train) %in% c(cls1, cls2)
    X_pair <- X_train[in_pair_all, , drop = FALSE]
    y_pair <- as.numeric(as.character(y_train[in_pair_all]) == cls1)
    
    for (f in seq_len(n_fold_all)) {
      tr_idx <- which(in_pair_all & fold_id != f)
      if (length(tr_idx) < 5) next
      y_tr <- as.numeric(as.character(y_train[tr_idx]) == cls1)
      if (length(unique(y_tr)) < 2) next
      dtr <- xgb.DMatrix(as.matrix(X_train[tr_idx, , drop = FALSE]), label = y_tr)
      m <- xgb.train(params = gb_params, data = dtr, nrounds = 200, verbose = 0)
      val_idx <- which(fold_id == f)
      oof_scores[val_idx, j] <- predict(m, xgb.DMatrix(as.matrix(X_train[val_idx, , drop = FALSE])))
    }
    
    dfull <- xgb.DMatrix(as.matrix(X_pair), label = y_pair)
    m_full <- xgb.train(params = gb_params, data = dfull, nrounds = 200, verbose = 0)
    test_scores[, j] <- predict(m_full, xgb.DMatrix(as.matrix(X_test)))
    pair_models[[pair_name]] <- m_full
  }
  
  for (j in seq_along(pairs_list)) {
    if (any(is.na(oof_scores[, j]))) {
      med <- median(oof_scores[, j], na.rm = TRUE); if (is.na(med)) med <- 0.5
      oof_scores[is.na(oof_scores[, j]), j] <- med
    }
  }
  colnames(oof_scores)  <- pair_names
  colnames(test_scores) <- pair_names
  
  X_train_meta <- cbind(as.data.frame(oof_scores),  X_train[, top_species, drop = FALSE])
  X_test_meta  <- cbind(as.data.frame(test_scores), X_test[,  top_species, drop = FALSE])
  all_feats    <- c(pair_names, top_species)
  X_train_meta$y <- y_train
  
  class_counts  <- table(y_train)
  class_weights <- sqrt(max(class_counts) / class_counts)
  names(class_weights) <- names(class_counts)
  for (tg in c("bacter", "fung")) {
    idx <- grep(tg, names(class_weights), ignore.case = TRUE)
    if (length(idx) > 0) class_weights[idx] <- class_weights[idx] * 1.6
  }
  
  seeds <- c(42, 123, 2024, 7, 99)
  rf_prob_train <- matrix(0, nrow(X_train_meta), length(classes))
  rf_prob_test  <- matrix(0, nrow(X_test_meta),  length(classes))
  colnames(rf_prob_train) <- classes
  colnames(rf_prob_test)  <- classes
  rf_last <- NULL
  
  for (s in seeds) {
    set.seed(s)
    rf_s <- randomForest(
      x = X_train_meta[, all_feats, drop = FALSE], y = y_train,
      ntree = 1500, mtry = 3, nodesize = 3, maxnodes = 60,
      classwt = class_weights, importance = TRUE)
    rf_prob_train <- rf_prob_train +
      predict(rf_s, X_train_meta[, all_feats, drop = FALSE], type = "prob")[, classes]
    rf_prob_test <- rf_prob_test +
      predict(rf_s, X_test_meta[, all_feats, drop = FALSE], type = "prob")[, classes]
    rf_last <- rf_s
  }
  rf_prob_train <- rf_prob_train / length(seeds)
  rf_prob_test  <- rf_prob_test  / length(seeds)
  rf_model <- rf_last
  
  rf_train_pred <- factor(classes[apply(rf_prob_train, 1, which.max)], levels = classes)
  rf_test_pred  <- factor(classes[apply(rf_prob_test,  1, which.max)], levels = classes)
  
  cm_train <- confusionMatrix(rf_train_pred, y_train)
  cm_test  <- confusionMatrix(rf_test_pred,  y_test)
  
  if (verbose) {
    cat("\n----- seed =", seed, "-----\n")
    cat("Train Acc:", round(cm_train$overall["Accuracy"], 4),
        " | Test Acc:", round(cm_test$overall["Accuracy"], 4), "\n")
    cat("Test 每类:\n")
    print(round(cm_test$byClass[, c("Sensitivity","Specificity","F1")], 3))
  }
  
  list(X_train = X_train, X_test = X_test, y_train = y_train, y_test = y_test,
       classes = classes, pairs_list = pairs_list, pair_names = pair_names,
       oof_scores = oof_scores, test_scores = test_scores,
       pair_models = pair_models,
       X_train_meta = X_train_meta, X_test_meta = X_test_meta,
       all_feats = all_feats, seeds = seeds, class_weights = class_weights,
       rf_model = rf_model,
       rf_prob_train = rf_prob_train, rf_prob_test = rf_prob_test,
       rf_train_pred = rf_train_pred, rf_test_pred = rf_test_pred,
       cm_train = cm_train, cm_test = cm_test, gb_params = gb_params)
}

cat("\n########## 主运行：seed = 42 ##########\n")
res <- fit_stacking(X_clr, y, top_species = top_species, seed = 42, verbose = TRUE)
list2env(res, envir = .GlobalEnv)

# =========================================================
# 两两 ROC（含 DeLong 95% CI，UI 美化版）
# =========================================================
pair_roc_data <- data.frame()
pair_box_data <- data.frame()
for (j in seq_along(pairs_list)) {
  cls1 <- pairs_list[[j]][1]; cls2 <- pairs_list[[j]][2]; pair_name <- pair_names[j]
  for (set_name in c("Train", "Test")) {
    if (set_name == "Train") { y_true <- y_train; score_raw <- oof_scores[, j] } else
    { y_true <- y_test;  score_raw <- test_scores[, j] }
    keep <- as.character(y_true) %in% c(cls1, cls2)
    y_sub <- as.character(y_true[keep]); sc_sub <- score_raw[keep]
    label <- as.numeric(y_sub == cls1)
    if (length(unique(label)) < 2) next
    r <- pROC::roc(label, sc_sub, levels = c(0, 1), direction = "<", quiet = TRUE)
    ci <- safe_ci(r); se <- safe_se(r)
    pair_roc_data <- rbind(pair_roc_data, data.frame(
      Pair = pair_name, Set = set_name,
      FDR = 1 - r$specificities, TPR = r$sensitivities,
      AUC = as.numeric(r$auc),
      AUC_L95 = ci[1], AUC_U95 = ci[3], AUC_SE = se))
    pair_box_data <- rbind(pair_box_data, data.frame(
      Pair = pair_name, Set = set_name, Score = sc_sub,
      True = ifelse(label == 1, "Positive", "Negative")))
  }
}

pair_plot_list <- list()
for (k in seq_along(pairs_list)) {
  pair <- pairs_list[[k]]; pair_name <- pair_names[k]
  roc_sub <- pair_roc_data[pair_roc_data$Pair == pair_name, ]
  box_sub <- pair_box_data[pair_box_data$Pair == pair_name, ]
  if (nrow(roc_sub) == 0 || nrow(box_sub) == 0) next
  box_sub$Set <- factor(box_sub$Set, levels = c("Train", "Test"))
  
  tr_row <- roc_sub[roc_sub$Set == "Train", , drop = FALSE]
  te_row <- roc_sub[roc_sub$Set == "Test",  , drop = FALSE]
  tr_row <- if (nrow(tr_row) > 0) tr_row[1, ] else NULL
  te_row <- if (nrow(te_row) > 0) te_row[1, ] else NULL
  
  legend_labels <- c()
  if (!is.null(tr_row)) {
    legend_labels["Train"] <- sprintf("Train: %.2f (95%% CI %.2f-%.2f)", tr_row$AUC, tr_row$AUC_L95, tr_row$AUC_U95)
  } else {
    legend_labels["Train"] <- "Train: NA"
  }
  if (!is.null(te_row)) {
    legend_labels["Test"] <- sprintf("Test: %.2f (95%% CI %.2f-%.2f)", te_row$AUC, te_row$AUC_L95, te_row$AUC_U95)
  } else {
    legend_labels["Test"] <- "Test: NA"
  }
  
  p_roc <- ggplot(roc_sub, aes(FDR, TPR, color = Set)) +
    geom_abline(linetype = "dashed", color = "gray70", linewidth = 0.6) +
    geom_line(linewidth = 1.5) +
    scale_color_manual(values = PAL_SET, labels = legend_labels) +
    labs(title = paste0(tools::toTitleCase(pair[1]), " vs ", tools::toTitleCase(pair[2])),
         x = "1 - Specificity", y = "TPR") +
    theme_top(15) +
    theme(
      legend.position = c(0.98, 0.02),
      legend.justification = c("right", "bottom"),
      legend.background = element_rect(fill = alpha("white", 0.8), color = "gray50", linewidth = 0.5),
      legend.title = element_blank(),
      legend.key.width = unit(1.2, "cm"),
      legend.text = element_text(size = 10)
    )
  
  p_box <- ggplot(box_sub, aes("Score", Score, fill = True)) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA, linewidth = 0.8, width = 0.6) +
    geom_jitter(aes(color = True), width = 0.15, height = 0, size = 2.5, alpha = 0.9) +
    facet_wrap(~Set, ncol = 1, scales = "free_x") +
    scale_fill_manual(values = PAL_PN) + scale_color_manual(values = PAL_PN) +
    labs(x = NULL, y = NULL) +
    theme_top(15) +
    theme(legend.position = "none",
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 10),
          axis.title.y = element_text(size = 12))
  
  pair_plot_list <- c(pair_plot_list,
                      list(p_roc + p_box + plot_layout(widths = c(3, 1))))
}

pair_plot_list <- Filter(Negate(is.null), pair_plot_list)
if (length(pair_plot_list) > 0)
  ggsave(file.path(plot_dir, "Pairwise_ROC_6panels.pdf"),
         wrap_plots(pair_plot_list, ncol = min(3, length(pair_plot_list))),
         width = 20, height = 12, dpi = 600)

# =========================================================
# 混淆矩阵
# =========================================================
acc_train <- round(cm_train$overall["Accuracy"] * 100, 0)
acc_test  <- round(cm_test$overall["Accuracy"]  * 100, 0)
make_cm_df <- function(cm, tt, lvls) {
  df <- as.data.frame(cm$table)
  df$Reference  <- factor(df$Reference, levels = lvls)
  df$Prediction <- factor(df$Prediction, levels = rev(lvls))
  df$Title <- tt
  df <- df %>% group_by(Title, Reference) %>%
    mutate(Pct = Freq / sum(Freq) * 100) %>% ungroup()
  df$Label <- paste0(df$Freq, "\n(", round(df$Pct, 1), "%)")
  df
}
tt_tr <- paste0("Train  |  Accuracy ", acc_train, "%")
tt_te <- paste0("Test  |  Accuracy ",  acc_test,  "%")
cm_all <- rbind(make_cm_df(cm_train, tt_tr, classes),
                make_cm_df(cm_test,  tt_te, classes))
cm_all$Title <- factor(cm_all$Title, levels = c(tt_tr, tt_te))
cm_all$TextColor <- ifelse(cm_all$Freq > 15, "white", "black")

p_cm <- ggplot(cm_all, aes(Reference, Prediction, fill = Freq)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = Label, color = TextColor), size = 5,
            fontface = "bold", lineheight = 0.8) +
  scale_fill_gradient(low = "#E0F2F1", high = "#004D40",
                      trans = "log1p", guide = "none") +
  scale_color_identity() + facet_wrap(~Title, ncol = 2) +
  labs(x = "Reference", y = "Prediction") + theme_top(16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(plot_dir, "Confusion_Matrix.pdf"), p_cm,
       width = 12, height = 7, dpi = 600)

# =========================================================
# ★ 每类完整分类指标（含 LR+ / LR-），供审稿人表格使用
#   Sensitivity / Specificity / PPV / NPV / LR+ / LR- / Balanced Acc / F1
# =========================================================
export_cm_metrics <- function(cm, set_name) {
  df <- as.data.frame(cm$byClass)
  df$Class <- rownames(df)
  df$Set   <- set_name
  
  df$PPV    <- df$`Pos Pred Value`
  df$NPV    <- df$`Neg Pred Value`
  df$LR_Pos <- safe_lr_pos(df$Sensitivity, df$Specificity)
  df$LR_Neg <- safe_lr_neg(df$Sensitivity, df$Specificity)
  
  df[, c("Set", "Class",
         "Sensitivity", "Specificity", "PPV", "NPV",
         "LR_Pos", "LR_Neg",
         "Balanced Accuracy", "F1",
         "Precision", "Recall")]
}

per_class_all <- rbind(
  export_cm_metrics(cm_train, "Train"),
  export_cm_metrics(cm_test,  "Test")
)
write.csv(per_class_all,
          file.path(outdir, "PerClass_metrics_train_test.csv"),
          row.names = FALSE)
cat("\n===== PerClass 指标（每类：Sens/Spec/PPV/NPV/LR+/LR-/BA/F1） =====\n")
print(per_class_all, row.names = FALSE)

# =========================================================
# OvR ROC（Test，含 DeLong 95% CI + LR + Balanced Accuracy）
# =========================================================
ovr_roc_data <- data.frame(); ovr_metrics <- data.frame()
for (cls in classes) {
  y_bin <- as.numeric(y_test == cls)
  r <- pROC::roc(y_bin, rf_prob_test[, cls], levels = c(0, 1),
                 direction = "<", quiet = TRUE)
  auc_val <- as.numeric(r$auc)
  ci <- safe_ci(r); auc_se <- safe_se(r)
  
  best <- pROC::coords(r, "best", best.method = "youden", transpose = FALSE,
                       ret = c("threshold","sensitivity","specificity","ppv","npv","accuracy"))
  sens <- best$sensitivity; spec <- best$specificity
  ppv  <- best$ppv;        npv  <- best$npv
  acc  <- best$accuracy;   thr  <- best$threshold
  f1 <- ifelse(is.na(ppv) || is.na(sens) || (ppv + sens) == 0, NA,
               2 * ppv * sens / (ppv + sens))
  ba     <- (sens + spec) / 2
  lr_pos <- safe_lr_pos(sens, spec)
  lr_neg <- safe_lr_neg(sens, spec)
  
  ovr_metrics <- rbind(ovr_metrics, data.frame(
    Class = cls, AUC = round(auc_val, 3),
    AUC_L95 = round(ci[1], 3), AUC_U95 = round(ci[3], 3), AUC_SE = round(auc_se, 3),
    Threshold = round(thr, 4),
    Sensitivity = round(sens, 3), Specificity = round(spec, 3),
    PPV = round(ppv, 3), NPV = round(npv, 3),
    LR_Pos = round(lr_pos, 3), LR_Neg = round(lr_neg, 3),
    Balanced_Accuracy = round(ba, 3),
    F1 = round(f1, 3), Accuracy = round(acc, 3)))
  
  ovr_roc_data <- rbind(ovr_roc_data, data.frame(
    Class = cls, FDR = 1 - r$specificities, TPR = r$sensitivities, AUC = auc_val))
}
write.csv(ovr_metrics, file.path(outdir, "OvR_metrics_test.csv"), row.names = FALSE)

ovr_metrics <- ovr_metrics %>%
  mutate(Label = sprintf("%s vs Rest (AUC=%.2f, CI: %.2f-%.2f)",
                         tools::toTitleCase(Class), AUC, AUC_L95, AUC_U95))
ovr_roc_data$Label <- factor(paste0(tools::toTitleCase(ovr_roc_data$Class), " vs Rest"),
                             levels = paste0(tools::toTitleCase(classes), " vs Rest"))
pal_ovr <- setNames(PAL[seq_along(classes)], paste0(tools::toTitleCase(classes), " vs Rest"))

p_ovr <- ggplot(ovr_roc_data, aes(FDR, TPR, color = Label)) +
  geom_abline(linetype = "dashed", color = "gray70", linewidth = 0.6) +
  geom_line(linewidth = 1.5) +
  scale_color_manual(values = pal_ovr, labels = ovr_metrics$Label) +
  scale_x_continuous(expand = c(0.02, 0.02)) +
  scale_y_continuous(expand = c(0.02, 0.02)) +
  labs(title = "One-vs-Rest ROC (Test Set)",
       x = "1 - Specificity", y = "Sensitivity", color = NULL) +
  theme_top(16) +
  theme(legend.position = "bottom",
        legend.key.width = unit(1.2, "cm"),
        legend.text = element_text(size = 11)) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))

ggsave(file.path(plot_dir, "OvR_ROC_4panels.pdf"), p_ovr,
       width = 10.5, height = 8.5, dpi = 600)

# =========================================================
# ★ OvO 参数导出（含 LR+ / LR-）
# =========================================================
pair_metrics <- data.frame()
for (j in seq_along(pairs_list)) {
  cls1 <- pairs_list[[j]][1]; cls2 <- pairs_list[[j]][2]; pair_name <- pair_names[j]
  for (set_name in c("Train", "Test")) {
    if (set_name == "Train") { y_true <- y_train; score_raw <- oof_scores[, j] } else
    { y_true <- y_test;  score_raw <- test_scores[, j] }
    keep   <- as.character(y_true) %in% c(cls1, cls2)
    y_sub  <- as.character(y_true[keep]); sc_sub <- score_raw[keep]
    label  <- as.numeric(y_sub == cls1)
    n_pos <- sum(label == 1); n_neg <- sum(label == 0)
    if (length(unique(label)) < 2) {
      pair_metrics <- rbind(pair_metrics, data.frame(
        Pair = pair_name, Set = set_name, N = length(label),
        N_pos = n_pos, N_neg = n_neg,
        AUC = NA, AUC_L95 = NA, AUC_U95 = NA, AUC_SE = NA, Threshold = NA,
        Sensitivity = NA, Specificity = NA, PPV = NA, NPV = NA,
        LR_Pos = NA, LR_Neg = NA,
        F1 = NA, Accuracy = NA, Balanced_Accuracy = NA, Youden = NA,
        Acc_at_0.5 = NA, F1_at_0.5 = NA, LogLoss = NA, Brier = NA))
      next
    }
    r <- pROC::roc(label, sc_sub, levels = c(0, 1), direction = "<", quiet = TRUE)
    auc_val <- as.numeric(r$auc); ci <- safe_ci(r); auc_se <- safe_se(r)
    
    best <- pROC::coords(r, "best", best.method = "youden", transpose = FALSE,
                         ret = c("threshold","sensitivity","specificity",
                                 "ppv","npv","accuracy"))
    sens <- best$sensitivity; spec <- best$specificity
    ppv  <- best$ppv;        npv  <- best$npv
    acc  <- best$accuracy;   thr  <- best$threshold
    f1 <- ifelse(is.na(ppv) || is.na(sens) || (ppv + sens) == 0, NA,
                 2 * ppv * sens / (ppv + sens))
    ba     <- (sens + spec) / 2
    youden <- sens + spec - 1
    lr_pos <- safe_lr_pos(sens, spec)
    lr_neg <- safe_lr_neg(sens, spec)
    
    pred05 <- as.integer(sc_sub >= 0.5)
    tp <- sum(pred05 == 1 & label == 1); fp <- sum(pred05 == 1 & label == 0)
    fn <- sum(pred05 == 0 & label == 1); tn <- sum(pred05 == 0 & label == 0)
    acc05 <- (tp + tn) / length(label)
    ppv05 <- ifelse((tp + fp) == 0, NA, tp / (tp + fp))
    sen05 <- ifelse((tp + fn) == 0, NA, tp / (tp + fn))
    f1_05 <- ifelse(is.na(ppv05) || is.na(sen05) || (ppv05 + sen05) == 0, NA,
                    2 * ppv05 * sen05 / (ppv05 + sen05))
    eps <- 1e-15; p_hat <- pmin(pmax(sc_sub, eps), 1 - eps)
    ll <- -mean(label * log(p_hat) + (1 - label) * log(1 - p_hat))
    brier <- mean((sc_sub - label)^2)
    
    pair_metrics <- rbind(pair_metrics, data.frame(
      Pair = pair_name, Set = set_name, N = length(label),
      N_pos = n_pos, N_neg = n_neg,
      AUC = round(auc_val, 4), AUC_L95 = round(ci[1], 4), AUC_U95 = round(ci[3], 4),
      AUC_SE = round(auc_se, 4), Threshold = round(thr, 4),
      Sensitivity = round(sens, 4), Specificity = round(spec, 4),
      PPV = round(ppv, 4), NPV = round(npv, 4),
      LR_Pos = round(lr_pos, 4), LR_Neg = round(lr_neg, 4),
      F1 = round(f1, 4),
      Accuracy = round(acc, 4), Balanced_Accuracy = round(ba, 4),
      Youden = round(youden, 4), Acc_at_0.5 = round(acc05, 4),
      F1_at_0.5 = round(f1_05, 4), LogLoss = round(ll, 4), Brier = round(brier, 4)))
  }
}
write.csv(pair_metrics, file.path(outdir, "OvO_metrics_train_test.csv"), row.names = FALSE)
pair_wide <- reshape(pair_metrics, idvar = "Pair", timevar = "Set",
                     direction = "wide", sep = "_")
write.csv(pair_wide, file.path(outdir, "OvO_metrics_wide.csv"), row.names = FALSE)
if (nrow(pair_roc_data) > 0)
  write.csv(pair_roc_data, file.path(outdir, "OvO_ROC_curve_points.csv"), row.names = FALSE)
if (nrow(pair_box_data) > 0)
  write.csv(pair_box_data, file.path(outdir, "OvO_score_distribution.csv"), row.names = FALSE)
write.csv(as.data.frame(oof_scores),  file.path(outdir, "OvO_oof_scores_train.csv"), row.names = FALSE)
write.csv(as.data.frame(test_scores), file.path(outdir, "OvO_test_scores.csv"),      row.names = FALSE)

cat("\n===== OvO Test 指标 =====\n")
print(pair_metrics[pair_metrics$Set == "Test", ], row.names = FALSE)
cat("\n===== OvO AUC Train vs Test =====\n")
print(reshape(pair_metrics[, c("Pair","Set","AUC")],
              idvar = "Pair", timevar = "Set", direction = "wide"), row.names = FALSE)

# =========================================================
# SHAP
# =========================================================
X_train_mat <- as.matrix(X_train)
for (k in seq_along(pairs_list)) {
  pair <- pairs_list[[k]]; cls1 <- pair[1]; cls2 <- pair[2]; pair_name <- pair_names[[k]]
  in_pair <- as.character(y_train) %in% c(cls1, cls2)
  X_sub <- X_train_mat[in_pair, , drop = FALSE]; n_feat <- ncol(X_sub)
  contrib <- predict(pair_models[[pair_name]], X_sub, predcontrib = TRUE)
  contrib <- contrib[, 1:n_feat, drop = FALSE]
  top_feat <- names(sort(colMeans(abs(contrib)), decreasing = TRUE))[1:min(15, n_feat)]
  long_df <- data.frame(
    feature = rep(top_feat, each = nrow(contrib)),
    shap    = as.vector(contrib[, top_feat, drop = FALSE]),
    value   = as.vector(scale(X_sub[, top_feat, drop = FALSE])))
  long_df$feature <- factor(long_df$feature, levels = rev(top_feat))
  p <- ggplot(long_df, aes(shap, feature, color = value)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray70", linewidth = 0.6) +
    ggbeeswarm::geom_quasirandom(size = 0.8, alpha = 0.6, width = 0.2) +
    scale_color_gradientn(colors = c("#377EB8", "#F7F7F7", "#E41A1C"),
                          name = "Feature\nvalue") +
    labs(title = paste0("SHAP: ", tools::toTitleCase(cls1), " vs ", tools::toTitleCase(cls2)),
         x = "SHAP value", y = NULL) +
    theme_top(15)
  ggsave(file.path(plot_dir, paste0("SHAP_GB_", cls1, "_vs_", cls2, ".pdf")),
         p, width = 9, height = 7, dpi = 600)
}

# =========================================================
# 特征重要性
# =========================================================
imp_rf <- as.data.frame(importance(rf_model)); imp_rf$Feature <- rownames(imp_rf)
imp_rf$Type <- ifelse(imp_rf$Feature %in% pair_names, "Pairwise score", "Species")
imp_rf <- imp_rf[order(-imp_rf$MeanDecreaseAccuracy), ]
imp_score <- imp_rf[imp_rf$Type == "Pairwise score", ]
p_imp_score <- ggplot(imp_score, aes(reorder(Feature, MeanDecreaseAccuracy),
                                     MeanDecreaseAccuracy)) +
  geom_col(fill = "#377EB8", width = 0.65) + coord_flip() +
  geom_text(aes(label = round(MeanDecreaseAccuracy, 1)),
            hjust = -0.2, size = 5, fontface = "bold") +
  labs(title = "Pairwise score importance", x = NULL, y = "Mean Decrease Accuracy") +
  theme_top(15) + expand_limits(y = max(imp_score$MeanDecreaseAccuracy) * 1.1)
ggsave(file.path(plot_dir, "RF_importance_6scores.pdf"), p_imp_score,
       width = 9, height = 5, dpi = 600)

imp_top <- head(imp_rf[imp_rf$Type == "Species", ], 20)
p_imp_sp <- ggplot(imp_top, aes(reorder(Feature, MeanDecreaseAccuracy), MeanDecreaseAccuracy)) +
  geom_col(fill = "#4DAF4A", width = 0.65) + coord_flip() +
  geom_text(aes(label = round(MeanDecreaseAccuracy, 1)),
            hjust = -0.2, size = 5, fontface = "bold") +
  labs(title = "Top 20 species importance", x = NULL, y = "Mean Decrease Accuracy") +
  theme_top(15) + expand_limits(y = max(imp_top$MeanDecreaseAccuracy) * 1.1)
ggsave(file.path(plot_dir, "RF_importance_species.pdf"), p_imp_sp,
       width = 9, height = 7, dpi = 600)

# =========================================================
# Permutation importance
# =========================================================
base_acc <- mean(rf_test_pred == y_test)
perm_imp <- data.frame(Feature = pair_names, Drop = NA_real_)
set.seed(42)
n_perm <- 15
for (k in seq_along(pair_names)) {
  feat <- pair_names[k]
  drops <- numeric(n_perm)
  for (i in seq_len(n_perm)) {
    X_test_perm <- X_test_meta
    X_test_perm[[feat]] <- sample(X_test_perm[[feat]])
    prob_perm <- matrix(0, nrow(X_test_perm), length(classes))
    colnames(prob_perm) <- classes
    for (s in seeds) {
      set.seed(s)
      rf_s <- randomForest(
        x = X_train_meta[, all_feats, drop = FALSE], y = y_train,
        ntree = 1500, mtry = 3, nodesize = 3, maxnodes = 60,
        classwt = class_weights, importance = FALSE)
      prob_perm <- prob_perm +
        predict(rf_s, X_test_perm[, all_feats, drop = FALSE], type = "prob")[, classes]
    }
    prob_perm <- prob_perm / length(seeds)
    pred_perm <- factor(classes[apply(prob_perm, 1, which.max)], levels = classes)
    drops[i] <- base_acc - mean(pred_perm == y_test)
  }
  perm_imp$Drop[k] <- mean(drops)
}
perm_imp <- perm_imp[order(-perm_imp$Drop), ]
p_rf_perm <- ggplot(perm_imp, aes(reorder(Feature, Drop), Drop)) +
  geom_col(fill = "#E41A1C", width = 0.65) + coord_flip() +
  geom_text(aes(label = round(Drop, 3)), hjust = -0.2, size = 5, fontface = "bold") +
  labs(title = "Permutation importance (test set)",
       subtitle = "Accuracy drop after shuffling each score",
       x = NULL, y = "Accuracy Drop") +
  theme_top(15) + expand_limits(y = max(perm_imp$Drop) * 1.1)
ggsave(file.path(plot_dir, "RF_permutation_importance.pdf"), p_rf_perm,
       width = 9, height = 5, dpi = 600)

# =========================================================
# 稳健性评估（新增 Macro_Balanced_Accuracy）
# =========================================================
if (RUN_ROBUSTNESS) {
  cat("\n########## 稳健性评估：", N_REP, "次重复划分 ##########\n")
  rob_rows <- list()
  rob_seeds <- 42 + seq_len(N_REP) - 1
  for (s in rob_seeds) {
    cat(sprintf("  [%d/%d] seed = %d ... ",
                which(rob_seeds == s), N_REP, s))
    r <- fit_stacking(X_clr, y, top_species = top_species, seed = s, verbose = FALSE)
    cm <- confusionMatrix(r$rf_test_pred, r$y_test)
    byc <- cm$byClass
    get_rec <- function(nm) {
      idx <- grep(nm, rownames(byc), ignore.case = TRUE)
      if (length(idx) == 0) NA_real_ else byc[idx[1], "Sensitivity"]
    }
    rob_rows[[length(rob_rows) + 1]] <- data.frame(
      Seed = s,
      Test_Accuracy = as.numeric(cm$overall["Accuracy"]),
      Test_Kappa    = as.numeric(cm$overall["Kappa"]),
      Macro_F1      = mean(byc[, "F1"], na.rm = TRUE),
      Macro_Recall  = mean(byc[, "Sensitivity"], na.rm = TRUE),
      Macro_Balanced_Accuracy = mean(byc[, "Balanced Accuracy"], na.rm = TRUE),
      Bact_Recall   = get_rec("bacter"),
      Fung_Recall   = get_rec("fung"),
      Viral_Recall  = get_rec("viral"))
    cat(sprintf("Acc=%.3f\n", rob_rows[[length(rob_rows)]]$Test_Accuracy))
    rm(r); gc()
  }
  robust_df <- bind_rows(rob_rows)
  write.csv(robust_df, file.path(outdir, "robustness_repeated_holdout.csv"), row.names = FALSE)
  summary_df <- robust_df %>%
    pivot_longer(-Seed, names_to = "Metric", values_to = "Value") %>%
    group_by(Metric) %>%
    summarise(Mean = mean(Value, na.rm = TRUE), SD = sd(Value, na.rm = TRUE),
              Min = min(Value, na.rm = TRUE), Max = max(Value, na.rm = TRUE),
              .groups = "drop")
  write.csv(summary_df, file.path(outdir, "robustness_summary.csv"), row.names = FALSE)
  cat("\n===== 稳健性汇总 =====\n")
  print(summary_df, row.names = FALSE)
}

# =========================================================
# 保存
# =========================================================
saveRDS(list(pair_models = pair_models, rf_model = rf_model,
             rf_prob_train = rf_prob_train, rf_prob_test = rf_prob_test,
             oof_scores = oof_scores, test_scores = test_scores,
             pair_names = pair_names, top_species = top_species,
             all_feats = all_feats, seeds = seeds,
             levels = levels(y_train), y_train = y_train, y_test = y_test),
        file.path(outdir, "model_bundle_stacking_hybrid.rds"))
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))

cat("\n=====================================\n")
cat("完成！\n")
cat("图片目录：", plot_dir, "\n")
cat("结果目录：", outdir, "\n")
cat("=====================================\n")
cat("【图片】\n"); print(list.files(plot_dir))
cat("\n【CSV】\n"); print(list.files(outdir, pattern = "\\.csv$"))
cat("=====================================\n")