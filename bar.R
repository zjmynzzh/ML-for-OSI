# ============================================
# 三组数据（BV / BF / FV）中 23 种算法的
# 平均 AUC ± SD 汇总表
# 输出：控制台表格 + 桌面 CSV
# ============================================

if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
library(dplyr)

# ============================================
# 步骤1：定位桌面目录
# ============================================

desktop_dir <- local({
  candidates <- c(
    file.path(Sys.getenv("USERPROFILE"), "Desktop"),
    file.path(Sys.getenv("USERPROFILE"), "桌面"),
    file.path(Sys.getenv("HOME"), "Desktop"),
    file.path(Sys.getenv("HOME"), "桌面")
  )
  candidates <- candidates[nzchar(candidates)]
  hit <- candidates[dir.exists(candidates)]
  if (length(hit) == 0) stop("未找到桌面目录，请手动修改 desktop_dir")
  hit[1]
})

cat("桌面目录：", desktop_dir, "\n\n")

# ============================================
# 步骤2：读取单个 CSV 并统一列名
# ============================================

read_one <- function(fname, dataset_name) {
  path <- file.path(desktop_dir, fname)
  if (!file.exists(path)) stop("找不到文件：", path)
  
  d <- read.csv(
    path,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    fileEncoding = "UTF-8"
  )
  
  names(d)[1:2] <- c("Model", "Training")
  d <- d[, c("Model", "Training")]
  
  d %>%
    mutate(
      Model    = trimws(as.character(Model)),
      Training = as.numeric(trimws(as.character(Training))),
      Dataset  = dataset_name
    ) %>%
    filter(!is.na(Model), Model != "", !is.na(Training))
}

# ============================================
# 步骤3：合并三个数据集
# ============================================

df_all <- bind_rows(
  read_one("BV.csv", "BV"),
  read_one("BF.csv", "BF"),
  read_one("FV.csv", "FV")
)

# ============================================
# 步骤4：按算法计算均值、SD、n
# ============================================

summary_df <- df_all %>%
  group_by(Model) %>%
  summarise(
    Mean = mean(Training),
    SD   = sd(Training),        # 样本标准差，n = 3
    n    = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(Mean)) %>%
  mutate(
    `Mean±SD` = sprintf("%.3f ± %.3f", Mean, SD)
  )

# ============================================
# 步骤5：控制台打印规整表格
# ============================================

cat("==================== 汇总表（按均值降序） ====================\n\n")

print_table <- summary_df %>%
  transmute(
    算法     = Model,
    均值     = sprintf("%.3f", Mean),
    标准差   = sprintf("%.3f", SD),
    `均值±SD` = `Mean±SD`,
    n        = n
  )

print(as.data.frame(print_table), row.names = FALSE)

cat("\n")

# ============================================
# 步骤6：导出 CSV（Excel 可直接打开）
# ============================================

csv_out <- file.path(desktop_dir, "Algorithm_AUC_mean_sd_table.csv")

# 用 UTF-8 BOM 便于 Excel 正确识别中文
con <- file(csv_out, open = "wb")
writeBin(charToRaw("\xEF\xBB\xBF"), con)  # BOM
close(con)

write.table(
  print_table,
  file = csv_out,
  sep = ",",
  row.names = FALSE,
  col.names = TRUE,
  quote = TRUE,
  fileEncoding = "UTF-8",
  append = TRUE
)

cat("表格已保存到：", csv_out, "\n")