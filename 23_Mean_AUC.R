# ============================================
# 棒棒糖图：模型性能对比（Training组）
# 数据：桌面 FV.csv
# 最大值在上，实心圆点，横线从0开始，右侧渐变色条
# 圆点右侧显示数值
# 只输出 PDF（cairo_pdf，避免 useDingbats 报错）
# ============================================

if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")

library(ggplot2)
library(dplyr)

# ============================================
# 步骤1：自动查找桌面上的 FV.csv
# ============================================

csv_path <- local({
  candidates <- c(
    file.path(Sys.getenv("USERPROFILE"), "Desktop", "FV.csv"),
    file.path(Sys.getenv("USERPROFILE"), "桌面", "FV.csv"),
    file.path(Sys.getenv("HOME"), "Desktop", "FV.csv"),
    file.path(Sys.getenv("HOME"), "桌面", "FV.csv")
  )
  candidates <- candidates[nzchar(candidates)]
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) file.choose() else hit[1]
})

cat("正在读取：", csv_path, "\n")

# ============================================
# 步骤2：读取并整理数据
# ============================================

model_data <- read.csv(
  csv_path,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8"
)

names(model_data)[1:2] <- c("Model", "Training")
model_data <- model_data[, c("Model", "Training")]

model_data <- model_data %>%
  mutate(
    Model = trimws(as.character(Model)),
    Training = as.numeric(trimws(as.character(Training)))
  ) %>%
  filter(!is.na(Model), Model != "", !is.na(Training))

# ============================================
# 步骤3：排序（最大值在上）
# ============================================

model_data <- model_data %>%
  arrange(desc(Training)) %>%
  mutate(Model = factor(Model, levels = rev(Model)))

# ============================================
# 步骤4：绘制棒棒糖图（带数值标签）
# ============================================

p <- ggplot(model_data, aes(x = Training, y = Model)) +
  # 连接线：从0开始
  geom_segment(
    aes(x = 0, xend = Training, y = Model, yend = Model),
    color = "gray70", linewidth = 0.8
  ) +
  # 实心圆点，颜色映射 Training
  geom_point(
    aes(color = Training),
    size = 5, stroke = 0, shape = 16
  ) +
  # 数值标签：标在圆点右侧
  geom_text(
    aes(x = Training + 0.025, label = sprintf("%.2f", Training)),
    hjust = 0,
    size = 3.5,
    fontface = "bold",
    color = "black"
  ) +
  # 颜色渐变 + 右侧色条图例
  scale_color_viridis_c(
    option = "D",
    direction = -1,
    begin = 0.2,
    end = 0.9,
    name = "AUC\n(Training)",
    guide = guide_colorbar(
      barwidth = 1.2,
      barheight = 12,
      title.position = "top",
      title.theme = element_text(size = 10, face = "bold"),
      label.theme = element_text(size = 9),
      frame.colour = "black",
      frame.linewidth = 0.5,
      ticks = TRUE
    )
  ) +
  # X轴：上限放宽到 1.1，给数值留空间
  scale_x_continuous(
    limits = c(0, 1.1),
    breaks = seq(0, 1.0, by = 0.2)
  ) +
  # 主题
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 11, color = "black"),
    axis.text.y = element_text(face = "bold"),
    legend.position = "right",
    legend.margin = margin(0, 0, 0, 10),
    plot.margin = margin(15, 15, 15, 15)
  ) +
  labs(x = "AUC (Training)")

print(p)

# ============================================
# 步骤5：只保存 PDF 到 FV.csv 所在目录
# ============================================

out_dir  <- dirname(csv_path)
pdf_file <- file.path(out_dir, "FV_Training_solid.pdf")

cairo_pdf(filename = pdf_file, width = 7, height = 8)
print(p)
dev.off()

cat("PDF已保存到：", pdf_file, "\n")