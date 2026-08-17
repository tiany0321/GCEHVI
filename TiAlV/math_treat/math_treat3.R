# =======================================================
# 🔍 多目标优化结果可视化程序（P1 文件夹）
# 第一行：美观箱线图（矩形表示 ± SE，胡须表示 ± SD）
# 第二行：迭代次数差值 PDF
# 第三行：OC 曲线 + HV 曲线（± SD 阴影）
# =======================================================

setwd("D:\\work\\paper\\CMOBO\\code\\TiAlV\\repeat")

library(tidyverse)
library(patchwork)
library(RColorBrewer)
library(purrr)

# =======================================================
# 全局颜色（与 OC / HV 完全一致）
# =======================================================

method_levels <- c("all_EHVI", "GEHVI", "DEHVI")

method_colors <- setNames(
  brewer.pal(3, "Set1"),
  method_levels
)

# =======================================================
# 工具函数：读取单个文件
# =======================================================

read_one <- function(filename, name, folder_path) {
  full_path <- file.path(folder_path, filename)
  if (!file.exists(full_path)) {
    warning(paste("⚠️ 文件不存在:", full_path))
    return(NULL)
  }
  dat <- read.csv(full_path, row.names = 1)
  vals <- as.numeric(unlist(dat, use.names = FALSE))
  tibble(method = name, value = vals)
}

# =======================================================
# 第一行：箱线图（矩形 ± SE，胡须 ± SD，胡须穿矩形部分隐藏）
# =======================================================

generate_box <- function(folder_path) {
  
  files <- c(
    "all_EHVI-5count.csv",
    "given_constrain_CEHVI-5count.csv",
    "3points_dist_CEHVI-5count.csv"
  )
  
  df <- map2(files, method_levels, ~read_one(.x, .y, folder_path)) %>%
    compact() %>%
    bind_rows() %>%
    mutate(method = factor(method, levels = method_levels))
  
  stat_df <- df %>%
    group_by(method) %>%
    summarise(
      mean = mean(value),
      sd   = sd(value),
      se   = sd / sqrt(n()),
      .groups = "drop"
    ) %>%
    mutate(
      # 胡须上下端点（避开矩形 ±SE）
      whisker_top    = mean + sd,
      whisker_bottom = mean - sd,
      rect_top       = mean + se * 2,
      rect_bottom    = mean - se * 2
    )
  
  ggplot(stat_df, aes(method, mean, fill = method)) +
    # ---------------------------
  # 矩形 ± SE
  # ---------------------------
  geom_rect(aes(
    xmin = as.numeric(method) - 0.25,
    xmax = as.numeric(method) + 0.25,
    ymin = rect_bottom,
    ymax = rect_top
  ),
  color = "black",
  linewidth = 0.9,
  alpha = 0.85
  ) +
    # ---------------------------
  # 上胡须（垂直线 + 顶部横线）
  # ---------------------------
  geom_segment(aes(
    x = as.numeric(method),
    xend = as.numeric(method),
    y = rect_top,
    yend = whisker_top
  ),
  linewidth = 1.1
  ) +
    geom_segment(aes(
      x = as.numeric(method) - 0.1,
      xend = as.numeric(method) + 0.1,
      y = whisker_top,
      yend = whisker_top
    ),
    linewidth = 1.1
    ) +
    # ---------------------------
  # 下胡须（垂直线 + 底部横线）
  # ---------------------------
  geom_segment(aes(
    x = as.numeric(method),
    xend = as.numeric(method),
    y = whisker_bottom,
    yend = rect_bottom
  ),
  linewidth = 1.1
  ) +
    geom_segment(aes(
      x = as.numeric(method) - 0.1,
      xend = as.numeric(method) + 0.1,
      y = whisker_bottom,
      yend = whisker_bottom
    ),
    linewidth = 1.1
    ) +
    # ---------------------------
  # 中心线表示均值
  # ---------------------------
  geom_segment(aes(
    x = as.numeric(method) - 0.25,
    xend = as.numeric(method) + 0.25,
    y = mean,
    yend = mean
  ),
  linewidth = 1.1
  ) +
    scale_fill_manual(values = method_colors) +
    theme_classic(base_size = 20) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      legend.position = "none"
    ) +
    labs(x = "Method", y = "Iteration count")
}

# =======================================================
# 第二行：PDF（迭代次数差值）
# =======================================================

generate_pdf_diff <- function(folder_path) {
  
  files <- c(
    "all_EHVI-5count.csv",
    "given_constrain_CEHVI-5count.csv",
    "3points_dist_CEHVI-5count.csv"
  )
  
  df <- map2(files, method_levels, ~read_one(.x, .y, folder_path)) %>%
    compact() %>%
    bind_rows()
  
  df_wide <- df %>%
    pivot_wider(names_from = method, values_from = value, values_fn = list) %>%
    unnest(cols = everything())
  
  df_diff <- df_wide %>%
    mutate(
      GEHVI_diff = GEHVI - all_EHVI,
      DEHVI_diff = DEHVI - all_EHVI
    ) %>%
    dplyr::select(GEHVI_diff, DEHVI_diff) %>%
    pivot_longer(everything(), names_to = "method", values_to = "diff") %>%
    mutate(
      method = recode(method,
                      GEHVI_diff = "GEHVI",
                      DEHVI_diff = "DEHVI"),
      method = factor(method, levels = method_levels[-1])
    )
  
  ggplot(df_diff, aes(diff, color = method, fill = method)) +
    geom_density(alpha = 0.3, linewidth = 1.2) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_color_manual(values = method_colors) +
    scale_fill_manual(values = method_colors) +
    theme_minimal(base_size = 20) +
    theme(
      legend.position = c(0.1, 0.9),  # 左上角
      legend.justification = c(0, 1),
      legend.background = element_rect(fill = "white", color = NA)
    ) +
    labs(
      x = "Iteration count difference",
      y = "Density",
      color = "Method",
      fill  = "Method"
    )
}

# =======================================================
# 第三行：OC & HV（± SD 阴影）
# =======================================================

generate_OC_HV <- function(folder_path) {
  
  OC <- read.csv("5ave.dis.OC.csv")
  HV <- read.csv("5ave.delta_hv.csv")
  
  colnames(OC)[1] <- "X"
  colnames(HV)[1] <- "X"
  
  df_OC <- map_df(method_levels, \(m)
                  tibble(
                    X = OC$X,
                    method = m,
                    mean = OC[[m]],
                    sd   = OC[[paste0(m, "_sd")]]
                  )
  )
  
  df_HV <- map_df(method_levels, \(m)
                  tibble(
                    X = HV$X,
                    method = m,
                    mean = HV[[m]],
                    sd   = HV[[paste0(m, "_sd")]]
                  )
  )
  
  p_OC <- ggplot(df_OC, aes(X, mean, color = method, fill = method)) +
    geom_ribbon(
      aes(ymin = mean - sd, ymax = mean + sd),
      alpha = 0.15,
      color = NA
    ) +
    geom_line(linewidth = 1.2) +
    scale_color_manual(values = method_colors) +
    scale_fill_manual(values = method_colors) +
    theme_minimal(base_size = 20) +
    theme(
      legend.position = c(0.98, 0.98),  # 右上角
      legend.justification = c(1, 1),
      legend.background = element_rect(fill = "white", color = NA)
    ) +
    labs(x = "Iteration", y = "OC")
  
  p_HV <- ggplot(df_HV, aes(X, mean, color = method, fill = method)) +
    geom_ribbon(
      aes(ymin = mean - sd, ymax = mean + sd),
      alpha = 0.15,
      color = NA
    ) +
    geom_line(linewidth = 1.2) +
    scale_color_manual(values = method_colors) +
    scale_fill_manual(values = method_colors) +
    theme_minimal(base_size = 20) +
    theme(
      legend.position = c(0.98, 0.02),  # 右下角
      legend.justification = c(1, 0),
      legend.background = element_rect(fill = "white", color = NA)
    ) +
    labs(x = "Iteration", y = "HV")
  
  list(p_OC, p_HV)
}

# =======================================================
# 组合输出
# =======================================================

p1 <- generate_box(getwd())
p2 <- generate_pdf_diff(getwd())
p_OC_HV <- generate_OC_HV(getwd())

final_plot <- wrap_plots(
  list(p1, p2, p_OC_HV[[1]], p_OC_HV[[2]]),
  ncol = 2
) +
  plot_annotation(tag_levels = "a")

print(final_plot)

ggsave(
  "iteration_OC_HV_box_PDF.png",
  final_plot,
  width = 14,
  height = 12,
  dpi = 300
)
