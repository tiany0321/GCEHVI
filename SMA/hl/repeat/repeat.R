# =======================================================
# HL 多目标优化结果可视化
#
# 三个图：
# a1：箱线图
# b1：PDF 差值
# c1：OC
# =======================================================


# =======================================================
# 1. 加载程序包
# =======================================================

library(ggplot2)
library(ggExtra)
library(tidyverse)
library(patchwork)
library(RColorBrewer)
library(purrr)
library(emoa)
library(cowplot)


# =======================================================
# 2. 文件路径
# =======================================================

folder_HL <- "C:\\Users\\12758\\Desktop\\GCEHVI\\hl\\repeat\\size0.25"


# =======================================================
# 3. 方法 & 颜色
# =======================================================

# 内部方法名称保持不变
method_levels <- c(
  "all_EHVI",
  "GEHVI",
  "DEHVI"
)


# 方法颜色
method_colors <- setNames(
  brewer.pal(
    3,
    "Set1"
  ),
  method_levels
)


# =======================================================
# 4. 显示名称
# =======================================================

method_labels <- c(
  all_EHVI = "FEHVI",
  GEHVI    = "FCEHVI",
  DEHVI    = "GCEHVI"
)


# 字体
tag_font <- "Arial"


# =======================================================
# 5. 工具函数：读取单个 CSV
# =======================================================

read_one <- function(
    filename,
    name,
    folder_path
){
  
  full_path <- file.path(
    folder_path,
    filename
  )
  
  
  # 检查文件是否存在
  if(!file.exists(full_path)){
    
    warning(
      "⚠️ 文件不存在: ",
      full_path
    )
    
    return(NULL)
  }
  
  
  # 读取数据
  dat <- read.csv(
    full_path,
    row.names = 1
  )
  
  
  # 转换为长格式
  tibble(
    method = name,
    value = as.numeric(
      unlist(dat)
    )
  )
}


# =======================================================
# 6. 图 A：箱线图
# =======================================================

generate_box <- function(
    folder_path
){
  
  # -----------------------------------------------------
  # 文件
  # -----------------------------------------------------
  
  files <- c(
    "all_EHVI-5count.csv",
    "given_constrain_CEHVI-5count.csv",
    "3points_dist_CEHVI-5count.csv"
  )
  
  
  # -----------------------------------------------------
  # 读取三个方法的数据
  # -----------------------------------------------------
  
  df <- map2(
    files,
    method_levels,
    ~read_one(
      .x,
      .y,
      folder_path
    )
  ) %>%
    
    compact() %>%
    
    bind_rows()
  
  
  # -----------------------------------------------------
  # 方法设置为因子
  # -----------------------------------------------------
  
  df <- df %>%
    mutate(
      method = factor(
        method,
        levels = method_levels
      )
    )
  
  
  # -----------------------------------------------------
  # 计算统计量
  # -----------------------------------------------------
  
  stat_df <- df %>%
    
    group_by(method) %>%
    
    summarise(
      
      mean = mean(value),
      
      sd = sd(value),
      
      se = sd / sqrt(n()),
      
      .groups = "drop"
    ) %>%
    
    mutate(
      
      whisker_top =
        mean + sd,
      
      whisker_bottom =
        mean - sd,
      
      rect_top =
        mean + 2 * se,
      
      rect_bottom =
        mean - 2 * se
    )
  
  
  # -----------------------------------------------------
  # 绘图
  # -----------------------------------------------------
  
  p <- ggplot(
    stat_df,
    aes(
      method,
      mean,
      fill = method
    )
  ) +
    
    
    # ---------------------------------------------------
  # 2 × SE 矩形
  # ---------------------------------------------------
  
  geom_rect(
    aes(
      xmin =
        as.numeric(method) - 0.25,
      
      xmax =
        as.numeric(method) + 0.25,
      
      ymin =
        rect_bottom,
      
      ymax =
        rect_top
    ),
    
    color = "black",
    
    alpha = 0.85
  ) +
    
    
    # ---------------------------------------------------
  # 上方误差线
  # ---------------------------------------------------
  
  geom_segment(
    aes(
      x =
        as.numeric(method),
      
      xend =
        as.numeric(method),
      
      y =
        rect_top,
      
      yend =
        whisker_top
    ),
    
    linewidth = 1.5
  ) +
    
    
    # ---------------------------------------------------
  # 上方横线
  # ---------------------------------------------------
  
  geom_segment(
    aes(
      x =
        as.numeric(method) - 0.1,
      
      xend =
        as.numeric(method) + 0.1,
      
      y =
        whisker_top,
      
      yend =
        whisker_top
    ),
    
    linewidth = 1.5
  ) +
    
    
    # ---------------------------------------------------
  # 下方误差线
  # ---------------------------------------------------
  
  geom_segment(
    aes(
      x =
        as.numeric(method),
      
      xend =
        as.numeric(method),
      
      y =
        whisker_bottom,
      
      yend =
        rect_bottom
    ),
    
    linewidth = 1.5
  ) +
    
    
    # ---------------------------------------------------
  # 下方横线
  # ---------------------------------------------------
  
  geom_segment(
    aes(
      x =
        as.numeric(method) - 0.1,
      
      xend =
        as.numeric(method) + 0.1,
      
      y =
        whisker_bottom,
      
      yend =
        whisker_bottom
    ),
    
    linewidth = 1.5
  ) +
    
    
    # ---------------------------------------------------
  # 平均值横线
  # ---------------------------------------------------
  
  geom_segment(
    aes(
      x =
        as.numeric(method) - 0.25,
      
      xend =
        as.numeric(method) + 0.25,
      
      y =
        mean,
      
      yend =
        mean
    ),
    
    linewidth = 1.5
  ) +
    
    
    # ---------------------------------------------------
  # 填充颜色
  # ---------------------------------------------------
  
  scale_fill_manual(
    values = method_colors,
    labels = method_labels
  ) +
    
    
    # ---------------------------------------------------
  # X 轴
  # ---------------------------------------------------
  
  scale_x_discrete(
    labels = method_labels
  ) +
    
    
    # ---------------------------------------------------
  # 主题
  # ---------------------------------------------------
  
  theme_classic(
    base_size = 45
  ) +
    
    
    theme(
      
      axis.text.x =
        element_text(
          angle = 15,
          hjust = 1,
          face = "bold"
        ),
      
      legend.position = "none"
    ) +
    
    
    # ---------------------------------------------------
  # 坐标轴标题
  # ---------------------------------------------------
  
  labs(
    x = "Method",
    y = "Iteration count"
  )
  
  
  return(p)
}


# =======================================================
# 7. 图 B：PDF 差值
# =======================================================

generate_pdf_diff <- function(
    folder_path
){
  
  # -----------------------------------------------------
  # 文件
  # -----------------------------------------------------
  
  files <- c(
    "all_EHVI-5count.csv",
    "given_constrain_CEHVI-5count.csv",
    "3points_dist_CEHVI-5count.csv"
  )
  
  
  # -----------------------------------------------------
  # 读取数据
  # -----------------------------------------------------
  
  df <- map2(
    files,
    method_levels,
    ~read_one(
      .x,
      .y,
      folder_path
    )
  ) %>%
    
    compact() %>%
    
    bind_rows()
  
  
  # -----------------------------------------------------
  # 转换为宽格式
  # -----------------------------------------------------
  
  df_wide <- df %>%
    
    pivot_wider(
      names_from = method,
      values_from = value,
      values_fn = list
    ) %>%
    
    unnest(
      cols = everything()
    )
  
  
  # -----------------------------------------------------
  # 计算相对于 FEHVI 的迭代次数差值
  #
  # FCEHVI - FEHVI
  # GCEHVI - FEHVI
  # -----------------------------------------------------
  
  diffs <- list(
    
    GEHVI =
      df_wide$GEHVI -
      df_wide$all_EHVI,
    
    DEHVI =
      df_wide$DEHVI -
      df_wide$all_EHVI
  )
  
  
  # -----------------------------------------------------
  # 转换成长格式
  # -----------------------------------------------------
  
  df_diff <- as_tibble(
    diffs
  ) %>%
    
    pivot_longer(
      everything(),
      names_to = "method",
      values_to = "diff"
    ) %>%
    
    mutate(
      method = factor(
        method,
        levels = c(
          "GEHVI",
          "DEHVI"
        )
      )
    )
  
  
  # -----------------------------------------------------
  # 绘图
  # -----------------------------------------------------
  
  p <- ggplot(
    df_diff,
    aes(
      diff,
      color = method,
      fill = method
    )
  ) +
    
    
    geom_density(
      alpha = 0.3,
      linewidth = 1.5
    ) +
    
    
    # 零差值参考线
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 1.2
    ) +
    
    
    # ---------------------------------------------------
  # 颜色
  # ---------------------------------------------------
  
  scale_color_manual(
    values = method_colors[
      c(
        "GEHVI",
        "DEHVI"
      )
    ],
    
    labels = c(
      
      GEHVI =
        expression(
          Delta[FCEHVI]
        ),
      
      DEHVI =
        expression(
          Delta[GCEHVI]
        )
    )
  ) +
    
    
    scale_fill_manual(
      values = method_colors[
        c(
          "GEHVI",
          "DEHVI"
        )
      ],
      
      labels = c(
        
        GEHVI =
          expression(
            Delta[FCEHVI]
          ),
        
        DEHVI =
          expression(
            Delta[GCEHVI]
          )
      )
    ) +
    
    
    # ---------------------------------------------------
  # 主题
  # ---------------------------------------------------
  
  theme_minimal(
    base_size = 45
  ) +
    
    
    theme(
      
      legend.position =
        c(
          0.7,
          0.99
        ),
      
      legend.justification =
        c(
          0,
          1
        ),
      
      legend.background =
        element_rect(
          fill = "white",
          color = NA
        )
    ) +
    
    
    # ---------------------------------------------------
  # 坐标轴
  # ---------------------------------------------------
  
  labs(
    x = "Iteration count difference",
    y = "Density",
    color = NULL,
    fill = NULL
  )
  
  
  return(p)
}


# =======================================================
# 8. 图 C：OC
# =======================================================

generate_OC <- function(
    folder_path
){
  
  # -----------------------------------------------------
  # 读取 OC 数据
  # -----------------------------------------------------
  
  OC <- read.csv(
    file.path(
      folder_path,
      "5ave.dis.OC.csv"
    )
  )
  
  
  # 第一列设置为 X
  colnames(OC)[1] <- "X"
  
  
  # -----------------------------------------------------
  # 转换成长格式
  # -----------------------------------------------------
  
  df_OC <- map_df(
    
    method_levels,
    
    \(m)
    
    tibble(
      
      X = OC$X,
      
      method = m,
      
      mean = OC[[m]],
      
      sd =
        OC[[paste0(
          m,
          "_sd"
        )]]
    )
  )
  
  
  # -----------------------------------------------------
  # 绘图
  # -----------------------------------------------------
  
  p_OC <- ggplot(
    df_OC,
    aes(
      X,
      mean,
      color = method,
      fill = method
    )
  ) +
    
    
    # ---------------------------------------------------
  # ± SD 阴影
  # ---------------------------------------------------
  
  geom_ribbon(
    aes(
      ymin = mean - sd,
      ymax = mean + sd
    ),
    
    alpha = 0.15,
    
    color = NA
  ) +
    
    
    # ---------------------------------------------------
  # 平均曲线
  # ---------------------------------------------------
  
  geom_line(
    linewidth = 1.5
  ) +
    
    
    # ---------------------------------------------------
  # 颜色
  # ---------------------------------------------------
  
  scale_color_manual(
    values = method_colors,
    labels = method_labels
  ) +
    
    
    scale_fill_manual(
      values = method_colors,
      labels = method_labels
    ) +
    
    
    # ---------------------------------------------------
  # 主题
  # ---------------------------------------------------
  
  theme_minimal(
    base_size = 45
  ) +
    
    
    theme(
      
      legend.position =
        c(
          0.995,
          0.995
        ),
      
      legend.justification =
        c(
          1,
          1
        ),
      
      legend.background =
        element_rect(
          fill = "white",
          color = NA
        )
    ) +
    
    
    # ---------------------------------------------------
  # 坐标轴
  # ---------------------------------------------------
  
  labs(
    x = "Iteration",
    y = "Distance",
    color = NULL,
    fill = NULL
  )
  
  
  return(p_OC)
}


# =======================================================
# 9. 给单个图添加 a1 / b1 / c1 标签
# =======================================================

add_tag <- function(
    p,
    tag
){
  
  ggdraw() +
    
    draw_plot(
      p,
      x = 0.04,
      y = 0,
      width = 0.96,
      height = 0.95
    ) +
    
    draw_label(
      tag,
      x = 0.005,
      y = 0.96,
      hjust = -0.5,
      vjust = 1,
      fontface = "bold",
      size = 44
    )
}


# =======================================================
# 10. 生成 HL 三个图
# =======================================================

p_a1 <- generate_box(
  folder_HL
)


p_b1 <- generate_pdf_diff(
  folder_HL
)


p_c1 <- generate_OC(
  folder_HL
)


# =======================================================
# 11. 添加图编号
# =======================================================

p_a1 <- add_tag(
  p_a1,
  "a2"
)


p_b1 <- add_tag(
  p_b1,
  "b2"
)


p_c1 <- add_tag(
  p_c1,
  "c2"
)


# =======================================================
# 12. 横向排列三个图
# =======================================================

p_all <- (
  p_a1 +
    p_b1 +
    p_c1
) +
  
  plot_layout(
    ncol = 3
  )


# =======================================================
# 13. 调整图与图之间的间距
# =======================================================

p_all <- p_all &
  
  theme(
    
    plot.margin =
      margin(
        15,
        25,
        15,
        25
      )
  )


# =======================================================
# 14. 显示最终大图
# =======================================================

p_all


# =======================================================
# 15. 保存最终图片
# =======================================================

ggsave(
  
  "C:\\Users\\12758\\Desktop\\GCEHVI\\hl\\repeat\\HL_repeat0.25.png",
  
  p_all,
  
  width = 36,
  
  height = 10,
  
  dpi = 300
)