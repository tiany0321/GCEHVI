library(ggplot2)
library(viridis)
library(patchwork)

# =========================
# 自定义绘图函数（增加 legend_pos 参数）
# =========================
plot_density <- function(data_all, data_train, data_new, tar,
                         color_limits = NULL, label_text = NULL,
                         legend_pos = c(0.85, 0.77)){
  ggplot() +
    ## 背景密度
    stat_density_2d(
      data = data_all,
      aes(x = es1, y = es2, fill = after_stat(density)),
      geom = "raster",
      contour = FALSE,
      show.legend = FALSE
    ) +
    scale_fill_viridis_c(
      option = "viridis",
      direction = 1,
      begin = 0.55,
      end = 0.95,
      alpha = 0.45
    ) +
    ## 轮廓线
    stat_density_2d(
      data = data_all,
      aes(x = es1, y = es2),
      color = "#66C2A5",
      bins = 8,
      linewidth = 0.5,
      alpha = 0.7
    ) +
    ## 训练点
    geom_point(
      data = data_train,
      aes(x = es1, y = es2),
      color = "grey60",
      size = 3.5,
      alpha = 0.8
    ) +
    ## 新采样点
    geom_point(
      data = data_new,
      aes(x = es1, y = es2, color = X),
      size = 5.5
    ) +
    scale_color_viridis_c(
      name = "Iteration",
      option = "viridis",
      limits = color_limits
    ) +
    ## 目标点
    geom_point(
      data = tar,
      aes(x = es1, y = es2),
      shape  = 8,
      size   = 8,
      stroke = 1.2,
      color  = "red"
    ) +
    ## 图内标注
    annotate(
      "text",
      x = -Inf, y = Inf,
      label = label_text,
      hjust = -0.02,
      vjust = 1.05,
      size = 10,
      fontface = "bold",
      family = "Arial"
    ) +
    labs(
      x = expression(y[1]),
      y = expression(y[2])
    ) +
    theme_minimal(base_family = "Arial") +
    theme(
      legend.position = legend_pos,  # ← 根据参数设置位置
      legend.background = element_rect(fill = alpha("white", 0.5), color = NA),
      legend.title = element_text(size = 18, family = "Arial"),
      legend.text  = element_text(size = 16, family = "Arial"),
      panel.grid = element_blank(),
      axis.title = element_text(size = 26, family = "Arial"),
      axis.text  = element_text(size = 22, family = "Arial")
    )
}

# =========================
# 文件夹和方法标签
# =========================
problems <- c("P1/specific","poloni/specific",  "Kursawe/specific")
methods  <- c("sel1", "sel2", "sel3")
labels   <- paste0(
  rep(c("a","b","c"), times = 3),
  rep(1:3, each = 3),
  ")"
)
plots <- list()

# =========================
# 按行生成子图并统一每行颜色坐标
# =========================
counter <- 1
for(p_idx in 1:length(problems)){
  p <- problems[p_idx]
  
  # 获取该问题下所有 methods 的 X 范围
  all_X <- c()
  for(m in methods){
    path <- paste0("D:/work/paper/CMOBO/code/", p, "/", m)
    data_new <- read.csv(file.path(path, "data.new.csv"))
    data_new$X <- 1:nrow(data_new)
    all_X <- c(all_X, data_new$X)
  }
  color_limits <- range(all_X)
  
  # 根据行数设置图例位置：第三行放左上角
  legend_pos <- if(p_idx == 3) c(0.15, 0.77) else c(0.85, 0.77)
  
  # 生成该行的三列图
  for(m in methods){
    path <- paste0("D:/work/paper/CMOBO/code/", p, "/", m)
    
    data_all   <- read.csv(file.path(path, "data.all.csv"))
    data_train <- read.csv(file.path(path, "data.train.csv"))
    data_new   <- read.csv(file.path(path, "data.new.csv"))
    tar        <- read.csv(file.path(path, "tar_point00.csv"))
    data_new$X <- 1:nrow(data_new)
    
    plots[[counter]] <- plot_density(
      data_all, data_train, data_new, tar,
      color_limits = color_limits,
      label_text   = labels[counter],
      legend_pos   = legend_pos
    )
    
    counter <- counter + 1
  }
}

# =========================
# 3×3 网格排列
# =========================
final_plot <- (plots[[1]] | plots[[2]] | plots[[3]]) /
  (plots[[4]] | plots[[5]] | plots[[6]]) /
  (plots[[7]] | plots[[8]] | plots[[9]])

print(final_plot)

# =========================
# 写出图片
# =========================
ggsave(
  filename = "math_specific.png",
  plot     = final_plot,
  width    = 50,
  height   = 42,
  units    = "cm",
  dpi      = 600,
  bg       = "white"
)
