library(ggplot2)
library(viridis)
library(patchwork)

# =========================
# 自定义绘图函数
# =========================
plot_density <- function(data_all, data_train, data_new, tar,
                         color_limits = NULL, label_text = NULL){
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
    
    ## 轮廓线（更柔和）
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
    
    ## 图内标注 a) b) c)
    annotate(
      "text",
      x = -Inf, y = Inf,
      label = label_text,
      hjust = -0.02,
      vjust = 1.05,
      size = 10,                 # 字体放大一倍
      fontface = "bold",
      family = "Arial"
    ) +
    
    ## 坐标轴
    labs(
      x = expression(y[1]),
      y = expression(y[2])
    ) +
    
    theme_minimal(base_family = "Arial") +
    theme(
      legend.position = c(0.85, 0.8),
      legend.background = element_rect(
        fill = alpha("white", 0.5), color = NA
      ),
      legend.title = element_text(size = 18, family = "Arial"),
      legend.text  = element_text(size = 16, family = "Arial"),
      panel.grid = element_blank(),
      axis.title = element_text(size = 26, family = "Arial"),
      axis.text  = element_text(size = 22, family = "Arial")
    )
}

# =========================
# 文件夹与标签
# =========================
folders <- c("sel1", "sel2", "sel3")
labels  <- c("a)", "b)", "c)")
plots   <- list()

# =========================
# 统一 Iteration 颜色坐标
# =========================
all_X <- c()
for(i in 1:3){
  path <- paste0("D:/work/paper/CMOBO/code/poloni/specific/", folders[i])
  data_new <- read.csv(file.path(path, "data.new.csv"))
  data_new$X <- 1:nrow(data_new)
  all_X <- c(all_X, data_new$X)
}
color_limits <- range(all_X)

# =========================
# 生成三个子图
# =========================
for(i in 1:3){
  path <- paste0("D:/work/paper/CMOBO/code/poloni/specific/", folders[i])
  
  data_all   <- read.csv(file.path(path, "data.all.csv"))
  data_train <- read.csv(file.path(path, "data.train.csv"))
  data_new   <- read.csv(file.path(path, "data.new.csv"))
  tar        <- read.csv(file.path(path, "tar_point00.csv"))
  data_new$X <- 1:nrow(data_new)
  
  plots[[i]] <- plot_density(
    data_all, data_train, data_new, tar,
    color_limits = color_limits,
    label_text   = labels[i]
  )
}

# =========================
# 一行三列
# =========================
final_plot <- plots[[1]] + plots[[2]] + plots[[3]] +
  plot_layout(nrow = 1)

print(final_plot)

# =========================
# 写出图片
# =========================
ggsave(
  filename = "density_sel1_sel2_sel3.png",
  plot     = final_plot,
  width    = 50,
  height   = 15,
  units    = "cm",
  dpi      = 600,
  bg       = "white"
)
