# =======================================================
# 完整大图：左图散点图 + 右图多目标优化结果
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
# 文件路径
# =======================================================
data_folder <- "D:\\work\\CMOBO\\manuscript-ty1\\Figure\\code\\Mg\\repeat"

# =======================================================
# 方法 & 颜色（内部名不变）
# =======================================================
method_levels <- c("all_EHVI","GEHVI","DEHVI")
method_colors <- setNames(brewer.pal(3,"Set1"), method_levels)

# =======================================================
# 显示名称（只影响图中文字）
# =======================================================
method_labels <- c(
  all_EHVI = "FEHVI",
  GEHVI    = "FCEHVI",
  DEHVI    = "GCEHVI"
)

# =======================================================
# 左图：散点图
# =======================================================
data <- read.csv(file.path(data_folder,"df_last.csv"))
data <- data[,c("EL","UTS")]
colnames(data) <- c("Elongation","UTS")

tar_point00 <- c(26.62, 309.73)
special_point <- data.frame(
  Elongation = tar_point00[1],
  UTS = tar_point00[2]
)

ref_es1 <- (max(data$Elongation) - min(data$Elongation)) * 0.5 + min(data$Elongation)
ref_es2 <- (max(data$UTS) - min(data$UTS)) * 0.5 + min(data$UTS)

remain_num <- which(
  data$Elongation < ref_es1 &
    data$UTS < ref_es2
)

# =======================================================
# 工具函数
# =======================================================
read_one <- function(filename, name, folder_path){
  full_path <- file.path(folder_path, filename)
  if(!file.exists(full_path)){
    warning("⚠️ 文件不存在: ", full_path)
    return(NULL)
  }
  dat <- read.csv(full_path, row.names=1)
  tibble(method=name, value=as.numeric(unlist(dat)))
}

# =======================================================
# 图 c：PDF 差值（去掉图例标题 Method）
# =======================================================
generate_pdf_diff <- function(folder_path){
  files <- c(
    "all_EHVI-5count.csv",
    "given_constrain_CEHVI-5count.csv",
    "3points_dist_CEHVI-5count.csv"
  )
  
  df <- map2(files, method_levels, ~read_one(.x,.y,folder_path)) %>%
    compact() %>% bind_rows()
  
  if(nrow(df)==0) return(NULL)
  
  df_wide <- df %>%
    pivot_wider(names_from=method, values_from=value, values_fn=list) %>%
    unnest(cols=everything())
  
  diffs <- list()
  if(all(c("GEHVI","all_EHVI") %in% colnames(df_wide)))
    diffs$GEHVI <- df_wide$GEHVI - df_wide$all_EHVI
  if(all(c("DEHVI","all_EHVI") %in% colnames(df_wide)))
    diffs$DEHVI <- df_wide$DEHVI - df_wide$all_EHVI
  
  df_diff <- as_tibble(diffs) %>%
    pivot_longer(everything(), names_to="method", values_to="diff") %>%
    mutate(method=factor(method, levels=names(diffs)))
  
  ggplot(df_diff, aes(diff, color=method, fill=method)) +
    geom_density(alpha=0.3, linewidth=1.5) +
    geom_vline(xintercept=0, linetype="dashed", linewidth=1.2) +
    scale_color_manual(
      values = method_colors[names(diffs)],
      labels = c(
        GEHVI = expression(Delta[FCEHVI]),
        DEHVI = expression(Delta[GCEHVI])
      )
    ) +
    scale_fill_manual(
      values = method_colors[names(diffs)],
      labels = c(
        GEHVI = expression(Delta[FCEHVI]),
        DEHVI = expression(Delta[GCEHVI])
      )
    ) +
    theme_minimal(base_size=45) +
    theme(
      legend.position=c(0.005,0.995),
      legend.justification=c(0,1),
      legend.background=element_rect(fill="white", color=NA)
    ) +
    labs(
      x="Iteration count difference",
      y="Density",
      color=NULL,
      fill=NULL
    )
}

# =======================================================
# 图 d / e：OC & HV（去掉图例标题 Method）
# =======================================================
generate_OC_HV <- function(folder_path){
  OC <- read.csv(file.path(folder_path,"5ave.dis.OC.csv"))
  HV <- read.csv(file.path(folder_path,"5ave.delta_hv.csv"))
  
  colnames(OC)[1] <- "X"
  colnames(HV)[1] <- "X"
  
  df_OC <- map_df(method_levels, \(m)
                  tibble(X=OC$X, method=m, mean=OC[[m]], sd=OC[[paste0(m,"_sd")]])
  )
  
  df_HV <- map_df(method_levels, \(m)
                  tibble(X=HV$X, method=m, mean=HV[[m]], sd=HV[[paste0(m,"_sd")]])
  )
  
  p_OC <- ggplot(df_OC, aes(X, mean, color=method, fill=method)) +
    geom_ribbon(aes(ymin=mean-sd, ymax=mean+sd), alpha=0.15, color=NA) +
    geom_line(linewidth=1.5) +
    scale_color_manual(values=method_colors, labels=method_labels) +
    scale_fill_manual(values=method_colors, labels=method_labels) +
    theme_minimal(base_size=45) +
    theme(
      legend.position=c(0.995,0.995),
      legend.justification=c(1,1),
      legend.background=element_rect(fill="white", color=NA)
    ) +
    labs(
      x="Iteration",
      y="OC",
      color=NULL,
      fill=NULL
    )
  
  p_HV <- ggplot(df_HV, aes(X, mean, color=method, fill=method)) +
    geom_ribbon(aes(ymin=mean-sd, ymax=mean+sd), alpha=0.15, color=NA) +
    geom_line(linewidth=1.5) +
    scale_color_manual(values=method_colors, labels=method_labels) +
    scale_fill_manual(values=method_colors, labels=method_labels) +
    theme_minimal(base_size=45) +
    theme(
      legend.position=c(0.995,0.005),
      legend.justification=c(1,0),
      legend.background=element_rect(fill="white", color=NA)
    ) +
    labs(
      x="Iteration",
      y="HV",
      color=NULL,
      fill=NULL
    )
  
  list(p_OC, p_HV)
}

# =======================================================
# 图 a：散点图
# =======================================================
p_base <- ggplot(data, aes(x=Elongation, y=UTS)) +
  geom_point(color=method_colors["GEHVI"], size=8) +
  geom_point(
    data=data[remain_num,],
    aes(x=Elongation,y=UTS),
    color=method_colors["DEHVI"],
    size=8
  ) +
  geom_point(
    data=special_point,
    aes(x=Elongation,y=UTS),
    color=method_colors["all_EHVI"],
    shape=17,
    size=13
  ) +
  labs(x="Elongation (%)", y="UTS (MPa)") +
  theme_minimal(base_size=45)

p_left_marginal <- ggMarginal(
  p_base,
  type="density",
  fill=method_colors["GEHVI"],
  alpha=0.3,
  color=method_colors["GEHVI"]
)

p_left_grob <- ggdraw() +
  draw_plot(p_left_marginal, x=0, y=0.05, width=1, height=0.9)

# =======================================================
# 图 b：箱线图（保持不变）
# =======================================================
generate_box <- function(folder_path){
  files <- c(
    "all_EHVI-5count.csv",
    "given_constrain_CEHVI-5count.csv",
    "3points_dist_CEHVI-5count.csv"
  )
  
  df <- map2(files, method_levels, ~read_one(.x,.y,folder_path)) %>%
    compact() %>% bind_rows()
  
  df <- df %>% mutate(method=factor(method, levels=method_levels))
  
  stat_df <- df %>%
    group_by(method) %>%
    summarise(
      mean=mean(value),
      sd=sd(value),
      se=sd/sqrt(n()),
      .groups="drop"
    ) %>%
    mutate(
      whisker_top=mean+sd,
      whisker_bottom=mean-sd,
      rect_top=mean+2*se,
      rect_bottom=mean-2*se
    )
  
  ggplot(stat_df, aes(method, mean, fill=method)) +
    geom_rect(
      aes(
        xmin=as.numeric(method)-0.25,
        xmax=as.numeric(method)+0.25,
        ymin=rect_bottom,
        ymax=rect_top
      ),
      color="black",
      alpha=0.85
    ) +
    geom_segment(aes(x=as.numeric(method), xend=as.numeric(method),
                     y=rect_top, yend=whisker_top), linewidth=1.5) +
    geom_segment(aes(x=as.numeric(method)-0.1, xend=as.numeric(method)+0.1,
                     y=whisker_top, yend=whisker_top), linewidth=1.5) +
    geom_segment(aes(x=as.numeric(method), xend=as.numeric(method),
                     y=whisker_bottom, yend=rect_bottom), linewidth=1.5) +
    geom_segment(aes(x=as.numeric(method)-0.1, xend=as.numeric(method)+0.1,
                     y=whisker_bottom, yend=whisker_bottom), linewidth=1.5) +
    geom_segment(aes(x=as.numeric(method)-0.25, xend=as.numeric(method)+0.25,
                     y=mean, yend=mean), linewidth=1.5) +
    scale_fill_manual(values=method_colors, labels=method_labels) +
    scale_x_discrete(labels=method_labels) +
    theme_classic(base_size=45) +
    theme(axis.text.x=element_text(angle=45, hjust=1, face="bold"),
          legend.position="none") +
    labs(x="Method", y="Iteration count")
}

# =======================================================
# 拼图 & 保存
# =======================================================
p_box <- generate_box(data_folder)
p_pdf <- generate_pdf_diff(data_folder)
p_OC_HV <- generate_OC_HV(data_folder)

right_plot <- wrap_plots(
  list(p_box, p_pdf, p_OC_HV[[1]], p_OC_HV[[2]]),
  ncol=2
)

tag_font <- "Arial"

p_left_tagged <- p_left_grob +
  labs(tag="a") +
  theme(
    plot.tag = element_text(size=60, face="bold", family=tag_font),
    plot.tag.position = c(0,1)
  )

final_plot <- wrap_plots(
  p_left_tagged,
  right_plot,
  ncol=2,
  widths=c(1,1.2)
) +
  plot_annotation(
    tag_levels="a",
    theme = theme(
      plot.tag = element_text(size=40, face="bold", family=tag_font)
    )
  )

print(final_plot)

ggsave(
  "D:\\work\\CMOBO\\manuscript-ty1\\Figure\\code\\Mg\\math_treat\\Mg.png",
 # "D:/work/paper/CMOBO/code/Mg/math_treat/Mg.png",
  final_plot,
  width=36,
  height=20,
  dpi=300
)
