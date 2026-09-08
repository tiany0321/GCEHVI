
# =======================================================
# Ti-6Al-4V
# 三个案例：
# sel1 | sel2 | sel3
#
# 原始物理目标：
#   Elongation：越大越好
#   UTS       ：越大越好
#
# Pareto 判断：
#   内部转换为：
#   -Elongation
#   -UTS
#
# 绘图：
#   恢复为原始正值
#
# Iteration 图例：
#   三个案例分别按照各自的最大 iteration 设置
#   不强制统一到 50
#   使用规则整数间隔，使颜色条刻度更加均匀
#
# =======================================================


# =======================================================
# 1. 加载程序包
# =======================================================

library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(patchwork)


# =======================================================
# 2. 三个案例路径
# =======================================================

data_folders <- c(
  
  "C:\\Users\\12758\\Desktop\\GCEHVI\\TiAlV\\specific\\more_iter\\sel1",
  
  "C:\\Users\\12758\\Desktop\\GCEHVI\\TiAlV\\specific\\more_iter\\sel2",
  
  "C:\\Users\\12758\\Desktop\\GCEHVI\\TiAlV\\specific\\more_iter\\sel3"
  
)


# =======================================================
# 3. 文件名称
# =======================================================

database_file_name <-
  "41467_2025_56267_MOESM3_ESM.csv"

train_file_name <-
  "data.train.csv"

new_file_name <-
  "data.new.csv"


# =======================================================
# 4. 目标点
# =======================================================

tar_points <- list(
  
  sel1 = c(
    -12.90,
    -1221
  ),
  
  sel2 = c(
    -12.90,
    -1221
  ),
  
  sel3 = c(
    -12.90,
    -1221
  )
  
)


# =======================================================
# 5. 原始数据库目标
# =======================================================

objective1 <- "Elongation"

objective2 <- "UTS"


# =======================================================
# 6. Pareto 判断
# =======================================================

is_pareto_min <- function(
    x,
    y
){
  
  n <- length(x)
  
  pareto <- rep(
    TRUE,
    n
  )
  
  for(
    i in seq_len(n)
  ){
    
    dominated <- (
      
      x <= x[i] &
        y <= y[i] &
        (
          x < x[i] |
            y < y[i]
        )
      
    )
    
    dominated[i] <- FALSE
    
    if(
      any(dominated)
    ){
      
      pareto[i] <- FALSE
      
    }
    
  }
  
  return(
    pareto
  )
  
}


# =======================================================
# 7. 自动识别目标列
# =======================================================

detect_objective_columns <- function(
    dat
){
  
  if(
    all(
      c(
        "es1",
        "es2"
      ) %in% colnames(dat)
    )
  ){
    
    return(
      c(
        "es1",
        "es2"
      )
    )
    
  }
  
  
  if(
    all(
      c(
        "Elongation",
        "UTS"
      ) %in% colnames(dat)
    )
  ){
    
    return(
      c(
        "Elongation",
        "UTS"
      )
    )
    
  }
  
  
  stop(
    paste0(
      
      "\n无法识别目标列。\n\n",
      
      "程序寻找：\n",
      
      "es1 + es2\n",
      
      "或者\n",
      
      "Elongation + UTS\n\n",
      
      "当前列名：\n",
      
      paste(
        colnames(dat),
        collapse = ", "
      )
      
    )
    
  )
  
}


# =======================================================
# 8. 读取一个案例
# =======================================================

read_case <- function(
    data_folder,
    tar_point00
){
  
  # =====================================================
  # 文件
  # =====================================================
  
  database_file <- file.path(
    data_folder,
    database_file_name
  )
  
  train_file <- file.path(
    data_folder,
    train_file_name
  )
  
  new_file <- file.path(
    data_folder,
    new_file_name
  )
  
  
  # =====================================================
  # 文件检查
  # =====================================================
  
  if(
    !file.exists(database_file)
  ){
    
    stop(
      paste0(
        "\n数据库不存在：\n",
        database_file
      )
    )
    
  }
  
  
  if(
    !file.exists(train_file)
  ){
    
    stop(
      paste0(
        "\n训练集不存在：\n",
        train_file
      )
    )
    
  }
  
  
  if(
    !file.exists(new_file)
  ){
    
    stop(
      paste0(
        "\n新数据不存在：\n",
        new_file
      )
    )
    
  }
  
  
  # =====================================================
  # 读取原始数据库
  # =====================================================
  
  data <- read.csv(
    database_file,
    check.names = FALSE
  )
  
  
  data <- data[
    ,
    -c(1:2),
    drop = FALSE
  ]
  
  
  # =====================================================
  # 检查原始目标
  # =====================================================
  
  if(
    !all(
      c(
        objective1,
        objective2
      ) %in% colnames(data)
    )
  ){
    
    stop(
      paste0(
        
        "\n原始数据库没有 Elongation / UTS。\n\n",
        
        "当前列名：\n",
        
        paste(
          colnames(data),
          collapse = ", "
        )
        
      )
    )
    
  }
  
  
  # =====================================================
  # 提取原始物理数据
  # =====================================================
  
  database_physical <- data[
    ,
    c(
      objective1,
      objective2
    ),
    drop = FALSE
  ]
  
  
  database_physical$Elongation <-
    as.numeric(
      database_physical$Elongation
    )
  
  
  database_physical$UTS <-
    as.numeric(
      database_physical$UTS
    )
  
  
  # =====================================================
  # 数据库范围
  # =====================================================
  
  elong_min <- min(
    database_physical$Elongation,
    na.rm = TRUE
  )
  
  elong_max <- max(
    database_physical$Elongation,
    na.rm = TRUE
  )
  
  
  uts_min <- min(
    database_physical$UTS,
    na.rm = TRUE
  )
  
  uts_max <- max(
    database_physical$UTS,
    na.rm = TRUE
  )
  
  
  # =====================================================
  # 读取训练集
  # =====================================================
  
  data.train <- read.csv(
    train_file,
    check.names = FALSE
  )
  
  
  train_columns <-
    detect_objective_columns(
      data.train
    )
  
  
  print(
    paste0(
      "\n训练集目标列：",
      train_columns[1],
      " / ",
      train_columns[2]
    )
  )
  
  
  # =====================================================
  # 提取训练集
  # =====================================================
  
  train_raw <- data.train[
    ,
    train_columns,
    drop = FALSE
  ]
  
  
  train_raw[[train_columns[1]]] <-
    as.numeric(
      train_raw[[train_columns[1]]]
    )
  
  
  train_raw[[train_columns[2]]] <-
    as.numeric(
      train_raw[[train_columns[2]]]
    )
  
  
  # =====================================================
  # 转换为物理正值
  # =====================================================
  
  if(
    identical(
      train_columns,
      c(
        "es1",
        "es2"
      )
    )
  ){
    
    train_physical <- data.frame(
      
      Elongation =
        -train_raw$es1,
      
      UTS =
        -train_raw$es2
      
    )
    
  } else {
    
    train_physical <- data.frame(
      
      Elongation =
        train_raw$Elongation,
      
      UTS =
        train_raw$UTS
      
    )
    
  }
  
  
  # =====================================================
  # 训练集归一化
  # =====================================================
  
  train_plot <- data.frame(
    
    Elongation = (
      
      train_physical$Elongation -
        elong_min
      
    ) /
      
      (
        elong_max -
          elong_min
      ),
    
    
    UTS = (
      
      train_physical$UTS -
        uts_min
      
    ) /
      
      (
        uts_max -
          uts_min
      )
    
  )
  
  
  # =====================================================
  # 保存物理坐标
  # =====================================================
  
  train_plot$Elongation_physical <-
    train_physical$Elongation
  
  
  train_plot$UTS_physical <-
    train_physical$UTS
  
  
  train_plot$iter <- 1
  
  train_plot$dataset <- "Training"
  
  
  # =====================================================
  # 读取新数据
  # =====================================================
  
  data.new <- read.csv(
    new_file,
    check.names = FALSE
  )
  
  
  new_columns <-
    detect_objective_columns(
      data.new
    )
  
  
  print(
    paste0(
      "新数据目标列：",
      new_columns[1],
      " / ",
      new_columns[2]
    )
  )
  
  
  # =====================================================
  # 提取新数据
  # =====================================================
  
  new_raw <- data.new[
    ,
    new_columns,
    drop = FALSE
  ]
  
  
  new_raw[[new_columns[1]]] <-
    as.numeric(
      new_raw[[new_columns[1]]]
    )
  
  
  new_raw[[new_columns[2]]] <-
    as.numeric(
      new_raw[[new_columns[2]]]
    )
  
  
  # =====================================================
  # 转换为原始物理正值
  # =====================================================
  
  if(
    identical(
      new_columns,
      c(
        "es1",
        "es2"
      )
    )
  ){
    
    new_physical <- data.frame(
      
      Elongation =
        -new_raw$es1,
      
      UTS =
        -new_raw$es2
      
    )
    
  } else {
    
    new_physical <- data.frame(
      
      Elongation =
        new_raw$Elongation,
      
      UTS =
        new_raw$UTS
      
    )
    
  }
  
  
  # =====================================================
  # 新数据
  # =====================================================
  
  new_plot <- data.frame(
    
    Elongation =
      new_physical$Elongation,
    
    UTS =
      new_physical$UTS
    
  )
  
  
  # =====================================================
  # iteration
  # =====================================================
  
  iteration_candidates <- c(
    
    "iteration",
    
    "Iteration",
    
    "iter",
    
    "Iter",
    
    "iteration_number",
    
    "Iteration_number"
    
  )
  
  
  iteration_col <-
    iteration_candidates[
      iteration_candidates %in%
        colnames(data.new)
    ]
  
  
  if(
    length(iteration_col) > 0
  ){
    
    iteration_col <- iteration_col[1]
    
    
    print(
      paste0(
        "Iteration 列：",
        iteration_col
      )
    )
    
    
    new_plot$iter <- as.numeric(
      data.new[[iteration_col]]
    )
    
  } else {
    
    print(
      "未检测到 iteration 列，按照数据顺序编号。"
    )
    
    
    new_plot$iter <- seq_len(
      nrow(new_plot)
    )
    
  }
  
  
  new_plot$dataset <- "New"
  
  
  # =====================================================
  # 目标点
  # =====================================================
  
  target_physical <- data.frame(
    
    Elongation =
      -tar_point00[1],
    
    UTS =
      -tar_point00[2]
    
  )
  
  
  # =====================================================
  # Pareto
  # =====================================================
  
  train_plot$Pareto <-
    is_pareto_min(
      
      -train_plot$Elongation_physical,
      
      -train_plot$UTS_physical
      
    )
  
  
  # =====================================================
  # Pareto 前沿
  # =====================================================
  
  train_pareto <- train_plot[
    train_plot$Pareto,
    ,
    drop = FALSE
  ]
  
  
  # =====================================================
  # 按 Elongation 排序
  # =====================================================
  
  train_pareto <- train_pareto[
    order(
      train_pareto$Elongation_physical
    ),
    ,
    drop = FALSE
  ]
  
  
  # =====================================================
  # 输出 Pareto
  # =====================================================
  
  cat(
    "\n==========================================\n"
  )
  
  cat(
    "案例：",
    data_folder,
    "\n"
  )
  
  cat(
    "Pareto 前沿点数量：",
    nrow(train_pareto),
    "\n"
  )
  
  cat(
    "最大 iteration：",
    max(
      new_plot$iter,
      na.rm = TRUE
    ),
    "\n"
  )
  
  cat(
    "==========================================\n"
  )
  
  
  return(
    
    list(
      
      train = train_plot,
      
      new = new_plot,
      
      pareto = train_pareto,
      
      target = target_physical
      
    )
    
  )
  
}


# =======================================================
# 9. 读取三个案例
# =======================================================

case1 <- read_case(
  data_folders[1],
  tar_points$sel1
)


case2 <- read_case(
  data_folders[2],
  tar_points$sel2
)


case3 <- read_case(
  data_folders[3],
  tar_points$sel3
)


# =======================================================
# 10. 三个图统一物理坐标范围
# =======================================================

all_x <- c(
  
  case1$train$Elongation_physical,
  case1$new$Elongation,
  case1$target$Elongation,
  
  case2$train$Elongation_physical,
  case2$new$Elongation,
  case2$target$Elongation,
  
  case3$train$Elongation_physical,
  case3$new$Elongation,
  case3$target$Elongation
  
)


all_y <- c(
  
  case1$train$UTS_physical,
  case1$new$UTS,
  case1$target$UTS,
  
  case2$train$UTS_physical,
  case2$new$UTS,
  case2$target$UTS,
  
  case3$train$UTS_physical,
  case3$new$UTS,
  case3$target$UTS
  
)


# =======================================================
# 11. 坐标范围
# =======================================================

x_range <- range(
  all_x,
  na.rm = TRUE
)

y_range <- range(
  all_y,
  na.rm = TRUE
)


x_margin <- diff(
  x_range
) * 0.06


y_margin <- diff(
  y_range
) * 0.06


common_xlim <- c(
  
  x_range[1] -
    x_margin,
  
  x_range[2] +
    x_margin
  
)


common_ylim <- c(
  
  y_range[1] -
    y_margin,
  
  y_range[2] +
    y_margin
  
)


# =======================================================
# 12. 为每个案例生成规则、均匀的 iteration 刻度
# =======================================================
#
# 核心原则：
#
#   1. 每个案例使用自己的最大 iteration
#   2. 刻度从 1 开始
#   3. 刻度尽可能等间距
#   4. 刻度为整数
#   5. 强制包含最大 iteration
#
# 例如：
#
#   max_iter = 16
#
#   得到：
#   1, 4, 7, 10, 13, 16
#
#   max_iter = 50
#
#   得到：
#   1, 11, 21, 31, 41, 50
#
# =======================================================

get_iteration_breaks <- function(
    max_iter,
    n_breaks = 6
){
  
  max_iter <- as.integer(
    max_iter
  )
  
  
  # -----------------------------------------------
  # 检查
  # -----------------------------------------------
  
  if(
    is.na(max_iter) ||
    max_iter < 1
  ){
    
    stop(
      "max_iter 必须是大于等于 1 的整数。"
    )
    
  }
  
  
  # -----------------------------------------------
  # iteration 较少时全部显示
  # -----------------------------------------------
  
  if(
    max_iter <= n_breaks
  ){
    
    return(
      seq(
        1,
        max_iter,
        by = 1
      )
    )
    
  }
  
  
  # -----------------------------------------------
  # 理想步长
  # -----------------------------------------------
  
  raw_step <- (
    max_iter - 1
  ) / (
    n_breaks - 1
  )
  
  
  # -----------------------------------------------
  # 选择最接近的整数步长
  # -----------------------------------------------
  
  step <- max(
    1,
    round(
      raw_step
    )
  )
  
  
  # -----------------------------------------------
  # 按整数步长生成
  # -----------------------------------------------
  
  breaks <- seq(
    1,
    max_iter,
    by = step
  )
  
  
  # -----------------------------------------------
  # 强制加入最大 iteration
  # -----------------------------------------------
  
  breaks <- unique(
    c(
      breaks,
      max_iter
    )
  )
  
  
  breaks <- sort(
    breaks
  )
  
  
  # -----------------------------------------------
  # 如果因为步长导致刻度太多
  # 则尽可能均匀抽取
  # -----------------------------------------------
  
  if(
    length(breaks) > n_breaks
  ){
    
    index <- round(
      seq(
        1,
        length(breaks),
        length.out = n_breaks
      )
    )
    
    breaks <- breaks[
      unique(index)
    ]
    
  }
  
  
  # -----------------------------------------------
  # 最终确保首尾存在
  # -----------------------------------------------
  
  breaks <- sort(
    unique(
      c(
        1,
        breaks,
        max_iter
      )
    )
  )
  
  
  return(
    as.integer(
      breaks
    )
  )
  
}


# =======================================================
# 13. 三个案例分别计算 iteration 范围
# =======================================================

max_iter1 <- max(
  case1$new$iter,
  na.rm = TRUE
)

max_iter2 <- max(
  case2$new$iter,
  na.rm = TRUE
)

max_iter3 <- max(
  case3$new$iter,
  na.rm = TRUE
)


iteration_breaks1 <- get_iteration_breaks(
  max_iter1,
  n_breaks = 6
)


iteration_breaks2 <- get_iteration_breaks(
  max_iter2,
  n_breaks = 6
)


iteration_breaks3 <- get_iteration_breaks(
  max_iter3,
  n_breaks = 6
)


# =======================================================
# 14. 输出 iteration 信息
# =======================================================

cat(
  "\n==========================================\n"
)

cat(
  "Iteration 图例设置\n"
)

cat(
  "==========================================\n"
)

cat(
  "sel1 最大 iteration：",
  max_iter1,
  "\n"
)

cat(
  "sel1 图例刻度：",
  paste(
    iteration_breaks1,
    collapse = ", "
  ),
  "\n\n"
)

cat(
  "sel2 最大 iteration：",
  max_iter2,
  "\n"
)

cat(
  "sel2 图例刻度：",
  paste(
    iteration_breaks2,
    collapse = ", "
  ),
  "\n\n"
)

cat(
  "sel3 最大 iteration：",
  max_iter3,
  "\n"
)

cat(
  "sel3 图例刻度：",
  paste(
    iteration_breaks3,
    collapse = ", "
  ),
  "\n"
)

cat(
  "==========================================\n\n"
)


# =======================================================
# 15. 单图函数
# =======================================================

generate_plot <- function(
    case_data,
    iteration_breaks,
    max_iteration
){
  
  train_obj <-
    case_data$train
  
  new_obj <-
    case_data$new
  
  train_pareto <-
    case_data$pareto
  
  target_obj <-
    case_data$target
  
  
  ggplot() +
    
    
    # ===================================================
  # 训练集普通点
  # ===================================================
  
  geom_point(
    
    data =
      train_obj[
        !train_obj$Pareto,
        ,
        drop = FALSE
      ],
    
    aes(
      x = Elongation_physical,
      y = UTS_physical
    ),
    
    shape = 1,
    
    size = 14.4,
    
    stroke = 1.2,
    
    color = "grey55"
    
  ) +
    
    
    # ===================================================
  # Pareto 前沿虚线
  # ===================================================
  
  geom_path(
    
    data =
      train_pareto,
    
    aes(
      x = Elongation_physical,
      y = UTS_physical
    ),
    
    color = "grey45",
    
    linewidth = 1.8,
    
    linetype = "dashed",
    
    lineend = "round",
    
    linejoin = "round"
    
  ) +
    
    
    # ===================================================
  # Pareto 实心点
  # ===================================================
  
  geom_point(
    
    data =
      train_pareto,
    
    aes(
      x = Elongation_physical,
      y = UTS_physical
    ),
    
    shape = 16,
    
    size = 16.5,
    
    color = "grey35"
    
  ) +
    
    
    # ===================================================
  # 新数据
  # ===================================================
  
  geom_point(
    
    data =
      new_obj,
    
    aes(
      x = Elongation,
      y = UTS,
      color = iter
    ),
    
    shape = 16,
    
    size = 14.4
    
  ) +
    
    
    # ===================================================
  # 目标点
  # 空心菱形
  # ===================================================
  
  geom_point(
    
    data =
      target_obj,
    
    aes(
      x = Elongation,
      y = UTS
    ),
    
    shape = 5,
    
    size = 15,
    
    stroke = 1.8,
    
    color = "black"
    
  ) +
    
    
    # ===================================================
  # Iteration 色阶
  # ===================================================
  
  scale_color_gradientn(
    
    colours =
      brewer.pal(
        9,
        "YlOrRd"
      ),
    
    name = "Iteration",
    
    breaks = iteration_breaks,
    
    labels = as.character(
      iteration_breaks
    ),
    
    limits = c(
      1,
      max_iteration
    ),
    
    guide = guide_colorbar(
      
      title.position = "top",
      
      title.hjust = 0.5,
      
      label.position = "bottom",
      
      ticks = TRUE,
      
      draw.ulim = TRUE,
      
      draw.llim = TRUE,
      
      # ---------------------------------------------
      # 颜色条长度
      # ---------------------------------------------
      
      barwidth =
        unit(
          7.5,
          "cm"
        ),
      
      barheight =
        unit(
          0.55,
          "cm"
        )
      
    )
    
  ) +
    
    
    # ===================================================
  # 物理坐标
  # ===================================================
  
  coord_cartesian(
    
    xlim = common_xlim,
    
    ylim = common_ylim,
    
    expand = FALSE,
    
    clip = "off"
    
  ) +
    
    
    # ===================================================
  # 坐标标签
  # ===================================================
  
  labs(
    
    x = "Elongation (%)",
    
    y = "UTS (MPa)"
    
  ) +
    
    
    # ===================================================
  # 主题
  # ===================================================
  
  theme_classic(
    
    base_size = 37.5
    
  ) +
    
    
    theme(
      
      # -------------------------------------------------
      # X 标题
      # -------------------------------------------------
      
      axis.title.x = element_text(
        
        size = 42,
        
        family = "Arial",
        
        margin =
          margin(
            t = 18
          )
        
      ),
      
      
      # -------------------------------------------------
      # Y 标题
      # -------------------------------------------------
      
      axis.title.y = element_text(
        
        size = 42,
        
        family = "Arial",
        
        margin =
          margin(
            r = 18
          )
        
      ),
      
      
      # -------------------------------------------------
      # 刻度文字
      # -------------------------------------------------
      
      axis.text.x = element_text(
        
        size = 33,
        
        family = "Arial"
        
      ),
      
      
      axis.text.y = element_text(
        
        size = 33,
        
        family = "Arial"
        
      ),
      
      
      # -------------------------------------------------
      # 坐标轴
      # -------------------------------------------------
      
      axis.line =
        element_blank(),
      
      panel.border =
        element_rect(
          
          color = "black",
          
          fill = NA,
          
          linewidth = 1.1
          
        ),
      
      
      # -------------------------------------------------
      # 内向刻度
      # -------------------------------------------------
      
      axis.ticks =
        element_line(
          
          color = "black",
          
          linewidth = 0.75
          
        ),
      
      
      axis.ticks.length =
        unit(
          
          -0.27,
          
          "cm"
          
        ),
      
      
      # -------------------------------------------------
      # 上、右不显示刻度
      # -------------------------------------------------
      
      axis.ticks.x.top =
        element_blank(),
      
      axis.ticks.y.right =
        element_blank(),
      
      
      # -------------------------------------------------
      # 图例
      # -------------------------------------------------
      
      legend.title =
        element_text(
          
          size = 33,
          
          family = "Arial"
          
        ),
      
      
      legend.text =
        element_text(
          
          size = 27,
          
          family = "Arial"
          
        ),
      
      
      legend.key.height =
        unit(
          1.2,
          "cm"
        ),
      
      
      legend.key.width =
        unit(
          1.2,
          "cm"
        ),
      
      
      legend.background =
        element_rect(
          
          fill = "white",
          
          color = NA
          
        ),
      
      
      legend.position =
        c(
          
          0.99,
          
          0.99
          
        ),
      
      
      legend.justification =
        c(
          
          1,
          
          1
          
        ),
      
      
      legend.direction =
        "horizontal",
      
      
      # -------------------------------------------------
      # 图边距
      # -------------------------------------------------
      
      plot.margin =
        margin(
          
          t = 10,
          
          r = 20,
          
          b = 20,
          
          l = 20
          
        )
      
    )
  
}


# =======================================================
# 16. 生成三个图
# =======================================================

p1 <- generate_plot(
  
  case1,
  
  iteration_breaks1,
  
  max_iter1
  
)


p2 <- generate_plot(
  
  case2,
  
  iteration_breaks2,
  
  max_iter2
  
)


p3 <- generate_plot(
  
  case3,
  
  iteration_breaks3,
  
  max_iter3
  
)


# =======================================================
# 17. 横向排列
# =======================================================

p_all <- (
  
  p1 |
    
    p2 |
    
    p3
  
) +
  
  plot_layout(
    
    nrow = 1,
    
    widths = c(
      1,
      1,
      1
    )
    
  )


# =======================================================
# 18. 三图间距
# =======================================================

p_all <- p_all &
  
  theme(
    
    plot.margin =
      margin(
        
        t = 10,
        
        r = 35,
        
        b = 20,
        
        l = 35
        
      )
    
  )


# =======================================================
# 19. 显示
# =======================================================

p_all


# =======================================================
# 20. 保存
# =======================================================

ggsave(
  
  "C:\\Users\\12758\\Desktop\\GCEHVI\\TiAlV\\specific\\more_iter\\TiAlV_specific.png",
  
  p_all,
  
  width = 36,
  
  height = 10,
  
  dpi = 300
  
)

