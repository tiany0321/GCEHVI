# =======================================================
# =======================================================
#              Mg alloy - LOOCV prediction
#              EL and UTS parity plots
#
#              No R2 / RMSE / MAE calculation
# =======================================================
# =======================================================


# =======================================================
# 1. Working directory and packages
# =======================================================

setwd(
  "C:\\Users\\12758\\Desktop\\GCEHVI\\Mg\\performance"
)


RNGkind("L'Ecuyer-CMRG")

set.seed(8)


library(DiceKriging)
library(ggplot2)
library(patchwork)


t1 <- proc.time()


# =======================================================
# 2. Parameters
# =======================================================

iter_num <- 80

repeat_num <- 100

threshold <- 0


# =======================================================
# 3. Read data
# =======================================================

data <- read.csv(
  "df_last.csv"
)


# -------------------------------------------------------
# Remove the 16th original column
# -------------------------------------------------------

data <- data[, -c(16)]


# -------------------------------------------------------
# Convert EL and UTS to minimization objectives
#
# Original:
#   EL  -> maximize
#   UTS -> maximize
#
# Optimization:
#   -EL
#   -UTS
# -------------------------------------------------------

data[, c(15, 16)] <- -data[, c(15, 16)]


# =======================================================
# 4. Save original objective values
# =======================================================

df <- data[, c("EL", "UTS")]


# =======================================================
# 5. Normalize objectives to [0,1]
# =======================================================

fn_scale <- function(
    scale_data,
    df
){
  
  
  es_norm <- (
    
    scale_data -
      min(df)
    
  ) / (
    
    max(df) -
      min(df)
    
  )
  
  
  return(
    es_norm
  )
  
}


data$EL <- fn_scale(
  
  data[, "EL"],
  
  df[, "EL"]
  
)


data$UTS <- fn_scale(
  
  data[, "UTS"],
  
  df[, "UTS"]
  
)


# =======================================================
# 6. Target point
# =======================================================

tar_point0 <- c(
  
  -26.62,
  -309.73
  
)


tar_point00 <- c(
  
  fn_scale(
    tar_point0[1],
    df[, "EL"]
  ),
  
  fn_scale(
    tar_point0[2],
    df[, "UTS"]
  )
  
)


cat(
  "\n===============================================\n"
)


cat(
  "Target point\n"
)


cat(
  "===============================================\n"
)


cat(
  "Original target EL  :",
  tar_point0[1],
  "\n"
)


cat(
  "Original target UTS :",
  tar_point0[2],
  "\n"
)


cat(
  "Normalized target es1:",
  tar_point00[1],
  "\n"
)


cat(
  "Normalized target es2:",
  tar_point00[2],
  "\n"
)


cat(
  "===============================================\n"
)


# =======================================================
# 7. Objective names and input variables
# =======================================================

scale <- TRUE


obj <- c(
  
  "es1",
  "es2"
  
)


keep.y <- obj


keep.x <- colnames(data)[
  
  c(
    1:5,
    7:14
  )
  
]


# =======================================================
# 8. Construct design and response matrices
# =======================================================

test.grid <- data[
  
  ,
  
  keep.x,
  
  drop = FALSE
  
]


response.grid <- data[
  
  ,
  
  keep.y,
  
  drop = FALSE
  
]


colnames(test.grid) <- keep.x


colnames(response.grid) <- keep.y


# =======================================================
# 9. Combine design and response data
# =======================================================

data.all <- cbind(
  
  test.grid,
  
  response.grid
  
)


order.data <- data.frame(
  
  order.num = seq_len(
    nrow(data.all)
  )
  
)


data.all <- cbind(
  
  data.all,
  
  order.data
  
)


# =======================================================
# 10. Pareto front
# =======================================================

nd_points <- nondominated_points(
  
  t(response.grid)
  
)


nd_df <- data.frame(
  
  f1 = nd_points[1, ],
  
  f2 = nd_points[2, ]
  
)


# =======================================================
# 11. Objective-space visualization
# =======================================================

plot(
  
  response.grid[, 1],
  
  response.grid[, 2],
  
  xlim = c(
    0,
    1
  ),
  
  ylim = c(
    0,
    1
  ),
  
  col = "grey80",
  
  pch = 16,
  
  xlab = "es1",
  
  ylab = "es2"
  
)


points(
  
  nd_df[, 1],
  
  nd_df[, 2],
  
  col = "red",
  
  pch = 16
  
)


# =======================================================
# 12. Construct training set
#
# IMPORTANT:
#
# No independent test set is used.
#
# Training set:
#
#   es1 > ref_es1
#   es2 > ref_es2
# =======================================================


ref_es1 <- (
  
  max(
    response.grid[, "es1"]
  ) -
    
    min(
      response.grid[, "es1"]
    )
  
) * 0.45


ref_es2 <- (
  
  max(
    response.grid[, "es2"]
  ) -
    
    min(
      response.grid[, "es2"]
    )
  
) * 0.45


remain_num <- which(
  
  response.grid[, "es1"] >
    ref_es1 &
    
    response.grid[, "es2"] >
    ref_es2
  
)


# -------------------------------------------------------
# Only selected samples are used as training samples
# -------------------------------------------------------

train_num <- remain_num


train.x <- test.grid[
  
  train_num,
  
  ,
  
  drop = FALSE
  
]


train.y <- response.grid[
  
  train_num,
  
  ,
  
  drop = FALSE
  
]


# =======================================================
# 13. Dataset information
# =======================================================

cat(
  "\n===============================================\n"
)


cat(
  "Dataset information\n"
)


cat(
  "===============================================\n"
)


cat(
  "Total samples    :",
  nrow(data.all),
  "\n"
)


cat(
  "Training samples :",
  nrow(train.x),
  "\n"
)


cat(
  "Training ratio   :",
  
  round(
    
    nrow(train.x) /
      nrow(data.all) *
      100,
    
    2
    
  ),
  
  "%\n"
)


cat(
  "Validation method: LOOCV\n"
)


cat(
  "===============================================\n"
)


# =======================================================
# 14. Training-set visualization
# =======================================================

plot(
  
  response.grid[, 1],
  
  response.grid[, 2],
  
  col = "grey80",
  
  pch = 16,
  
  xlab = "es1",
  
  ylab = "es2"
  
)


points(
  
  response.grid[train_num, 1],
  
  response.grid[train_num, 2],
  
  col = "#2166AC",
  
  pch = 16
  
)


# =======================================================
# 15. LOOCV
#
# Only the selected training set is used.
#
# No R2 / RMSE / MAE calculation.
# =======================================================

n_train <- nrow(
  train.x
)


loocv_pred_es1 <- rep(
  
  NA_real_,
  
  n_train
  
)


loocv_pred_es2 <- rep(
  
  NA_real_,
  
  n_train
  
)


cat(
  "\n===============================================\n"
)


cat(
  "Starting LOOCV\n"
)


cat(
  "===============================================\n"
)


for(
  
  i in seq_len(n_train)
  
){
  
  
  cat(
    
    "\rLOOCV: ",
    
    i,
    
    " / ",
    
    n_train,
    
    sep = ""
    
  )
  
  
  # -----------------------------------------------------
  # Leave one sample out
  # -----------------------------------------------------
  
  train_index <- setdiff(
    
    seq_len(n_train),
    
    i
    
  )
  
  
  x_cv <- train.x[
    
    train_index,
    
    ,
    
    drop = FALSE
    
  ]
  
  
  y_cv_es1 <- train.y[
    
    train_index,
    
    "es1"
    
  ]
  
  
  y_cv_es2 <- train.y[
    
    train_index,
    
    "es2"
    
  ]
  
  
  x_valid <- train.x[
    
    i,
    
    ,
    
    drop = FALSE
    
  ]
  
  
  # =====================================================
  # EL: es1
  # =====================================================
  
  set.seed(243)
  
  
  mf1_cv <- km(
    
    formula = ~.,
    
    design = x_cv,
    
    response = y_cv_es1,
    
    nugget.estim = TRUE,
    
    covtype = "exp",
    
    control = list(
      
      trace = FALSE
      
    )
    
  )
  
  
  pred_es1 <- predict(
    
    mf1_cv,
    
    newdata = x_valid,
    
    type = "UK"
    
  )
  
  
  loocv_pred_es1[i] <- as.numeric(
    
    pred_es1$mean
    
  )
  
  
  # =====================================================
  # UTS: es2
  # =====================================================
  
  set.seed(243)
  
  
  mf2_cv <- km(
    
    formula = ~.,
    
    design = x_cv,
    
    response = y_cv_es2,
    
    nugget.estim = TRUE,
    
    covtype = "exp",
    
    control = list(
      
      trace = FALSE
      
    )
    
  )
  
  
  pred_es2 <- predict(
    
    mf2_cv,
    
    newdata = x_valid,
    
    type = "UK"
    
  )
  
  
  loocv_pred_es2[i] <- as.numeric(
    
    pred_es2$mean
    
  )
  
  
}


cat(
  "\nLOOCV finished.\n"
)


# =======================================================
# 16. LOOCV results
# =======================================================

loocv_result <- data.frame(
  
  true_es1 =
    train.y[, "es1"],
  
  pred_es1 =
    loocv_pred_es1,
  
  true_es2 =
    train.y[, "es2"],
  
  pred_es2 =
    loocv_pred_es2
  
)


# =======================================================
# 17. Save normalized LOOCV prediction results
# =======================================================

write.csv(
  
  loocv_result,
  
  "LOOCV_prediction_results_normalized.csv",
  
  row.names = FALSE
  
)


# =======================================================
# 18. Final EL model
# =======================================================

cat(
  "\n===============================================\n"
)


cat(
  "Building final EL model\n"
)


cat(
  "===============================================\n"
)


set.seed(243)


mf1_final <- km(
  
  formula = ~.,
  
  design = train.x,
  
  response = train.y[, "es1"],
  
  nugget.estim = TRUE,
  
  covtype = "exp",
  
  control = list(
    
    trace = FALSE
    
  )
  
)


# =======================================================
# 19. Final UTS model
# =======================================================

cat(
  "\n===============================================\n"
)


cat(
  "Building final UTS model\n"
)


cat(
  "===============================================\n"
)


set.seed(243)


mf2_final <- km(
  
  formula = ~.,
  
  design = train.x,
  
  response = train.y[, "es2"],
  
  nugget.estim = TRUE,
  
  covtype = "exp",
  
  control = list(
    
    trace = FALSE
    
  )
  
)


# =======================================================
# 20. Inverse normalization
#
# IMPORTANT:
#
# df contains the negative original objectives:
#
#     df$EL  = -original EL
#     df$UTS = -original UTS
#
# Therefore:
#
#     normalized
#          ↓
#     negative physical value
#          ↓
#     multiply by -1
#
# Final values:
#
#     EL  (%)
#     UTS (MPa)
# =======================================================

fn_inverse_scale <- function(
    
  norm_data,
  
  df
  
){
  
  
  physical_negative <-
    
    norm_data *
    
    (
      max(df) -
        min(df)
    ) +
    
    min(df)
  
  
  return(
    
    physical_negative
    
  )
  
}


# =======================================================
# 21. Recover real EL and UTS values
# =======================================================

loocv_result$true_EL_real <-
  
  -fn_inverse_scale(
    
    loocv_result$true_es1,
    
    df[, "EL"]
    
  )


loocv_result$pred_EL_real <-
  
  -fn_inverse_scale(
    
    loocv_result$pred_es1,
    
    df[, "EL"]
    
  )


loocv_result$true_UTS_real <-
  
  -fn_inverse_scale(
    
    loocv_result$true_es2,
    
    df[, "UTS"]
    
  )


loocv_result$pred_UTS_real <-
  
  -fn_inverse_scale(
    
    loocv_result$pred_es2,
    
    df[, "UTS"]
    
  )


# =======================================================
# 22. Check inverse-normalized values
# =======================================================

cat(
  "\n===============================================\n"
)


cat(
  "Inverse-normalization check\n"
)


cat(
  "===============================================\n"
)


cat(
  "\nEL range in original dataset:\n"
)


cat(
  "Min EL =",
  
  min(
    loocv_result$true_EL_real,
    na.rm = TRUE
  ),
  
  "\n"
)


cat(
  "Max EL =",
  
  max(
    loocv_result$true_EL_real,
    na.rm = TRUE
  ),
  
  "\n"
)


cat(
  "\nUTS range in original dataset:\n"
)


cat(
  "Min UTS =",
  
  min(
    loocv_result$true_UTS_real,
    na.rm = TRUE
  ),
  
  "\n"
)


cat(
  "Max UTS =",
  
  max(
    loocv_result$true_UTS_real,
    na.rm = TRUE
  ),
  
  "\n"
)


cat(
  "\n===============================================\n"
)


# =======================================================
# 23. Save physical-space prediction results
# =======================================================

write.csv(
  
  loocv_result,
  
  "LOOCV_prediction_results_real.csv",
  
  row.names = FALSE
  
)


# =======================================================
# 24. Parity plot function
# =======================================================

plot_parity <- function(
    
  true,
  
  pred,
  
  xlab,
  
  ylab,
  
  panel_label,
  
  point_color
  
){
  
  
  d <- data.frame(
    
    True = true,
    
    Predicted = pred
    
  )
  
  
  r <- range(
    
    c(
      d$True,
      d$Predicted
    ),
    
    na.rm = TRUE
    
  )
  
  
  d_range <- diff(
    r
  )
  
  
  if(
    d_range == 0
  ){
    
    d_range <- 1
    
  }
  
  
  margin <- 0.06 *
    d_range
  
  
  xmin <- r[1] -
    margin
  
  
  xmax <- r[2] +
    margin
  
  
  # -----------------------------------------------------
  # 45-degree reference line
  # -----------------------------------------------------
  
  line_df <- data.frame(
    
    x = c(
      xmin,
      xmax
    ),
    
    y = c(
      xmin,
      xmax
    )
    
  )
  
  
  ggplot(
    
    d,
    
    aes(
      
      x = True,
      
      y = Predicted
      
    )
    
  ) +
    
    
    geom_line(
      
      data = line_df,
      
      aes(
        
        x = x,
        
        y = y
        
      ),
      
      inherit.aes = FALSE,
      
      linewidth = 1.2,
      
      linetype = "dashed",
      
      color = "#8B0000"
      
    ) +
    
    
    geom_point(
      
      size = 4.2,
      
      shape = 16,
      
      color = point_color
      
    ) +
    
    
    coord_fixed(
      
      xlim = c(
        xmin,
        xmax
      ),
      
      ylim = c(
        xmin,
        xmax
      ),
      
      expand = FALSE
      
    ) +
    
    
    labs(
      
      x = xlab,
      
      y = ylab
      
    ) +
    
    
    annotate(
      
      "text",
      
      x = xmin +
        0.04 *
        (
          xmax -
            xmin
        ),
      
      y = xmax -
        0.04 *
        (
          xmax -
            xmin
        ),
      
      label = panel_label,
      
      hjust = 0,
      
      vjust = 1,
      
      size = 8,
      
      fontface = "bold"
      
    ) +
    
    
    theme_classic(
      
      base_size = 21
      
    ) +
    
    
    theme(
      
      panel.background =
        element_blank(),
      
      plot.background =
        element_blank(),
      
      
      # -------------------------------------------------
      # Four-sided border
      # -------------------------------------------------
      
      panel.border =
        element_rect(
          
          colour = "black",
          
          fill = NA,
          
          linewidth = 1.0
          
        ),
      
      
      # -------------------------------------------------
      # Avoid duplicated axis lines
      # -------------------------------------------------
      
      axis.line =
        element_blank(),
      
      
      panel.grid.major =
        element_blank(),
      
      panel.grid.minor =
        element_blank(),
      
      
      # -------------------------------------------------
      # Inward ticks
      # -------------------------------------------------
      
      axis.ticks =
        element_line(
          
          colour = "black",
          
          linewidth = 0.85
          
        ),
      
      
      axis.ticks.length =
        unit(
          
          -0.18,
          
          "cm"
          
        ),
      
      
      # -------------------------------------------------
      # No ticks on top and right
      # -------------------------------------------------
      
      axis.ticks.x.top =
        element_blank(),
      
      axis.ticks.y.right =
        element_blank(),
      
      
      axis.text =
        element_text(
          
          colour = "black",
          
          size = 18
          
        ),
      
      
      axis.title =
        element_text(
          
          colour = "black",
          
          size = 22
          
        ),
      
      
      plot.margin =
        margin(
          
          5,
          
          5,
          
          5,
          
          5
          
        )
      
    ) +
    
    
    # ---------------------------------------------------
  # Top axis without labels or ticks
  # ---------------------------------------------------
  
  scale_x_continuous(
    
    sec.axis = dup_axis(
      
      name = NULL,
      
      labels = NULL
      
    )
    
  ) +
    
    
    # ---------------------------------------------------
  # Right axis without labels or ticks
  # ---------------------------------------------------
  
  scale_y_continuous(
    
    sec.axis = dup_axis(
      
      name = NULL,
      
      labels = NULL
      
    )
    
  )
  
}


# =======================================================
# 25. EL parity plot
#
# Plot real EL (%), NOT normalized es1
# =======================================================

p_EL <- plot_parity(
  
  true =
    loocv_result$true_EL_real,
  
  pred =
    loocv_result$pred_EL_real,
  
  xlab =
    "True EL (%)",
  
  ylab =
    "Predicted EL (%)",
  
  panel_label =
    "(a)",
  
  point_color =
    "#2166AC"
  
)


# =======================================================
# 26. UTS parity plot
#
# Plot real UTS (MPa), NOT normalized es2
# =======================================================

p_UTS <- plot_parity(
  
  true =
    loocv_result$true_UTS_real,
  
  pred =
    loocv_result$pred_UTS_real,
  
  xlab =
    "True UTS (MPa)",
  
  ylab =
    "Predicted UTS (MPa)",
  
  panel_label =
    "(b)",
  
  point_color =
    "#E69F00"
  
)


# =======================================================
# 27. Combine EL and UTS
#
# Final figure:
#
#       (a) EL       (b) UTS
#
# =======================================================

p_final <-
  
  p_EL +
  
  p_UTS +
  
  plot_layout(
    
    ncol = 2
    
  )


# =======================================================
# 28. Display final figure
# =======================================================

print(
  p_final
)


# =======================================================
# 29. Save combined figure
# =======================================================

ggsave(
  
  "Mg_EL_UTS_LOOCV_real.png",
  
  p_final,
  
  width = 11.5,
  
  height = 5.8,
  
  dpi = 600,
  
  bg = "white"
  
)






# =======================================================
# 32. End
# =======================================================

t2 <- proc.time()


cat(
  "\n===============================================\n"
)


cat(
  "Program finished.\n"
)


cat(
  "===============================================\n"
)


print(
  t2 - t1
)