

setwd(
  "C:\\Users\\12758\\Desktop\\GCEHVI\\hl\\performance"
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

tar_point00 <- c(
  
  25.40,
  -23.60
  
)


iter_num <- 80

repeat_num <- 100

threshold <- 0

scale <- TRUE


obj <- c(
  
  "es1",
  "es2"
  
)


keep.y <- obj


# =======================================================
# 3. Read data
# =======================================================

data0 <- read.csv(
  
  "dataall_train.csv"
  
)


# =======================================================
# 4. Variables
# =======================================================

keep.heat <- c(
  
  "Lambda",
  "VEN",
  "delta.x",
  "ou"
  
)


keep.hysteresis <- c(
  
  "delta.r",
  "VEN",
  "Lambda"
  
)


keep.name <- c(
  
  "Lambda",
  "delta.x",
  "VEN",
  "ou",
  "delta.r"
  
)


keep.x <- keep.name


# =======================================================
# 5. Construct data
# =======================================================

test.grid <- data0[
  
  ,
  
  keep.x,
  
  drop = FALSE
  
]


response.grid <- data0[
  
  ,
  
  keep.y,
  
  drop = FALSE
  
]


colnames(test.grid) <- keep.x


colnames(response.grid) <- keep.y


data <- cbind(
  
  test.grid,
  
  response.grid
  
)


colnames(data) <- c(
  
  keep.x,
  
  keep.y
  
)


# =======================================================
# 6. Construct training set
#
# IMPORTANT:
#
# No independent test set is used.
#
# Training criterion:
#
#   es1 > 35
#   es2 > -30
# =======================================================


remain_num <- which(
  
  response.grid[, "es1"] > 35 &
    
    response.grid[, "es2"] > -30
  
)


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
# 7. Dataset information
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
  nrow(data),
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
      nrow(data) *
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
# 8. Noise parameters
#
# These are model noise parameters, NOT evaluation
# error metrics, so they are retained.
# =======================================================

noiselevel_hysteresis <- 0.03 * (
  
  max(
    response.grid[, "es1"]
  ) -
    
    min(
      response.grid[, "es1"]
    )
  
)


noiselevel_heat <- 0.07 * (
  
  max(
    response.grid[, "es2"]
  ) -
    
    min(
      response.grid[, "es2"]
    )
  
)


cat(
  "\n===============================================\n"
)


cat(
  "Noise parameters\n"
)


cat(
  "===============================================\n"
)


cat(
  "Hysteresis (es1) noise SD :",
  noiselevel_hysteresis,
  "\n"
)


cat(
  "Latent heat (es2) noise SD:",
  noiselevel_heat,
  "\n"
)


cat(
  "===============================================\n"
)


# =======================================================
# 9. Kriging formulas
# =======================================================

formula1 <- es1 ~ delta.r + VEN + Lambda


formula2 <- es2 ~ Lambda + VEN + delta.x + ou


# =======================================================
# 10. Training-set visualization
# =======================================================

plot(
  
  response.grid[, 1],
  
  response.grid[, 2],
  
  col = "grey80",
  
  pch = 16,
  
  xlab = "Hysteresis",
  
  ylab = "Latent heat"
  
)


points(
  
  response.grid[
    train_num,
    1
  ],
  
  response.grid[
    train_num,
    2
  ],
  
  col = "#2166AC",
  
  pch = 16
  
)


# =======================================================
# 11. LOOCV
#
# Only the selected training set is used.
#
# No R2 / RMSE / MAE / MSE calculation.
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
  # Hysteresis: es1
  # =====================================================
  
  noise_es1_cv <- rep(
    
    noiselevel_hysteresis^2,
    
    length(train_index)
    
  )
  
  
  set.seed(243)
  
  
  mf1_cv <- km(
    
    formula = formula1,
    
    design = x_cv[
      
      ,
      
      c(
        "delta.r",
        "VEN",
        "Lambda"
      ),
      
      drop = FALSE
      
    ],
    
    response = y_cv_es1,
    
    noise.var = noise_es1_cv,
    
    covtype = "exp",
    
    control = list(
      
      trace = FALSE
      
    )
    
  )
  
  
  pred_es1 <- predict(
    
    mf1_cv,
    
    newdata = x_valid[
      
      ,
      
      c(
        "delta.r",
        "VEN",
        "Lambda"
      ),
      
      drop = FALSE
      
    ],
    
    type = "UK"
    
  )
  
  
  loocv_pred_es1[i] <- as.numeric(
    
    pred_es1$mean
    
  )
  
  
  # =====================================================
  # Latent heat: es2
  # =====================================================
  
  noise_es2_cv <- rep(
    
    noiselevel_heat^2,
    
    length(train_index)
    
  )
  
  
  set.seed(243)
  
  
  mf2_cv <- km(
    
    formula = formula2,
    
    design = x_cv[
      
      ,
      
      c(
        "Lambda",
        "VEN",
        "delta.x",
        "ou"
      ),
      
      drop = FALSE
      
    ],
    
    response = y_cv_es2,
    
    noise.var = noise_es2_cv,
    
    covtype = "exp",
    
    control = list(
      
      trace = FALSE
      
    )
    
  )
  
  
  pred_es2 <- predict(
    
    mf2_cv,
    
    newdata = x_valid[
      
      ,
      
      c(
        "Lambda",
        "VEN",
        "delta.x",
        "ou"
      ),
      
      drop = FALSE
      
    ],
    
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
# 12. LOOCV results
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
# 13. Save LOOCV prediction results
# =======================================================

write.csv(
  
  loocv_result,
  
  "LOOCV_prediction_results.csv",
  
  row.names = FALSE
  
)


# =======================================================
# 14. Final hysteresis model
# =======================================================

cat(
  "\n===============================================\n"
)


cat(
  "Building final hysteresis model\n"
)


cat(
  "===============================================\n"
)


noise_es1_final <- rep(
  
  noiselevel_hysteresis^2,
  
  nrow(train.x)
  
)


set.seed(243)


mf1_final <- km(
  
  formula = formula1,
  
  design = train.x[
    
    ,
    
    c(
      "delta.r",
      "VEN",
      "Lambda"
    ),
    
    drop = FALSE
    
  ],
  
  response = train.y[, "es1"],
  
  noise.var = noise_es1_final,
  
  covtype = "exp",
  
  control = list(
    
    trace = FALSE
    
  )
  
)


# =======================================================
# 15. Final latent heat model
# =======================================================

cat(
  "\n===============================================\n"
)


cat(
  "Building final latent heat model\n"
)


cat(
  "===============================================\n"
)


noise_es2_final <- rep(
  
  noiselevel_heat^2,
  
  nrow(train.x)
  
)


set.seed(243)


mf2_final <- km(
  
  formula = formula2,
  
  design = train.x[
    
    ,
    
    c(
      "Lambda",
      "VEN",
      "delta.x",
      "ou"
    ),
    
    drop = FALSE
    
  ],
  
  response = train.y[, "es2"],
  
  noise.var = noise_es2_final,
  
  covtype = "exp",
  
  control = list(
    
    trace = FALSE
    
  )
  
)


# =======================================================
# 16. Parity plot function
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
  
  
  # -----------------------------------------------------
  # Determine plotting range
  # -----------------------------------------------------
  
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
  
  
  # -----------------------------------------------------
  # Plot
  # -----------------------------------------------------
  
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
# 17. Hysteresis parity plot
# =======================================================

p_hysteresis <- plot_parity(
  
  true =
    loocv_result$true_es1,
  
  pred =
    loocv_result$pred_es1,
  
  xlab =
    "True hysteresis (MPa)",
  
  ylab =
    "Predicted hysteresis (MPa)",
  
  panel_label =
    "(a)",
  
  point_color =
    "#2166AC"
  
)


# =======================================================
# 18. Latent heat parity plot
# =======================================================

p_heat <- plot_parity(
  
  true =
    loocv_result$true_es2,
  
  pred =
    loocv_result$pred_es2,
  
  xlab =
    "True latent heat",
  
  ylab =
    "Predicted latent heat",
  
  panel_label =
    "(b)",
  
  point_color =
    "#E69F00"
  
)


# =======================================================
# 19. Combine two performance plots
#
# Final figure:
#
#       (a) Hysteresis | (b) Latent heat
#
# =======================================================

p_final <-
  
  p_hysteresis +
  
  p_heat +
  
  plot_layout(
    
    ncol = 2
    
  )


# =======================================================
# 20. Display final figure
# =======================================================

print(
  p_final
)


# =======================================================
# 21. Save combined figure
# =======================================================

ggsave(
  
  "Hysteresis_LatentHeat_LOOCV.png",
  
  p_final,
  
  width = 11.5,
  
  height = 5.8,
  
  dpi = 600,
  
  bg = "white"
  
)






# =======================================================
# 24. End
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