library(ggplot2)
library(ggExtra)

# ============================
# 原始数据，不标准化
# ============================

data <- read.csv("41467_2025_56267_MOESM3_ESM.csv")
data <- data[,-c(1:2)]

# 注意你之前对第6、7列做了取负号，这里保留
data[,6:7] <- data[,6:7]
colnames(data)[6:7] <- c("Elongation", "UTS")  # 原始单位

# 设置参考目标点（原始单位）
tar_point00 <- c(16.5, 1190)  # 和你之前一样
special_point <- data.frame(
  Elongation = tar_point00[1],
  UTS        = tar_point00[2]
)

# ============================
# 提取非劣点集 (Pareto front)
# ============================

library(emoa) # 如果没有安装请 install.packages("emoa")
fmat <- data[, c("Elongation", "UTS")]
nd_points <- nondominated_points(t(fmat)) 
nd_df <- data.frame(
  Elongation = nd_points[1, ],
  UTS        = nd_points[2, ]
)

# ============================
# 基础散点图
# ============================

p <- ggplot(data, aes(x = Elongation, y = UTS)) +
  geom_point(color = "blue", alpha = 0.6) +        # 所有点
  geom_point(data = nd_df, aes(x = Elongation, y = UTS),
             color = "red", size = 2) +           # Pareto 点
  geom_point(data = special_point, aes(x = Elongation, y = UTS),
             color = "orange", shape = 17, size = 4) + # 目标点
  labs(
       x = "Elongation (%)",
       y = "UTS (MPa)") +
  theme_minimal(base_size = 15)

# ============================
# 添加边缘高斯分布
# ============================

p_marginal <- ggMarginal(
  p,
  type = "density",
  fill = "skyblue",
  color = "black",
  alpha = 0.6
)

# 显示图
p_marginal

# ============================
# 计算比目标点小的频率
# ============================

freqen_x <- mean(data$Elongation < tar_point00[1])
freqen_y <- mean(data$UTS < tar_point00[2])

cat("比目标值还小的频率:", freqen_x, freqen_y, "\n")

# ============================
# 参考直线斜率 (tan θ)
# ============================

ref00 <- c(max(data$Elongation), max(data$UTS))
tan_theta <- (ref00[2] - tar_point00[2]) / (ref00[1] - tar_point00[1])
tan_theta
