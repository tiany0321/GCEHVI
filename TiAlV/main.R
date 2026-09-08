setwd("C:\\Users\\12758\\Desktop\\GCEHVI\\TiAlV")

RNGkind("L'Ecuyer-CMRG")
set.seed(8)


t1=proc.time()

#############参考点设置
#
iter_num=50#50
repeat_num=100#100
threshold=0


data=read.csv("41467_2025_56267_MOESM3_ESM.csv")
data=data[,-c(1:2)]
data[,c(6,7)]=-data[,c(6,7)]
data[,c(7)]=data[,c(7)]
df=data[,c("Elongation","UTS")]
fn_scale=function(scale_data,df){
  es1_norm <- (scale_data - min(df)) / (max(df) - min(df))
  return(es1_norm)
}
data$Elongation<-fn_scale(data[,"Elongation"],df[,"Elongation"])
data$UTS<-fn_scale(data[,"UTS"],df[,"UTS"])
#############参考点设置
tar_point00=c(-12.90, -1221)#c(-8.90, -1247)#c(-16.5,  -1190)
tar_point00=c(fn_scale(tar_point00[1],df[,"Elongation"]),fn_scale(tar_point00[2],df[,"UTS"]))

colnames(data)[c(6,7)]=c("es1" ,"es2")
scale=TRUE
obj=c("es1","es2")
keep.x=colnames(data)[c(1:5)]
keep.y=obj

size=0.1

#########虚拟空间构造
test.grid <- data[,keep.x]
colnames(test.grid)=keep.x
n_appr <- round(size*nrow(test.grid))


response.grid<- data[,keep.y]
colnames(response.grid)=keep.y



n_appr <- round(size*nrow(test.grid))


# 3. 提取非劣点集
nd_points <- nondominated_points(t(response.grid)) # emoa 要输入目标在列
nd_df <- data.frame(f1 = nd_points[1, ], f2 = nd_points[2, ])



response.grid <-data[,keep.y] 
colnames(response.grid)=keep.y
plot(response.grid[,1],response.grid[,2],xlim=c(0,1),ylim=c(0,1))
par(new=T)
plot(nd_df[,1],nd_df [,2] ,xlim=c(0,1),ylim=c(0,1),col="red")

extend.per_fa = 0.2  # 扩大 20%
range1 = max(response.grid[,1]) - min(response.grid[,1])
range2 = max(response.grid[,2]) - min(response.grid[,2])

ref00=c(max(response.grid[,1]),max(response.grid[,2]))
ref_min=c(min(response.grid[,1])-extend.per_fa * range1, 
          min(response.grid[,2])- extend.per_fa * range2)

data.all=cbind(test.grid,response.grid)

order.data=as.data.frame(rep(1:nrow(data.all)))
colnames(order.data)="order.num"
data.all=cbind(data.all,order.data)

#########（只有在已知虚拟空间时可以得到）虚拟空间所有点的真实值、真实pareto前沿、和目标最佳值

test.response.grid=response.grid

remain_num =which(test.response.grid[,"es1"]>0.4&test.response.grid[,"es2"]>0.4)

write.csv(remain_num,"remain_num.csv")