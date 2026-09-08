setwd("C:\\Users\\12758\\Desktop\\GCEHVI\\Mg")

RNGkind("L'Ecuyer-CMRG")
set.seed(8)

source("fn_cbomo.R")
t1=proc.time()

#############参考点设置

iter_num=100#50
repeat_num=100#100
threshold=0


data=read.csv("df_last.csv")
#data=data[which(data[,"homo_T"]==0),]
data=data[,-c(16)]
data[,c(15,16)]=-data[,c(15,16)]

df=data[,c("EL","UTS")]
fn_scale=function(scale_data,df){
  es1_norm <- (scale_data - min(df)) / (max(df) - min(df))
  return(es1_norm)
}
data$EL<-fn_scale(data[,"EL"],df[,"EL"])
data$UTS<-fn_scale(data[,"UTS"],df[,"UTS"])
#############参考点设置
tar_point0=c(-26.62, -309.73)#c(-29.00,-290.00)#c(-22.30, -325.10)#c(-26.62, -309.73)
tar_point00=c(fn_scale(tar_point0[1],df[,"EL"]),fn_scale(tar_point0[2],df[,"UTS"]))

colnames(data)[c(15,16)]=c("es1" ,"es2")
scale=TRUE
obj=c("es1","es2")
keep.x=colnames(data)[c(1:5,7:14)]
keep.y=obj

size=0.1

#########虚拟空间构造
test.grid <- data[,keep.x]
colnames(test.grid)=keep.x
n_appr <- round(size*nrow(test.grid))

#write.csv(test.grid,"test.grid.csv")
response.grid<- data[,keep.y]
colnames(response.grid)=keep.y



n_appr <- round(size*nrow(test.grid))

#write.csv(test.grid,"test.grid.csv")
#response.grid<- df
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
#write.csv(response.grid,"response.grid.csv")
data.all=cbind(test.grid,response.grid)

order.data=as.data.frame(rep(1:nrow(data.all)))
colnames(order.data)="order.num"
data.all=cbind(data.all,order.data)

#########（只有在已知虚拟空间时可以得到）虚拟空间所有点的真实值、真实pareto前沿、和目标最佳值
#test.response.grid <- as.data.frame(t(apply(test.grid, 1, f_name)))
#colnames(test.response.grid)=obj
# test.response.grid = fn_response(obj_target=obj_opt,test.response.grid)
test.response.grid=response.grid
pareto_vir=pareto.front(data=test.response.grid, obj)$front
pareto_vir_num=as.numeric(row.names(pareto_vir))
#remain_num <- setdiff(c(1:nrow(data.all)), pareto_vir_num)
ref_es1=(max(test.response.grid[,"es1"])-min(test.response.grid[,"es1"]))*0.45
ref_es2=(max(test.response.grid[,"es2"])-min(test.response.grid[,"es2"]))*0.45
remain_num =which(test.response.grid[,"es1"]>ref_es1&test.response.grid[,"es2"]>ref_es2)
#remain_num <- setdiff(c(1:nrow(data.all)), remain_region_num)
write.csv(remain_num,"remain_num.csv")

plot(
  response.grid[,1], response.grid[,2],
  col = "grey80", pch = 16,
  xlab = "es1", ylab = "es2"
)
points(response.grid[remain_num,1], response.grid[remain_num,2],
       col = "red", pch = 16)
