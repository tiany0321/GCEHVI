setwd("C:\\Users\\12758\\Desktop\\GCEHVI\\TiAlV\\specific\\more_iter\\sel2")

RNGkind("L'Ecuyer-CMRG")
set.seed(8)

source("fn_cbomo.R")
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
tar=t(tar_point00)
colnames(tar)=c("es1","es2")
tar_n=-tar
write.csv(tar_n,"tar_point00.csv")
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
pareto_vir=pareto.front(data=test.response.grid, obj)$front
pareto_vir_num=as.numeric(row.names(pareto_vir))
remain_num =which(test.response.grid[,"es1"]>0.4&test.response.grid[,"es2"]>0.4)



plot(
  response.grid[,1], response.grid[,2],
  col = "grey80", pch = 16,
  xlab = "es1", ylab = "es2"
)
points(response.grid[remain_num,1], response.grid[remain_num,2],
       col = "red", pch = 16)

##############determine the tar_point00
vir_cos_theta=cos_theta_calculate(vir_point0=pareto_vir,ref00,tar_point00,scale)
vir_opt_data=pareto_vir[which(vir_cos_theta[]==max(vir_cos_theta[])),]
vir_target=vir_opt_data
vir_target

set.seed(11*n_appr)
sample.num0<-vector()
for(ii in 1:1000){
  s <-  sample(remain_num,n_appr,replace=F,prob=NULL) 
  sample.num0<-as.data.frame(rbind(sample.num0,s))
}

#######################结果数据框准备

sel_all=c("all_EHVI","given_constrain_CEHVI","3points_dist_CEHVI","all_EI_EI","EI_EI","3points_dist_EI_EI")
ave.dis=matrix(,iter_num+1,2*length(sel_all))
ave.dis=as.data.frame(ave.dis)
colnames(ave.dis)=c("all_EHVI","all_EHVI_sd","GEHVI","GEHVI_sd","DEHVI","DEHVI_sd","all_EI_EI","all_EI_EI_sd","GEI_EI","GEI_EI_sd","DEI_EI","DEI_EI_sd")

ave.distance=matrix(,iter_num+1,2*length(sel_all))
ave.distance=as.data.frame(ave.distance)
colnames(ave.distance)=c("all_EHVI","all_EHVI_sd","GEHVI","GEHVI_sd","DEHVI","DEHVI_sd","all_EI_EI","all_EI_EI_sd","GEI_EI","GEI_EI_sd","DEI_EI","DEI_EI_sd")

ave.delta_hv=matrix(,iter_num+1,2*length(sel_all))
ave.delta_hv=as.data.frame(ave.delta_hv)
colnames(ave.delta_hv)=c("all_EHVI","all_EHVI_sd","GEHVI","GEHVI_sd","DEHVI","DEHVI_sd","all_EI_EI","all_EI_EI_sd","GEI_EI","GEI_EI_sd","DEI_EI","DEI_EI_sd")

ave.count=matrix(,1,2*length(sel_all))
ave.count=as.data.frame(ave.count)
colnames(ave.count)=c("all_EHVI","all_EHVI_sd","GEHVI","GEHVI_sd","DEHVI","DEHVI_sd","all_EI_EI","all_EI_EI_sd","GEI_EI","GEI_EI_sd","DEI_EI","DEI_EI_sd")

ave.delta_hv_given=matrix(,iter_num+1,2*length(sel_all))
ave.delta_hv_given=as.data.frame(ave.delta_hv_given)
colnames(ave.dis)=c("all_EHVI","all_EHVI_sd","GEHVI","GEHVI_sd","DEHVI","DEHVI_sd","all_EI_EI","all_EI_EI_sd","GEI_EI","GEI_EI_sd","DEI_EI","DEI_EI_sd")

colnames(ave.count)=c("all_EHVI","all_EHVI_sd","GEHVI","GEHVI_sd","DEHVI","DEHVI_sd","all_EI_EI","all_EI_EI_sd","GEI_EI","GEI_EI_sd","DEI_EI","DEI_EI_sd")


density=matrix(,repeat_num,length(sel_all))
density=as.data.frame(density)
colnames(density)=sel_all

testerror=matrix(,repeat_num,2)
testerror=as.data.frame(testerror)
colnames(testerror)=c("obj1","obj2")
##################不同的优化算法迭代
t_data=cbind(0,1)

sel=c(2)

#for(selector in sel){
  selector = sel
  if(sel_all[selector]=="given_constrain_CEHVI"){
    crit="given_constrain"
  }
  if(sel_all[selector]=="3points_dist_CEHVI"){
    crit="3points_dist"
  }  
  if(sel_all[selector]=="3points_angle_CEHVI"){
    crit="3points"
  } 
  if(sel_all[selector]=="all_EHVI"){
    crit="all_pareto"
  } 
  
  if(sel_all[selector]=="EI_EI"){
    crit="given_constrain"
  }
  if(sel_all[selector]=="3points_dist_EI_EI"){
    crit="3points_dist"
  }
  if(sel_all[selector]=="3points_angle_EI_EI"){
    crit="3points"
  }
  if(sel_all[selector]=="all_EI_EI"){
    crit="all_pareto"
  } 
  
  if(sel_all[selector]=="given_constrain_CEMI"){
    crit="given_constrain"}
  if(sel_all[selector]=="3points_dist_CEMI"){
    crit="3points_dist"
  }
  if(sel_all[selector]=="3points_angle_CEMI"){
    crit="3points"
  }
  
  if(sel_all[selector]=="given_constrain_CSMS"){
    crit="given_constrain"
  }
  if(sel_all[selector]=="3points_dist_CSMS"){
    crit="3points_dist"
  }
 ############
  if(sel_all[selector]=="SUR"){
    crit="given_constrain"
    
  }
  if(sel_all[selector]=="USeMO"){
    crit="given_constrain"
    
  }
  if(sel_all[selector]=="EHVI_sigmoid"){
    crit="given_constrain"
    
  }
  if(sel_all[selector]=="NSGA-II"){
    crit="given_constrain"
  }

  #######################每一种优化算法的结果数据框准备
  distance=matrix(0,iter_num+1,repeat_num)
  distance=as.data.frame(distance)
  dis=matrix(0,iter_num+1,repeat_num)
  dis=as.data.frame(dis)
  delta_hv=matrix(,iter_num+1,repeat_num)
  delta_hv=as.data.frame(delta_hv)
  count=matrix(iter_num,1,repeat_num)
  count=as.data.frame(count)
  
  delta_hv_given=matrix(,iter_num+1,repeat_num)
  delta_hv_given=as.data.frame(delta_hv_given)
  #####################不同的训练集训练
  #for(seed_num in 1:repeat_num){
  seed_num=44
    trial_seed <- 100000 +selector * 1000 +seed_num
    set.seed(trial_seed)
    
    sample.num<-sample.num0[,1:n_appr]
    
    split=sample.num[seed_num,]
    split=unlist(split)
    training.data <- data.all[split,]
    design.grid0 <- data.frame(training.data[,keep.x]) 
    response.grid0 <-  data.frame(training.data[,keep.y]) 
    response.grid0_noscale=-df[split,]
    write.csv(response.grid0_noscale,"data.train.csv") 
    
    
    test=data.all[-split,]
    test.grid0=test[,keep.x]
    test.response.grid0=test[,keep.y]
    
   
    ################model construction
    
    mf1_0 <- km(~., design = design.grid0, response = response.grid0[,1], noise.var=rep(10^(-5),nrow(design.grid0)),covtype = "exp")
    mf2_0 <- km(~., design = design.grid0, response = response.grid0[,2], noise.var=rep(10^(-5),nrow(design.grid0)),covtype = "exp")
    model0 = list(mf1_0, mf2_0)
    

    
    
    #############现有训练集的pareto前沿
    train_pareto0=pareto.front(data=response.grid0, obj)$front
    if(crit=="given_constrain"){
      pareto_P <- train_pareto0
    }
    
     pareto_P0=fn_pareto_P(crit,vir_target,train_data=response.grid0,train_pareto0,ref00,tar_point00,scale)
    ############现有训练集pareto前沿与参考点之间的超体积
    hv <- dominated_hypervolume(points=t(train_pareto0), ref=ref00)
    delta_hv[1,seed_num]=hv
    ###################基于可行域的超体积
    pareto_P0_min=fn_pareto_P_min(crit,ref_min,vir_target,train_data=response.grid0,train_pareto0,ref00,tar_point00,scale)
    hv_given <- dominated_hypervolume(points=t(pareto_P0_min), ref=ref00)
    delta_hv_given[1,seed_num]=hv_given
    #################（需要知道真实的目标点）离虚拟空间中目标点的距离及OC
   
    distances <- apply(train_pareto0, 1, function(x) sqrt((x[1] - vir_opt_data[1])^2 + (x[2] - vir_opt_data[2])^2))
    distance[1,seed_num]=min(unlist(distances))
    dis[1,seed_num]= distance[1,seed_num]
    
    ##############################多目标优化算法对虚拟空间的预测与评估
    iter=0
    train_pareto0=as.data.frame(train_pareto0)
    beta=0.5
    multi_obj_calculate0=fn_multi_obj(crit=NULL,vir_target,train_data=response.grid0,selector_num=sel_all[selector],test.grid0,train_pareto0,ref00,tar_point00,scale, model0)
    multi_obj_score0=multi_obj_calculate0$score
    model.result0=multi_obj_calculate0[,c("mean1","mean2")]
    
    ##########推荐的点
    multi_obj_num0=as.data.frame(multi_obj_score0)
    multi_obj_order0=as.data.frame(as.numeric(rownames(test.grid0)))
    colnames(multi_obj_order0)="order"
    multi_obj_num0=cbind(multi_obj_num0,multi_obj_order0)
    choose0=multi_obj_num0[order(-multi_obj_num0[,1]),]
    choose_p0=choose0[1,"order"]
    
  
     fn_plot(selector,test.response.grid,response.grid=response.grid0,model.result=model.result0,ref00,tar_point00,train_pareto=train_pareto0,pareto_P=pareto_P0,choose_p=choose_p0)
    
    ###############数据框准备与start iteration#########################    
    choose_p=choose_p0
    choose_data=as.data.frame(matrix(,500,1))
    choose_data[1,]=choose_p
    split0=split
    
    for(iter in 1:iter_num){
      ###############首先判断现有训练集中的点与目标点的距离是否在阈值内，阈值内说明找到了目标值，停止迭代
      if(distance[iter,seed_num]<=threshold){
        count[1,seed_num]=iter-1
        delta_hv[(iter+1):nrow(delta_hv),seed_num]=hv
        delta_hv_given[(iter+1):nrow(delta_hv_given),seed_num]=hv_given
        break
      }
      distance[iter+1,seed_num]=((vir_opt_data[1]-response.grid[choose_p,1])^2+(vir_opt_data[2]-response.grid[choose_p,2])^2)^0.5
      dis[iter+1,seed_num]=min(distance[1:(iter+1),seed_num])
      
      
      
      new.design.grid=test.grid[choose_p,]
      colnames(new.design.grid)=keep.x
      design.grid0 <- rbind(design.grid0,new.design.grid)
      ###########
      new.response.grid=response.grid[choose_p,]
      colnames(new.response.grid)=keep.y
      response.grid0 <- rbind(response.grid0,new.response.grid)
      ##############
      split0=c(split0,choose_p)
      test=data.all[-split0,]
      test.grid0=test[,keep.x]
      test.response.grid0=test[,keep.y]
      ###########
      ################model construction
      mf1 <- km(~., design = design.grid0, response = response.grid0[,1], noise.var=rep(10^(-5),nrow(design.grid0)),covtype = "exp")
      mf2 <- km(~., design = design.grid0, response = response.grid0[,2], noise.var=rep(10^(-5),nrow(design.grid0)),covtype = "exp")
      model = list(mf1, mf2)
      
      
      
      train_pareto=pareto.front(data=response.grid0, obj)$front
      pareto_P=fn_pareto_P(crit,vir_target,train_data=response.grid0,train_pareto,ref00,tar_point00,scale)

      hv <- dominated_hypervolume(points=t(train_pareto), ref=ref00)
      delta_hv[iter+1,seed_num]=hv
      
      ###################基于可行域的超体积
      pareto_P_min=fn_pareto_P_min(crit,ref_min,vir_target,train_data=response.grid0,train_pareto,ref00,tar_point00,scale)

      hv_given <- dominated_hypervolume(points=t(pareto_P_min), ref=ref00)
 
      delta_hv_given[iter+1,seed_num]=hv_given
      ###########新算法
      train_pareto=as.data.frame(train_pareto)
      multi_obj_calculate=fn_multi_obj(crit,vir_target,train_data=response.grid0,sel_all[selector],test.grid0,train_pareto,ref00,tar_point00,scale, model)
      multi_obj_score=multi_obj_calculate$score
      
      ##########新pareto前沿推荐的点
      multi_obj_num=as.data.frame(multi_obj_score)
      multi_obj_order=as.data.frame(as.numeric(rownames(test.grid0)))
      colnames(multi_obj_order)="order"
      multi_obj_num=cbind(multi_obj_num,multi_obj_order)

      choose=multi_obj_num[order(-multi_obj_num[,1]),]
      choose_p=choose[1,"order"]
      choose_data[iter+1,]=choose_p
       fn_plot(selector,test.response.grid,response.grid,model.result,ref00,tar_point00,train_pareto,pareto_P,choose_p)
      
    }
    iter_response=response.grid0[-c(1:n_appr),]
    iter_response_order=iter_response
    row.names(iter_response_order)=1:nrow(iter_response_order)
    iter_response_noscale=-df[as.numeric(rownames(iter_response)),]
    write.csv(iter_response_noscale,"data.new.csv")
  


t2=proc.time()
t=t2-t1
print(paste0('执行时间：',t[3][[1]],'秒'))
write.csv(t[3][[1]],"total_time.csv")