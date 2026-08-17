
##################library
library(GPareto)
library(DiceDesign)
library(DiceKriging)
library(dplyr)
library(RColorBrewer)
library(MASS)
library(mco)              # NSGA-II 多目标优化
library(emoa)             # Hypervolume 计算
library(randtoolbox)      # 替代 lhs::sobol
library(lhs)              # 用于 Sobol 序列生成
#library(DiceOptim)
display.brewer.all(type = "seq")

fn_dens=function(iter_response,tar_point00){
  x=iter_response[,1] 
  y=iter_response[,2]
  dens <- kde2d(x, y, n = 100)
  tar_point00_x=tar_point00[1]
  tar_point00_y=tar_point00[2]
  
  nearest_index_x <- which.min(abs(dens$x - tar_point00_x))
  nearest_index_y <- which.min(abs(dens$y - tar_point00_y))
  
  density_at_point <- dens$z[nearest_index_x, nearest_index_y]
  return(density_at_point)
}

######both_min
#####################plot function
fn_plot=function(selector,test.response.grid,response.grid,model.result,ref00,tar_point00,train_pareto,pareto_P,choose_p){
  
  x_lim=c(-20,-4)
  y_lim=c(-10,25)
  plot(test.response.grid[,"es1"],test.response.grid[,"es2"],cex=1.5, col=brewer.pal(9, "YlOrRd")[3],pch=16,xlab="Objective1",ylab="Objective2",ylim=y_lim,xlim=x_lim, lwd=3, tck=0.01,font=2,font.lab=2,cex.axis=1.2,cex.lab=1.2)
  
  par(new=T)
  plot(response.grid0[,"es1"],response.grid0[,"es2"], cex=1.5,col=brewer.pal(9, "Greens")[5],pch=16,xlab="Objective1",ylab="Objective2",ylim=y_lim,xlim=x_lim, lwd=3, tck=0.01,font=2,font.lab=2,cex.axis=1.2,cex.lab=1.2)
  
  #####reference point
  par(new=T)
  points(ref00[1],ref00[2],cex=1.5,pch=15,col="red",xlab="Objective1",ylab="Objective2",ylim=y_lim,xlim=x_lim, lwd=3, tck=0.01,font=2,font.lab=2,cex.axis=1.2,cex.lab=1.2)
  ######ideal point
  par(new=T)
  points(tar_point00[1],tar_point00[2],cex=1.5,pch=17,col="red",xlab="Objective1",ylab="Objective2",ylim=y_lim,xlim=x_lim, lwd=3, tck=0.01,font=2,font.lab=2,cex.axis=1.2,cex.lab=1.2)
  
  ###########pareto for training data
  par(new=T)
  plot(train_pareto[,"es1"],train_pareto[,"es2"],cex=1.5, col=brewer.pal(9, "YlOrRd")[6],type="p",pch=16,xlab="Objective1",ylab="Objective2",ylim=y_lim,xlim=x_lim, lwd=3, tck=0.01,font=2,font.lab=2,cex.axis=1.2,cex.lab=1.2)
  par(new=T)
  plot(train_pareto[,"es1"],train_pareto[,"es2"],cex=1.5,type="s",pch=16,col="red",xlab="Objective1",ylab="Objective2",ylim=y_lim,xlim=x_lim, lwd=3, tck=0.01,font=2,font.lab=2,cex.axis=1.2,cex.lab=1.2)
  #QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQq
  ###########current best trade-off options
  pareto_P=pareto_P[rev(1:nrow(pareto_P)), ]
  par(new=T)
  plot(pareto_P[,"es1"],pareto_P[,"es2"],cex=1.5, col="cyan",type="p",pch=16,xlab="Objective1",ylab="Objective2",ylim=y_lim,xlim=x_lim, lwd=3, tck=0.01,font=2,font.lab=2,cex.axis=1.2,cex.lab=1.2)
  par(new=T)
  plot(pareto_P[,"es1"],pareto_P[,"es2"],cex=1.5,type="s",pch=16,col="cyan",xlab="Objective1",ylab="Objective2",ylim=y_lim,xlim=x_lim, lwd=3, tck=0.01,font=2,font.lab=2,cex.axis=1.2,cex.lab=1.2)
  #QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQq
  ##########constrain
  par(new=T)
  points(test.response.grid[choose_p,1],test.response.grid[choose_p,2],cex=1.5,pch=17,col="blue",xlab="Objective1",ylab="Objective2",ylim=y_lim,xlim=x_lim, lwd=3, tck=0.01,font=2,font.lab=2,cex.axis=1.2,cex.lab=1.2)
  
  text(x = -20, y = 25, labels = paste(selector,as.character(iter),sep="-iter"),cex = 1,font=1.7)
  #par(new=T)
 # arrows(ref00[1], ref00[2], tar_point00[1], tar_point00[2], col = "red")
}

#######pareto front function
pareto.front=function(data=response, obj){
  response_obj=as.data.frame(data[,obj])
  
  response_obj=response_obj[order(response_obj$es1,response_obj$es2,decreasing=FALSE),]
  num=which(!duplicated(cummin(response_obj$es2)))
  front = response_obj[num,]
  front_num=as.numeric(rownames(front))
  
  return(list(front=front,front_num=front_num))
}

###########decide the coordinates of the options which will construct the constrains
fn_pareto_P <- function(crit=NULL,vir_target,train_data,train_pareto,ref00,tar_point00,scale){
  
  if(crit=="3points"){
    ###########decide the current best trade-off option by maximizing cos_theta
    cos_theta=cos_theta_calculate(vir_point0=train_pareto,ref00,tar_point00,scale)
    opt_data=train_pareto[which(cos_theta[]==max(cos_theta[])),]
    
    ###########range the optimal options in training data according to cos_theta
    train_pareto_th=cbind(train_pareto,cos_theta)
    train_pareto_th_order=train_pareto_th[order(-train_pareto_th[,"cos_theta"]),]
    ###########choose three current best trade-off options
    point2=train_pareto_th_order[which(train_pareto_th_order[,"cos_theta"]==max(train_pareto_th_order[,"cos_theta"])),]
    train_pareto_th_order$abs=train_pareto_th_order[,"es1"]-point2[,"es1"]
    point3_c=train_pareto_th_order[which(train_pareto_th_order[,"abs"]<0),]
    point3=point3_c[1,]
    point1_c=train_pareto_th_order[which(train_pareto_th_order[,"abs"]>0),]
    point1=point1_c[1,]
    
    keep=c("es1","es2")
    pareto_P0=t(as.data.frame(c(point1[,"es1"],-Inf)))
    pareto_P1=point1[,keep]
    pareto_P2=point2[,keep]
    pareto_P3=point3[,keep]
    pareto_P4=t(as.data.frame(c(-Inf,point3[,"es2"])))
    
    colnames(pareto_P0)=keep
    colnames(pareto_P4)=keep
    ###########return the coordinates of the options which will construct the constrains
    pareto_P=as.data.frame(rbind(pareto_P0,pareto_P1,pareto_P2,pareto_P3,pareto_P4))
    if(is.na(pareto_P[2,1])){
      pareto_P[1,1]=pareto_P[3,1]
      pareto_P=pareto_P[-2,]
    }
    if(is.na(pareto_P[4,2])){
      pareto_P[5,2]=pareto_P[3,2]
      pareto_P=pareto_P[-4,]
    }
  }
  
  if(crit!="3points"){
    
    if(crit=="dist"){
      dist=unlist(apply(train_data,1,fn_o_dist,vir_target))
      opt_data=train_data[which(dist[]==min(dist[])),]
      
      pareto_P1=as.data.frame(c(opt_data[1],-Inf)) # pareto_P1=as.data.frame(c(opt_data[1],tar_point00[2]))
      pareto_P2=as.data.frame(c(opt_data[1],opt_data[2]))
      pareto_P3=as.data.frame(c(-Inf,opt_data[2]))#  pareto_P3=as.data.frame(c(tar_point00[1],opt_data[2]))
      colnames(pareto_P1)=colnames(opt_data)
      colnames(pareto_P3)=colnames(opt_data)
      pareto_P=as.data.frame(rbind(pareto_P1,pareto_P2,pareto_P3))
    }
    
    if(crit=="angle"){
      cos_theta=cos_theta_calculate(vir_point0=train_pareto,ref00,tar_point00,scale)
      opt_data=train_pareto[which(cos_theta[]==max(cos_theta[])),]
      
      pareto_P1=as.data.frame(c(opt_data[1],-Inf)) # pareto_P1=as.data.frame(c(opt_data[1],tar_point00[2]))
      pareto_P2=as.data.frame(c(opt_data[1],opt_data[2]))
      pareto_P3=as.data.frame(c(-Inf,opt_data[2]))#  pareto_P3=as.data.frame(c(tar_point00[1],opt_data[2]))
      colnames(pareto_P1)=colnames(opt_data)
      colnames(pareto_P3)=colnames(opt_data)
      pareto_P=as.data.frame(rbind(pareto_P1,pareto_P2,pareto_P3))
    }

    if(crit=="3points_dist"){
      ###########decide the current best trade-off option by minimizing the distance

      dist=t(as.data.frame((apply(train_pareto,1,fn_o_dist,vir_target))))
      colnames(dist)="distance"
      rownames(dist)=rep(1:nrow(dist))
      ###########choose three current best trade-off options
      sorted_indices <- sort(dist[,1], index.return = TRUE)$ix
      if(length(sorted_indices)>=3){
      opt_data=train_pareto[sorted_indices[1:3],]
      }else{
        opt_data=train_pareto[sorted_indices[1:length(sorted_indices)],]
      }
      ###########range the optimal options in training data
      opt_data_order=opt_data[order(-opt_data[,"es1"]),]
      
      keep=c("es1","es2")
      pareto_P0=t(as.data.frame(c(opt_data_order[1,"es1"],-Inf)))
      pareto_P4=t(as.data.frame(c(-Inf,opt_data_order[nrow(opt_data_order),"es2"])))
      colnames(pareto_P0)=keep
      colnames(pareto_P4)=keep
      

      ###########return the coordinates of the options which will construct the constrains
      pareto_P=as.data.frame(rbind(pareto_P0, opt_data_order,pareto_P4))
    }
    
    if(crit=="given_constrain"){

        
        condition <- train_pareto$es1 < given_point[1] & train_pareto$es2 < given_point[2]
        row_indices <- which(condition)
        
        if(length(row_indices)>=1){
        opt_data=train_pareto[row_indices,]
        ###########range the optimal options in training data
        opt_data_order=opt_data[order(-opt_data[,"es1"]),]
        
        keep=c("es1","es2")
        pareto_P0=t(as.data.frame(c(given_point[1],-Inf)))
        pareto_P1=t(as.data.frame(c(given_point[1],opt_data_order[1,"es2"])))
        pareto_P3=t(as.data.frame(c(opt_data_order[nrow(opt_data_order),"es1"],given_point[2])))
        pareto_P4=t(as.data.frame(c(-Inf,given_point[2])))
        colnames(pareto_P0)=keep
        colnames(pareto_P1)=keep
        colnames(pareto_P3)=keep
        colnames(pareto_P4)=keep
        
     
      ###########return the coordinates of the options which will construct the constrains
      pareto_P=as.data.frame(rbind(pareto_P0,pareto_P1,opt_data_order,pareto_P3,pareto_P4))
        }else{
          pareto_P0=t(as.data.frame(c(given_point[1],-Inf)))
          pareto_P1=t(as.data.frame(given_point))
          pareto_P4=t(as.data.frame(c(-Inf,given_point[2])))
          keep=c("es1","es2")
          colnames(pareto_P0)=keep
          colnames(pareto_P1)=keep
          colnames(pareto_P4)=keep
          pareto_P=as.data.frame(rbind(pareto_P0,pareto_P1,pareto_P4))
        }
          
        
    }
  }
  if(crit=="all_pareto"){
    pareto_P=train_pareto
  }
  
  
  return(pareto_P)
}
fn_pareto_P_min <- function(crit=NULL,ref_min,vir_target,train_data,train_pareto,ref00,tar_point00,scale){
  
  if(crit=="3points"){
    ###########decide the current best trade-off option by maximizing cos_theta
    cos_theta=cos_theta_calculate(vir_point0=train_pareto,ref00,tar_point00,scale)
    opt_data=train_pareto[which(cos_theta[]==max(cos_theta[])),]
    
    ###########range the optimal options in training data according to cos_theta
    train_pareto_th=cbind(train_pareto,cos_theta)
    train_pareto_th_order=train_pareto_th[order(-train_pareto_th[,"cos_theta"]),]
    ###########choose three current best trade-off options
    point2=train_pareto_th_order[which(train_pareto_th_order[,"cos_theta"]==max(train_pareto_th_order[,"cos_theta"])),]
    train_pareto_th_order$abs=train_pareto_th_order[,"es1"]-point2[,"es1"]
    point3_c=train_pareto_th_order[which(train_pareto_th_order[,"abs"]<0),]
    point3=point3_c[1,]
    point1_c=train_pareto_th_order[which(train_pareto_th_order[,"abs"]>0),]
    point1=point1_c[1,]
    
    keep=c("es1","es2")
    pareto_P0=t(as.data.frame(c(point1[,"es1"],-Inf)))
    pareto_P1=point1[,keep]
    pareto_P2=point2[,keep]
    pareto_P3=point3[,keep]
    pareto_P4=t(as.data.frame(c(-Inf,point3[,"es2"])))
    
    colnames(pareto_P0)=keep
    colnames(pareto_P4)=keep
    ###########return the coordinates of the options which will construct the constrains
    pareto_P=as.data.frame(rbind(pareto_P0,pareto_P1,pareto_P2,pareto_P3,pareto_P4))
    if(is.na(pareto_P[2,1])){
      pareto_P[1,1]=pareto_P[3,1]
      pareto_P=pareto_P[-2,]
    }
    if(is.na(pareto_P[4,2])){
      pareto_P[5,2]=pareto_P[3,2]
      pareto_P=pareto_P[-4,]
    }
  }
  
  if(crit!="3points"){
    
    if(crit=="dist"){
      dist=unlist(apply(train_data,1,fn_o_dist,vir_target))
      opt_data=train_data[which(dist[]==min(dist[])),]
      
      pareto_P1=as.data.frame(c(opt_data[1],-Inf)) # pareto_P1=as.data.frame(c(opt_data[1],tar_point00[2]))
      pareto_P2=as.data.frame(c(opt_data[1],opt_data[2]))
      pareto_P3=as.data.frame(c(-Inf,opt_data[2]))#  pareto_P3=as.data.frame(c(tar_point00[1],opt_data[2]))
      colnames(pareto_P1)=colnames(opt_data)
      colnames(pareto_P3)=colnames(opt_data)
      pareto_P=as.data.frame(rbind(pareto_P1,pareto_P2,pareto_P3))
    }
    
    if(crit=="angle"){
      cos_theta=cos_theta_calculate(vir_point0=train_pareto,ref00,tar_point00,scale)
      opt_data=train_pareto[which(cos_theta[]==max(cos_theta[])),]
      
      pareto_P1=as.data.frame(c(opt_data[1],-Inf)) # pareto_P1=as.data.frame(c(opt_data[1],tar_point00[2]))
      pareto_P2=as.data.frame(c(opt_data[1],opt_data[2]))
      pareto_P3=as.data.frame(c(-Inf,opt_data[2]))#  pareto_P3=as.data.frame(c(tar_point00[1],opt_data[2]))
      colnames(pareto_P1)=colnames(opt_data)
      colnames(pareto_P3)=colnames(opt_data)
      pareto_P=as.data.frame(rbind(pareto_P1,pareto_P2,pareto_P3))
    }
    
    if(crit=="3points_dist"){
      ###########decide the current best trade-off option by minimizing the distance
      
      dist=t(as.data.frame((apply(train_pareto,1,fn_o_dist,vir_target))))
      colnames(dist)="distance"
      rownames(dist)=rep(1:nrow(dist))
      ###########choose three current best trade-off options
      sorted_indices <- sort(dist[,1], index.return = TRUE)$ix
      if(length(sorted_indices)>=3){
        opt_data=train_pareto[sorted_indices[1:3],]
      }else{
        opt_data=train_pareto[sorted_indices[1:length(sorted_indices)],]
      }
      ###########range the optimal options in training data
      opt_data_order=opt_data[order(-opt_data[,"es1"]),]
      
      keep=c("es1","es2")
      pareto_P0=t(as.data.frame(c(opt_data_order[1,"es1"],-Inf)))
      pareto_P4=t(as.data.frame(c(-Inf,opt_data_order[nrow(opt_data_order),"es2"])))
      colnames(pareto_P0)=keep
      colnames(pareto_P4)=keep
      
      
      ###########return the coordinates of the options which will construct the constrains
      pareto_P=as.data.frame(rbind(pareto_P0, opt_data_order,pareto_P4))
    }
    
    if(crit=="given_constrain"){
      
      
      condition <- train_pareto$es1 < given_point[1] & train_pareto$es2 < given_point[2]
      row_indices <- which(condition)
      
      if(length(row_indices)>=1){
        opt_data=train_pareto[row_indices,]
        ###########range the optimal options in training data
        opt_data_order=opt_data[order(-opt_data[,"es1"]),]
        
        keep=c("es1","es2")
        pareto_P0=t(as.data.frame(c(given_point[1],-Inf)))
        pareto_P1=t(as.data.frame(c(given_point[1],opt_data_order[1,"es2"])))
        pareto_P3=t(as.data.frame(c(opt_data_order[nrow(opt_data_order),"es1"],given_point[2])))
        pareto_P4=t(as.data.frame(c(-Inf,given_point[2])))
        colnames(pareto_P0)=keep
        colnames(pareto_P1)=keep
        colnames(pareto_P3)=keep
        colnames(pareto_P4)=keep
        
        
        ###########return the coordinates of the options which will construct the constrains
        pareto_P=as.data.frame(rbind(pareto_P0,pareto_P1,opt_data_order,pareto_P3,pareto_P4))
      }else{
        pareto_P0=t(as.data.frame(c(given_point[1],-Inf)))
        pareto_P1=t(as.data.frame(given_point))
        pareto_P4=t(as.data.frame(c(-Inf,given_point[2])))
        keep=c("es1","es2")
        colnames(pareto_P0)=keep
        colnames(pareto_P1)=keep
        colnames(pareto_P4)=keep
        pareto_P=as.data.frame(rbind(pareto_P0,pareto_P1,pareto_P4))
      }
      
      
    }
    
    if(crit=="all_pareto"){
      pareto_P=train_pareto
    }
  }
  
  if(pareto_P[1,"es2"]==-Inf){
    pareto_P[1,"es2"]=ref_min[2]
  }
  if(pareto_P[nrow(pareto_P),"es1"]==-Inf){
    pareto_P[nrow(pareto_P),"es1"]=ref_min[1]
  }

  
  return(pareto_P)
}
###########decide the coordinates of the options of EI
set_point_3points <- function(crit,vir_target,train_data,train_pareto,ref00,tar_point00,scale){
  
  if(crit=="3points"){
    ###########decide the current best trade-off option by maximizing cos_theta
    cos_theta=cos_theta_calculate(vir_point0=train_pareto,ref00,tar_point00,scale)
    opt_data=train_pareto[which(cos_theta[]==max(cos_theta[])),]
    
    ###########range the optimal options in training data according to cos_theta
    train_pareto_th=cbind(train_pareto,cos_theta)
    train_pareto_th_order=train_pareto_th[order(-train_pareto_th[,"cos_theta"]),]
    ###########choose three current best trade-off options
    point2=train_pareto_th_order[which(train_pareto_th_order[,"cos_theta"]==max(train_pareto_th_order[,"cos_theta"])),]
    train_pareto_th_order$abs=train_pareto_th_order[,"es1"]-point2[,"es1"]
    point3_c=train_pareto_th_order[which(train_pareto_th_order[,"abs"]<0),]
    point3=point3_c[1,]
    point1_c=train_pareto_th_order[which(train_pareto_th_order[,"abs"]>0),]
    point1=point1_c[1,]
    
    keep=c("es1","es2")
    pareto_P0=t(as.data.frame(c(point1[,"es1"],-Inf)))
    pareto_P1=point1[,keep]
    pareto_P2=point2[,keep]
    pareto_P3=point3[,keep]
    pareto_P4=t(as.data.frame(c(-Inf,point3[,"es2"])))
    
    colnames(pareto_P0)=keep
    colnames(pareto_P4)=keep
    ###########return the coordinates of the options which will construct the constrains
    pareto_P=as.data.frame(rbind(pareto_P0,pareto_P1,pareto_P2,pareto_P3,pareto_P4))
    if(is.na(pareto_P[2,1])){
      pareto_P[1,1]=pareto_P[3,1]
      pareto_P=pareto_P[-2,]
    }
    set_point_3points=t(as.data.frame(c(point1[,"es1"],point3[,"es2"])))
    colnames(set_point_3points)=c("es1","es2")
  }
  if(crit=="3points_dist"){
    ###########decide the current best trade-off option by minimizing the distance
    
    dist=t(as.data.frame((apply(train_pareto,1,fn_o_dist,vir_target))))
    colnames(dist)="distance"
    rownames(dist)=rep(1:nrow(dist))
    ###########choose three current best trade-off options
    sorted_indices <- sort(dist[,1], index.return = TRUE)$ix
    if(length(sorted_indices)>=3){
      opt_data=train_pareto[sorted_indices[1:3],]
    }else{
      opt_data=train_pareto[sorted_indices[1:length(sorted_indices)],]
    }
    ###########range the optimal options in training data
    opt_data_order=opt_data[order(-opt_data[,"es1"]),]
    
    keep=c("es1","es2")
    pareto_P0=t(as.data.frame(c(opt_data_order[1,"es1"],-Inf)))
    pareto_P4=t(as.data.frame(c(-Inf,opt_data_order[nrow(opt_data_order),"es2"])))
    colnames(pareto_P0)=keep
    colnames(pareto_P4)=keep
    
    
    ###########return the coordinates of the options which will construct the constrains
    pareto_P=as.data.frame(rbind(pareto_P0, opt_data_order,pareto_P4))
    
    set_point_3points=t(as.data.frame(c(opt_data_order[1,"es1"],opt_data_order[nrow(opt_data_order),"es2"])))
    colnames(set_point_3points)=c("es1","es2")
  }
  
  return(set_point_3points)
}

###########cos_theta calculation function
cos_theta_calculate=function(vir_point0,ref00,tar_point00,scale){
  
  ref0=as.data.frame(t(as.data.frame(ref00)))
  colnames(ref0)=colnames(vir_point0)
  tar_point0=as.data.frame(t(as.data.frame(tar_point00)))
  colnames(tar_point0)=colnames(vir_point0)
  
  all=rbind(ref0,tar_point0,vir_point0)
  
  #######scale
  if(scale==TRUE){
    all[,1]=scale(all[,1])
    all[,2]=scale(all[,2])
  }
  ref=unlist(all[1,])
  tar_point=unlist(all[2,])
  vir_point=all[-c(1,2),]
  
  a1=(ref[1]-vir_point[,1])*(ref[1]-tar_point[1])+(ref[2]-vir_point[,2])*(ref[2]-tar_point[2])
  a2=(((ref[1]-vir_point[,1])^2+(ref[2]-vir_point[,2])^2)*((ref[1]-tar_point[1])^2+(ref[2]-tar_point[2])^2))^0.5
  cos_theta=a1/a2
  # angle=as.data.frame(EHI_grid*cos_theta)
  return(cos_theta)
}

###############SMS##########################################################

################################CSMS
crit_CSMS <- function(x, model, ref_min, paretoFront, critcontrol=NULL, type="UK")
{
  if(paretoFront[1,"es2"]==-Inf){
    paretoFront[1,"es2"]=ref_min[2]
  }
  if(paretoFront[nrow(paretoFront),"es1"]==-Inf){
    paretoFront[nrow(paretoFront),"es1"]=ref_min[1]
  }
  #paretoFront
  n.obj <- length(model)
  d <- model[[1]]@d
  x.new <- matrix(x, 1, d)
  
  distp <- 0  # penalty if too close in the checkPredict sense
  
  if(checkPredict(x.new, model, type = type, distance = critcontrol$distance, threshold = critcontrol$threshold)){
    # return(0) Not compatible with penalty with SMS
    distp <- 1 # may be changed
  }#else{
  
  refPoint  <- critcontrol$refPoint
  currentHV <- critcontrol$currentHV
  epsilon   <- critcontrol$epsilon
  gain      <- critcontrol$gain
  nsteps.remaining <- critcontrol$nsteps.remaining
  
  if(is.null(paretoFront) || is.null(refPoint)) {
    observations <- Reduce(cbind, lapply(model, slot, "y"))
  }
  if(is.null(paretoFront)) paretoFront <- t(nondominated_points(t(observations)))
  if (is.null(refPoint)){
    if(is.null(critcontrol$extendper)) critcontrol$extendper <- 0.2
    # refPoint    <- matrix(apply(paretoFront, 2, max) + 1, 1, n.obj)
    PF_range <- apply(paretoFront, 2, range)
    refPoint <- matrix(PF_range[2,] + pmax(1, (PF_range[2,] - PF_range[1,]) * critcontrol$extendper), 1, n.obj)
    cat("No refPoint provided, ", signif(refPoint, 3), "used \n")
  }
  n.pareto <- nrow(paretoFront)
  
  if (is.null(currentHV)) currentHV <- dominated_hypervolume(points=t(paretoFront), ref=refPoint)
  if (is.null(nsteps.remaining)) nsteps.remaining <- 1
  # if (is.null(gain))  gain <- -qnorm( 0.5*(0.5^(1/n.obj)) )
  if (is.null(gain))  gain <- -qnorm( 0.5*(0.5^(1/n.obj)) )
  
  if (is.null(epsilon)) {
    if (n.pareto < 2){
      epsilon <- rep(0, n.obj)
    } else {
      spread <- apply(paretoFront,2,max) - apply(paretoFront,2,min)
      c <- 1 - (1 / (2^n.obj) )
      epsilon <- spread / (n.pareto + c * (nsteps.remaining-1))
    }
  }
  
  # mu    <- rep(NaN, n.obj)
  # sigma <- rep(NaN, n.obj)
  # for (i in 1:n.obj){    
  #   pred     <- predict(object=model[[i]], newdata=x.new, type=type, checkNames = FALSE, light.return = TRUE, cov.compute = FALSE)
  #   mu[i]    <- pred$mean
  #   sigma[i] <- pred$sd
  # }
  pred <- predict_kms(model, newdata=x.new, type=type, checkNames = FALSE, light.return = TRUE, cov.compute = FALSE)
  mu <- as.numeric(pred$mean)
  sigma <- as.numeric(pred$sd)
  ##################new addition  
  pre_mean=t(as.data.frame(mu))
  colnames(pre_mean)=c("mean1","mean2")
  pre_sd=t(as.data.frame(sigma))
  colnames(pre_sd)=c("sd1","sd2")
  pre=cbind(pre_mean,pre_sd)
  #####################  
  potSol <- mu + gain*sigma
  penalty <- distp
  for (j in 1:n.pareto){
    # assign penalty to all epsilon-dominated solutions
    if (min(paretoFront[j,] <= potSol + epsilon)){
      p <- -1 + prod(1 + pmax(potSol - paretoFront[j,], rep(0, n.obj)))
      penalty <- max(penalty, p)
    }
  }
  if (penalty == 0){
    # non epsilon-dominated solution
    potFront <- rbind(paretoFront, potSol)
    mypoints <- t(nondominated_points(t(potFront)))
    if (is.null(nrow(mypoints))) mypoints <- matrix(mypoints, 1, n.obj)
    myhv <- dominated_hypervolume(points=t(mypoints), ref=refPoint)
    f    <- currentHV - myhv
  } else{
    f <- penalty
  }
  #}
  penal=as.data.frame(-f)
  colnames(penal)="score"
  pre_fin=cbind(pre,penal)
  return(pre_fin)
}
fn_CSMSEGO_grid<-function(crit=NULL,ref_min,vir_target,train_data,test.grid,train_pareto,ref00,tar_point00,scale,opt_algri=crit_CSMS, model){
  pareto_P=fn_pareto_P(crit,vir_target,train_data,train_pareto,ref00,tar_point00,scale)
  CSMSEGO_grid <- apply(test.grid, 1, opt_algri, model,ref_min=ref_min,
                        paretoFront = pareto_P, critcontrol =  list(refPoint = ref00))
  
  CSMSEGO_grid=data.frame(matrix(unlist(CSMSEGO_grid), byrow = T, nrow = nrow(test.grid)))
  colnames(CSMSEGO_grid)=c("mean1","mean2","sd1","sd2","score")
  #CSMSEGO_grid=as.data.frame(CSMSEGO_grid)
  return(CSMSEGO_grid)
}
########
################################CSMS_error
crit_CSMS00 <- function(x, model, paretoFront=NULL,pareto_P=NULL, critcontrol=NULL, type="UK")
{
  n.obj <- length(model)
  d <- model[[1]]@d
  x.new <- matrix(x, 1, d)
  
  distp <- 0  # penalty if too close in the checkPredict sense
  
  if(checkPredict(x.new, model, type = type, distance = critcontrol$distance, threshold = critcontrol$threshold)){
    # return(0) Not compatible with penalty with SMS
    distp <- 1 # may be changed
  }#else{
  
  refPoint  <- critcontrol$refPoint
  currentHV <- critcontrol$currentHV
  epsilon   <- critcontrol$epsilon
  gain      <- critcontrol$gain
  nsteps.remaining <- critcontrol$nsteps.remaining
  
  if(is.null(paretoFront) || is.null(refPoint)) {
    observations <- Reduce(cbind, lapply(model, slot, "y"))
  }
  if(is.null(paretoFront)) paretoFront <- t(nondominated_points(t(observations)))
  if (is.null(refPoint)){
    if(is.null(critcontrol$extendper)) critcontrol$extendper <- 0.2
    # refPoint    <- matrix(apply(paretoFront, 2, max) + 1, 1, n.obj)
    PF_range <- apply(paretoFront, 2, range)
    refPoint <- matrix(PF_range[2,] + pmax(1, (PF_range[2,] - PF_range[1,]) * critcontrol$extendper), 1, n.obj)
    cat("No refPoint provided, ", signif(refPoint, 3), "used \n")
  }
  n.pareto <- nrow(paretoFront)
  
  if (is.null(currentHV)) currentHV <- dominated_hypervolume(points=t(paretoFront), ref=refPoint)
  if (is.null(nsteps.remaining)) nsteps.remaining <- 1
  # if (is.null(gain))  gain <- -qnorm( 0.5*(0.5^(1/n.obj)) )
  if (is.null(gain))  gain <- -qnorm( 0.5*(0.5^(1/n.obj)) )
  
  if (is.null(epsilon)) {
    if (n.pareto < 2){
      epsilon <- rep(0, n.obj)
    } else {
      spread <- apply(paretoFront,2,max) - apply(paretoFront,2,min)
      c <- 1 - (1 / (2^n.obj) )
      epsilon <- spread / (n.pareto + c * (nsteps.remaining-1))
    }
  }
  
  # mu    <- rep(NaN, n.obj)
  # sigma <- rep(NaN, n.obj)
  # for (i in 1:n.obj){    
  #   pred     <- predict(object=model[[i]], newdata=x.new, type=type, checkNames = FALSE, light.return = TRUE, cov.compute = FALSE)
  #   mu[i]    <- pred$mean
  #   sigma[i] <- pred$sd
  # }
  pred <- predict_kms(model, newdata=x.new, type=type, checkNames = FALSE, light.return = TRUE, cov.compute = FALSE)
  mu <- as.numeric(pred$mean)
  sigma <- as.numeric(pred$sd)
  ##################new addition  
  pre_mean=t(as.data.frame(mu))
  colnames(pre_mean)=c("mean1","mean2")
  pre_sd=t(as.data.frame(sigma))
  colnames(pre_sd)=c("sd1","sd2")
  pre=cbind(pre_mean,pre_sd)
  #####################  
  potSol <- mu + gain*sigma
  penalty <- distp
  for (j in 1:n.pareto){
    # assign penalty to all epsilon-dominated solutions
    if (min(paretoFront[j,] <= potSol + epsilon)){
      p <- -1 + prod(1 + pmax(potSol - paretoFront[j,], rep(0, n.obj)))
      penalty <- max(penalty, p)
    }
  }
  if (penalty == 0){
    # non epsilon-dominated solution
    potFront <- rbind(pareto_P, potSol)
    mypoints <- t(nondominated_points(t(potFront)))
    if (is.null(nrow(mypoints))) mypoints <- matrix(mypoints, 1, n.obj)
    myhv <- dominated_hypervolume(points=t(mypoints), ref=refPoint)
    f    <- currentHV - myhv
  } else{
    f <- penalty
  }
  #}
  penal=as.data.frame(-f)
  colnames(penal)="score"
  pre_fin=cbind(pre,penal)
  return(pre_fin)
}
fn_CSMSEGO_grid00<-function(crit=NULL,vir_target,train_data,test.grid,train_pareto,ref00,tar_point00,scale,opt_algri=crit_CSMS, model){
  pareto_P=fn_pareto_P(crit,vir_target,train_data,train_pareto,ref00,tar_point00,scale)
  CSMSEGO_grid <- apply(test.grid, 1, opt_algri, model,
                        paretoFront = train_pareto, pareto_P=pareto_P, critcontrol =  list(refPoint = ref00))
  
  CSMSEGO_grid=data.frame(matrix(unlist(CSMSEGO_grid), byrow = T, nrow = nrow(test.grid)))
  colnames(CSMSEGO_grid)=c("mean1","mean2","sd1","sd2","score")
  #CSMSEGO_grid=as.data.frame(CSMSEGO_grid)
  return(CSMSEGO_grid)
}
####################SMS
fn_crit_SMS=function (x, model, paretoFront = NULL, critcontrol = NULL, type = "UK") 
{
  n.obj <- length(model)
  d <- model[[1]]@d
  x.new <- matrix(x, 1, d)
  distp <- 0
  if (checkPredict(x.new, model, type = type, distance = critcontrol$distance, 
                   threshold = critcontrol$threshold)) {
    distp <- 1
  }
  refPoint <- critcontrol$refPoint
  currentHV <- critcontrol$currentHV
  epsilon <- critcontrol$epsilon
  gain <- critcontrol$gain
  nsteps.remaining <- critcontrol$nsteps.remaining
  if (is.null(paretoFront) || is.null(refPoint)) {
    observations <- Reduce(cbind, lapply(model, slot, "y"))
  }
  if (is.null(paretoFront)) 
    paretoFront <- t(nondominated_points(t(observations)))
  if (is.null(refPoint)) {
    if (is.null(critcontrol$extendper)) 
      critcontrol$extendper <- 0.2
    PF_range <- apply(paretoFront, 2, range)
    refPoint <- matrix(PF_range[2, ] + pmax(1, (PF_range[2, 
    ] - PF_range[1, ]) * critcontrol$extendper), 1, n.obj)
    cat("No refPoint provided, ", signif(refPoint, 3), "used \n")
  }
  n.pareto <- nrow(paretoFront)
  if (is.null(currentHV)) 
    currentHV <- dominated_hypervolume(points = t(paretoFront), 
                                       ref = refPoint)
  if (is.null(nsteps.remaining)) 
    nsteps.remaining <- 1
  if (is.null(gain)) 
    gain <- -qnorm(0.5 * (0.5^(1/n.obj)))
  if (is.null(epsilon)) {
    if (n.pareto < 2) {
      epsilon <- rep(0, n.obj)
    }
    else {
      spread <- apply(paretoFront, 2, max) - apply(paretoFront, 
                                                   2, min)
      c <- 1 - (1/(2^n.obj))
      epsilon <- spread/(n.pareto + c * (nsteps.remaining - 
                                           1))
    }
  }
  pred <- predict_kms(model, newdata = x.new, type = type, 
                      checkNames = FALSE, light.return = TRUE, cov.compute = FALSE)
  mu <- as.numeric(pred$mean)
  sigma <- as.numeric(pred$sd)
  ##################new addition  
  pre_mean=t(as.data.frame(mu))
  colnames(pre_mean)=c("mean1","mean2")
  pre_sd=t(as.data.frame(sigma))
  colnames(pre_sd)=c("sd1","sd2")
  pre=cbind(pre_mean,pre_sd)
  ###########
  potSol <- mu - gain * sigma
  penalty <- distp
  for (j in 1:n.pareto) {
    if (min(paretoFront[j, ] <= potSol + epsilon)) {
      p <- -1 + prod(1 + pmax(potSol - paretoFront[j, ], 
                              rep(0, n.obj)))
      penalty <- max(penalty, p)
    }
  }
  if (penalty == 0) {
    potFront <- rbind(paretoFront, potSol)
    mypoints <- t(nondominated_points(t(potFront)))
    if (is.null(nrow(mypoints))) 
      mypoints <- matrix(mypoints, 1, n.obj)
    myhv <- dominated_hypervolume(points = t(mypoints), ref = refPoint)
    f <- currentHV - myhv
  }
  else {
    f <- penalty
  }
  penal=as.data.frame(-f)
  colnames(penal)="score"
  pre_fin=cbind(pre,penal)
  return(pre_fin)
  
}
fn_SMS<-function(test.grid,train_pareto,ref00,tar_point00,scale,opt_algri=fn_crit_SMS, model){
  
  SMS <- apply(test.grid, 1, opt_algri, model,
               paretoFront = train_pareto, critcontrol =  list(refPoint = ref00))
  
  SMS=data.frame(matrix(unlist(SMS), byrow = T, nrow = nrow(test.grid)))
  colnames(SMS)=c("mean1","mean2","sd1","sd2","score")
  #CSMSEGO_grid=as.data.frame(CSMSEGO_grid)
  return(SMS)
}
#####################different epsilon SMS
fn.epsilon=function(paretoFront = NULL){
  ###########ordering the points in pareto front according to theta
  cos_theta_pareto=cos_theta_calculate(paretoFront,ref00,tar_point00,scale)
  train_pareto_th=cbind(paretoFront,cos_theta_pareto)
  train_pareto_th_order=train_pareto_th[order(-train_pareto_th[,"cos_theta_pareto"]),]
  ########calculate epsilon according to weight
  spread <- apply(paretoFront,2,max) - apply(paretoFront,2,min)
  n.pareto=nrow(paretoFront)
  epsilon <- spread / (n.pareto)
  num_epsilon=matrix(0:(n.pareto-1),nrow=n.pareto,ncol=1)
  epsilon_x1=spread[1] / (n.pareto-1*num_epsilon)
  epsilon_x2=spread[2] / (n.pareto-1*num_epsilon)
  train_pareto_th_all=cbind(train_pareto_th_order,epsilon_x1,epsilon_x2)
  train_pareto_th_all=train_pareto_th_all[order(train_pareto_th_all[,"es1"]),]
  epsilon.name=c("epsilon_x1","epsilon_x2")
  epsilon=train_pareto_th_all[,epsilon.name]
  return(epsilon)
}

fn_epsilon_CSMS=function (x, model, paretoFront = NULL, critcontrol = NULL, type = "UK") 
{
  n.obj <- length(model)
  d <- model[[1]]@d
  x.new <- matrix(x, 1, d)
  distp <- 0
  if (checkPredict(x.new, model, type = type, distance = critcontrol$distance, 
                   threshold = critcontrol$threshold)) {
    distp <- 1
  }
  refPoint <- critcontrol$refPoint
  currentHV <- critcontrol$currentHV
  #epsilon <- critcontrol$epsilon
  gain <- critcontrol$gain
  nsteps.remaining <- critcontrol$nsteps.remaining
  
  if (is.null(paretoFront) || is.null(refPoint)) {
    observations <- Reduce(cbind, lapply(model, slot, "y"))
  }
  if (is.null(paretoFront)) 
    paretoFront <- t(nondominated_points(t(observations)))
  if (is.null(refPoint)) {
    if (is.null(critcontrol$extendper)) 
      critcontrol$extendper <- 0.2
    PF_range <- apply(paretoFront, 2, range)
    refPoint <- matrix(PF_range[2, ] + pmax(1, (PF_range[2, 
    ] - PF_range[1, ]) * critcontrol$extendper), 1, n.obj)
    cat("No refPoint provided, ", signif(refPoint, 3), "used \n")
  }
  n.pareto <- nrow(paretoFront)
  if (is.null(currentHV)) 
    currentHV <- dominated_hypervolume(points = t(paretoFront), 
                                       ref = refPoint)
  if (is.null(nsteps.remaining)) 
    nsteps.remaining <- 1
  if (is.null(gain)) 
    gain <- -qnorm(0.5 * (0.5^(1/n.obj)))
  #if (is.null(epsilon)) {
  #  if (n.pareto < 2) {
  #    epsilon <- rep(0, n.obj)
  #  }
  #  else {
  #   spread <- apply(paretoFront, 2, max) - apply(paretoFront,  2, min)
  #    c <- 1 - (1/(2^n.obj))
  #   epsilon <- spread/(n.pareto + c * (nsteps.remaining -  1))
  # }
  #}
  pred <- predict_kms(model, newdata = x.new, type = type, 
                      checkNames = FALSE, light.return = TRUE, cov.compute = FALSE)
  mu <- as.numeric(pred$mean)
  sigma <- as.numeric(pred$sd)
  ##################new addition  
  pre_mean=t(as.data.frame(mu))
  colnames(pre_mean)=c("mean1","mean2")
  pre_sd=t(as.data.frame(sigma))
  colnames(pre_sd)=c("sd1","sd2")
  pre=cbind(pre_mean,pre_sd)
  epsilon=fn.epsilon(paretoFront)
  ###########
  potSol <- mu - gain * sigma
  penalty <- distp
  for (j in 1:n.pareto) {
    if (min(paretoFront[j, ] <= potSol + epsilon[j,])) {
      p <- -1 + prod(1 + pmax(potSol - paretoFront[j, ], 
                              rep(0, n.obj)))
      penalty <- max(penalty, p)
    }
  }
  if (penalty == 0) {
    potFront <- rbind(paretoFront, potSol)
    mypoints <- t(nondominated_points(t(potFront)))
    if (is.null(nrow(mypoints))) 
      mypoints <- matrix(mypoints, 1, n.obj)
    myhv <- dominated_hypervolume(points = t(mypoints), ref = refPoint)
    f <- currentHV - myhv
  }
  else {
    f <- penalty
  }
  penal=as.data.frame(-f)
  colnames(penal)="score"
  pre_fin=cbind(pre,penal)
  return(pre_fin)
  
}
fn_SMS_epsilon<-function(test.grid,train_pareto,ref00,tar_point00,scale,opt_algri=fn_epsilon_CSMS, model){
  
  SMS_epsilon <- apply(test.grid, 1, opt_algri, model,
                       paretoFront = train_pareto, critcontrol =  list(refPoint = ref00))
  
  SMS_epsilon=data.frame(matrix(unlist(SMS_epsilon), byrow = T, nrow = nrow(test.grid)))
  colnames(SMS_epsilon)=c("mean1","mean2","sd1","sd2","score")
  #CSMSEGO_grid=as.data.frame(CSMSEGO_grid)
  return(SMS_epsilon)
}

##########################SUR
fn_crit_SUR=function (x, model, paretoFront = NULL, critcontrol = NULL, type = "UK") 
{
  X.new <- matrix(x, nrow = 1, ncol = model[[1]]@d)
  if (checkPredict(X.new, model, type = type, distance = critcontrol$distance, 
                   threshold = critcontrol$threshold)) {
    crit <- -1
  }
  else {
    n.obj <- length(model)
    if (is.null(critcontrol$integration.points)) {
      d <- model[[1]]@d
      integration.param <- integration_design_optim(critcontrol, 
                                                    d, critcontrol$lower, critcontrol$upper, model = model)
      integration.points <- as.matrix(integration.param$integration.points)
      integration.weights <- integration.param$integration.weights
    }
    if (is.null(critcontrol$mn.X) || is.null(critcontrol$sn.X) || 
        is.null(critcontrol$precalc.data)) {
      precalc.data <- vector("list", n.obj)
      mn.X <- sn.X <- matrix(0, n.obj, nrow(integration.points))
      p.tst <- predict_kms(model, newdata = integration.points, 
                           checkNames = FALSE, type = type, cov.compute = FALSE, 
                           light.return = TRUE)
      mn.X <- p.tst$mean
      sn.X <- p.tst$sd
      wrapped_precomputeUpdateData <- function(model, integration.points) {
        if (!is(model, "km")) 
          return(NULL)
        else return(precomputeUpdateData(model, integration.points))
      }
      precalc.data <- lapply(model, FUN = wrapped_precomputeUpdateData, 
                             integration.points = integration.points)
    }
    else {
      integration.points <- critcontrol$integration.points
      integration.weights <- critcontrol$integration.weights
      mn.X <- critcontrol$mn.X
      sn.X <- critcontrol$sn.X
      precalc.data <- critcontrol$precalc.data
    }
    n.integration.points <- nrow(integration.points)
    if (is.null(integration.weights)) {
      integration.weights <- rep(1/n.integration.points, 
                                 n.integration.points)
    }
    if (is.null(paretoFront)) {
      observations <- Reduce(cbind, lapply(model, slot, 
                                           "y"))
      paretoFront <- t(nondominated_points(t(observations)))
    }
    I <- which(sn.X[1, ] < sqrt(model[[1]]@covariance@sd2)/10000)
    if (length(I) > 0) {
      integration.points <- integration.points[-I, , drop = FALSE]
      mn.X <- mn.X[, -I, drop = FALSE]
      sn.X <- sn.X[, -I, drop = FALSE]
      for (i in 1:n.obj) {
        precalc.data[[i]]$Kinv.c.olddata <- precalc.data[[i]]$Kinv.c.olddata[, 
                                                                             -I, drop = FALSE]
        precalc.data[[i]]$first.member <- precalc.data[[i]]$first.member[-I]
      }
      if (length(integration.weights) > 1) {
        integration.weights <- integration.weights[-I]
      }
    }
    if (length(integration.weights) < 1) {
      stop("Unable to compute the SUR criterion (crit_SUR): all the integration points have a too small kriging variance")
      pred.xnew <- NULL
      crit <- 0
    }
    else {
      n.integration.points <- nrow(integration.points)
      integration.weights <- integration.weights/sum(integration.weights)
      if (n.obj == 2) {
        if (is.unsorted(paretoFront[, 1])) 
          paretoFront <- paretoFront[sort(paretoFront[, 
                                                      1], index.return = TRUE)[[2]], , drop = FALSE]
      }
      else {
        if (is.unsorted(-paretoFront[, 3])) 
          paretoFront <- paretoFront[sort(paretoFront[, 
                                                      3], decreasing = TRUE, index.return = TRUE)[[2]], 
                                     , drop = FALSE]
      }
      n.pareto <- nrow(paretoFront)
      mn.xnew <- sn.xnew <- rep(0, n.obj)
      kn <- rho <- eta <- r <- matrix(0, nrow = n.obj, 
                                      ncol = n.integration.points)
      pred.xnew <- vector("list", n.obj)
      for (i in 1:n.obj) {
        krig <- predict(object = model[[i]], newdata = data.frame(x = (X.new)), 
                        type = type, se.compute = TRUE, cov.compute = FALSE, 
                        checkNames = FALSE)
        mn.xnew[i] <- krig$mean
        sn.xnew[i] <- krig$sd
        
        ##################new addition  
        pre_mean=t(as.data.frame(mn.xnew))
        colnames(pre_mean)=c("mean1","mean2")
        pre_sd=t(as.data.frame(sn.xnew))
        colnames(pre_sd)=c("sd1","sd2")
        pre=cbind(pre_mean,pre_sd)
        ########################

        pred.xnew[[i]] <- krig
        if (krig$sd != 0) {
          kn[i, ] = computeQuickKrigcov2(model[[i]], 
                                         integration.points = integration.points, 
                                         X.new = (X.new), precalc.data = precalc.data[[i]], 
                                         F.newdata = krig$F.newdata, c.newdata = krig$c)
          rho[i, ] <- kn[i, ]/(sn.xnew[i] * sn.X[i, ])
          denom <- sqrt(pmax(1e-12, sn.xnew[i]^2 + sn.X[i, 
          ]^2 - 2 * kn[i, ]))
          eta[i, ] <- (mn.xnew[i] - mn.X[i, ])/denom
          r[i, ] <- (kn[i, ] - sn.xnew[i]^2)/sn.xnew[i]/denom
        }
      }
      phi.eta <- pnorm(eta)
      objx.bar <- (paretoFront[, 1] - mn.xnew[1])/sn.xnew[1]
      objx.tilde <- (matrix(rep(paretoFront[, 1], n.integration.points), 
                            ncol = n.integration.points) - t(matrix(rep(mn.X[1, 
                            ], n.pareto), ncol = n.pareto)))/t(matrix(rep(sn.X[1, 
                            ], n.pareto), ncol = n.pareto))
      phi.x.bar <- pnorm(objx.bar)
      phi.x.tilde <- pnorm(objx.tilde)
      phi.eta.x <- pnorm(eta[1, ])
      phi2.x.x <- phi2.x.eta <- matrix(0, (n.pareto), n.integration.points)
      for (i in 1:(n.pareto)) {
        phi2.x.x[i, ] <- pbivnorm(rep(objx.bar[i], n.integration.points), 
                                  objx.tilde[i, ], rho[1, ])
        phi2.x.eta[i, ] <- pbivnorm(rep(objx.bar[i], 
                                        n.integration.points), eta[1, ], r[1, ])
      }
      if (sn.xnew[2] != 0) {
        objy.bar <- (paretoFront[, 2] - mn.xnew[2])/sn.xnew[2]
        objy.tilde <- (matrix(rep(paretoFront[, 2], n.integration.points), 
                              ncol = n.integration.points) - t(matrix(rep(mn.X[2, 
                              ], n.pareto), ncol = n.pareto)))/t(matrix(rep(sn.X[2, 
                              ], n.pareto), ncol = n.pareto))
        phi.y.bar <- pnorm(objy.bar)
        phi.y.tilde <- pnorm(objy.tilde)
        phi.eta.y <- pnorm(eta[2, ])
        phi2.y.y <- phi2.y.eta <- matrix(0, (n.pareto), 
                                         n.integration.points)
        for (i in 1:(n.pareto)) {
          phi2.y.y[i, ] <- pbivnorm(rep(objy.bar[i], 
                                        n.integration.points), objy.tilde[i, ], rho[2, 
                                        ])
          phi2.y.eta[i, ] <- pbivnorm(rep(objy.bar[i], 
                                          n.integration.points), eta[2, ], r[2, ])
        }
      }
      else {
        phi.y.bar <- as.numeric(paretoFront[, 2] > mn.xnew[2])
        phi.y.tilde <- 1 * (matrix(rep(paretoFront[, 
                                                   2], n.integration.points), ncol = n.integration.points) > 
                              t(matrix(rep(mn.X[2, ], n.pareto), ncol = n.pareto)))
        phi.eta.y <- as.numeric(mn.xnew[2] > mn.X[2, 
        ])
        phi2.y.y <- phi2.y.eta <- matrix(0, (n.pareto), 
                                         n.integration.points)
        for (i in 1:(n.pareto)) {
          phi2.y.y[i, ] <- phi.y.bar[i] * phi.y.tilde[i, 
          ]
          phi2.y.eta[i, ] <- phi.y.bar[i] * phi.eta.y
        }
      }
      if (n.obj == 3) {
        if (sn.xnew[3] != 0) {
          objz.bar <- (paretoFront[, 3] - mn.xnew[3])/sn.xnew[3]
          objz.tilde <- (matrix(rep(paretoFront[, 3], 
                                    n.integration.points), ncol = n.integration.points) - 
                           t(matrix(rep(mn.X[3, ], n.pareto), ncol = n.pareto)))/t(matrix(rep(sn.X[3, 
                           ], n.pareto), ncol = n.pareto))
          phi.z.bar <- pnorm(objz.bar)
          phi.z.tilde <- pnorm(objz.tilde)
          phi.eta.z <- pnorm(eta[3, ])
          phi2.z.z <- phi2.z.eta <- matrix(0, (n.pareto), 
                                           n.integration.points)
          for (i in 1:(n.pareto)) {
            phi2.z.z[i, ] <- pbivnorm(rep(objz.bar[i], 
                                          n.integration.points), objz.tilde[i, ], 
                                      rho[3, ])
            phi2.z.eta[i, ] <- pbivnorm(rep(objz.bar[i], 
                                            n.integration.points), eta[3, ], r[3, ])
          }
        }
        else {
          phi.z.bar <- as.numeric(paretoFront[, 3] > 
                                    mn.xnew[3])
          phi.z.tilde <- 1 * (matrix(rep(paretoFront[, 
                                                     3], n.integration.points), ncol = n.integration.points) > 
                                t(matrix(rep(mn.X[3, ], n.pareto), ncol = n.pareto)))
          phi.eta.z <- as.numeric(mn.xnew[3] > mn.X[3, 
          ])
          phi2.z.z <- phi2.z.eta <- matrix(0, (n.pareto), 
                                           n.integration.points)
          for (i in 1:(n.pareto)) {
            phi2.z.z[i, ] <- phi.z.bar[i] * phi.z.tilde[i, 
            ]
            phi2.z.eta[i, ] <- phi.z.bar[i] * phi.eta.z
          }
        }
      }
      if (n.obj == 2) {
        res <- EEV.2D.computation(phi.x.bar, phi.x.tilde, 
                                  phi.eta.x, phi.y.bar, phi.y.tilde, phi.eta.y, 
                                  phi2.x.x, phi2.x.eta, phi2.y.y, phi2.y.eta)
        piold <- colSums(res[[2]])
        pinew <- colSums(res[[3]])
        pinew[is.na(pinew)] <- piold[is.na(pinew)]
        newpn <- sum(pinew * integration.weights)
        oldpn <- sum(piold * integration.weights)
        crit <- oldpn - newpn
      }
      else {
        pn <- pold <- rep(0, n.integration.points)
        for (i in 1:(n.pareto)) {
          nondominated.sub <- which(!is_dominated((nondominated_points(t(paretoFront[(i:n.pareto), 
                                                                                     1:2, drop = FALSE])))))
          nondominated.sub <- i - 1 + nondominated.sub[sort(paretoFront[(i:n.pareto)[nondominated.sub], 
                                                                        1], index.return = TRUE)[[2]]]
          n.pareto.sub <- length(nondominated.sub)
          res <- EEV.2D.computation(phi.x.bar[nondominated.sub], 
                                    phi.x.tilde[nondominated.sub, , drop = FALSE], 
                                    phi.eta.x, phi.y.bar[nondominated.sub], phi.y.tilde[nondominated.sub, 
                                                                                        , drop = FALSE], phi.eta.y, phi2.x.x[nondominated.sub, 
                                                                                                                             , drop = FALSE], phi2.x.eta[nondominated.sub, 
                                                                                                                                                         , drop = FALSE], phi2.y.y[nondominated.sub, 
                                                                                                                                                                                   , drop = FALSE], phi2.y.eta[nondominated.sub, 
                                                                                                                                                                                                               , drop = FALSE])
          pijold <- colSums(res[[2]])
          pijnew <- colSums(res[[3]])
          if (i == 1) {
            p.joint.1 <- 1 - phi2.z.z[i, ] - phi.eta.z + 
              phi2.z.eta[i, ]
            p.joint.2 <- phi2.z.z[i, ] - phi.z.tilde[i, 
            ] + phi.eta.z - phi2.z.eta[i, ]
          }
          else {
            p.joint.1 <- phi2.z.z[i - 1, ] - phi2.z.z[i, 
            ] - phi2.z.eta[i - 1, ] + phi2.z.eta[i, 
            ]
            p.joint.2 <- phi.z.tilde[i - 1, ] - phi.z.tilde[i, 
            ] - phi2.z.z[i - 1, ] + phi2.z.z[i, ] + 
              phi2.z.eta[i - 1, ] - phi2.z.eta[i, ]
          }
          pold <- pold + pijold * (p.joint.1 + p.joint.2)
          pn <- pn + pijnew * (p.joint.1) + pijold * 
            (p.joint.2)
        }
        pold <- pold + phi.z.tilde[n.pareto, ]
        pa <- phi.z.tilde[n.pareto, ] - phi2.z.z[n.pareto, 
        ] + phi2.z.eta[n.pareto, ]
        px <- phi.eta.x
        py <- phi.eta.y
        pb <- (phi2.z.z[n.pareto, ] - phi2.z.eta[n.pareto, 
        ]) * (px + py - px * py)
        pn <- pn + (pa + pb)
        pn[is.na(pn)] <- pold[is.na(pn)]
        newpn <- sum(pn * integration.weights)
        oldpn <- sum(pold * integration.weights)
        crit <- oldpn - newpn
      }
    }
    if (crit < 1e-32) 
      crit <- prob.of.non.domination(paretoFront = paretoFront, 
                                     model, X.new, predictions = pred.xnew) - 1
    
    crit_result=cbind(pre,crit)
    return(crit_result)
  }
}
fn_SUR<-function(test.grid,train_pareto,ref00,tar_point00,scale,opt_algri=fn_crit_SUR, model){
  
  SUR <- apply(test.grid, 1, opt_algri, model,
               paretoFront = train_pareto, critcontrol =  list(refPoint = ref00))
  
  SUR=data.frame(matrix(unlist(SUR), byrow = T, nrow = nrow(test.grid)))
  colnames(SUR)=c("mean1","mean2","sd1","sd2","score")
  #CSMSEGO_grid=as.data.frame(CSMSEGO_grid)
  return(SUR)
}
#########################EMI
fn_crit_EMI=function (x, model, paretoFront = NULL, critcontrol = list(nb.samp = 50, 
                                                           seed = 42), type = "UK") 
{
  nobj <- length(model)
  if (nobj < 2) {
    cat("Incorrect Number of objectives \n")
    return(NA)
  }
  if (is.null(critcontrol$nb.samp)) 
    critcontrol$nb.samp <- 50
  if (is.null(critcontrol$seed)) 
    critcontrol$seed <- 42

  if (nobj == 2) 
    if (max(abs(model[[1]]@y)) <= 2 && max(abs(model[[2]]@y)) <= 
        2) 
      return(EMI_2d(x = x, model = model, critcontrol = critcontrol, 
                    type = type, paretoFront = paretoFront))
  if (is.null(critcontrol)) {
    critcontrol <- list()
  }
  critcontrol$type <- "maximin"
  
 
  crit_re=SAA_mEI(x = x, model = model, critcontrol = critcontrol, 
                 type = "UK", paretoFront = paretoFront)
 # crit_result=cbind(pre,crit_re)
  return(crit_re)
}
fn_EMI<-function(test.grid,train_pareto,ref00,tar_point00,scale,opt_algri=fn_crit_EMI, model){
  
  EMI <- apply(test.grid, 1, opt_algri, model,
               paretoFront = train_pareto, critcontrol =  list(refPoint = ref00))
  
  EMI=data.frame(matrix(unlist(EMI), byrow = T, nrow = nrow(test.grid)))
  #colnames(EMI)=c("mean1","mean2","sd1","sd2","score")
  #CSMSEGO_grid=as.data.frame(CSMSEGO_grid)
  return(EMI)
}
SAA_mEI <- function(x, model, critcontrol, type = "UK", paretoFront = NULL) {
  # x 可以是向量或矩阵
  if (is.null(dim(x))) x <- matrix(x, nrow = 1)
  
  nobj <- length(model)
  npts <- nrow(x)
  ei_mat <- matrix(NA, nrow = npts, ncol = 1)
  
  # 多目标 EI（基于 Monte Carlo 模拟）
  nb.samp <- critcontrol$nb.samp %||% 50  # 使用 critcontrol$nb.samp，默认50
  set.seed(critcontrol$seed %||% 42)
  
  for (i in 1:npts) {
    xnew <- matrix(x[i, ], nrow = 1)
    samples <- matrix(NA, nrow = nb.samp, ncol = nobj)
    
    for (j in 1:nobj) {
      pred <- predict(model[[j]], newdata = xnew, type = type)
      samples[, j] <- rnorm(nb.samp, mean = pred$mean, sd = pred$sd)
    }
    
    # 比较每个 sample 是否支配当前 Pareto Front（这里只是近似 EI）
    if (!is.null(paretoFront)) {
      dominated <- apply(samples, 1, function(samp) {
        any(apply(paretoFront, 1, function(p) all(p <= samp) && any(p < samp)))
      })
      ei_mat[i] <- mean(!dominated)
    } else {
      ei_mat[i] <- mean(apply(samples, 1, function(s) sum(s)))
    }
  }
  
  return(ei_mat)
}
`%||%` <- function(a, b) {
  if (!is.null(a)) a else b
}
#########################CEMI

fn_CEMI<-function(crit=NULL,test.grid,train_pareto,ref00,tar_point00,scale,opt_algri=fn_crit_CEMI, model){
  pareto_P=fn_pareto_P(crit,vir_target,train_data,train_pareto,ref00,tar_point00,scale)
  EMI <- apply(test.grid, 1, opt_algri, model,
               paretoFront =  pareto_P, critcontrol =  list(refPoint = ref00))
  
  EMI=data.frame(matrix(unlist(EMI), byrow = T, nrow = nrow(test.grid)))
  #colnames(EMI)=c("mean1","mean2","sd1","sd2","score")
  #CSMSEGO_grid=as.data.frame(CSMSEGO_grid)
  return(EMI)
}

fn_crit_CEMI=function (x, model, paretoFront = NULL, critcontrol = list(nb.samp = 50, 
                                                           seed = 42), type = "UK") 
{
  nobj <- length(model)
  if (nobj < 2) {
    cat("Incorrect Number of objectives \n")
    return(NA)
  }
  if (is.null(critcontrol$nb.samp)) 
    critcontrol$nb.samp <- 50
  if (is.null(critcontrol$seed)) 
    critcontrol$seed <- 42
  if (nobj == 2) 
    if (max(abs(model[[1]]@y)) <= 2 && max(abs(model[[2]]@y)) <= 
        2) 
      return(EMI_2d(x = x, model = model, critcontrol = critcontrol, 
                    type = type, paretoFront = paretoFront))
  if (is.null(critcontrol)) {
    critcontrol <- list()
  }
  critcontrol$type <- "maximin"
  return(SAA_mEI(x = x, model = model, critcontrol = critcontrol, 
                 type = "UK", paretoFront = paretoFront))
}
####################psi function
fn_psi_ab=function(a_val,b_val,mean_val,sd_val){
  nor_val = (b_val-mean_val )/sd_val
  z = nor_val
  psi_ab= sd_val*dnorm(z)+(a_val-mean_val)*pnorm(z)   
  return(psi_ab)
}

###############Expected Hypervolume Improvement(EHVI)##########################################################

#########Expected Hypervolume Improvement(EHVI)
fn_EHVI=function(ref_poi=ref00,pareto_Data=pareto_Data,model.result){
  
  pareto_Data=pareto_Data[order(-pareto_Data[,"es1"]),]
  y_0=t(as.data.frame(c(ref_poi[1],-Inf)))
  colnames(y_0)=c("es1","es2")
  y_i1=t(as.data.frame(c(-Inf,ref_poi[2])))
  colnames(y_i1)=c("es1","es2")   
  pareto_Data=rbind(y_0,pareto_Data,y_i1)
  i_max=nrow(pareto_Data)-1
  EHVI=rep()
  
  for (i in 1:i_max){
    z_1=(pareto_Data[i+1,"es1"]-model.result["mean1"] )/model.result["sd1"]
    iterm1=(pareto_Data[i,"es1"]-pareto_Data[i+1,"es1"])*pnorm(z_1)*fn_psi_ab(pareto_Data[i+1,"es2"],pareto_Data[i+1,"es2"],model.result["mean2"],model.result["sd2"])
    iterm2=(fn_psi_ab(pareto_Data[i,"es1"],pareto_Data[i,"es1"],model.result["mean1"],model.result["sd1"])-fn_psi_ab(pareto_Data[i,"es1"],pareto_Data[i+1,"es1"],model.result["mean1"],model.result["sd1"]))*fn_psi_ab(pareto_Data[i+1,"es2"],pareto_Data[i+1,"es2"],model.result["mean2"],model.result["sd2"])
    EHVI[i]=iterm1+iterm2
  }
  EHVI=as.data.frame(EHVI)
  EHVI=na.omit(EHVI)
  EHVI_sum=sum(EHVI)
  return(EHVI_sum)
}
#########EHVI evaluation function
fn_EHI=function(test.grid,train_pareto,ref_poi=ref00,tar_point00,scale,opt_algri=fn_EHVI,model){
  #km prediction
  pred <- predict_kms(model, newdata=test.grid, type="UK", checkNames = FALSE, light.return = TRUE, cov.compute = FALSE)
  
  mu <- as.data.frame(t(as.data.frame(pred$mean)))
  colnames(mu)=c("mean1","mean2")
  
  sigma <- as.data.frame(t(as.data.frame(pred$sd)))
  colnames(sigma)=c("sd1","sd2")
  ##################new addition 
  model.result=cbind(mu,sigma)
  #########EHVI evaluation
  EHI_result=apply(model.result,1,opt_algri,ref_poi=ref_poi,pareto_Data=train_pareto)
  EHI_result=as.data.frame(EHI_result)
  EHI=cbind(model.result,EHI_result)
  colnames(EHI)=c("mean1","mean2","sd1","sd2","score")
  return(EHI)
}

################################################################################################################


####################################constrained EHVI###########################################################
fn_1P_EHI=function(crit,vir_target,test.grid,train_data,train_pareto,ref00,tar_point00,scale,opt_algri=fn_EHVI,model){
  pred <- predict_kms(model, newdata=test.grid, type="UK", checkNames = FALSE, light.return = TRUE, cov.compute = FALSE)
  mu <- as.data.frame(t(as.data.frame(pred$mean)))
  sigma <- as.data.frame(t(as.data.frame(pred$sd)))
  ##################new addition  
  #pre_mean=t(as.data.frame(mu))
  colnames(mu)=c("mean1","mean2")
  # pre_sd=t(as.data.frame(sigma))
  colnames(sigma )=c("sd1","sd2")
  model.result=cbind(mu,sigma)
  pareto_P=fn_pareto_P(crit,vir_target,train_data,train_pareto,ref00,tar_point00,scale)
  
  EHI_result=apply(model.result,1,opt_algri,ref_poi=ref00,pareto_Data=pareto_P)
  EHI_result=as.data.frame(EHI_result)
  EHI=cbind(model.result,EHI_result)
  colnames(EHI)=c("mean1","mean2","sd1","sd2","score")
  return(EHI)
}
####################constrained EHVI
fn_P_EHVI=function(ref_poi=ref00,pareto_Data=pareto_Data,model.result){
  
  pareto_Data=pareto_Data[order(-pareto_Data[,"es1"]),]
  #if(nrow(pareto_Data)==3){
  #  pareto_Data=pareto_Data
 # }else{
  pareto_Data=pareto_Data[c(1,3,5),]#}
  colnames(pareto_Data)=c("es1","es2")  
  i_max=nrow(pareto_Data)-1
  EHVI=rep()
  
  for (i in 1:i_max){
    z_1=(pareto_Data[i+1,"es1"]-model.result["mean1"] )/model.result["sd1"]
    iterm1=(pareto_Data[i,"es1"]-pareto_Data[i+1,"es1"])*pnorm(z_1)*fn_psi_ab(pareto_Data[i+1,"es2"],pareto_Data[i+1,"es2"],model.result["mean2"],model.result["sd2"])
    iterm2=(fn_psi_ab(pareto_Data[i,"es1"],pareto_Data[i,"es1"],model.result["mean1"],model.result["sd1"])-fn_psi_ab(pareto_Data[i,"es1"],pareto_Data[i+1,"es1"],model.result["mean1"],model.result["sd1"]))*fn_psi_ab(pareto_Data[i+1,"es2"],pareto_Data[i+1,"es2"],model.result["mean2"],model.result["sd2"])
    EHVI[i]=iterm1+iterm2
  }
  EHVI=as.data.frame(EHVI)
  EHVI=na.omit(EHVI)
  EHVI_sum=sum(EHVI)
  return(EHVI_sum)
}

#########constrained EHVI evaluation function
fn_3points_EHI=function(crit,vir_target,test.grid,train_data,train_pareto,ref00,tar_point00,scale,opt_algri=fn_P_EHVI,model){
  pred <- predict_kms(model, newdata=test.grid,type="UK", checkNames = FALSE, light.return = TRUE, cov.compute = FALSE)
  mu <- as.data.frame(t(as.data.frame(pred$mean))) 
  colnames(mu)=c("mean1","mean2")
  sigma <- as.data.frame(t(as.data.frame(pred$sd)))
  colnames(sigma )=c("sd1","sd2")
  
  ##################new addition  
  model.result=cbind(mu,sigma)
  #the coordinates of the options which will construct the constrains
  pareto_P=fn_pareto_P(crit,vir_target,train_data,train_pareto,ref00,tar_point00,scale)
  
  EHI_result=apply(model.result,1,opt_algri,ref_poi=ref00,pareto_Data=pareto_P)
  EHI_result=as.data.frame(EHI_result)
  EHI=cbind(model.result,EHI_result)
  colnames(EHI)=c("mean1","mean2","sd1","sd2","score")
  return(EHI)
}
################################################################################################################


###############minimize the distance to the ideal point##########################################################

#########distance calculation function
fn_dist=function(model.result.option,tar_point00){
  dist_cal=((model.result.option["mean1"]-tar_point00[1])^2+(model.result.option["mean2"]-tar_point00[2])^2)^0.5
  return(dist_cal)
}
fn_o_dist=function(train_data,tar_point00){
  dist_cal=((train_data["es1"]-tar_point00[1])^2+(train_data["es2"]-tar_point00[2])^2)^0.5
  return(dist_cal)
}
#########mean_dist evaluation function
fn_mean_dist=function(vir_target,test.grid,train_data,train_pareto,ref00,tar_point00,scale,opt_algri=fn_dist,model){
  pred <- predict_kms(model, newdata=test.grid, type="UK", checkNames = FALSE, light.return = TRUE, cov.compute = FALSE)
  
  mu <- as.data.frame(t(as.data.frame(pred$mean)))
  sigma <- as.data.frame(t(as.data.frame(pred$sd)))
  
  colnames(mu)=c("mean1","mean2")
  colnames(sigma )=c("sd1","sd2")
  ##################new addition  
  model.result=cbind(mu,sigma)
  
  dist_result=unlist(apply(model.result,1,opt_algri,tar_point00))
  dist_result=as.data.frame(-dist_result)
  dist=cbind(model.result,dist_result)
  colnames(dist)=c("mean1","mean2","sd1","sd2","score")
  return(dist)
}
##############EI*EI set_point
fn.ego.ei = function(set_point,model.result)
{
 # ego = (min(train_data[,"es"]) - model.result[,"mean"])/model.result[,"sd"]
  ego = (set_point[,"es"] - model.result[,"mean"])/model.result[,"sd"]
  z = ego
  ei.ego = model.result[,"sd"]*z*pnorm(z) + model.result[,"sd"]*dnorm(z) 
  ei.ego = data.frame(ei.ego)
  return (ei.ego)
}
fn_EIEI=function(vir_target,test.grid,set_point,train_pareto,ref00,tar_point00,scale,opt_algri=fn.ego.ei,model){
  pred <- predict_kms(model, newdata=test.grid, type="UK", checkNames = FALSE, light.return = TRUE, cov.compute = FALSE)
  
  mu <- as.data.frame(t(as.data.frame(pred$mean)))
  sigma <- as.data.frame(t(as.data.frame(pred$sd)))
  
  colnames(mu)=c("mean1","mean2")
  colnames(sigma )=c("sd1","sd2")

  ##################new addition  
  model.result=cbind(mu,sigma)
  model.result1=model.result[,c("mean1","sd1")]
  colnames(model.result1)=c("mean","sd")

  set_point1=as.data.frame(set_point[,"es1"])
  colnames(set_point1)="es"
  EI1= fn.ego.ei(set_point1,model.result1)
  EI1[which(EI1[,]=="NaN"),]=0
  EI1=as.data.frame(EI1)
  colnames(EI1)="EI1"
  
  model.result2=model.result[,c("mean2","sd2")]
  colnames(model.result2)=c("mean","sd")
  set_point2=as.data.frame(set_point[,"es2"])
  colnames(set_point2)="es"
  EI2= fn.ego.ei(set_point2,model.result2)
  EI2[which(EI2[,]=="NaN"),]=0
  EI2=as.data.frame(EI2)
  colnames(EI2)="EI2"
  
  EI_EI=EI1*EI2
  
  ###############uncertainty
  
  ub1=model.result[,"mean1"] + beta * model.result[,"sd1"] 
  lb1=model.result[,"mean1"] - beta * model.result[,"sd1"] 
  
  
  ub2=model.result[,"mean2"] + beta * model.result[,"sd2"] 
  lb2=model.result[,"mean2"] - beta * model.result[,"sd2"] 
  
  uncertainty=(ub1 - lb1) * (ub2 - lb2)
  
  EI_result=cbind(model.result,EI_EI,EI1,EI2,uncertainty)
   colnames(EI_result)=c("mean1","mean2","sd1","sd2","score","EI1","EI2","uncertainty")
  return(EI_result)
}
################################################################################################################

#############USeMO##################


fn_USeMO=function(EI_result,test.grid){
  
  # 提取要优化的 EI1 和 EI2
  acq_matrix <- EI_result[, c("EI1", "EI2")]
  
  obj_fn <- function(x) {
    if (is.null(dim(x))) {
      # x 是向量
      idx <- round(x[1])
      idx <- min(max(idx, 1), nrow(acq_matrix))
     # print("1")
     # print(matrix(-as.numeric(acq_matrix[idx, ]), nrow = 1))
      return(matrix(-as.numeric(acq_matrix[idx, ]), nrow = 1))
    } else {
      # x 是矩阵（每行一个个体）
      res <- apply(x, 1, function(row) {
        idx <- round(row[1])
        idx <- min(max(idx, 1), nrow(acq_matrix))
      #  print("2")
        return(-as.numeric(acq_matrix[idx, ]))  # 返回长度为 2 的向量
      })
      return(t(res))  # 转置回来，得到 n行 × 2列
    }
  }
  nsga_result <- nsga2(
    fn = obj_fn,
    idim = 1,
    odim = 2,
    lower.bounds = 1,
    upper.bounds = as.numeric(nrow(acq_matrix)),
    popsize = 20,
    generations = 50,
    vectorized = F  # ← 关键区别在这里
  )
test_EI_result=cbind(test.grid,EI_result)
selected_idx <- unique(round(nsga_result$par))
test_EI_result[,"uncertainty"] <- ifelse(1:nrow(test_EI_result) %in% selected_idx, test_EI_result[,"uncertainty"], 0)
test_EI_result_score=as.data.frame(test_EI_result[,"uncertainty"]) 
colnames(test_EI_result_score)=c("score")
USeMO=cbind(EI_result[,c("mean1","mean2","sd1","sd2")],test_EI_result_score)

return(USeMO)
}
#selected_points <- test_EI_result[selected_idx, , drop = FALSE]

#top_index <- order(selected_points[,"uncertainty"], decreasing = TRUE)[1:batch_size]

#################qEHVI+
compute_qehvi_optim <- function(x, models, ref_point, pareto_front, q, n_mc = 128, d) {
  X_cand <- matrix(x, nrow = q, byrow = TRUE)
  n_obj <- length(models)
  
  samples <- lapply(1:n_obj, function(i) {
    pred <- predict(models[[i]], newdata = X_cand, type = "UK", checkNames = FALSE)
    matrix(rep(pred$mean, each = n_mc), nrow = q) +
      matrix(rnorm(q * n_mc, sd = sqrt(models[[i]]@covariance@sd2)), nrow = q)
  })
  
  samples <- array(unlist(samples), dim = c(q, n_mc, n_obj))
  
  hv_improvements <- sapply(1:n_mc, function(i) {
    sampled_obj <- samples[, i, ]
    front_aug <- cbind(pareto_front, t(sampled_obj))
    emoa::dominated_hypervolume(front_aug, ref = ref_point) - 
      emoa::dominated_hypervolume(pareto_front, ref = ref_point)
  })
  
  return(mean(hv_improvements))
}
# ---- 随机组合并计算 qEHVI ----
qehvi=function(n_points,q=3,n_samples=10,X_cand,ref_point,models,pf,d)
{
combinations <- combn(n_points, q, simplify = FALSE)
selected_combinations <- sample(combinations, n_samples)

X_list <- lapply(selected_combinations, function(idx) X_cand[idx, , drop = FALSE])

qehvi_values <- sapply(X_list, function(x_batch) {
  compute_qehvi_optim(x = x_batch, models = models, ref_point = ref_point,
                      pareto_front = pf, q = q, n_mc, d = d)
})

best_idx <- which.max(qehvi_values)
best_combination <- selected_combinations[[best_idx]]
best_value <- qehvi_values[best_idx]
return(best_combination)
}
############constrained EHVI


# ---- 使用 compute_qehvi_optim ----
compute_qehvi_optim <- function(test.grid, models, ref_point, pareto_front, n_mc,
                                constraints, eta) {
  
  constraints <- list(
    function(X_cand) { predict(models[[1]], newdata = X_cand, type = "UK")$mean- given_point[1] },  # 强度 > 400 → 400 - y1 < 0
    function(X_cand) {predict(models[[2]] , newdata = X_cand, type = "UK")$mean- given_point[2] }   # 塑性 > 0.1 → 0.1 - y2 < 0
  )
  
  
  X_cand <- test.grid#t(as.data.frame(x))
  # X_cand <- (x)
  q <- nrow(X_cand)
  n_obj <- length(models)
  
  # means_list <- list()
  #sds_list <- list()
  
  model.result <- matrix(NA, nrow = 1, ncol = 4)
  colnames(model.result) <- c("mean1", "mean2", "sd1", "sd2")
  
  samples <- lapply(1:n_obj, function(i) {
    pred <- predict(models[[i]], newdata = X_cand, type = "UK", checkNames = FALSE)
    
    model.result[1, i] <<- pred$mean   # 注意：<<-
    model.result[1, i + 2] <<- pred$sd # 注意：sd放在第3和第4列
    
    matrix(rep(pred$mean, each = n_mc), nrow = q) +
      matrix(rnorm(q * n_mc, sd = rep(pred$sd, each = n_mc)), nrow = q)
  })
  
  samples <- array(unlist(samples), dim = c(q, n_mc, n_obj))
  
  if (length(constraints) > 0) {
    constraint_probs <- sapply(constraints, function(con_fun) {
      pred <- con_fun(X_cand)
      sigmoid <- function(z) 1 / (1 + exp(-eta * z))
      sigmoid(-pred)
    })
    feasibility <- if (is.matrix(constraint_probs)) {
      apply(constraint_probs, 1, prod)
    } else {
      constraint_probs
    }
  } else {
    feasibility <- rep(1, q)
  }
  
  hv_improvements <- sapply(1:n_mc, function(i) {
    sampled_obj <- samples[, i, ]
    front_aug <- rbind(pareto_front, sampled_obj)
    emoa::dominated_hypervolume(t(front_aug), ref = ref_point) -
      emoa::dominated_hypervolume(t(pareto_front), ref = ref_point)
  })
  
  mean_hvi <- mean(hv_improvements)
  final_score <- mean_hvi * mean(feasibility)
  
  result_df0 <- data.frame(
    
    feasibility = t(as.data.frame(feasibility)),
    constrained_ehvi = rep(final_score, q)
  )
  result_df <- cbind(model.result,result_df0)
  return(result_df)
}



fn_EHI_sigmoid <- function(test.grid, compute_qehvi_optim,
                           models, ref_point, pareto_front, n_mc,
                           constraints, eta) {
  ehvi_result0 <- lapply(1:nrow(test.grid), function(i) {
    compute_qehvi_optim(test.grid[i, , drop = FALSE],
                        models, ref_point, pareto_front, n_mc,
                        constraints, eta)
  })
  ehvi_result=do.call(rbind, ehvi_result0)
  rownames (ehvi_result)=1:nrow(ehvi_result)
  colnames(ehvi_result)[which(colnames(ehvi_result)=="constrained_ehvi")]="score"
  return(ehvi_result)
}

#############NSGA-II
make_obj_fn <- function(test.grid,models, given_point) {
  x_df <- test.grid#as.data.frame(matrix(test.grid, nrow = 1))
  colnames(x_df) <- colnames(test.grid)
  #x_df=test.grid
  # 预测目标值
  y1_pred <- predict(models[[1]], newdata = x_df, type = "UK")$mean
  y2_pred <- predict(models[[2]], newdata = x_df, type = "UK")$mean
  
  # 目标函数值（需最小化）
  obj1 <- y1_pred
  obj2 <- y2_pred
  
  # 约束函数：目标是 y1 < 400, y2 < 0.1
  g1 <- pmax(0, y1_pred - given_point[1])
  g2 <-pmax(0, y2_pred - given_point[2])
  g_violation <- g1 + g2
  res_matrix=cbind(obj1, obj2, g_violation)
  return(res_matrix)
}

# -------- 多目标遗传算法 NSGA-II --------
fn_NSGA_II <- function(test.grid, models) {
  
  res_matrix=make_obj_fn(test.grid,models, given_point)
  obj_fn <- function(x) {
    idx <- round(x[1])
    return(res_matrix[idx, ])
  }
  res <- nsga2(fn = obj_fn,
               idim = 1,
               odim = 3,
               lower.bounds = c(1),
               upper.bounds = as.numeric(nrow(test.grid)),
               popsize = 40,
               generations = 50)
# -------- 提取可行解 --------
feasible_idx <- which(res$value[, 3] <= 1e-6)
selected_idx <-unique(round(res$par[feasible_idx,]))
#pareto_obj <- res$value[feasible_idx, 1:2]  # mean1, mean2
# -------- 构建输出数据框 --------

  mean1 = res_matrix[,1]
  mean2 = res_matrix[,2]
  score = mean1*mean2*0  # 简单评分：目标加和
df_result <- data.frame(mean1,mean2,score)

set.seed(seed_num *(iter+1))
id=sample(seq_along(selected_idx), 1)
df_result[selected_idx[id],"score"]=1

return(df_result)
}
############main function for bi-objective optimization #####################################################

fn_multi_obj=function(crit,vir_target,train_data,selector_num,test.grid,train_pareto,ref00,tar_point00,scale,model){
  #############CEHVI
  if(selector_num=="3points_dist_CEHVI"){
    opt_algri=fn_P_EHVI#fn_EHVI
    score_calculate= fn_3points_EHI(crit="3points_dist",vir_target,test.grid,train_data,train_pareto,ref00,tar_point00,scale,opt_algri,model)
  }
  if(selector_num=="given_constrain_CEHVI"){
    opt_algri=fn_P_EHVI#fn_EHVI
    score_calculate= fn_3points_EHI(crit="given_constrain",vir_target,test.grid,train_data,train_pareto,ref00,tar_point00,scale,opt_algri,model)
  }
  if(selector_num=="3points_angle_CEHVI"){
    opt_algri=fn_P_EHVI#fn_EHVI
    score_calculate= fn_3points_EHI(crit="3points",vir_target,test.grid,train_data,train_pareto,ref00,tar_point00,scale,opt_algri,model)
  }
  if(selector_num=="all_EHVI"){
    opt_algri=fn_P_EHVI
    score_calculate= fn_3points_EHI(crit="all_pareto",vir_target,test.grid,train_data,train_pareto,ref00,tar_point00,scale,opt_algri,model)
  }
  #####EI*EI
  if(selector_num=="EI_EI"){
    set_point=as.data.frame(t(given_point))
    colnames(set_point)=c("es1","es2")
    opt_algri=fn.ego.ei
    score_calculate=fn_EIEI(vir_target,test.grid,set_point,train_pareto,ref00,tar_point00,scale,opt_algri=fn.ego.ei,model)
  }
  if(selector_num=="3points_dist_EI_EI"){
    set_point0=set_point_3points(crit="3points_dist",vir_target,train_data,train_pareto,ref00,tar_point00,scale)
    opt_algri=fn.ego.ei
    score_calculate=fn_EIEI(vir_target,test.grid,set_point0,train_pareto,ref00,tar_point00,scale,opt_algri=fn.ego.ei,model)
  }
  if(selector_num=="3points_angle_EI_EI"){
    set_point0=set_point_3points(crit="3points",vir_target,train_data,train_pareto,ref00,tar_point00,scale)
    opt_algri=fn.ego.ei
    score_calculate=fn_EIEI(vir_target,test.grid,set_point0,train_pareto,ref00,tar_point00,scale,opt_algri=fn.ego.ei,model)
  }
  if(selector_num=="all_EI_EI"){
    set_point0=t(as.data.frame(ref00))
    colnames(set_point0)=c("es1","es2")
    opt_algri=fn.ego.ei
    score_calculate=fn_EIEI(vir_target,test.grid,set_point0,train_pareto,ref00,tar_point00,scale,opt_algri=fn.ego.ei,model)
  }
  ######EMI
  if(selector_num=="given_constrain_CEMI"){
    ###############prediction
    pred <- predict_kms(model, newdata = test.grid, type = "UK", 
                        checkNames = FALSE, light.return = TRUE, cov.compute = FALSE)
    pre_mean <- t(as.data.frame(pred$mean))
    pre_sd <- t(as.data.frame(pred$sd))
    ##################new addition  
    # pre_mean=t(as.data.frame(mu))
    colnames(pre_mean)=c("mean1","mean2")
    # pre_sd=t(as.data.frame(sigma))
    colnames(pre_sd)=c("sd1","sd2")
    pre=cbind(pre_mean,pre_sd)
    
    opt_algri=fn_crit_CEMI
    score_calculate0=fn_CEMI(crit="given_constrain",test.grid,train_pareto,ref00,tar_point00,scale,opt_algri, model)
    colnames(score_calculate0)="score"
    score_calculate=cbind(pre,score_calculate0)
    
  }
  if(selector_num=="3points_dist_CEMI"){
    ###############prediction
    pred <- predict_kms(model, newdata = test.grid, type = "UK", 
                        checkNames = FALSE, light.return = TRUE, cov.compute = FALSE)
    pre_mean <- t(as.data.frame(pred$mean))
    pre_sd <- t(as.data.frame(pred$sd))
    ##################new addition  
    # pre_mean=t(as.data.frame(mu))
    colnames(pre_mean)=c("mean1","mean2")
    # pre_sd=t(as.data.frame(sigma))
    colnames(pre_sd)=c("sd1","sd2")
    pre=cbind(pre_mean,pre_sd)
    
    opt_algri=fn_crit_CEMI
    score_calculate0=fn_CEMI(crit="3points_dist",test.grid,train_pareto,ref00,tar_point00,scale,opt_algri, model)
    colnames(score_calculate0)="score"
    score_calculate=cbind(pre,score_calculate0)
    
  }
  if(selector_num=="3points_angle_CEMI"){
    ###############prediction
    pred <- predict_kms(model, newdata = test.grid, type = "UK", 
                        checkNames = FALSE, light.return = TRUE, cov.compute = FALSE)
    pre_mean <- t(as.data.frame(pred$mean))
    pre_sd <- t(as.data.frame(pred$sd))
    ##################new addition  
    # pre_mean=t(as.data.frame(mu))
    colnames(pre_mean)=c("mean1","mean2")
    # pre_sd=t(as.data.frame(sigma))
    colnames(pre_sd)=c("sd1","sd2")
    pre=cbind(pre_mean,pre_sd)
    
    opt_algri=fn_crit_CEMI
    score_calculate0=fn_CEMI(crit="3points",test.grid,train_pareto,ref00,tar_point00,scale,opt_algri, model)
    colnames(score_calculate0)="score"
    score_calculate=cbind(pre,score_calculate0)
    
  }
  
  #################SMS
  if(selector_num=="3points_dist_CSMS"){
    opt_algri=crit_CSMS
    score_calculate=fn_CSMSEGO_grid(crit="3points_dist",ref_min,vir_target,train_data,test.grid,train_pareto,ref00,tar_point00,scale,opt_algri, model)
  }
  if(selector_num=="given_constrain_CSMS"){
    opt_algri=crit_CSMS
    score_calculate=fn_CSMSEGO_grid(crit="given_constrain",ref_min,vir_target,train_data,test.grid,train_pareto,ref00,tar_point00,scale,opt_algri, model)
  }
  
  
  
  if(selector_num=="SMS"){
    opt_algri=fn_crit_SMS
    score_calculate=fn_SMS(test.grid,train_pareto,ref00,tar_point00,scale,opt_algri, model)
  }
  
  ###### EHVI
  if(selector_num=="EHVI"){
    opt_algri=fn_EHVI
    score_calculate=fn_EHI(test.grid,train_pareto,ref00,tar_point00,scale,opt_algri,model)
  }

  ######constrained EHVI
  if(selector_num=="1point_CEHVI"){
    opt_algri=fn_EHVI
    score_calculate= fn_1P_EHI(crit="angle",vir_target,test.grid,train_data,train_pareto,ref00,tar_point00,scale,opt_algri,model)
  }
  ######constrained EHVI
  if(selector_num=="1point_dis_CEHVI"){
    opt_algri=fn_EHVI
    score_calculate= fn_1P_EHI(crit="dist",vir_target,test.grid,train_data,train_pareto,ref00,tar_point00,scale,opt_algri,model)
  }
  #####minimize the distance to the ideal point
  if(selector_num=="mean_dist"){
    opt_algri=fn_dist
    score_calculate=fn_mean_dist(vir_target,test.grid,train_data,train_pareto,ref00,tar_point00,scale,opt_algri=fn_dist,model)
  }

  ######SUR
  if(selector_num=="SUR"){
    opt_algri=fn_crit_SUR
    score_calculate=fn_SUR(test.grid,train_pareto,ref00,tar_point00,scale,opt_algri, model)
  }
  ######EMI
  if(selector_num=="EMI"){
    ###############prediction
    pred <- predict_kms(model, newdata = test.grid, type = "UK", 
                      checkNames = FALSE, light.return = TRUE, cov.compute = FALSE)
    pre_mean <- t(as.data.frame(pred$mean))
    pre_sd <- t(as.data.frame(pred$sd))
    ##################new addition  
    # pre_mean=t(as.data.frame(mu))
     colnames(pre_mean)=c("mean1","mean2")
   # pre_sd=t(as.data.frame(sigma))
     colnames(pre_sd)=c("sd1","sd2")
     pre=cbind(pre_mean,pre_sd)
    
    opt_algri=fn_crit_EMI
    score_calculate0=fn_EMI(test.grid,train_pareto,ref00,tar_point00,scale,opt_algri, model)
    colnames(score_calculate0)="score"
    score_calculate=cbind(pre,score_calculate0)
    
     }


  # if(selector_num=="SMS_epsilon"){
  #   opt_algri=fn_epsilon_SMS_self
  #   SMS_calculate=fn_SMS_epsilon(test.grid,train_pareto,ref00,tar_point00,scale,opt_algri, model)
  # }
  #####USeMO
  if(selector_num=="USeMO"){
    set_point=as.data.frame(t(given_point))
    colnames(set_point)=c("es1","es2")
    opt_algri=fn.ego.ei
    EI_result=fn_EIEI(vir_target,test.grid,set_point,train_pareto,ref00,tar_point00,scale,opt_algri=fn.ego.ei,model)
    score_calculate=fn_USeMO(EI_result,test.grid)
  
  }
  ###### EHVI_sigmoid
  if(selector_num=="EHVI_sigmoid"){
    
    score_calculate=fn_EHI_sigmoid(test.grid, compute_qehvi_optim,models=model,ref_point=ref00, pareto_front=train_pareto, n_mc = 10,
                                            constraints, eta = 10)
  }
  
  if(selector_num=="NSGA-II"){
    score_calculate= fn_NSGA_II(test.grid, models=model)
  }
  return(score_calculate)
}

#######################################################################################################