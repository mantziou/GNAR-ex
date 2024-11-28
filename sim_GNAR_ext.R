########################################################
############## SIMULATIONS ############################
########################################################
library(igraph)
library(dplyr)
library(intervals)
library(forecast)
library(reshape2)
library(ggplot2)

########### ########### ##### ##### ##### ##### ##### ##### 
##### ##### ########### REPETITIONS ########### ##### ##### 
##### ##### ##### ##### ##### ##### ########### ########### 

# Regime 1: GNAR(1,[1])
alpha_par_1 <- .2
beta_par_1 <- list(c(0.2))
gamma_par_1 <- c(.3)
delta_par_1 <- list(c(0.2))

# Regime 2: GNAR(1,[2])
alpha_par_2 <- .2
beta_par_2 <- list(c(-0.2,0.1))
gamma_par_2 <-  c(.1)
delta_par_2 <- list(c(0.05,-0.2))

# Regime 3: GNAR(2,[2,2])
alpha_par_3 <- c(-.1,.3)
beta_par_3 <- list(c(.1,-.2),c(-.02,.03))
gamma_par_3 <- c(.01,.01)
delta_par_3 <- list(c(.02,-.01),c(-.02,.01))


set.seed(30)
seeds <- sample(1:500,50,replace=FALSE)
#SBM parameters
probmat <- matrix(c(.7,.2,.1,.7),nrow = 2)


alphaOrder_reg <- c(1,1,2)
betaOrder_reg <- list(c(1),c(2),rep(2,2))
gammaOrder_reg <- c(1,1,2)
deltaOrder_reg <- list(c(1),c(2),c(2,2))

num_param <- c(4,6,12) 

for (reg in 1:3){
  for (mod in c("grg","er","sbm")){
    print(c("regime ",reg, " model ",mod))
    cov_mat <- matrix(nrow=num_param[reg],ncol = 50)
    rmse_mat <- matrix(nrow=num_param[reg],ncol = 50)
    for (i in 1:50){
      print(c("seed ",i))
      set.seed(seeds[i])
      
      if (mod=="sbm"){# SBM model
        net <- sample_sbm(20,probmat,c(10,10),directed = TRUE)
      }
      if (mod=="er"){# ER model
        net <- erdos.renyi.game(20,p.or.m = 168,type = "gnm",directed = TRUE) # p=.4 fix the number of edges in graph ("gnm") rather than give prob of edge ("gnp")
      }
      if (mod=="grg"){# Random Geometric Graph
        lpvs2 <- sample_sphere_surface(dim = 2, n = 20,radius = .7)#.35 for 0.1 density, .7 for density 0.4
        net <- sample_dot_product(lpvs2,directed = TRUE)
      }
      
      nedges <- ecount(net)
      nnodes <- vcount(net)
      n_edges_nodes <- nedges+nnodes
      data_edges <- get.edgelist(net)
      data_nodes <- V(net)
      E(net)$name <- c(1:nedges)
      V(net)$name <- c(1:nnodes)#c((nedges+1):(nnodes+nedges))
      
      simdata <- gnar_x_sim(n=200, net, alphaParams=get(paste("alpha_par_",reg,sep = "")), betaParams=get(paste("beta_par_",reg,sep = "")),
                            gammaParams=get(paste("gamma_par_",reg,sep = "")), deltaParams=get(paste("delta_par_",reg,sep = "")), 
                            sigma=1, meann=0, nedges,nnodes,n_edges_nodes,
                            data_edges,data_nodes,noise="n",vt=3,rcor=0.5,pay_now = FALSE)
      
      rownames(simdata) <- c(paste(data_edges[,1],data_edges[,2]),seq(1,nnodes))
      
      fit <- gnar_x_fit(simdata,nnodes,nedges,n_edges_nodes,data_edges,data_nodes,alphaOrder_reg[reg],betaOrder_reg[[reg]],
                        gammaOrder_reg[reg],deltaOrder_reg[[reg]],
                        net,lead_lag_mat,globalalpha=TRUE,lead_lag_weights=FALSE,pay_now = FALSE,lag_0_sep = FALSE)
      #print(fit$mod$coefficients)
      
      comp_truevec <- c()
      if (pay_now==FALSE){
        ordstart <- 1
      }else{
        ordstart <- 0
      }
      for (ord in ordstart:max(alphaOrder_reg[reg],gammaOrder_reg[reg])){
        if (ord==0){
          if (deltaOrder_reg[[reg]][ord+1]>0){
            comp_truevec <- c(comp_truevec,get(paste("gamma_par_",reg,sep = ""))[ord+1],get(paste("delta_par_",reg,sep = ""))[[ord+1]])
          }else{
            comp_truevec <- c(comp_truevec,get(paste("gamma_par_",reg,sep = ""))[ord+1])
          }
        }else{
          if (betaOrder_reg[[reg]][ord]>0){
            comp_truevec <- c(comp_truevec,get(paste("alpha_par_",reg,sep = ""))[ord],get(paste("beta_par_",reg,sep = ""))[[ord]])
          }else{
            comp_truevec <- c(comp_truevec,get(paste("alpha_par_",reg,sep = ""))[ord])
          }
          if (pay_now==FALSE){
            ordstag <- ord
          }else{
            ordstag <- ord+1
          }
          if (deltaOrder_reg[[reg]][ordstag]>0){
            comp_truevec <- c(comp_truevec,get(paste("gamma_par_",reg,sep = ""))[ordstag],get(paste("delta_par_",reg,sep = ""))[[ordstag]])
          }else{
            comp_truevec <- c(comp_truevec,get(paste("gamma_par_",reg,sep = ""))[ordstag])
          }
        }
      }
      
      if (lag_0_sep==FALSE){
        for(tr in 1:length(comp_truevec)){
          cov_mat[tr,i] <- between(comp_truevec[tr],confint(fit$mod)[tr,1],confint(fit$mod)[tr,2])
        }
        rmse_mat[,i] <- (fit$mod$coefficients-comp_truevec)^2
      }else{
        indtr <- 1
        for(tr in 1:length(comp_truevec)){
          if (tr<=length(fit$modlag0$coefficients)){
            cov_mat[tr,i] <- between(comp_truevec[tr],confint(fit$modlag0)[tr,1],confint(fit$modlag0)[tr,2])
          }else{
            cov_mat[tr,i] <- between(comp_truevec[tr],confint(fit$mod)[indtr,1],confint(fit$mod)[indtr,2])
            indtr <- indtr+1
          }
        }
        rmse_mat[,i] <- (c(fit$modlag0$coefficients,fit$mod$coefficients)-comp_truevec)^2
      }
    }
    if (lag_0_sep==FALSE){
      rownames(cov_mat) <- rownames(confint(fit$mod,level=.95))
      rownames(rmse_mat) <- names(fit$mod$coefficients)
    }else{
      rownames(cov_mat) <- c(rownames(confint(fit$modlag0,level=.95)),rownames(confint(fit$mod,level=.95)))
      rownames(rmse_mat) <- c(names(fit$modlag0$coefficients),names(fit$mod$coefficients))
    }
    assign(paste("coverage_reg_",reg,"_mod_",mod,sep = ""),apply(cov_mat,1,sum)/50)
    assign(paste("rmse_reg_",reg,"_mod_",mod,sep = ""),apply(rmse_mat,1,function(x) sqrt(mean(x))))
  }
}

########### ########### ########### ########### ########### ########### ########### ########### ########### 
########### PREDICTIVE PERFORMANCE: REPETITIONS OF SBM,ER,SBM GRAPHS WITH 3RD SIM REGIME ###########
########### ########### ########### ########### ########### ########### ########### ########### ########### 

set.seed(12) #used for 0.4 density, and 0.23 density
seeds <- sample(seq(1,1000),50)# used for 0.4 density, and 0.23 density

#SBM parameters
probmat <- matrix(c(.7,.2,.1,.7),nrow = 2) #for 0.4 density

rmse_mat <- rmse_mat_all <-  matrix(ncol = 12,nrow = 50)
colnames(rmse_mat_all) <- colnames(rmse_mat) <- c("gnarnei_er","gnar_er","hlag_er","arima_er",
                                                  "gnarnei_sbm","gnar_sbm","hlag_sbm","arima_sbm",
                                                  "gnarnei_grg","gnar_grg","hlag_grg","arima_grg")

for (mod in c("grg","er","sbm")){
  for (i in 1:50){
    print(i)
    set.seed(seeds[i])
    if (mod=="sbm"){# SBM model
      net <- sample_sbm(20,probmat,c(10,10),directed = TRUE)
    }
    if (mod=="er"){# ER model
      net <- erdos.renyi.game(20,p.or.m = 168,type = "gnm",directed = TRUE) # p=.4 fix the number of edges in graph ("gnm") rather than give prob of edge ("gnp"), 168 edge gnm for 0.4 density
    }
    if (mod=="grg"){# Random Geometric Graph
      lpvs2 <- sample_sphere_surface(dim = 2, n = 20,radius = .7)#.35 for 0.1 density, .7 for density 0.4, 0.55 for density 0.24
      net <- sample_dot_product(lpvs2,directed = TRUE)
    }
    
    nedges <- ecount(net)
    nnodes <- vcount(net)
    n_edges_nodes <- nedges+nnodes
    data_edges <- get.edgelist(net)
    data_nodes <- V(net)
    E(net)$name <- c(1:nedges)
    V(net)$name <- c(1:nnodes)#c((nedges+1):(nnodes+nedges))
    
    simdata <- gnar_x_sim(n=200, net, alphaParams=alpha_par_3, betaParams=beta_par_3,
                          gammaParams=gamma_par_3, deltaParams=delta_par_3, 
                          sigma=1, meann=0, nedges,nnodes,n_edges_nodes,
                          data_edges,data_nodes,noise="n",vt=3,rcor=0.5,pay_now = FALSE)
    
    rownames(simdata) <- c(paste(data_edges[,1],data_edges[,2]),seq(1,nnodes))
    
    datasim_train <- simdata[,1:199]
    
    # Fit GNAR(2,[2,2])
    fit_train <- gnar_x_fit(datasim_train,nnodes,nedges,n_edges_nodes,data_edges,data_nodes,alphaOrder_reg[3],betaOrder_reg[[3]],
                            gammaOrder_reg[3],deltaOrder_reg[[3]],
                            net,lead_lag_mat,globalalpha=TRUE,lead_lag_weights=FALSE,pay_now = FALSE, lag_0_sep = FALSE)
    fit_predict <- gnar_x_predict(fit_train,datasim_train,alphaOrder_reg[3],betaOrder_reg[[3]],
                                  gammaOrder_reg[3],deltaOrder_reg[[3]],n_edges_nodes,nnodes,nedges,1,fit_train$wei_mat,fit_train$data_loc_mat,
                                  pay_now=FALSE,pay_data_now=NULL) 
    rmse_mat[i,paste("gnarnei_",mod,sep = "")] <- sqrt(mean((fit_predict[(nedges+1):n_edges_nodes,alphaOrder_reg[3]+1]-simdata[(nedges+1):n_edges_nodes,200])^2))
    rmse_mat_all[i,paste("gnarnei_",mod,sep = "")] <- sqrt(mean((fit_predict[,alphaOrder_reg[3]+1]-simdata[,200])^2))
    
    print(c("seed: ",i," gnar nei done, model ", mod))
    
    # Fit GNAR(2,[0,0],2,[0,0])
    fit_train_nonei <- gnar_x_fit(datasim_train,nnodes,nedges,n_edges_nodes,data_edges,data_nodes,alphaOrder_reg[3],c(0,0),
                                  gammaOrder_reg[3],c(0,0),net,lead_lag_mat,globalalpha=TRUE,lead_lag_weights=FALSE,pay_now = FALSE,lag_0_sep = FALSE)

    fit_predict_nonei <- gnar_x_predict(fit_train_nonei,datasim_train,alphaOrder_reg[3],c(0,0),
                                        gammaOrder_reg[3],c(0,0),n_edges_nodes,nnodes,nedges,1,fit_train_nonei$wei_mat,fit_train_nonei$data_loc_mat,
                                        pay_now=FALSE,pay_data_now=NULL)
    rmse_mat[i,paste("gnar_",mod,sep = "")] <- sqrt(mean((fit_predict_nonei[(nedges+1):n_edges_nodes,alphaOrder_reg[3]+1]-simdata[(nedges+1):n_edges_nodes,200])^2))
    rmse_mat_all[i,paste("gnar_",mod,sep = "")] <- sqrt(mean((fit_predict_nonei[,alphaOrder_reg[3]+1]-simdata[,200])^2))
    
    print(c("seed: ",i," gnar nonei done, model ", mod))
    
    simtrainvar <- t(datasim_train)
    
    # Autoarima
    simple_arima <- apply(simtrainvar, 2, function(x){forecast(auto.arima(x),h=1)$mean})
    rmse_mat[i,paste("arima_",mod,sep = "")] <- sqrt(mean((simple_arima[(nedges+1):n_edges_nodes]-simdata[(nedges+1):n_edges_nodes,200])^2))
    rmse_mat_all[i,paste("arima_",mod,sep = "")] <- sqrt(mean((simple_arima-simdata[,200])^2))
    print(c("seed: ",i," ar done"))
    
    # Fit bigvar  
    Model2 <- constructModel(simtrainvar,p=2,struct='HLAGC',gran=c(50,10),verbose=FALSE) # HLAG 
    results2 <- cv.BigVAR(Model2)
    varfor2 <- predict(results2,n.ahead=1)
    rmse_mat[i,paste("hlag_",mod,sep = "")] <- sqrt(mean((varfor2[(nedges+1):n_edges_nodes,]-simdata[(nedges+1):n_edges_nodes,200])^2))
    rmse_mat_all[i,paste("hlag_",mod,sep = "")] <- sqrt(mean((varfor2-simdata[,200])^2))
    
    print(c("seed: ",i," var done"))
  }
}

########################################################################
################## VISUALISATION OF RESULTS ##################
########################################################################

rmse_df <- as.data.frame(rmse_mat_all[,c(1,2,4,3)])
rmse_df <- as.data.frame(rmse_mat_all[,c(5,6,8,7)])
rmse_df <- as.data.frame(rmse_mat_all[,c(9,10,12,11)])

colnames(rmse_df) <- rep(c("GNAR-ex(2,[2,2])","GNAR-ex(2,[0,0])","ARIMA","HLAG"),1)
rmse_melt <- melt(rmse_df)
rmse_melt$variable <- as.factor(rmse_melt$variable)
colnames(rmse_melt) <- c("m","RMSE")
p1 <- ggplot(rmse_melt, aes(x=m, y=RMSE)) +
  geom_boxplot(fill="red")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))+xlab(" ")+ggtitle("ER")+theme(plot.title = element_text(hjust = 0.5))
p2 <- ggplot(rmse_melt, aes(x=m, y=RMSE)) +
  geom_boxplot(fill="yellowgreen")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))+xlab(" ")+ggtitle("SBM")+theme(plot.title = element_text(hjust = 0.5))+ylab("")
p3 <- ggplot(rmse_melt, aes(x=m, y=RMSE)) +
  geom_boxplot(fill="cornflowerblue")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))+xlab(" ")+ggtitle("RDP")+
  theme(plot.title = element_text(hjust = 0.5))+ylab("")
library(patchwork)
p1+p2+p3


########### ########### ########### ########### ########### ########### ########### ########### ########### 
########### MODEL AVERAGING PREDICTIVE PERFORMANCE: REPETITIONS OF SBM,ER,SBM GRAPHS WITH 3rd SIM REGIME ###########
########### ########### ########### ########### ########### ########### ########### ########### ########### 



set.seed(12) #used for 0.4 density, and 0.23 density
seeds <- sample(seq(1,1000),50)# used for 0.4 density, and 0.23 density

#SBM parameters
probmat <- matrix(c(.7,.2,.1,.7),nrow = 2) #for 0.4 density

errormat <- matrix(ncol = 9,nrow = 50)

colnames(errormat) <- c("aver_pred_rmse_nodeonly_er","aver_pred_rmse_nodeonly_sbm","aver_pred_rmse_nodeonly_grg",
                        "aver_pred_rmse_all_er","aver_pred_rmse_all_sbm","aver_pred_rmse_all_grg",
                        "aver_totalnode_pred_relerr_er","aver_totalnode_pred_relerr_sbm","aver_totalnode_pred_relerr_grg")

for (mod in c("grg","er","sbm")){#
  pred_incl <- c()
  for (i in 1:50){
    print(i)
    set.seed(seeds[i])
    if (mod=="sbm"){# SBM model
      net <- sample_sbm(20,probmat,c(10,10),directed = TRUE)
    }
    if (mod=="er"){# ER model
      net <- erdos.renyi.game(20,p.or.m = 168,type = "gnm",directed = TRUE) # p=.4 fix the number of edges in graph ("gnm") rather than give prob of edge ("gnp"), 168 edge gnm for 0.4 density
    }
    if (mod=="grg"){# Random Geometric Graph
      lpvs2 <- sample_sphere_surface(dim = 2, n = 20,radius = .7)#.35 for 0.1 density, .7 for density 0.4, 0.55 for density 0.24
      net <- sample_dot_product(lpvs2,directed = TRUE)
    }
    
    nedges <- ecount(net)
    nnodes <- vcount(net)
    n_edges_nodes <- nedges+nnodes
    data_edges <- get.edgelist(net)
    data_nodes <- V(net)
    E(net)$name <- c(1:nedges)
    V(net)$name <- c(1:nnodes)#c((nedges+1):(nnodes+nedges))
    
    simdata <- gnar_x_sim(n=200, net, alphaParams=alpha_par_3, betaParams=beta_par_3,
                          gammaParams=gamma_par_3, deltaParams=delta_par_3, 
                          sigma=1, meann=0, nedges,nnodes,n_edges_nodes,
                          data_edges,data_nodes,noise="n",vt=3,rcor=0.5,pay_now = FALSE)
    
    rownames(simdata) <- c(paste(data_edges[,1],data_edges[,2]),seq(1,nnodes))
    
    datasim_train <- simdata[,1:199]
    
    # Fit GNAR-ex 1 stage nei and various lag=1,...,9
    fit_pred_mat <- matrix(nrow = n_edges_nodes,ncol = 9)
    upper <- lower <- matrix(nrow = nnodes,ncol = 9)
    unionint <- lapply(1:length(1:50), function(x) vector("list",nnodes))
    for (lagit in 1:9){
      fit_train <- gnar_x_fit(datasim_train,nnodes,nedges,n_edges_nodes,data_edges,data_nodes,lagit,rep(1,lagit),
                              lagit,rep(1,lagit),
                              net,lead_lag_mat,globalalpha=TRUE,lead_lag_weights=FALSE,pay_now = FALSE, lag_0_sep = FALSE)
      fit_predict <- gnar_x_predict(fit_train,datasim_train,lagit,rep(1,lagit),
                                    lagit,rep(1,lagit),n_edges_nodes,nnodes,nedges,1,fit_train$wei_mat,fit_train$data_loc_mat,
                                    pay_now=FALSE,pay_data_now=NULL) 
      fit_pred_mat[,lagit] <- fit_predict[,lagit+1]

      # Prediction intervals 
      resmat <- matrix(fit_train$mod$residuals,  ncol=n_edges_nodes, byrow=FALSE)
      covresmat <- (t(resmat) %*% resmat)/ncol(datasim_train) # covariance matrix of residuals
      predsd_ver <- sqrt(diag(covresmat))
      
      # Get confidence intervals 
      upper[,lagit] <- tail(fit_predict[,lagit+1],n=nnodes)+1.96*tail(as.vector(predsd_ver),n=nnodes)
      lower[,lagit]  <- tail(fit_predict[,lagit+1],n=nnodes)-1.96*tail(as.vector(predsd_ver),n=nnodes)
      
    }
    # rmse of MA forecast (node and all)
    average_pred <- apply(fit_pred_mat,1,mean)
    errormat[i,paste("aver_pred_rmse_nodeonly_",mod,sep = "")] <- sqrt(mean((average_pred[(nedges+1):n_edges_nodes]-simdata[(nedges+1):n_edges_nodes,200])^2))
    errormat[i,paste("aver_pred_rmse_all_",mod,sep = "")] <- sqrt(mean((average_pred-simdata[,200])^2))
    
    # relative error of MA total node forecast
    average_totalnode_pred <- mean(apply(fit_pred_mat[(nedges+1):n_edges_nodes,], 2, sum))
    totalnode_true <- sum(simdata[(nedges+1):n_edges_nodes,200])
    errormat[i,paste("aver_totalnode_pred_relerr_",mod,sep = "")] <- abs((average_totalnode_pred-totalnode_true))/totalnode_true
    
    # get union of intervals per node
    for(inode in 1:nnodes){
      unionint[[i]][[inode]] <- interval_union(Intervals(t(rbind(lower[inode,],upper[inode,]))))
    }
    unionintmat <- do.call("rbind",unionint[[i]])
    pred_incl[i] <- length(which(between(simdata[(nedges+1):n_edges_nodes,200],unionintmat[,1],
                                                                unionintmat[,2])))/nnodes
    print(c("seed: ",i," gnar nei done, model ", mod))
    
  }
  assign(paste("pred_incl_",mod,sep = ""),pred_incl)
}

########################################################################
################## VISUALISATION OF RESULTS MA REG 3 ##################
########################################################################

rmse_df <- as.data.frame(cbind(errormat[,1],rmse_mat[,4])) #nodeonly
rmse_df <- as.data.frame(cbind(errormat[,2],rmse_mat[,8])) #nodeonly
rmse_df <- as.data.frame(cbind(errormat[,3],rmse_mat[,12])) #nodeonly

rmse_df <- as.data.frame(cbind(errormat[,4],rmse_mat_all[,4])) #all
rmse_df <- as.data.frame(cbind(errormat[,5],rmse_mat_all[,8])) #all
rmse_df <- as.data.frame(cbind(errormat[,6],rmse_mat_all[,12])) #all


colnames(rmse_df) <- rep(c("MA GNAR-ex ","ARIMA"),1)
rmse_melt <- melt(rmse_df)
rmse_melt$variable <- as.factor(rmse_melt$variable)
colnames(rmse_melt) <- c("m","RMSE")
p1 <- ggplot(rmse_melt, aes(x=m, y=RMSE)) +
  geom_boxplot(fill="red")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))+xlab(" ")+ggtitle("ER")+theme(plot.title = element_text(hjust = 0.5))
p2 <- ggplot(rmse_melt, aes(x=m, y=RMSE)) +
  geom_boxplot(fill="yellowgreen")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))+xlab(" ")+ggtitle("SBM")+theme(plot.title = element_text(hjust = 0.5))+ylab("")
p3 <- ggplot(rmse_melt, aes(x=m, y=RMSE)) +
  geom_boxplot(fill="cornflowerblue")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))+xlab(" ")+ggtitle("RDP")+
  theme(plot.title = element_text(hjust = 0.5))+ylab("")
library(patchwork)
p1+p2+p3




########### ########### ########### ########### ########### ########### ########### ########### ########### 
########### MODEL AVERAGING PREDICTIVE PERFORMANCE: REPETITIONS OF SBM,ER,SBM GRAPHS WITH 1st SIM REGIME ###########
########### ########### ########### ########### ########### ########### ########### ########### ########### 



set.seed(12) #used for 0.4 density, and 0.23 density
seeds <- sample(seq(1,1000),50)# used for 0.4 density, and 0.23 density

#SBM parameters
probmat <- matrix(c(.7,.2,.1,.7),nrow = 2) #for 0.4 density

errormat <- matrix(ncol = 9,nrow = 50)

colnames(errormat) <- c("aver_pred_rmse_nodeonly_er","aver_pred_rmse_nodeonly_sbm","aver_pred_rmse_nodeonly_grg",
                        "aver_pred_rmse_all_er","aver_pred_rmse_all_sbm","aver_pred_rmse_all_grg",
                        "aver_totalnode_pred_relerr_er","aver_totalnode_pred_relerr_sbm","aver_totalnode_pred_relerr_grg")

for (mod in c("grg","er","sbm")){
  pred_incl <- c()
  for (i in 1:50){
    print(i)
    set.seed(seeds[i])
    if (mod=="sbm"){# SBM model
      net <- sample_sbm(20,probmat,c(10,10),directed = TRUE)
    }
    if (mod=="er"){# ER model
      net <- erdos.renyi.game(20,p.or.m = 168,type = "gnm",directed = TRUE) # p=.4 fix the number of edges in graph ("gnm") rather than give prob of edge ("gnp"), 168 edge gnm for 0.4 density
    }
    if (mod=="grg"){# Random Geometric Graph
      lpvs2 <- sample_sphere_surface(dim = 2, n = 20,radius = .7)#.35 for 0.1 density, .7 for density 0.4, 0.55 for density 0.24
      net <- sample_dot_product(lpvs2,directed = TRUE)
    }
    
    nedges <- ecount(net)
    nnodes <- vcount(net)
    n_edges_nodes <- nedges+nnodes
    data_edges <- get.edgelist(net)
    data_nodes <- V(net)
    E(net)$name <- c(1:nedges)
    V(net)$name <- c(1:nnodes)#c((nedges+1):(nnodes+nedges))
    
    simdata <- gnar_x_sim(n=200, net, alphaParams=alpha_par_1, betaParams=beta_par_1,
                          gammaParams=gamma_par_1, deltaParams=delta_par_1, 
                          sigma=1, meann=0, nedges,nnodes,n_edges_nodes,
                          data_edges,data_nodes,noise="n",vt=3,rcor=0.5,pay_now = FALSE)
    
    rownames(simdata) <- c(paste(data_edges[,1],data_edges[,2]),seq(1,nnodes))
    
    datasim_train <- simdata[,1:199]
    
    # Fit GNAR-ex 1 stage nei and various lag=1,...,9
    fit_pred_mat <- matrix(nrow = n_edges_nodes,ncol = 9)
    upper <- lower <- matrix(nrow = nnodes,ncol = 9)
    unionint <- lapply(1:length(1:50), function(x) vector("list",nnodes))
    for (lagit in 1:9){
      fit_train <- gnar_x_fit(datasim_train,nnodes,nedges,n_edges_nodes,data_edges,data_nodes,lagit,rep(1,lagit),
                              lagit,rep(1,lagit),
                              net,lead_lag_mat,globalalpha=TRUE,lead_lag_weights=FALSE,pay_now = FALSE, lag_0_sep = FALSE)
      fit_predict <- gnar_x_predict(fit_train,datasim_train,lagit,rep(1,lagit),
                                    lagit,rep(1,lagit),n_edges_nodes,nnodes,nedges,1,fit_train$wei_mat,fit_train$data_loc_mat,
                                    pay_now=FALSE,pay_data_now=NULL) # pay_data_now=simdata[1:nedges,200]
      fit_pred_mat[,lagit] <- fit_predict[,lagit+1]
      
      # Prediction intervals 
      resmat <- matrix(fit_train$mod$residuals,  ncol=n_edges_nodes, byrow=FALSE)
      covresmat <- (t(resmat) %*% resmat)/ncol(datasim_train) # covariance matrix of residuals
      predsd_ver <- sqrt(diag(covresmat))
      
      # Get confidence intervals 
      upper[,lagit] <- tail(fit_predict[,lagit+1],n=nnodes)+1.96*tail(as.vector(predsd_ver),n=nnodes)
      lower[,lagit]  <- tail(fit_predict[,lagit+1],n=nnodes)-1.96*tail(as.vector(predsd_ver),n=nnodes)
      
    }
    # rmse of MA forecast (node and all)
    average_pred <- apply(fit_pred_mat,1,mean)
    errormat[i,paste("aver_pred_rmse_nodeonly_",mod,sep = "")] <- sqrt(mean((average_pred[(nedges+1):n_edges_nodes]-simdata[(nedges+1):n_edges_nodes,200])^2))
    errormat[i,paste("aver_pred_rmse_all_",mod,sep = "")] <- sqrt(mean((average_pred-simdata[,200])^2))
    
    # relative error of MA total node forecast
    average_totalnode_pred <- mean(apply(fit_pred_mat[(nedges+1):n_edges_nodes,], 2, sum))
    totalnode_true <- sum(simdata[(nedges+1):n_edges_nodes,200])
    errormat[i,paste("aver_totalnode_pred_relerr_",mod,sep = "")] <- abs((average_totalnode_pred-totalnode_true))/totalnode_true
    
    # get union of intervals per node
    for(inode in 1:nnodes){
      unionint[[i]][[inode]] <- interval_union(Intervals(t(rbind(lower[inode,],upper[inode,]))))
    }
    unionintmat <- do.call("rbind",unionint[[i]])
    pred_incl[i] <- length(which(between(simdata[(nedges+1):n_edges_nodes,200],unionintmat[,1],
                                         unionintmat[,2])))/nnodes
    print(c("seed: ",i," gnar nei done, model ", mod))
    
  }
  assign(paste("pred_incl_",mod,sep = ""),pred_incl)
}

############### ############### ############### ############### ############### 
############### PREDICTION PERFORMANCE OF AUTO.ARIMA FOR SIM REG 1 ############### 
############### ############### ############### ############### ############### 
rmse_mat <- rmse_mat_all <-  matrix(ncol = 3,nrow = 50)
colnames(rmse_mat_all) <- colnames(rmse_mat) <- c("arima_er","arima_sbm","arima_grg")
for (mod in c("grg","er","sbm")){
  for (i in 1:50){
    print(i)
    set.seed(seeds[i])
    if (mod=="sbm"){# SBM model
      net <- sample_sbm(20,probmat,c(10,10),directed = TRUE)
    }
    if (mod=="er"){# ER model
      net <- erdos.renyi.game(20,p.or.m = 168,type = "gnm",directed = TRUE) # p=.4 fix the number of edges in graph ("gnm") rather than give prob of edge ("gnp"), 168 edge gnm for 0.4 density
    }
    if (mod=="grg"){# Random Geometric Graph
      lpvs2 <- sample_sphere_surface(dim = 2, n = 20,radius = .7)#.35 for 0.1 density, .7 for density 0.4, 0.55 for density 0.24
      net <- sample_dot_product(lpvs2,directed = TRUE)
    }
    
    nedges <- ecount(net)
    nnodes <- vcount(net)
    n_edges_nodes <- nedges+nnodes
    data_edges <- get.edgelist(net)
    data_nodes <- V(net)
    E(net)$name <- c(1:nedges)
    V(net)$name <- c(1:nnodes)#c((nedges+1):(nnodes+nedges))
    
    simdata <- gnar_x_sim(n=200, net, alphaParams=alpha_par_1, betaParams=beta_par_1,
                          gammaParams=gamma_par_1, deltaParams=delta_par_1, 
                          sigma=1, meann=0, nedges,nnodes,n_edges_nodes,
                          data_edges,data_nodes,noise="n",vt=3,rcor=0.5,pay_now = FALSE)
    
    rownames(simdata) <- c(paste(data_edges[,1],data_edges[,2]),seq(1,nnodes))
    
    datasim_train <- simdata[,1:199]
    
    simtrainvar <- t(datasim_train)
    
    
    # Autoarima
    simple_arima <- apply(simtrainvar, 2, function(x){forecast(auto.arima(x),h=1)$mean})
    rmse_mat[i,paste("arima_",mod,sep = "")] <- sqrt(mean((simple_arima[(nedges+1):n_edges_nodes]-simdata[(nedges+1):n_edges_nodes,200])^2))
    rmse_mat_all[i,paste("arima_",mod,sep = "")] <- sqrt(mean((simple_arima-simdata[,200])^2))
    print(c("seed: ",i," ar done"))
    
  }
}

########################################################################
################## VISUALISATION OF RESULTS MA REG 1 ##################
########################################################################

rmse_df <- as.data.frame(cbind(errormat[,1],rmse_mat[,1])) #nodeonly
rmse_df <- as.data.frame(cbind(errormat[,2],rmse_mat[,2])) #nodeonly
rmse_df <- as.data.frame(cbind(errormat[,3],rmse_mat[,3])) #nodeonly

rmse_df <- as.data.frame(cbind(errormat[,4],rmse_mat_all[,1])) #all
rmse_df <- as.data.frame(cbind(errormat[,5],rmse_mat_all[,2])) #all
rmse_df <- as.data.frame(cbind(errormat[,6],rmse_mat_all[,3])) #all

colnames(rmse_df) <- rep(c("MA GNAR-ex","ARIMA"),1)
rmse_melt <- melt(rmse_df)
rmse_melt$variable <- as.factor(rmse_melt$variable)
colnames(rmse_melt) <- c("m","RMSE")
p1 <- ggplot(rmse_melt, aes(x=m, y=RMSE)) +
  geom_boxplot(fill="red")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))+xlab(" ")+ggtitle("ER")+theme(plot.title = element_text(hjust = 0.5))
p2 <- ggplot(rmse_melt, aes(x=m, y=RMSE)) +
  geom_boxplot(fill="yellowgreen")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))+xlab(" ")+ggtitle("SBM")+theme(plot.title = element_text(hjust = 0.5))+ylab("")
p3 <- ggplot(rmse_melt, aes(x=m, y=RMSE)) +
  geom_boxplot(fill="cornflowerblue")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))+xlab(" ")+ggtitle("RDP")+
  theme(plot.title = element_text(hjust = 0.5))+ylab("")
library(patchwork)
p1+p2+p3


