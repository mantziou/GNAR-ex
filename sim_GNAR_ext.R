########################################################
############## SIMULATIONS ############################
########################################################
library(igraph)
library(dplyr)
library(intervals)
library(forecast)
library(reshape2)
library(ggplot2)

dir_path <- "/Users/u2371456/Documents/R_directory_Turing_ONS/collab"

##### ##### ##### 
##### TEST ##### 
##### ##### ##### 
set.seed(15)
net <- erdos.renyi.game(20,0.3,directed = TRUE)
edgelist_net <- get.edgelist(net)
nedges <- ecount(net)
nnodes <- vcount(net)
n_edges_nodes <- nedges+nnodes
data_edges <- get.edgelist(net)
data_nodes <- V(net)
E(net)$name <- c(1:nedges)
V(net)$name <- c(1:nnodes)#c((nedges+1):(nnodes+nedges))
lead_lag_weights <- FALSE
alphaParams <- c(-.1,.3)
betaParams <- list(c(.1,-.2),c(-.02,.03))
gammaParams <- c(.01,.01,.03)
deltaParams <- list(c(.02,-.01),c(-.02,.01),c(.03,.1))
testsimdata <- gnar_x_sim(n=200, net, alphaParams, betaParams, gammaParams, deltaParams, sigma=1, meann=0, nedges,nnodes,n_edges_nodes,
                          data_edges,data_nodes,noise="n",vt=3,rcor=0.5)
rownames(testsimdata) <- c(paste(edgelist_net[,1],edgelist_net[,2]),seq(1,nnodes))

alphaOrder <- length(alphaParams)
betaOrder <- sapply(betaParams,length)
gammaOrder <- length(gammaParams)-1
deltaOrder <- sapply(deltaParams,length)

testsimdata_res <- gnar_x_fit(testsimdata,nnodes,nedges,n_edges_nodes,data_edges,data_nodes,alphaOrder,betaOrder,gammaOrder,deltaOrder,
                              net,lead_lag_mat,globalalpha=TRUE,lead_lag_weights=FALSE)

########### ########### ##### ##### ##### ##### ##### ##### 
##### ##### ########### REPETITIONS ########### ##### ##### 
##### ##### ##### ##### ##### ##### ########### ########### 

# Regime 1: GNAR(1,[1],1,[1])
alpha_par_1 <- .2
beta_par_1 <- list(c(0.2))
gamma_par_1 <- c(.3,.03)
delta_par_1 <- list(c(0.2),c(0.05))
# alternative for lag>=1
alpha_par_1 <- .2
beta_par_1 <- list(c(0.2))
gamma_par_1 <- c(.3)
delta_par_1 <- list(c(0.2))

# Regime 2: GNAR(1,[2],1,[2])
alpha_par_2 <- .2
beta_par_2 <- list(c(-0.2,0.1))
gamma_par_2 <-  c(.1,0.02)
delta_par_2 <- list(c(0.05,-0.2),c(0.1,-0.02))
# alternative for lag>=1
alpha_par_2 <- .2
beta_par_2 <- list(c(-0.2,0.1))
gamma_par_2 <-  c(.1)
delta_par_2 <- list(c(0.05,-0.2))

# # Regime 3: GNAR(3,[1,1,1],3,[1,1,1])
# alpha_par_3 <- c(.2,.02,-0.1)
# beta_par_3 <- list(c(0.02),c(0.05),c(-.2))
# gamma_par_3 <- c(.02,-.1,.01,0.1)
# delta_par_3 <- list(c(0.2),c(0.01),c(-.02),c(.01))

# Regime 3: GNAR(2,[2,2],2,[2,2])
alpha_par_3 <- c(-.1,.3)
beta_par_3 <- list(c(.1,-.2),c(-.02,.03))
gamma_par_3 <- c(.01,.01,.03)
delta_par_3 <- list(c(.02,-.01),c(-.02,.01),c(.03,.1))
# alternative for lag>=1
alpha_par_3 <- c(-.1,.3)
beta_par_3 <- list(c(.1,-.2),c(-.02,.03))
gamma_par_3 <- c(.01,.01)
delta_par_3 <- list(c(.02,-.01),c(-.02,.01))

# Regime 4: GNAR(2,[3,3],2,[3,3])
alpha_par_4 <- c(-.1,.1)
beta_par_4 <- list(c(.1,-.2,0.05),c(-.02,.03,-0.05))
gamma_par_4 <- c(.1,.01)
delta_par_4 <- list(c(.02,-.03,0.1),c(-.02,.01,.03))

# Regime 4: GNAR(3,[2,2,2],3,[2,2,2])
alpha_par_4 <- c(-.05,.05,.1)
beta_par_4 <- list(c(.05,-.2),c(-.02,.03),c(.1,0.05))
gamma_par_4 <- c(.1,.01,.05)
delta_par_4 <- list(c(.02,-.03),c(-.02,.01),c(.03,0.07))

# # Regime 5: GNAR(3,[2,0,0])
# alpha_par_5 <- c(.2,.4,-0.6)
# beta_par_5 <- list(c(0.3,.4),c(0),c(0))


set.seed(30)
seeds <- sample(1:500,50,replace=FALSE)
#SBM parameters
probmat <- matrix(c(.7,.2,.1,.7),nrow = 2)


alphaOrder_reg <- c(1,1,2,3)
betaOrder_reg <- list(c(1),c(2),rep(2,2),rep(2,3))
gammaOrder_reg <- c(1,1,2,3)
deltaOrder_reg <- list(c(1),c(2),c(2,2),rep(2,3))

num_param <- c(6,9,15) # for lag>=0
num_param <- c(4,6,12,18) # for lag>=1

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

# save(coverage_reg_1_mod_er,coverage_reg_1_mod_grg,coverage_reg_1_mod_sbm,coverage_reg_2_mod_er,coverage_reg_2_mod_grg,coverage_reg_2_mod_sbm,
#      coverage_reg_3_mod_er,coverage_reg_3_mod_grg,coverage_reg_3_mod_sbm,file = paste(dir_path,"/results/coverage_sim_netden0.4.RData",sep=""))
# save(rmse_reg_1_mod_er,rmse_reg_1_mod_grg,rmse_reg_1_mod_sbm,rmse_reg_2_mod_er,rmse_reg_2_mod_grg,rmse_reg_2_mod_sbm,
#      rmse_reg_3_mod_er,rmse_reg_3_mod_grg,rmse_reg_3_mod_sbm, file = paste(dir_path,"/results/rmse_estimation_sim_netden0.4.RData",sep = ""))

########### ########### ########### ########### ########### ########### ########### ########### ########### 
########### PREDICTIVE PERFORMANCE: REPETITIONS OF SBM,ER,SBM GRAPHS WITH 4TH SIM REGIME ###########
########### ########### ########### ########### ########### ########### ########### ########### ########### 

# Regime 3: GNAR(2,[2,2],2,[2,2])
#alpha_par_3 <- c(-.1,.3)
#beta_par_3 <- list(c(.1,-.2),c(-.02,.03))
#gamma_par_3 <- c(.01,.03)
#delta_par_3 <- list(c(-.02,.01),c(.03,.1))

set.seed(12) #used for 0.4 density, and 0.23 density
seeds <- sample(seq(1,1000),50)# used for 0.4 density, and 0.23 density
# set.seed(30)
# seeds <- sample(1:500,50,replace=FALSE)
#SBM parameters
probmat <- matrix(c(.7,.2,.1,.7),nrow = 2) #for 0.4 density
probmat <- matrix(c(.7,.1,.1,.7),nrow = 2)
#probmat <- matrix(c(.5,.005,.01,.5),nrow = 2) #for 0.24 density
#probmat <- matrix(c(.2,.02,.02,.2),nrow = 2)
rmse_mat <- rmse_mat_all <-  matrix(ncol = 15,nrow = 50)
colnames(rmse_mat_all) <- colnames(rmse_mat) <- c("gnarnei_er","gnar_er","ar_er","hlag_er","arima_er",
                                                  "gnarnei_sbm","gnar_sbm","ar_sbm","hlag_sbm","arima_sbm",
                                                  "gnarnei_grg","gnar_grg","ar_grg","hlag_grg","arima_grg")

for (mod in c("grg","er","sbm")){#,"er","sbm"
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
    
    # Fit GNAR(2,[2,2],2,[2,2])
    fit_train <- gnar_x_fit(datasim_train,nnodes,nedges,n_edges_nodes,data_edges,data_nodes,alphaOrder_reg[3],betaOrder_reg[[3]],
                            gammaOrder_reg[3],deltaOrder_reg[[3]],
                            net,lead_lag_mat,globalalpha=TRUE,lead_lag_weights=FALSE,pay_now = FALSE, lag_0_sep = FALSE)
    fit_predict <- gnar_x_predict(fit_train,datasim_train,alphaOrder_reg[3],betaOrder_reg[[3]],
                                  gammaOrder_reg[3],deltaOrder_reg[[3]],n_edges_nodes,nnodes,nedges,1,fit_train$wei_mat,fit_train$data_loc_mat,
                                  pay_now=FALSE,pay_data_now=NULL) # pay_data_now=simdata[1:nedges,200]
    rmse_mat[i,paste("gnarnei_",mod,sep = "")] <- sqrt(mean((fit_predict[(nedges+1):n_edges_nodes,alphaOrder_reg[3]+1]-simdata[(nedges+1):n_edges_nodes,200])^2))
    rmse_mat_all[i,paste("gnarnei_",mod,sep = "")] <- sqrt(mean((fit_predict[,alphaOrder_reg[3]+1]-simdata[,200])^2))
    
    print(c("seed: ",i," gnar nei done, model ", mod))
    
    # Fit GNAR(2,[0,0],2,[0,0])
    fit_train_nonei <- gnar_x_fit(datasim_train,nnodes,nedges,n_edges_nodes,data_edges,data_nodes,alphaOrder_reg[3],c(0,0),
                                  gammaOrder_reg[3],c(0,0),net,lead_lag_mat,globalalpha=TRUE,lead_lag_weights=FALSE,pay_now = FALSE,lag_0_sep = FALSE)
    print("done")
    fit_predict_nonei <- gnar_x_predict(fit_train_nonei,datasim_train,alphaOrder_reg[3],c(0,0),
                                        gammaOrder_reg[3],c(0,0),n_edges_nodes,nnodes,nedges,1,fit_train_nonei$wei_mat,fit_train_nonei$data_loc_mat,
                                        pay_now=FALSE,pay_data_now=NULL)
    rmse_mat[i,paste("gnar_",mod,sep = "")] <- sqrt(mean((fit_predict_nonei[(nedges+1):n_edges_nodes,alphaOrder_reg[3]+1]-simdata[(nedges+1):n_edges_nodes,200])^2))
    rmse_mat_all[i,paste("gnar_",mod,sep = "")] <- sqrt(mean((fit_predict_nonei[,alphaOrder_reg[3]+1]-simdata[,200])^2))
    
    print(c("seed: ",i," gnar nonei done, model ", mod))
    
    simtrainvar <- t(datasim_train)
    
    # # Fit VAR(1)
    # varforecast <- predict(restrict(VAR(simtrainvar,p=1,type = "none")),n.ahead=1)
    # getfcst <- function(x){return(x[,1])}
    # varfor <- unlist(lapply(varforecast$fcst, getfcst))
    # rmse_mat[i,paste("var_",mod,sep = "")] <- sqrt(mean((varfor[(nedges+1):n_edges_nodes]-simdata[(nedges+1):n_edges_nodes,200])^2))
    # rmse_mat_all[i,paste("gnar_",mod,sep = "")] <- sqrt(mean((varfor-simdata[,200])^2))
    # print(c("seed: ",i," var done"))
    
    # Fit simple AR max lag 2
    simple_ar <- apply(simtrainvar, 2, function(x){forecast(auto.arima(x,d=0,D=0,max.p = 2,max.q = 0,max.P = 0,max.Q = 0,stationary = TRUE,seasonal = FALSE,
                                                                       ic="bic",allowmean = FALSE,allowdrift = FALSE,trace = FALSE),h=1)$mean})
    rmse_mat[i,paste("ar_",mod,sep = "")] <- sqrt(mean((simple_ar[(nedges+1):n_edges_nodes]-simdata[(nedges+1):n_edges_nodes,200])^2))
    rmse_mat_all[i,paste("ar_",mod,sep = "")] <- sqrt(mean((simple_ar-simdata[,200])^2))
    print(c("seed: ",i," ar done"))
    
    # Autoarima
    simple_arima <- apply(simtrainvar, 2, function(x){forecast(auto.arima(x),h=1)$mean})
    rmse_mat[i,paste("arima_",mod,sep = "")] <- sqrt(mean((simple_arima[(nedges+1):n_edges_nodes]-simdata[(nedges+1):n_edges_nodes,200])^2))
    rmse_mat_all[i,paste("arima_",mod,sep = "")] <- sqrt(mean((simple_arima-simdata[,200])^2))
    print(c("seed: ",i," ar done"))
    
    # Fit bigvar  
    #Model1 <- constructModel(simtrainvar,p=3,struct='Basic',gran=c(50,10),verbose=FALSE) # VARX
    Model2 <- constructModel(simtrainvar,p=2,struct='HLAGC',gran=c(50,10),verbose=FALSE) # HLAG (most recent)
    #results <- cv.BigVAR(Model1)
    results2 <- cv.BigVAR(Model2)
    #varfor <- predict(results,n.ahead=1)
    varfor2 <- predict(results2,n.ahead=1)
    #rmse_mat[i,paste("varx_",j,sep = "")] <- sqrt(mean((varfor-ts[,200])^2))
    rmse_mat[i,paste("hlag_",mod,sep = "")] <- sqrt(mean((varfor2[(nedges+1):n_edges_nodes,]-simdata[(nedges+1):n_edges_nodes,200])^2))
    rmse_mat_all[i,paste("hlag_",mod,sep = "")] <- sqrt(mean((varfor2-simdata[,200])^2))
    
    print(c("seed: ",i," var done"))
  }
}

#saveRDS(rmse_mat,file = "/Users/amantziou/Documents/R_directory_Turing_ONS/collab/results/rmse_netden_0.4.rds",version = 2)
#saveRDS(rmse_mat,file = "/Users/amantziou/Documents/R_directory_Turing_ONS/collab/results/rmse_netden_0.4_lag1.rds",version = 2) # lag1 means for lag >=1 in model (i.e. not lag=0 for payments)
#saveRDS(rmse_mat,file = "/Users/amantziou/Documents/R_directory_Turing_ONS/collab/results/rmse_nodes_netden_0.4_lag1_reg2.rds",version = 2) # lag1 means for lag >=1 in model (i.e. not lag=0 for payments)
#saveRDS(rmse_mat,file = "/Users/amantziou/Documents/R_directory_Turing_ONS/collab/results/rmse_payandnodes_netden_0.4_lag1_reg2.rds",version = 2) # lag1 means for lag >=1 in model (i.e. not lag=0 for payments)
rmse_mat <- readRDS(paste(dir_path,"/results/rmse_netden_0.4.rds",sep = ""))

#save(rmse_mat,rmse_mat_all,file = paste(dir_path,"/results/rmse_netden_0.4_lag1_reg3.RData",sep="")) # lag1 means for lag >=1 in model (i.e. not lag=0 for payments), regime 3 from table with regimes for lag>=1
load(paste(dir_path,"/results/rmse_netden_0.4_lag1_reg3.RData",sep=""))

rmse_df <- as.data.frame(rmse_mat_all[,c(1,2,5,4)])
rmse_df <- as.data.frame(rmse_mat_all[,c(6,7,10,9)])
rmse_df <- as.data.frame(rmse_mat_all[,c(11,12,15,14)])

#rmse_df <- as.data.frame(rmse_mat)
colnames(rmse_df) <- rep(c("GNAR(2,[2,2],2,[2,2])","GNAR(2,[0,0],2,[0,0])","AR","HLAG","ARIMA"),1)
colnames(rmse_df) <- rep(c("GNAR-ex(2,[2,2])","GNAR-ex(2,[0,0])","ARIMA","HLAG"),1)
colnames(rmse_df) <- rep(c("GNAR(1,[2],1,[2])","GNAR(1,[0],1,[0])","AR","HLAG"),1)
rmse_melt <- melt(rmse_df)
#rmse_melt["model"] <- c(rep("GRG",200))#,rep("SBM",200))
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

# test the new data frame fun
xnewdata <- newdat(fit_train,datasim_train,alphaOrder_reg[2],betaOrder_reg[[2]],
                   gammaOrder_reg[2],deltaOrder_reg[[2]],n_edges_nodes,nnodes,nedges,1,
                   fit_train$wei_mat,fit_train$data_loc_mat)

apply(xnewdata, 1,function(x) t(matrix(rep(x,3),nrow = 6)))

xnewdata2 <- c()
for (i in 1:n_edges_nodes){
  xnewdata2 <- rbind(xnewdata2,t(matrix(rep(xnewdata[i,],198),nrow = 6)))
}

predfunres <- predict(fit_mod$mod,newdata=as.data.frame(xnewdata2),interval="prediction")
predfunres2 <- predict(fit_mod$mod,newdata=as.data.frame(xnewdata2),interval="predict")
# matrix of residuals
resmat <- (matrix(fit_mod$mod$residuals,  ncol=n_edges_nodes, byrow=FALSE))
all(t(resmat) %*% resmat==crossprod(resmat))
(t(resmat) %*% resmat)/199 # covariance matrix of residuals
apply((t(resmat) %*% resmat)/199,1,which.max)==seq(1:182) # diagonal has highest variance

########### ########### ########### ########### ########### ########### ########### ########### ########### 
########### MODEL AVERAGING PREDICTIVE PERFORMANCE: REPETITIONS OF SBM,ER,SBM GRAPHS WITH 3rd SIM REGIME ###########
########### ########### ########### ########### ########### ########### ########### ########### ########### 



set.seed(12) #used for 0.4 density, and 0.23 density
seeds <- sample(seq(1,1000),50)# used for 0.4 density, and 0.23 density
# set.seed(30)
# seeds <- sample(1:500,50,replace=FALSE)
#SBM parameters
probmat <- matrix(c(.7,.2,.1,.7),nrow = 2) #for 0.4 density
#probmat <- matrix(c(.5,.005,.01,.5),nrow = 2) #for 0.24 density
#probmat <- matrix(c(.2,.02,.02,.2),nrow = 2)
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
                                    pay_now=FALSE,pay_data_now=NULL) # pay_data_now=simdata[1:nedges,200]
      fit_pred_mat[,lagit] <- fit_predict[,lagit+1]

      # Prediction intervals 
      resmat <- matrix(fit_train$mod$residuals,  ncol=n_edges_nodes, byrow=FALSE)
      covresmat <- (t(resmat) %*% resmat)/ncol(datasim_train) # covariance matrix of residuals
      predsd_ver <- sqrt(diag(covresmat))
      
      # Get confidence intervals on pre-processed data scale
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
#saveRDS(errormat,file = paste(dir_path,"/results/ma_sim_netden_0.4_reg3.rds",sep = ""),version = 2) 
errormat <- readRDS(paste(dir_path,"/results/ma_sim_netden_0.4_reg3.rds",sep = ""))
#save(pred_incl_er,pred_incl_grg,pred_incl_sbm,file = paste(dir_path,"/results/ma_sim_netden_0.4_reg3_intervals.RData",sep = ""))

#save(rmse_mat,rmse_mat_all,file = paste(dir_path,"/results/rmse_netden_0.4_lag1_reg3.RData",sep="")) # lag1 means for lag >=1 in model (i.e. not lag=0 for payments), regime 3 from table with regimes for lag>=1

rmse_df <- as.data.frame(cbind(errormat[,1],rmse_mat[,5])) #nodeonly
rmse_df <- as.data.frame(cbind(errormat[,2],rmse_mat[,10])) #nodeonly
rmse_df <- as.data.frame(cbind(errormat[,3],rmse_mat[,15])) #nodeonly

rmse_df <- as.data.frame(cbind(errormat[,4],rmse_mat_all[,5])) #all
rmse_df <- as.data.frame(cbind(errormat[,5],rmse_mat_all[,5])) #all
rmse_df <- as.data.frame(cbind(errormat[,6],rmse_mat_all[,5])) #all

#rmse_df <- as.data.frame(rmse_mat)
#colnames(rmse_df) <- rep(c("RMSE MA nodes only","RMSE ARIMA"),1)#,"REL ERROR MA TOTAL NODE"
#colnames(rmse_df) <- rep(c("RMSE MA nodes only","RMSE MA"),1)#,"REL ERROR MA TOTAL NODE"
colnames(rmse_df) <- rep(c("MA GNAR-ex ","ARIMA"),1)#,"REL ERROR MA TOTAL NODE"
rmse_melt <- melt(rmse_df)
#rmse_melt["model"] <- c(rep("GRG",200))#,rep("SBM",200))
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
# set.seed(30)
# seeds <- sample(1:500,50,replace=FALSE)
#SBM parameters
probmat <- matrix(c(.7,.2,.1,.7),nrow = 2) #for 0.4 density
#probmat <- matrix(c(.5,.005,.01,.5),nrow = 2) #for 0.24 density
#probmat <- matrix(c(.2,.02,.02,.2),nrow = 2)
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
      
      # Get confidence intervals on pre-processed data scale
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
#save(errormat,pred_incl_er,pred_incl_grg,pred_incl_sbm,file = paste(dir_path,"/results/ma_sim_netden_0.4_reg1_intervals_errormat.RData",sep = ""))
load(paste(dir_path,"/results/ma_sim_netden_0.4_reg1_intervals_errormat.RData",sep = ""))
load(paste(dir_path,"/results/rmse_netden_0.4_lag1_reg1_autoarimaonly.RData",sep=""))

rmse_df <- as.data.frame(cbind(errormat[,1],rmse_mat[,1])) #nodeonly
rmse_df <- as.data.frame(cbind(errormat[,2],rmse_mat[,2])) #nodeonly
rmse_df <- as.data.frame(cbind(errormat[,3],rmse_mat[,3])) #nodeonly

rmse_df <- as.data.frame(cbind(errormat[,4],rmse_mat_all[,1])) #all
rmse_df <- as.data.frame(cbind(errormat[,5],rmse_mat_all[,2])) #all
rmse_df <- as.data.frame(cbind(errormat[,6],rmse_mat_all[,3])) #all

#rmse_df <- as.data.frame(rmse_mat)
#colnames(rmse_df) <- rep(c("RMSE MA nodes only","RMSE ARIMA"),1)#,"REL ERROR MA TOTAL NODE"
#colnames(rmse_df) <- rep(c("RMSE MA nodes only","RMSE MA"),1)#,"REL ERROR MA TOTAL NODE"
colnames(rmse_df) <- rep(c("MA GNAR-ex","ARIMA"),1)#,"REL ERROR MA TOTAL NODE"
rmse_melt <- melt(rmse_df)
#rmse_melt["model"] <- c(rep("GRG",200))#,rep("SBM",200))
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

############### ############### ############### ############### ############### 
############### PREDICTION PERFORMANCE OF AUTO.ARIMA FOR SIM REG 1 ############### 
############### ############### ############### ############### ############### 
rmse_mat <- rmse_mat_all <-  matrix(ncol = 3,nrow = 50)
colnames(rmse_mat_all) <- colnames(rmse_mat) <- c("arima_er","arima_sbm","arima_grg")
for (mod in c("grg","er","sbm")){#,"er","sbm"
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
#save(rmse_mat,rmse_mat_all,file = paste(dir_path,"/results/rmse_netden_0.4_lag1_reg1_autoarimaonly.RData",sep=""))
