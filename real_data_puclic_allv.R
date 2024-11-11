# PUBLIC DATA VERSION

library(dplyr)
library(zoo) # for yearmonth objects
library(stats)
library(stringr)
library(igraph)
library(reshape2)
library(ggplot2)
library(intervals)

dir_path <- "/Users/u2371456/Documents/R_directory_Turing_ONS/collab"

# load gdp all versions
load("/Users/amantziou/Documents/R_directory_Turing_ONS/collab/data/public_updated30apr_allversions/industry_panel_Pluto_ONS_monthly_CPA_public.RData") # 40 unique industries
load(paste(dir_path,"/data/public_updated13may_allversions/industry_panel_Pluto_ONS_monthly_CPA_public.RData",sep = "")) # 40 unique industries

# load updated version of Payments for CPA codes
load(paste(dir_path,"/data/public_updated18apr/IOT_flow_of_goods_Pluto_monthly_CPA_public.RData",sep = ""))

prepr_gva <- function(df,rev_mon_year){
  # df: data frame with data for all versions/revsions
  # rev_mon_year: month/year of revision in form e.g. "mar2024"
  
  # get column for GVA of corresp revision month/year
  gva <- paste("GVA sa_",rev_mon_year,sep="")
  gdp_data <- dplyr::select(df, time, code, gva) #all_of(gva)
  
  # remove NAs from gdp_data 
  gdp_data <- gdp_data[complete.cases(gdp_data),]
  
  no_ind <- length(unique(gdp_data$code))
  no_ts <- length(unique(gdp_data[,"time"])) # 321 from Jan 1997 until Sep 2023
  
  
  gdp_data_str <- matrix(,nrow = no_ind,ncol = no_ts+1)
  for(i in 1:no_ind){
    if(table(gdp_data[seq(i,nrow(gdp_data),no_ind),'code'])==no_ts){ 
      gdp_data_str[i,1] <- gdp_data[seq(i,nrow(gdp_data),no_ind),'code'][1]
      gdp_data_str[i,2:ncol(gdp_data_str)] <- gdp_data[seq(i,nrow(gdp_data),no_ind),gva]
    }else{
      warning(c("non similar codes for row ",i)) 
      gdp_data_str[i,1] <- gdp_data[seq(i,nrow(gdp_data),no_ind),'code'][1]
      gdp_data_str[i,2:ncol(gdp_data_str)] <- c(gdp_data[seq(i,nrow(gdp_data),no_ind),gva],rep(0,abs(table(gdp_data[seq(i,nrow(gdp_data),no_ind),'code'])-no_ts)))
    }
  }
  
  colnames(gdp_data_str) <- c("CPA code",format(gdp_data[seq(1,nrow(gdp_data),no_ind),]$time,"%Y-%m"))
  gdp_data_str_df <- as.data.frame(gdp_data_str)
  
  return(gdp_data_str_df)
}

prepr_pay <- function(df){
  # df: dataframe of payments data public
  
  # remove counts of transactions
  newd0 <- df[,-which(grepl("_cnt_", colnames(df), fixed = TRUE))]
  
  # public data: remove NAs industry codes in From/To columns
  newd0 <- newd0[-union(which(is.na(newd0$from)), which(is.na(newd0$to))),]
  
  return(newd0)  
}

prepr_comb <- function(df_gva,df_pay,date_start,date_end,date_end_2,all_nas_row=FALSE,atleast1_na=TRUE,thres){
  # df_gva: dataframe with GDP
  # df_pay: dataframe with payments
  # date_start: start date in payments in form "year-month"
  # date_end: end date in gdp in form "year-month"
  # date_end_2: end date in gdp in form "sep_amt_2023"
  # all_nas_row: remove rows with all NAs
  # atleast1_na: remove rows with at least one NA
  # thres: threshold number of 0 to keep per row in payments 
  
  # codes not in common
  notin_codes_to <- setdiff(unique(df_gva$`CPA code`),unique(df_pay$to))
  notin_codes_from <- setdiff(unique(df_gva$`CPA code`),unique(df_pay$from))
  notin_codes_all <- unique(union(notin_codes_from,notin_codes_to))
  
  # if there are, remove codes not in common
  if (length(notin_codes_all)!=0){
    if (length(which(df_gva$`CPA code` %in% notin_codes_all))!=0){
      df_gva <- df_gva[-which(df_gva$`CPA code` %in% notin_codes_all),] 
    }
    if ((length(which(df_pay$from %in% notin_codes_all))+length(which(df_pay$from %in% notin_codes_all)))!=0){
      df_pay <- df_pay[-which(df_pay$from %in% notin_codes_all),] 
      df_pay <- df_pay[-which(df_pay$to %in% notin_codes_all),] 
    }
  }
  
  
  # keep the columns/yearmonths of GDP corresponding to yearmonth columns in Payments data
  starttime <- which(names(df_gva)==date_start) # start time of payments ts ("2015-08")
  endtime <- which(names(df_gva)==date_end) # end time of payments ts is Oct 2023 but GDP data up to including Sep 2023
  df_gva <- df_gva[,c(1,starttime:endtime)]
  
  # remove from payments data the time stamps for which no GDP (i.e. after Oct 2023)
  df_pay <- df_pay[,1:which(names(df_pay)==date_end_2)]
  
  # align names of timestamps between GDP ts and PAY ts
  colnames(df_pay)[-c(1,2)] <- colnames(df_gva[-c(1)])
  
  ########### NAs handling ########### 
  
  # for public data: check which all NAs in rows
  if (all_nas_row){
    row_nas <- c()
    for (i in 1:nrow(df_pay)){
      if (all(is.na(df_pay[i,3:ncol(df_pay)]))){
        row_nas <- c(row_nas,i)
      }
    }
  
    # for public data: remove rows with all NAs
    df_pay <- df_pay[-row_nas,]
    # assign 0 to NAs for the rest
    df_pay[is.na(df_pay)] <- 0
    # remove ts with more than thres zeros
    df_pay <- df_pay[which(rowSums(df_pay[,-c(1,2)]==0)<thres),]
  }else if (atleast1_na){
    # remove all rows with at least one NA
    df_pay <- df_pay[complete.cases(df_pay[,3:ncol(df_pay)]),]
  }
  
  
  # remaining CPA codes after removal of transactions
  CPA_codes_pub <- unique(union(unique(df_pay$to),unique(df_pay$from)))
  
  # check which is not anymore and remove from GDP ts
  if (!all(c(df_gva$`CPA code`%in%CPA_codes_pub,CPA_codes_pub%in%df_gva$`CPA code`))){
    df_gva <- df_gva[-which(!(df_gva$`CPA code`%in%CPA_codes_pub)),]
  }
  
  # create numbering for CPA codes
  CPA_node <- cbind(sort(CPA_codes_pub),seq(1:length(CPA_codes_pub)))
  df_pay$to <- CPA_node[match(df_pay$to,CPA_node[,1]),2]
  df_pay$from <- CPA_node[match(df_pay$from,CPA_node[,1]),2]
  df_gva[,1] <- CPA_node[match(df_gva[,1],CPA_node[,1]),2]
  
  ########## create data frame for model input #########
  
  ts_pay_gdp <- rbind(df_pay[,-c(1,2)],df_gva[,-c(1)])
  ts_pay_gdp <- apply(ts_pay_gdp,2,as.numeric)
  rownames(ts_pay_gdp) <- c(paste(df_pay[,1],df_pay[,2]),df_gva[,1])
  
  ######### create graph #########
  
  g <- graph_from_edgelist(as.matrix(df_pay[,1:2]),directed = TRUE)
  #ecount(simplify(g,remove.loops = FALSE))==ecount(g) # check if non multiple edges
  nnodes <- vcount(g)
  nedges <- ecount(g)
  n_edges_nodes <- nnodes+nedges
  data_edges <- df_pay[,c(1,2)]
  data_nodes <- df_gva[,1]
  
  return(list(df_pay_2=df_pay,df_gva_2=df_gva,ts_pay_gdp_2=ts_pay_gdp,CPA_node_2=CPA_node,graph_2=g,nnodes_2=nnodes,
              nedges_2=nedges,
         n_edges_nodes_2=n_edges_nodes,data_edges_2=data_edges,data_nodes_2=data_nodes))
}




prepr_stl <- function(data,start_train_mon,end_train_mon,start_train_year,end_train_year,nts){
  # data: combined edge/nodal time series
  # start_train_mon: start month of time series training sample
  # end_train_mon: end month of time series training sample
  # start_train_year: start year of time series training sample
  # end_train_year: end year of time series training sample
  # nts: number of time stamps to predict ahead
  
  stl_list <- apply(data,1,function(x) stl(ts(as.vector(x),start = c(start_train_year,start_train_mon), end=c(end_train_year,end_train_mon),frequency = 12),s.window = "periodic"))
  res <- t(sapply(stl_list, function(x) x$time.series[,3]))
  seas <- t(sapply(stl_list, function(x) x$time.series[,1]))
  trend <- t(sapply(stl_list, function(x) x$time.series[,2]))
  colnames(res) <- colnames(data)[1:ncol(res)]
  rownames(res) <- rownames(data)
  # get seas component for next time stamp/ time stamps
  seas_end_mon <- ifelse(end_train_mon==12,1,as.numeric(end_train_mon)+nts)
  seas_comp_nts <- t(sapply(stl_list, function(x) x$time.series[,1][seas_end_mon]))
  # get pred with ar for trend
  x_trend <- seq(1,ncol(trend))
  fitlist <- apply(trend, 1, function(x) lm(x~x_trend))
  pred_trend <- sapply(fitlist, function(x) predict(x,data.frame(x_trend=ncol(res)+nts),type="response"))
  return(list(data_train_stl=res,seasonal=seas,trend=trend,seas_comp_nts=seas_comp_nts,pred_trend=pred_trend))
}

growthrates <- function(ts_data){ # (ending-starting)/starting
  new_ts <- apply(ts_data,1,diff)
  new_ts <- t(new_ts)
  new_ts <- new_ts/ts_data[,1:(ncol(ts_data)-1)]
  new_ts <- cbind(rep(1,nrow(new_ts)),new_ts)
  colnames(new_ts) <- colnames(ts_data)
  return(new_ts)
}

paydata <- prepr_pay(d0)
rev_mon_year <- c("dec2021","mar2022","jun2022","sep2022","dec2022","mar2023","jun2023","sep2023","dec2023","mar2024")
date_end_2_list <- c("jun_amt_2021","sep_amt_2021","dec_amt_2021","jun_amt_2022","jun_amt_2022","sep_amt_2022","dec_amt_2022","mar_amt_2023","jun_amt_2023",
                     "sep_amt_2023")
gdpdata <- prepr_gva(dt,rev_mon_year[1])
prep_data <- prepr_comb(gdpdata,paydata,"2016-01",tail(colnames(gdpdata),n=1),date_end_2_list[1],all_nas_row=FALSE,atleast1_na=TRUE,NULL)
data <- prep_data$ts_pay_gdp
stl_list <- apply(data,1,function(x) stl(ts(as.vector(x),start = c(2016,1), end=c(2021,6),frequency = 12),s.window = "periodic"))


for (rev in 1:10){
  gdpdata <- prepr_gva(dt,rev_mon_year[rev])
  #print(tail(colnames(gdpdata),n=1))
  assign(paste("prep_data_rev_",rev,sep=""),prepr_comb(gdpdata,paydata,"2016-01",tail(colnames(gdpdata),n=1),date_end_2_list[rev],all_nas_row=FALSE,atleast1_na=TRUE,NULL))
  # gdpdata2 <- prepr_gva(dt,rev_mon_year[rev+1])
  # #print(tail(colnames(gdpdata),n=1))
  # prep_data2 <- prepr_comb(gdpdata2,paydata,"2016-01",tail(colnames(gdpdata2),n=1),date_end_2_list[rev+1],all_nas_row=FALSE,atleast1_na=TRUE,NULL)
  # #print(ecount(prep_data$graph))
  # print(all(E(prep_data$graph)==E(prep_data2$graph)))
  assign(paste("graph_rev_",rev,sep=""),get(paste("prep_data_rev_",rev,sep = ""))$graph)
  assign(paste("edg_rev_",rev,sep=""),get(paste("prep_data_rev_",rev,sep = ""))$data_edges)
}

save(prep_data_rev_1,prep_data_rev_2,prep_data_rev_3,prep_data_rev_4,prep_data_rev_5,prep_data_rev_6,
     prep_data_rev_7,prep_data_rev_8,prep_data_rev_9,prep_data_rev_10,file="/Users/amantziou/Documents/R_directory_Turing_ONS/collab/data/data_preprocessing/data_allver.RData")

# # Rev 1,2,3:same graph
# # Rev 4,5,6,7: same graph
# # Rev 8,9,10: same graph

# find common graph among all versions: intersection of edges
graph_allv <- graph_rev_1 %s% graph_rev_2 %s% graph_rev_3 %s% graph_rev_4 %s% graph_rev_5 %s% graph_rev_6 %s%
              graph_rev_7 %s% graph_rev_8 %s% graph_rev_9 %s% graph_rev_10

# # union of edges not in common among all
# edge_rem <- E(graph_rev_3)[!(E(graph_rev_3) %in% E(graph_allv))]
# edge_pos_rem <- which(!(E(graph_rev_3) %in% E(graph_allv)))

# ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
# # reform datasets to all share the same network structure #
# ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
# for (i in 1:10){
#   if (!all(E(get(paste("graph_rev_",i,sep = ""))) %in% E(graph_allv))){ # check if there are not in common edges
#     edge_pos_rem <- which(!(E(get(paste("graph_rev_",i,sep = ""))) %in% E(graph_allv)))
#     edge_rem <- E(get(paste("graph_rev_",i,sep = "")))[!(E(get(paste("graph_rev_",i,sep = ""))) %in% E(graph_allv))]
#     if (all(matrix(as.numeric(unlist(strsplit(as_ids(edge_rem),split="\\|"))),ncol = 2,byrow = TRUE)==get(paste("prep_data_rev_",i,sep=""))$data_edges[edge_pos_rem,])){ # confirm that positions of edges in graph same as position in edge list, thus in df
#       df_pay2 <- get(paste("prep_data_rev_",i,sep = ""))$df_pay[-edge_pos_rem,]
#       data_edges2 <-  get(paste("prep_data_rev_",i,sep = ""))$data_edges[-edge_pos_rem,]
#       ts_pay_gdp2 <-  get(paste("prep_data_rev_",i,sep = ""))$ts_pay_gdp[-edge_pos_rem,]
#       #graph2 <-  graph_from_edgelist(as.matrix(data_edges2),directed = TRUE)
#       graph2 <-  graph_from_data_frame(as.data.frame(data_edges2),directed = TRUE,vertices = get(paste("prep_data_rev_",i,sep = ""))$data_nodes)
#       nedges2 <-  ecount(graph2)
#       n_edges_nodes2 <-  ecount(graph2)+vcount(graph2)
#       assign(paste("prep_data_rev_",i,"uni",sep = ""),list(df_pay=df_pay2,df_gva=get(paste("prep_data_rev_",i,sep = ""))$df_gva,ts_pay_gdp=ts_pay_gdp2,CPA_node=get(paste("prep_data_rev_",i,sep = ""))$CPA_node,graph=graph2,nnodes=get(paste("prep_data_rev_",i,sep = ""))$nnodes,nedges=nedges2,
#            n_edges_nodes=n_edges_nodes2,data_edges=data_edges2,data_nodes=get(paste("prep_data_rev_",i,sep = ""))$data_nodes))
#     }else{
#       warning("Not same positioning between edges in graph and edges in dataframe")
#     }
#   }else{
#     assign(paste("prep_data_rev_",i,"uni",sep = ""),get(paste("prep_data_rev_",i,sep = "")))
#     print(c("Version ",i," does not require changes"))
#   }
# }
# 
# # save(prep_data_rev_1uni,prep_data_rev_2uni,prep_data_rev_3uni,prep_data_rev_4uni,prep_data_rev_5uni,prep_data_rev_6uni,
# #      prep_data_rev_7uni,prep_data_rev_8uni,prep_data_rev_9uni,prep_data_rev_10uni,file="/Users/amantziou/Documents/R_directory_Turing_ONS/collab/data/data_preprocessing/data_allver_alignednet.RData")
# load("/Users/amantziou/Documents/R_directory_Turing_ONS/collab/data/data_preprocessing/data_allver_alignednet.RData")

ind_stag <- 1
nodes_only <- FALSE
gr <- TRUE
rmse_mat_real <- lapply(1:length(1:9),function(x) matrix(nrow=9,ncol = 4))
rmse_mat_real_nodes <- lapply(1:length(1:9),function(x) matrix(nrow=9,ncol = 4))
bicres4 <- bicres5 <-  aicres <- lapply(1:length(1:9),function(x) matrix(nrow=9,ncol = 4))
pred_ver <- lapply(1:length(1:9), function(x) vector("list",9))
predsd_ver <- lapply(1:length(1:9), function(x) vector("list",9))
covresmat_list <- resmat_list <- lapply(1:length(1:9), function(x) vector("list",9))
pred_trend_list <- lapply(1:9, function(x) list())
seas_comp_list <- lapply(1:9, function(x) list())
data_vec_list <- lapply(1:9, function(x) list())
upper_pre <- lower_pre <- upper_or <- lower_or <- lapply(1:length(1:9), function(x) vector("list",9))
pred_incl_pre <- pred_incl_or <- lapply(1:length(1:9),function(x) matrix(nrow = 9,ncol = 4))
for (rev in 9){
  # data 
  # gdpdata <- prepr_gva(dt,rev_mon_year[rev])
  # #print(tail(colnames(gdpdata),n=1))
  # prep_data <- prepr_comb(gdpdata,paydata,"2016-01",tail(colnames(gdpdata),n=1),date_end_2_list[rev],all_nas_row=FALSE,atleast1_na=TRUE,NULL)
  # data <- prep_data$ts_pay_gdp[-edge_pos_rem,]
  #prep_data <- get(paste("sparse_",rev,"_thres0.3",sep = ""))
  #prep_data <- get(paste("prep_data_rev_",rev,"uni",sep = ""))
  prep_data <- get(paste("sparse_",rev,"_nodedgerem_thres0.4",sep = ""))
  data <- prep_data$ts_pay_gdp_2[,1:48] # up to 48 for pre-covid
  # data <- log(data)
  # 
  # # stationary
  # data <- stationary_ts(data)
  # 
  # growth rates
  data <- growthrates(data)
  
  # stl
  res_stl <- prepr_stl(data,str_sub(colnames(data)[1],6,7),str_sub(tail(colnames(data),n=1),6,7),str_sub(colnames(data)[1],1,4),str_sub(tail(colnames(data),n=1),1,4),nts=1)
  data_train_stl <- res_stl$data_train_stl
  pred_trend_list[[rev]] <- res_stl$pred_trend
  seas_comp_list[[rev]] <- res_stl$seas_comp_nts
  simtrainvar <- t(data_train_stl)
  
  
  for (lagi in 1:9){
    ind_stag <- 1
    for (stag in 0:3){
      print(c("data rev",rev,"Lag iter",lagi," stag iter ",stag))
      alphaOrder <- lagi
      betaOrder <- rep(stag,lagi)
      gammaOrder <- lagi
      deltaOrder <- rep(stag,lagi)

      # fit GNAR-x
      fit_train <- gnar_x_fit(data_train_stl,prep_data$nnodes_2,prep_data$nedges_2,prep_data$n_edges_nodes_2,prep_data$data_edges_2,prep_data$data_nodes_2,
                              alphaOrder,betaOrder,gammaOrder,deltaOrder,prep_data$graph_2,lead_lag_mat,globalalpha=TRUE,lead_lag_weights=FALSE,pay_now = FALSE,lag_0_sep = FALSE
                              )#nodes_only = nodes_only)
      #print(summary(fit_train$mod)$coefficients[,4] < 0.05)
      # Predict
      pred <- gnar_x_predict(fit_train,data_train_stl,alphaOrder,betaOrder,
                             gammaOrder,deltaOrder,prep_data$n_edges_nodes_2,prep_data$nnodes_2,prep_data$nedges_2,1,fit_train$wei_mat,fit_train$data_loc_mat,
                             pay_now=FALSE,pay_data_now=NULL)

      pred_adj <- pred[,lagi+1]+res_stl$pred_trend+res_stl$seas_comp_nts

      if (gr==TRUE){
        pred_adj <- prep_data$ts_pay_gdp_2[,ncol(data)]*(1+pred_adj)
      }
      
      if (nodes_only==TRUE){
        pred_adj <- pred_adj[(prep_data$nedges_2+1):prep_data$n_edges_nodes_2]
      }
      
      pred_ver[[rev]][[lagi]][[ind_stag]] <- pred_adj
      data2 <- get(paste("sparse_",rev+1,"_nodedgerem_thres0.4",sep = ""))$ts_pay_gdp_2
      # check that two consecutive data revisions, same graph
      ifelse(all(rownames(data2)==rownames(data)),print("same graph"),warning(paste("graph not the same between rev ",rev," and rev ",rev+1,sep="")))

      if (rev!=4){
        data_vec <- data2[,ncol(data)+1]
        if (gr==TRUE){
          data_vec_pre <- (prep_data$ts_pay_gdp_2[,ncol(data)]-data_vec)/prep_data$ts_pay_gdp_2[,ncol(data)]
        }
      }else{
        data2 <- get(paste("sparse_",rev+2,"_nodedgerem_thres0.4",sep = ""))$ts_pay_gdp_2
        data_vec <- data2[,ncol(data)+1]
        if (gr==TRUE){
          data_vec_pre <- (prep_data$ts_pay_gdp_2[,ncol(data)]-data_vec)/prep_data$ts_pay_gdp_2[,ncol(data)]
        }
      }
      
      if(nodes_only==TRUE){
        data_vec <- data_vec[(prep_data$nedges_2+1):prep_data$n_edges_nodes_2]
      }
      data_vec_list[[rev]] <- data_vec
      rmse_mat_real[[rev]][lagi,ind_stag] <- sqrt(mean((pred_adj-data_vec)^2))
      rmse_mat_real_nodes[[rev]][lagi,ind_stag] <- sqrt(mean((tail(as.vector(pred_adj),n=prep_data$nnodes_2)-tail(data_vec,n=prep_data$nnodes_2))^2))
      
      if (nodes_only){
        bicres4[[rev]][lagi,ind_stag] <- gnarxbic(fit_train,ncol(data_train_stl),globalalpha=TRUE,makepd = TRUE,nodes_only = TRUE,no_nodes = prep_data$nnodes_2)
        bicres5[[rev]][lagi,ind_stag] <- gnarxbic(fit_train,ncol(data_train_stl),globalalpha=TRUE,makepd = FALSE,nodes_only = TRUE,no_nodes = prep_data$nnodes_2)
      }else{
        bicres4[[rev]][lagi,ind_stag] <- gnarxbic(fit_train,ncol(data_train_stl),globalalpha=TRUE,makepd = TRUE)
        bicres5[[rev]][lagi,ind_stag] <- gnarxbic(fit_train,ncol(data_train_stl),globalalpha=TRUE,makepd = FALSE)
      }
      
      #aicres[[rev]][lagi,ind_stag] <- gnarxaic(fit_train,ncol(data_train_stl),globalalpha=TRUE)
      # get standard error of prediction to calculate prediction interval
      
      if (nodes_only==TRUE){
        resmat <- matrix(fit_train$mod$residuals,  ncol=prep_data$nnodes_2, byrow=FALSE)
      }else{
        resmat <- matrix(fit_train$mod$residuals,  ncol=prep_data$n_edges_nodes_2, byrow=FALSE)
      }
      
      covresmat <- (t(resmat) %*% resmat)/ncol(data_train_stl) # covariance matrix of residuals
      covresmat_list[[rev]][[lagi]][[ind_stag]] <- covresmat
      resmat_list[[rev]][[lagi]][[ind_stag]] <- resmat
      predsd_ver[[rev]][[lagi]][[ind_stag]] <- sqrt(diag(covresmat))
      
      # Get confidence intervals on pre-processed data scale
      upper_pre[[rev]][[lagi]][[ind_stag]] <- tail(pred[,lagi+1],n=prep_data$nnodes_2)+1.96*tail(as.vector(predsd_ver[[rev]][[lagi]][[ind_stag]]),n=prep_data$nnodes_2)
      lower_pre[[rev]][[lagi]][[ind_stag]]  <- tail(pred[,lagi+1],n=prep_data$nnodes_2)-1.96*tail(as.vector(predsd_ver[[rev]][[lagi]][[ind_stag]]),n=prep_data$nnodes_2)
      pred_incl_pre[[rev]][lagi,ind_stag] <- length(which(between(tail(data_vec_pre,n=prep_data$nnodes_2),lower_pre[[rev]][[lagi]][[ind_stag]],
                                                    upper_pre[[rev]][[lagi]][[ind_stag]])))/prep_data$nnodes_2
      
      upper <- upper_pre[[rev]][[lagi]][[ind_stag]]+tail(res_stl$pred_trend,n=prep_data$nnodes_2)+res_stl$seas_comp_nts[,(prep_data$nedges_2+1):prep_data$n_edges_nodes_2]
      lower <- lower_pre[[rev]][[lagi]][[ind_stag]]+tail(res_stl$pred_trend,n=prep_data$nnodes_2)+res_stl$seas_comp_nts[,(prep_data$nedges_2+1):prep_data$n_edges_nodes_2]
      
      if (gr==TRUE){
        upper <- tail(prep_data$ts_pay_gdp_2[,ncol(prep_data$ts_pay_gdp_2)],n=prep_data$nnodes_2)*(1+upper)
        lower <- tail(prep_data$ts_pay_gdp_2[,ncol(prep_data$ts_pay_gdp_2)],n=prep_data$nnodes_2)*(1+lower)
      }
      
      upper_or[[rev]][[lagi]][[ind_stag]] <- upper
      lower_or[[rev]][[lagi]][[ind_stag]] <- lower
      pred_incl_or[[rev]][lagi,ind_stag] <- length(which(between(tail(data_vec,n=prep_data$nnodes_2),lower,upper)))/prep_data$nnodes_2
      ind_stag <- ind_stag+1
    }
    simple_ar <- apply(simtrainvar, 2, function(x){forecast(auto.arima(x,d=0,D=0,max.p = lagi,max.q = 0,max.P = 0,max.Q = 0,stationary = TRUE,seasonal = FALSE,
                                                                       ic="bic",allowmean = FALSE,allowdrift = FALSE,trace = FALSE),h=1)$mean})
    pred_adj_ar <- simple_ar+pred_trend_list[[rev]]+seas_comp_list[[rev]]
    if (gr==TRUE){
      pred_adj_ar <- prep_data$ts_pay_gdp_2[,ncol(data)]*(1+pred_adj_ar)
    }
    if(nodes_only==TRUE){
      pred_adj_ar <- pred_adj_ar[(prep_data$nedges_2+1):prep_data$n_edges_nodes_2]
    }
    pred_ver[[rev]][[lagi]][[5]] <- pred_adj_ar
  }
}

fit_arima <- lapply(1:9, function(x) list())
ord <- lapply(1:9, function(x) list())
pred_arima <- lapply(1:9, function(x) list())
for (rev in 1:9){
  print(rev)
  # data 
  # gdpdata <- prepr_gva(dt,rev_mon_year[rev])
  # #print(tail(colnames(gdpdata),n=1))
  # prep_data <- prepr_comb(gdpdata,paydata,"2016-01",tail(colnames(gdpdata),n=1),date_end_2_list[rev],all_nas_row=FALSE,atleast1_na=TRUE,NULL)
  # data <- prep_data$ts_pay_gdp[-edge_pos_rem,]
  #prep_data <- get(paste("prep_data_rev_",rev,"uni",sep = ""))
  prep_data <- get(paste("sparse_",rev,"_nodedgerem_thres0.4",sep = ""))
  data <- prep_data$ts_pay_gdp_2
  # data <- log(data)
  # 
  # # stationary
  # data <- stationary_ts(data)
  # 
  # growth rates
  #data <- growthrates(data)
  
  # # stl
  # res_stl <- prepr_stl(data,str_sub(colnames(data)[1],6,7),str_sub(tail(colnames(data),n=1),6,7),str_sub(colnames(data)[1],1,4),str_sub(tail(colnames(data),n=1),1,4),nts=1)
  # data_train_stl <- res_stl$data_train_stl
  # pred_trend_list[[rev]] <- res_stl$pred_trend
  # seas_comp_list[[rev]] <- res_stl$seas_comp_nts
  # simtrainvar <- t(data_train_stl)
  simtrainvar <- t(data)
  
  #fit_arima[[rev]] <- apply(simtrainvar,2,function(x){auto.arima(x)})
  #fit_arima[[rev]] <- apply(simtrainvar, 2, function(x){auto.arima(x,d=1,D=0,max.p = 0,max.q = 0,max.P = 0,max.Q = 0,max.d = 1,max.D = 1)})
  fit_arima[[rev]] <- apply(simtrainvar, 2, function(x){Arima(x,order=c(0,1,0))})
  ord[[rev]] <- do.call(rbind,lapply(fit_arima[[rev]],function(x) arimaorder(x)))

  pred_arima[[rev]] <- do.call(rbind,lapply(fit_arima[[rev]], function(x) forecast(x,h=1)$mean))
  # if (gr==TRUE){
  #   pred_arima[[rev]] <- prep_data$ts_pay_gdp_2[,ncol(prep_data$ts_pay_gdp_2)]*(1+pred_arima[[rev]])
  # }
}

# save(rmse_mat_real,
#      rmse_mat_real_nodes ,
#      pred_ver, 
#      predsd_ver,
#      pred_trend_list,
#      seas_comp_list,
#      data_vec_list,file = "/Users/amantziou/Documents/R_directory_Turing_ONS/collab/results/res_public_updated13may_allversions.RData")

# save(rmse_mat_real,
#      rmse_mat_real_nodes ,
#      bicres4,
#      pred_ver,
#      predsd_ver,
#      covresmat_list,
#      resmat_list,
#      pred_trend_list,
#      seas_comp_list,
#      data_vec_list,file = "/Users/amantziou/Documents/R_directory_Turing_ONS/collab/results/res_public_updated13may_allversions_stlonly.RData")
load("/Users/amantziou/Documents/R_directory_Turing_ONS/collab/results/res_public_updated13may_allversions_stlonly.RData")


# save(rmse_mat_real,
#      rmse_mat_real_nodes ,
#      bicres4,
#      pred_ver,
#      predsd_ver,
#      covresmat_list,
#      resmat_list,
#      pred_trend_list,
#      seas_comp_list,
#      data_vec_list,file = "/Users/amantziou/Documents/R_directory_Turing_ONS/collab/results/res_public_updated13may_allversions_growthstl.RData")


# # below revised results for grwoth rates and stl and evaluating predictions on original scale
# save(rmse_mat_real,
#      rmse_mat_real_nodes ,
#      bicres4,
#      pred_ver,
#      predsd_ver,
#      covresmat_list,
#      resmat_list,
#      pred_trend_list,
#      seas_comp_list,
#      data_vec_list,file = "/Users/amantziou/Documents/R_directory_Turing_ONS/collab/results/res_public_updated13may_allversions_growthstl_rev.RData")
load("/Users/amantziou/Documents/R_directory_Turing_ONS/collab/results/res_public_updated13may_allversions_growthstl_rev.RData")

# below results from auto.arima on raw data (no stl, no grwoth rates)
# save(fit_arima,ord,pred_arima,file = "/Users/amantziou/Documents/R_directory_Turing_ONS/collab/results/res_public_updated13may_allversions_arima_rawdata.RData")
load("/Users/amantziou/Documents/R_directory_Turing_ONS/collab/results/res_public_updated13may_allversions_arima_rawdata.RData")

# REVISED data preprocessing and revised results for growth and stl, evaulating prediction on original scale 
# save(rmse_mat_real,
#      rmse_mat_real_nodes ,
#      bicres4,
#      pred_ver,
#      predsd_ver,
#      covresmat_list,
#      resmat_list,
#      pred_trend_list,
#      seas_comp_list,
#      data_vec_list,file = "/Users/amantziou/Documents/R_directory_Turing_ONS/collab/results/revised_preprocessing/res_public_updated13may_allversions_growthstl_rev_prepr_rev.RData")

# REVISED data preprocessing and revised results results from auto.arima on raw data (no stl, no grwoth rates)
# NOTE: OBVIOUSLY FOR NODES SAME RESULTS AS THIS REVISION INVOLVES CHANGES IN NETWORK STRUCTURE NOT ON NODES, ARIMA FITTED INDIVIDUALLY USING NO NETWORK INFO
# save(fit_arima,ord,pred_arima,file = "/Users/amantziou/Documents/R_directory_Turing_ONS/collab/results/revised_preprocessing/res_public_updated13may_allversions_arima_rawdata_prepr_rev.RData")
load( "/Users/amantziou/Documents/R_directory_Turing_ONS/collab/results/revised_preprocessing/res_public_updated13may_allversions_arima_rawdata_prepr_rev.RData")

# # below  results for SPARSE network (pearson 0.4), grwoth rates and stl and evaluating predictions on original scale
# save(rmse_mat_real,
#      rmse_mat_real_nodes ,
#      bicres4,
#      pred_ver,
#      predsd_ver,
#      covresmat_list,
#      resmat_list,
#      pred_trend_list,
#      seas_comp_list,
#      data_vec_list,file = "/Users/amantziou/Documents/R_directory_Turing_ONS/collab/results/res_public_updated13may_allversions_growthstl_sparse_0.4pearson.RData")

# # below  results for SPARSE network (pearson 0.4), grwoth rates and stl and evaluating predictions on original scale, SUBCASE OF GNAR-ex (ONLY NODAL EQUATION 1)
# save(rmse_mat_real,
#      rmse_mat_real_nodes ,
#      bicres4,bicres5,
#      pred_ver,
#      predsd_ver,
#      covresmat_list,
#      resmat_list,
#      pred_trend_list,
#      seas_comp_list,
#      data_vec_list,file = "/Users/amantziou/Documents/R_directory_Turing_ONS/collab/results/res_public_updated13may_allversions_growthstl_sparse_0.4pearson_eq1nodes.RData")

# # below  results for SPARSE network (pearson 0.3), grwoth rates and stl and evaluating predictions on original scale
# save(rmse_mat_real,
#      rmse_mat_real_nodes ,
#      bicres4,bicres5,
#      pred_ver,
#      predsd_ver,
#      covresmat_list,
#      resmat_list,
#      pred_trend_list,
#      seas_comp_list,
#      data_vec_list,file = "/Users/amantziou/Documents/R_directory_Turing_ONS/collab/results/res_public_updated13may_allversions_growthstl_sparse_0.3pearson.RData")


# # below results from auto.arima on raw data (no stl, no grwoth rates) FOR SPARSE network (pearson 0.4) i.e. 37 NODES (node 36 excluded)
# save(fit_arima,ord,pred_arima,file = "/Users/amantziou/Documents/R_directory_Turing_ONS/collab/results/res_public_updated13may_allversions_arima_rawdata_sparse_0.4pearson.RData")

########################################################################################################
############# OBTAIN PROPORTION OF TRUE VALUES IN PREDICTIVE INTERVALS AMONG NODES #############
########################################################################################################

upper <- lower <- lapply(1:length(1:9), function(x) vector("list",9))
pred_incl <- lapply(1:length(1:9),function(x) matrix(nrow = 9,ncol = 4))
for (rev in 1:9){
  gdp_sum <- sum(tail(data_vec_list[[rev]],32))
  for (i in 1:9){
    for (j in 1:4){
      upper[[rev]][[i]][[j]] <- tail(as.vector(pred_ver[[rev]][[i]][[j]]),n=32)+1.96*tail(as.vector(predsd_ver[[rev]][[i]][[j]]),n=32)
      lower[[rev]][[i]][[j]]  <- tail(as.vector(pred_ver[[rev]][[i]][[j]]),n=32)-1.96*tail(as.vector(predsd_ver[[rev]][[i]][[j]]),n=32)
      pred_incl[[rev]][i,j] <- length(which(between(tail(data_vec_list[[rev]],n=32),tail(as.vector(pred_ver[[rev]][[i]][[j]]),n=32)-1.96*tail(as.vector(predsd_ver[[rev]][[i]][[j]]),n=32),
                                                    tail(as.vector(pred_ver[[rev]][[i]][[j]]),n=32)+1.96*tail(as.vector(predsd_ver[[rev]][[i]][[j]]),n=32))))/32
    }
  }
}

lapply(pred_incl,function(x) apply(x,1,which.max))

# # get the true gdp/payments values from subsequent release
# data_vec_list <- lapply(1:9, function(x) list())
# for (rev in 1:9){
#   # data 
#   data <- get(paste("prep_data_rev_",rev,"uni",sep = ""))$ts_pay_gdp  
#   
#   data2 <- get(paste("prep_data_rev_",rev+1,"uni",sep = ""))$ts_pay_gdp
#   # check that two consecutive data revisions, same graph
#   ifelse(all(rownames(data2)==rownames(data)),print("same graph"),warning(paste("graph not the same between rev ",rev," and rev ",rev+1,sep="")))
#   
#   data_vec_list[[rev]] <- data2[,ncol(data)+1]
# }

# ind_rev <- 1
# for (rev in 8:9){
#   # data
#   gdpdata <- prepr_gva(dt,rev_mon_year[rev])
#   #print(tail(colnames(gdpdata),n=1))
#   prep_data <- prepr_comb(gdpdata,paydata,"2016-01",tail(colnames(gdpdata),n=1),date_end_2_list[rev],all_nas_row=FALSE,atleast1_na=TRUE,NULL)
#   data <- prep_data$ts_pay_gdp
#   # stl
#   res_stl <- prepr_stl(data,str_sub(colnames(data)[1],6,7),str_sub(tail(colnames(data),n=1),6,7),str_sub(colnames(data)[1],1,4),str_sub(tail(colnames(data),n=1),1,4),nts=1)
#   data_train_stl <- res_stl$data_train_stl
#   simtrainvar <- t(data_train_stl)
#   for (i in 1:9){
#     simple_ar <- apply(simtrainvar, 2, function(x){forecast(auto.arima(x,d=0,D=0,max.p = i,max.q = 0,max.P = 0,max.Q = 0,stationary = TRUE,seasonal = FALSE,
#                                                                        ic="bic",allowmean = FALSE,allowdrift = FALSE,trace = FALSE),h=1)$mean})
#     pred_ver[[ind_rev]][[i]][[5]] <- simple_ar+pred_trend_list[[ind_rev]]+seas_comp_list[[ind_rev]]
#   }
#   ind_rev <- ind_rev+1
# }

########################################################################################
############# PLOT MEAN ABSOLUTE PERCENTAGE ERROR OF NODAL GDP PER RELEASE #############
########################################################################################

mape_nodes <- lapply(1:length(1:9),function(x) matrix(nrow=9,ncol = 5))

for (rev in 1:9){
  for (lagi in 1:9){
    for (stag in 1:5){
      predval <- tail(as.vector(pred_ver[[rev]][[lagi]][[stag]]),n=38)
      datval <- tail(as.vector(data_vec_list[[rev]]),n=38)
      mape_nodes[[rev]][lagi,stag] <- 100*(mean(abs((datval-predval))/datval))
    }
  }
}

##############################################################################
############# PLOT RELATIVE ERROR OF TOTAL GDP PER RELEASE #############
##############################################################################

gdp_sum_mat <- lapply(1:length(1:9),function(x) matrix(nrow = 9,ncol = 4))
gdp_sum_rel_err <- lapply(1:length(1:9),function(x) matrix(nrow = 9,ncol = 4))
for (rev in 1:9){
  gdp_sum <- sum(tail(data_vec_list[[rev]],32))
  for (i in 1:9){
    for (j in 1:4){
      gdp_sum_mat[[rev]][i,j] <- sum(tail(as.vector(pred_ver[[rev]][[i]][[j]]),n=32)) 
      gdp_sum_rel_err[[rev]][i,j] <- abs((gdp_sum_mat[[rev]][i,j]-gdp_sum))/gdp_sum
    }
  }
}

gdp_sum_arima_rel_err <- lapply(1:9,function(x)list())
for (rev in 1:9){
  gdp_sum <- sum(tail(data_vec_list[[rev]],32))
  gdp_sum_arima <- sum(tail(as.vector(pred_arima[[rev]]),n=32)) 
  gdp_sum_arima_rel_err[[rev]] <- abs((gdp_sum_arima-gdp_sum))/gdp_sum
}


titlenames <- c("Release Dec 2021","Release Mar 2022","Release Jun 2022","Release Sep 2022","Release Dec 2022","Release Mar 2023",
                "Release Jun 2023","Release Sep 2023","Release Dec 2023")
for (rev in 1:9){
  rmse_df <- data.frame(gdp_sum_rel_err[[rev]])
  #rmse_df <- cbind(rmse_df,rmse_ar)
  colnames(rmse_df) <- c("GNAR-ex 0","GNAR-ex 1","GNAR-ex 2","GNAR-ex 3")#,"AR")
  
  rmse_df_melt <- melt(rmse_df)
  rmse_df_melt["lag"] <- rep(1:9,4)
  colnames(rmse_df_melt) <- c("Model","error","lag")
  if(rev %in% c(1:8)){
    assign(paste("p_",rev,sep = ""),ggplot(rmse_df_melt, aes(x = lag, y = error, colour =Model, group = Model)) +  geom_line() + geom_point(size=1)+
             ylab("error")+
             theme(axis.text.x = element_text(size=12),axis.text.y = element_text(size=12),
                   axis.title.x = element_text(size=12),axis.title.y = element_text(size=12),legend.text = element_text(size=12),legend.title = element_blank(),
                   plot.title = element_text(hjust = 0.5),legend.position = "none")+
             scale_x_continuous("lag", labels = as.character(1:9), breaks = 1:9)+ylab("relative error")+ggtitle(titlenames[rev]))
    #ggtitle(paste("Release ",substr(colnames(dt)[9:18],start = 8,stop = 14)[rev],sep = ""))
  }else{
    assign(paste("p_",rev,sep = ""),ggplot(rmse_df_melt, aes(x = lag, y = error, colour =Model, group = Model)) +  geom_line() + geom_point(size=1)+
             ylab("error")+
             theme(axis.text.x = element_text(size=12),axis.text.y = element_text(size=12),
                   axis.title.x = element_text(size=12),axis.title.y = element_text(size=12),legend.text = element_text(size=12),legend.title = element_blank(),
                   plot.title = element_text(hjust = 0.5),legend.position = "bottom")+
             scale_x_continuous("lag", labels = as.character(1:9), breaks = 1:9)+ylab("relative error")+ggtitle(titlenames[rev]))
    #ggtitle(paste("Release ",substr(colnames(dt)[9:18],start = 8,stop = 14)[rev],sep = ""))
  }
}
library(patchwork)
p_1+p_2+p_3
p_4+p_5+p_6
p_7+p_8+p_9
library(gridExtra)
grid.arrange(p_1,p_2,p_3,p_4,p_5,p_6,p_7,p_8,p_9,nrow=3)
library(ggpubr)
ggarrange(p_1,p_2,p_3,p_4,p_5,p_6,p_7,p_8,p_9, ncol=3, nrow=3, common.legend = TRUE, legend="bottom")

#####################################################################################################################
############# HEATMAP OF CORRELATIONS BETWEEN TOTAL GDP FORECAST, AND RELATIVE ERROR OF TOTAL GDP #############
#####################################################################################################################

cor_mat_totalgdp <- lapply(1:9,function(x)matrix(rep(1,81),nrow = 9,ncol = 9))
for (rev in 1:9){
  for (i in 1:8){
    for (j in (i+1):9){
      cor_mat_totalgdp[[rev]][i,j] <- cor_mat_totalgdp[[rev]][j,i] <- cor(gdp_sum_mat[[rev]][i,],gdp_sum_mat[[rev]][j,],method = "pearson")
    }
  }
}
cor_mat_totalgdp_melt <- lapply(cor_mat_totalgdp, melt)

for (rev in 1:9){
  print(ggplot(as.data.frame(cor_mat_totalgdp_melt[[rev]]), aes(as.character(Var1), as.character(Var2), fill= value)) + 
    geom_tile() +
    labs(title = paste("Release ",substr(colnames(dt)[9:18],start = 8,stop = 14)[rev],sep = ""),
         x = "lag",
         y = "lag"))
}
######################################################
################ MODEL AVERAGING ################
######################################################

ma_rel_err <- lapply(1:9,function(x)list())
for (rev in 1:9){
  gdp_sum <- sum(tail(data_vec_list[[rev]],32))
  ma <- mean(gdp_sum_mat[[rev]][,2])
  ma_rel_err[[rev]] <- abs((ma-gdp_sum))/gdp_sum
}

ma_allmod <- cbind(round(unlist(ma_rel_err),5),round(unlist(gdp_sum_arima_rel_err),5))
colnames(ma_allmod) <- c("MA all GNAR-ex stage 1","ARIMA")
rownames(ma_allmod) <- paste("Release ",substr(colnames(dt)[9:18],start = 8,stop = 14)[1:9],sep = "")
xtable(ma_allmod)

# model averaging with weights depending on BIC
library(Rmpfr)
weights_bic <- lapply(1:9, function(x)list())
for (rev in 1:9){
  expnum <- mpfr(bicres4[[rev]]*(-1/2),precBits = 106)
  w_den <- sum(exp(expnum))# calculate the weights denominator
  weights_bic[[rev]] <- exp(expnum)/w_den
}

ma_rel_err_bicw <- lapply(1:9,function(x)list())
for (rev in 1:9){
  gdp_sum <- sum(tail(data_vec_list[[rev]],32))
  ma <- sum((weights_bic[[rev]]*gdp_sum_mat[[rev]][,-c(5)])[,1])
  ma_rel_err_bicw[[rev]] <- abs((ma-gdp_sum)/gdp_sum)
}

# union forecast interval from MA 

lower_or_mat <- lapply(lower_or, function(x) do.call(rbind,x))
upper_or_mat <- lapply(upper_or, function(x) do.call(rbind,x))

lower_or_mat <- lapply(lower_or_mat,function(x) do.call(rbind,x[,2]))
upper_or_mat <- lapply(upper_or_mat,function(x) do.call(rbind,x[,2]))

unionint <- lapply(1:9,function(x) list())
pred_incl <- list()
for (rev in 1:9){
  for(inode in 1:32){
    unionint[[rev]][[inode]] <- interval_union(Intervals(t(rbind(lower_or_mat[[rev]][,inode],upper_or_mat[[rev]][,inode]))))
  }
  unionintmat <- do.call("rbind",unionint[[rev]])
  pred_incl[rev] <- length(which(between(tail(data_vec_list[[rev]],n=32),unionintmat[,1],
                                       unionintmat[,2])))/32
}


#### model averaging industry-level #### 

pred_ver_ind <- lapply(1:length(1:9),function(x) matrix(nrow = 9,ncol = 32))
for (rev in 1:9){
  for (lagi in 1:9){
    pred_ver_ind[[rev]][lagi,] <- tail(as.vector(pred_ver[[rev]][[lagi]][[2]]),n=32)
  }
}

ma_rel_err <- lapply(1:9,function(x)list())
for (rev in 1:9){
  gdp_sum <- tail(data_vec_list[[rev]],32)
  ma <- apply(pred_ver_ind[[rev]],2,mean)
  ma_rel_err[[rev]] <- abs((ma-gdp_sum))/gdp_sum
}

load(paste(dir_path,"/data/public_updated20apr/sector_names_public.RData",sep = ""))
CPA_ind_name <- conc_public[,3:4]
CPA_ind_name <- unique(CPA_ind_name)
codes_notin <- setdiff(unique(conc_public[,3]),sparse_1_nodedgerem_thres0.4$CPA_node_2[,1])
CPA_ind_name <- CPA_ind_name[-which(CPA_ind_name[,1] %in% codes_notin),]
all(CPA_ind_name[,1]==sparse_1_nodedgerem_thres0.4$CPA_node_2[,1])

errordf <- lapply(1:9,function(x) list())
for (rev in 1:9){
  errordf[[rev]] <- data.frame(
    x=CPA_ind_name[,2],
    y=ma_rel_err[[rev]]
  )
}

# Horizontal version
for (rev in 1:9){
    assign(paste("p_",rev,sep = ""),ggplot(errordf[[rev]], aes(x=x, y=y)) +
    geom_segment( aes(x=x, xend=x, y=0, yend=y), color="skyblue") +
    geom_point( color="blue", size=4, alpha=0.6) +
    theme_light() +
    coord_flip() +
    theme(
      panel.grid.major.y = element_blank(),
      panel.border = element_blank(),
      axis.ticks.y = element_blank()
    )+ylab(" ")+xlab("relative error")+ggtitle(paste("Release ",substr(colnames(dt)[9:18],start = 8,stop = 14)[rev],sep = "")))
}

library(gridExtra)
grid.arrange(p_1,p_2,nrow=2)

# plot by industry, mean and sd of relative error across releases
ma_rel_err_conc <- Reduce(cbind,ma_rel_err)
ma_rel_err_summ_df <- data.frame(mean=apply(ma_rel_err_conc,1,mean),sd=apply(ma_rel_err_conc,1,sd),ind=CPA_ind_name[,2])
p<- ggplot(ma_rel_err_summ_df, aes(x=mean, y=ind)) + 
  geom_line() +
  geom_point(color="blue")+
  geom_errorbar(aes(xmin=mean-sd, xmax=mean+sd), width=.4,color="skyblue")+theme_classic()+ylab(" ")+xlab("relative error")


ma_rel_err_precov <- c(0.20626879, 0.20438890 ,0.15296496 ,0.11185910, 0.21471886 ,0.14856731, 0.16876466, 0.16843101, 0.11119359 ,0.10511052 ,
                       0.23325073, 0.16704475 ,0.16224809, 0.15753224 ,0.17476478, 0.16385308, 0.19501818 ,0.17487006, 0.12959336, 0.07425401 ,
                       0.12794837, 0.16052101, 0.18409527, 0.18228139, 0.15983526 ,0.13093793 ,0.21594079, 0.10705130, 0.13692580, 0.15230211 ,
                       0.17340790, 0.25320257 )
names(ma_rel_err_precov) <- c( 1,2,3,4,5,6,7,11,12,13 ,
                              14,15,16,17,18,19,20,21,22,23 ,
                              24,25,26,27,28,29,30,31,32,34 ,
                              37,38)
ma_rel_err_summ_df <- data.frame(precov=ma_rel_err_precov/sqrt(48),mean=apply(ma_rel_err_conc,1,mean),sd=apply(ma_rel_err_conc,1,sd),ind=CPA_ind_name[,2])
ma_rel_err_summ_df_melt <- melt(ma_rel_err_summ_df,id=c('sd','ind'))
p<- ggplot(ma_rel_err_summ_df_melt, aes(x=value, y=ind,group = variable)) + 
  #geom_line() +
  geom_point(aes(color=variable))+
  #geom_errorbar(stat="summary")
  geom_errorbar(aes(xmin=ma_rel_err_summ_df_melt[which(ma_rel_err_summ_df_melt$variable=='mean'),'value']-sd, xmax=ma_rel_err_summ_df_melt[which(ma_rel_err_summ_df_melt$variable=='mean'),'value']+sd), width=.4,color="skyblue")+theme_classic()+ylab(" ")+xlab("relative error")+
  scale_colour_discrete(name="Value",labels = c("Pre covid-19", "Mean"))

# MAYBE USE BOXPLOTS FOR DISTRIBUTION OF RELL ERROR ACROSS RELEASES

# take top 3 worst perfroming industries from each release
top3_release <- lapply(ma_rel_err,function(x) x[order(x,decreasing = TRUE)][1:3])
releases <- factor(rep(substr(colnames(dt)[9:18],start = 8,stop = 14)[1:9],each=3))
top3_release_df <- data.frame(ind=names(unlist(top3_release)),err=unlist(top3_release),release= releases)
top3_codes <- sparse_1_nodedgerem_thres0.4$CPA_node_2[match(top3_release_df$ind,sparse_1_nodedgerem_thres0.4$CPA_node_2[,2]),1]
top3_ind <- CPA_ind_name[match(top3_codes,CPA_ind_name[,1]),2]
top3_release_df$ind <- factor(paste(top3_ind,seq(1:27)),levels = paste(top3_ind,seq(1:27)))
ggplot(top3_release_df, aes(x=ind, y=err,colour=release)) +
  geom_segment( aes(x=ind, xend=ind, y=0, yend=err)) +
  geom_point( color="blue", size=4, alpha=0.6) +
  theme_light() +
  coord_flip() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  )+xlab(" ")+ylab("relative error")+
  scale_x_discrete(labels=top3_ind)+
  scale_colour_discrete(name="Release",breaks=c(as.character(rev(unique(top3_release_df$release)))))# c("Dec 2023","Sep 2023","Jun 2023","Mar 2023","Dec 2022","Sep 2022","Jun 2022","Mar 2022","Dec 2021")

# ordering of industries by error, for each release
lapply(errordf,function(x) x[order(x$y, decreasing = TRUE),])

##################################################################
########### GATHER IN TABLE BEST MODEL ERRORS VS ARIMA ###########
##################################################################
#lapply(gdp_sum_rel_err,function(x) x[,-5]) #without AR Model
bestmod <- t(sapply(gdp_sum_rel_err, function(x) which(x==min(x),arr.ind = TRUE)))
minbestmod <- sapply(gdp_sum_rel_err, function(x) min(x))
bestcomb <- cbind(bestmod,format(round(minbestmod, 5), nsmall = 2))
bestcomb <- apply(bestcomb,2,function(x) as.numeric(x))
bestcomb[,2] <- bestcomb[,2]-1
bestcomb <- cbind(bestcomb,round(unlist(gdp_sum_arima_rel_err),5))
colnames(bestcomb) <- c("GNAR-lag","GNAR-stage","GNAR-error","ARIMA error")
rownames(bestcomb) <- paste("Release ",substr(colnames(dt)[9:18],start = 8,stop = 14)[1:9],sep = "")

library(xtable)
xtable(bestcomb)
########### ########### ########### 
########### BIC TABLES ########### 
########### ########### ########### 
library(xtable) # library to extract tex code for table from matrix in R

table(sapply(bicres4,function(x) apply(x,1,which.min))) # best stag among lags, among releases
bicres4 <- lapply(bicres4,function(x) {colnames(x) <- paste("stage_",c(0:3),sep = "");rownames(x) <- paste("lag_",c(1:9),sep = "");x})
lapply(bicres4, function(x) xtable(x))
# best lag-stage among releases
modspec_rel <- t(sapply(bicres4,function(x) which(x==min(x),arr.ind = TRUE)))
rownames(modspec_rel) <-paste("Release ",substr(colnames(dt)[9:18],start = 8,stop = 14)[1:9],sep = "")
colnames(modspec_rel) <- c("lag","stage")
modspec_rel[,2] <-modspec_rel[,2]-1 
xtable(modspec_rel)

####################################################################### 
################ FIND CORELLATIONS BETWEEN RESIDUALS ########### 
##################################################################

# do that for best performing model, for each release (using relative error Total GDP)


# release 1
covresmat_rel <- covresmat_list[[1]][[8]][[4]]
covresmat_rel <- covresmat_list[[7]][[bestmod[7,1]]][[bestmod[7,2]]]
hist(covresmat_rel[upper.tri(covresmat_rel)])
pos_uptr <- which(upper.tri(matrix(nrow=871,ncol=871)),arr.ind = TRUE)
lowerb <- quantile(covresmat_rel[upper.tri(covresmat_rel)],.025)
upperb <- quantile(covresmat_rel[upper.tri(covresmat_rel)],.975)
outlierspos <- which((covresmat_rel[upper.tri(covresmat_rel)]<lowerb | 
                       covresmat_rel[upper.tri(covresmat_rel)]>upperb))

covpos <- pos_uptr[outlierspos,] # corresponding entries in covariance matrix
all(names(sort(table(outlierspos[,1])))==names(sort(table(outlierspos[,2]))))

pay_gdp_highcov <- covpos[which((covpos[,1]<=833 & covpos[,2]>833)),]

removepaym_pos <- unique(pay_gdp_highcov[,1]) # pos of payments time series to remove


data_edges2 <-  get(paste("prep_data_rev_",1,"uni",sep = ""))$data_edges[-removepaym_pos,]

# not in common nodes (after removing edges, removing also some nodes)
removegdp_pos <- setdiff(c(1:38),unique(c(unique(data_edges2[,2]),unique(data_edges2[,1]))))

df_pay2 <- get(paste("prep_data_rev_",1,"uni",sep = ""))$df_pay[-removepaym_pos,]
if (length(removegdp_pos)==0){
  df_gva2 <- get(paste("prep_data_rev_",1,"uni",sep = ""))$df_gva
  CPA_node2 <- get(paste("prep_data_rev_",1,"uni",sep = ""))$CPA_node
  data_nodes2 <- get(paste("prep_data_rev_",1,"uni",sep = ""))$data_nodes
}else{
  df_gva2 <- get(paste("prep_data_rev_",1,"uni",sep = ""))$df_gva[-removegdp_pos,]  
  CPA_node2 <- get(paste("prep_data_rev_",1,"uni",sep = ""))$CPA_node[-removegdp_pos,]
  data_nodes2 <- get(paste("prep_data_rev_",1,"uni",sep = ""))$data_nodes[-removegdp_pos]
}

ts_pay_gdp2 <-  get(paste("prep_data_rev_",1,"uni",sep = ""))$ts_pay_gdp[-c(removepaym_pos,removegdp_pos+833),]
graph2 <-  graph_from_edgelist(as.matrix(data_edges2),directed = TRUE)
nedges2 <-  ecount(graph2)
nnodes2 <- vcount(graph2)
n_edges_nodes2 <-  ecount(graph2)+vcount(graph2)

prep_data_rev_1uni_sparse <- list(df_pay=df_pay2,df_gva=df_gva2,ts_pay_gdp=ts_pay_gdp2,CPA_node=CPA_node2,graph=graph2,nnodes=nnodes2,nedges=nedges2,
                                                     n_edges_nodes=n_edges_nodes2,data_edges=data_edges2,data_nodes=data_nodes2)



########################################################################################
plot(diag(covresmat_list[[1]][[8]][[4]]),type = "l") # plot variance of residuals (non-constant, big bumps)


########################################################################################
############## residuals of this regime/lag-stag specification ############## 
########################################################################################
norm_test_bestmod_all <- norm_test_bestmod_nodes <- c()
for (rev in 1:9){
  res_mean_mat <- apply(resmat_list[[rev]][[bestmod[rev,1]]][[bestmod[rev,2]]],1,mean)
  if (nodes_only){
    res_mean_mat_nodes <- apply(resmat_list[[rev]][[bestmod[rev,1]]][[bestmod[rev,2]]],1,mean)#[,126:162]
  }else{
    res_mean_mat_nodes <- apply(resmat_list[[rev]][[bestmod[rev,1]]][[bestmod[rev,2]]][,tail(seq(n_edges_nodes3),n=nnodes3)],1,mean)#[,126:162]
    res_mean_mat_pay <- apply(resmat_list[[rev]][[bestmod[rev,1]]][[bestmod[rev,2]]][,head(seq(n_edges_nodes3),n=nedges3)],1,mean)#[,1:125]
  }
  
  plot(res_mean_mat,type = "l") 
  qqnorm(res_mean_mat,cex.lab=1.2,main=" ",cex.axis=1.2)
  qqline(res_mean_mat,col="red")
  norm_test_bestmod_all[rev] <- shapiro.test(res_mean_mat)[[2]]
  norm_test_bestmod_nodes[rev] <- shapiro.test(res_mean_mat_nodes)[[2]]
  #print(shapiro.test(res_mean_mat_pay)[[2]])
  
}
shap_test <- cbind(paste("Release ",substr(colnames(dt)[9:18],start = 8,stop = 14)[1:9],sep = ""),round(norm_test_bestmod_all,digits = 4))
colnames(shap_test) <- c("","p-value")
xtable(as.matrix(shap_test))
#res_mean_mat <- apply(resmat,1,mean)

# SHAPIRO-WILK NORMALITY TEST HO: THERE IS NO DIFFERENCE BETWEEN NORMAL DISTRIB AND DISTR OF DATA


############################
# PLOT GROWTH RATES GDP
############################
for (rev in 9){
  # data 
  # gdpdata <- prepr_gva(dt,rev_mon_year[rev])
  # #print(tail(colnames(gdpdata),n=1))
  # prep_data <- prepr_comb(gdpdata,paydata,"2016-01",tail(colnames(gdpdata),n=1),date_end_2_list[rev],all_nas_row=FALSE,atleast1_na=TRUE,NULL)
  # data <- prep_data$ts_pay_gdp[-edge_pos_rem,]
  prep_data <- get(paste("prep_data_rev_",rev,"uni",sep = ""))
  data <- prep_data$ts_pay_gdp
  # data <- log(data)
  # 
  # # stationary
  # data <- stationary_ts(data)
  # 
  # growth rates
  data <- growthrates(data)
  datadf <- data.frame(t(data),year=c(1:ncol(data)))
  colnames(datadf)[834:871] <- prep_data$CPA_node[,1]
  data_final <- melt(datadf[,c(866:(866+5),872)], id.vars = "year")
  
  # Plot the final data
  print(ggplot(data_final,                            
         aes(x = year,
             y = value,
             col = variable)) + geom_line())
  
}

nodesextr <- c(7,12,15,13,29,26,27,34,37,35,38)

CPA_node2[nodesextr,]
testtest <- which(CPA_ind_name[,1] %in% CPA_node2[nodesextr,1])
CPA_ind_name[testtest,]

################################################################
######## CORELLATIONS BETWEEN GDP AND PAYMENTS TS ########
################################################################

cormatrev <- lapply(1:9, function(x) matrix(,ncol=38,nrow = 833))
paymrem <- lapply(1:9,function(x) list())
for (rev in 1:9){
  prep_data <- get(paste("prep_data_rev_",rev,"uni",sep = ""))
  data <- prep_data$ts_pay_gdp
  # data <- log(data)
  # 
  # # stationary
  # data <- stationary_ts(data)
  # 
  # growth rates
  data <- growthrates(data)
  for (i in 1:833){
    ind_col <- 1
    for (j in 834:871){
      cormatrev[[rev]][i,ind_col] <- cor(data[i,],data[j,],method = "pearson")
      ind_col <- ind_col+1
    }
  }
  hist(cormatrev[[rev]])
  lowerb <- quantile(cormatrev[[rev]],.025)
  upperb <- quantile(cormatrev[[rev]],.975)
  abline(v=upperb,col="red")
  abline(v=lowerb,col="red")
  #thres <- quantile(cormatrev[[rev]],.95)
  thres <- .3
  outlierspos <- which((cormatrev[[rev]]<(-thres))|(cormatrev[[rev]]>thres),arr.ind = TRUE)
  paymrem[[rev]] <- unique(outlierspos[,1])
}


removepaym_pos <- unique(unlist(paymrem))
################################################
# graph related changes common for all versions
################################################
data_edges3 <-  get(paste("prep_data_rev_1uni",sep = ""))$data_edges[-removepaym_pos,]
# not in common nodes (after removing edges, removing also some nodes)
removegdp_pos <- as.numeric(setdiff(get(paste("prep_data_rev_1uni",sep = ""))$data_nodes,unique(c(unique(data_edges3[,2]),unique(data_edges3[,1])))))
if (length(removegdp_pos)!=0){
  CPA_node3 <- get(paste("prep_data_rev_1uni",sep = ""))$CPA_node[-removegdp_pos,]
  data_nodes3 <- get(paste("prep_data_rev_1uni",sep = ""))$data_nodes[-removegdp_pos]
}else{
  CPA_node3 <- get(paste("prep_data_rev_1uni",sep = ""))$CPA_node
  data_nodes3 <- get(paste("prep_data_rev_1uni",sep = ""))$data_nodes
}
graph3 <-  graph_from_data_frame(as.data.frame(data_edges3),directed = TRUE,vertices = data_nodes3)
nedges3 <-  ecount(graph3)
nnodes3 <- vcount(graph3)
n_edges_nodes3 <-  ecount(graph3)+vcount(graph3)
for (rev in 1:10){
  
  
    
  df_pay3 <- get(paste("prep_data_rev_",rev,"uni",sep = ""))$df_pay[-removepaym_pos,]
  if (length(removegdp_pos)==0){
    df_gva3 <- get(paste("prep_data_rev_",rev,"uni",sep = ""))$df_gva
    ts_pay_gdp3 <-  get(paste("prep_data_rev_",rev,"uni",sep = ""))$ts_pay_gdp[-removepaym_pos,]
  }else{
    df_gva3 <- get(paste("prep_data_rev_",rev,"uni",sep = ""))$df_gva[-removegdp_pos,]  
    ts_pay_gdp3 <-  get(paste("prep_data_rev_",rev,"uni",sep = ""))$ts_pay_gdp[-c(removepaym_pos,removegdp_pos+833),]
  }
  
  assign(paste("sparse_",rev,"_thres0.3",sep = ""),list(df_pay=df_pay3,df_gva=df_gva3,ts_pay_gdp=ts_pay_gdp3,CPA_node=CPA_node3,graph=graph3,nnodes=nnodes3,nedges=nedges3,
                                    n_edges_nodes=n_edges_nodes3,data_edges=data_edges3,data_nodes=data_nodes3))
  
}
#save(sparse_1,sparse_2,sparse_3,sparse_4,sparse_5,sparse_6,sparse_7,sparse_8,sparse_9,sparse_10,file="/Users/amantziou/Documents/R_directory_Turing_ONS/collab/data/data_preprocessing/sparse_data_allver.RData") # sparse network from Pearson's 0.4 threshold 
#save(sparse_1_thres0.7,sparse_2_thres0.7,sparse_3_thres0.7,sparse_4_thres0.7,sparse_5_thres0.7,sparse_6_thres0.7,sparse_7_thres0.7,sparse_8_thres0.7,sparse_9_thres0.7,sparse_10_thres0.7,file="/Users/amantziou/Documents/R_directory_Turing_ONS/collab/data/data_preprocessing/sparse_data_allver_pears0.7.RData") # sparse network from Pearson's 0.7 threshold 
#save(sparse_1_thres0.3,sparse_2_thres0.3,sparse_3_thres0.3,sparse_4_thres0.3,sparse_5_thres0.3,sparse_6_thres0.3,sparse_7_thres0.3,sparse_8_thres0.3,sparse_9_thres0.3,sparse_10_thres0.3,file="/Users/amantziou/Documents/R_directory_Turing_ONS/collab/data/data_preprocessing/sparse_data_allver_pears0.3.RData") # sparse network from Pearson's within [-0.3,0.3] threshold 
for (rev in 1:9){
  if (all(E(get(paste("sparse_",rev,sep = ""))$graph) %in% E(get(paste("sparse_",rev+1,sep = ""))$graph)) &
      all(E(get(paste("sparse_",rev+1,sep = ""))$graph) %in% E(get(paste("sparse_",rev,sep = ""))$graph))){
    print(c(rev,rev+1,"same"))
  }
}

################################################################# 
############# TABLE WITH INDUSTRY NAMES OF NODES REMOVED ############# 
################################################################# 
xtable(CPA_ind_name[which(CPA_ind_name[,1] %in% (get(paste("prep_data_rev_1uni",sep = ""))$CPA_node[removegdp_pos,1])),])

###########################################################
############ PLOT SPARSIFIED GRAPH #####################
###########################################################
plot(graph3,vertex.size=15,edge.arrow.size=0.3,vertex.color="lightblue",vertex.label.cex=.8,edge.width=2)

##########################################################################################
################## TOTAL GDP RESCALED BACK TO INDEX VALUE AND COMPARED WITH TOTAL ##################
################## GDP INDEX FOR DIFFERENT RELEASES ANDREW DATA ##################
##########################################################################################

rescale_dates <- c("2019-01","2019-07","2019-07","2019-10","2019-01","2019-04","2019-07")
date_ind <- 1
totalgdp_ind <- c()
for (i in c(3,seq(4,9))){
  prep_data <- get(paste("sparse_",i,"_nodedgerem_thres0.4",sep = ""))
  totalgdp_ind[date_ind] <- mean(gdp_sum_mat[[i]][,2])/sum(tail(prep_data$ts_pay_gdp_2[,which(colnames(prep_data$ts_pay_gdp_2)==rescale_dates[date_ind])],n = 32))
  date_ind <- date_ind+1
}
library(readxl)
newdata <- read_excel("/Users/u2371456/Downloads/Tail data for each release since August 2023.xlsx")
forecast_dates <- c("2022JAN","2022JUL","2022JUL","2022OCT","2023JAN","2023APR","2023JUL")
relerr <- list()
relerr_df <- list()
for (i in 1:7){
  true <- newdata$`Total GVA`[which(newdata$`Time period`==forecast_dates[i])]
  relerr[[i]] <- abs((totalgdp_ind[i]*100-true))/true
  relerr_df[[i]] <- data.frame(relerror=rev(relerr[[i]]),publication=factor(rev(newdata$Publication[which(newdata$`Time period`==forecast_dates[1])]),levels=rev(newdata$Publication[which(newdata$`Time period`==forecast_dates[1])])))
}

for(i in 1:7){
  print(ggplot(relerr_df[[i]],aes(x=as.factor(publication),y=relerror))+
          geom_bar(stat = "identity", fill="steelblue")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=12),axis.title.x = element_blank())+
          geom_text(aes(label=round(relerror,digits = 7)), hjust=1, color="white", size=6,angle=90)+ggtitle(forecast_dates[i])
        )
}

