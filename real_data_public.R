# PUBLIC DATA VERSION

library(dplyr)
library(zoo) # for yearmonth objects

load("/Users/amantziou/Documents/R_directory_Turing_ONS/collab/data/public_updated20apr/industry_panel_Pluto_ONS_monthly_CPA_public.RData") # 40 unique industries


gdp_data <- dplyr::select(dt, time, code, "GVA sa")

no_ind <- length(unique(gdp_data$code))
no_ts <- 321 # from Jan 1997 until Sep 2023

# remove NAs from gdp_data (for sic3 case)
gdp_data <- gdp_data[complete.cases(gdp_data),]


gdp_data_str <- matrix(,nrow = no_ind,ncol = length(gdp_data[seq(1,nrow(gdp_data),no_ind),]$time)+1)
for(i in 1:no_ind){
  if(table(gdp_data[seq(i,nrow(gdp_data),no_ind),'code'])==no_ts){ 
    gdp_data_str[i,1] <- gdp_data[seq(i,nrow(gdp_data),no_ind),'code'][1]
    gdp_data_str[i,2:ncol(gdp_data_str)] <- gdp_data[seq(i,nrow(gdp_data),no_ind),'GVA sa']
  }else{
    print(c("non similar codes for row ",i)) 
    gdp_data_str[i,1] <- gdp_data[seq(i,nrow(gdp_data),no_ind),'code'][1]
    gdp_data_str[i,2:ncol(gdp_data_str)] <- c(gdp_data[seq(i,nrow(gdp_data),no_ind),'GVA sa'],rep(0,abs(table(gdp_data[seq(i,nrow(gdp_data),no_ind),'code'])-no_ts)))
  }
}

colnames(gdp_data_str) <- c("CPA code",format(gdp_data[seq(1,nrow(gdp_data),no_ind),]$time,"%Y-%m"))
gdp_data_str_df <- as.data.frame(gdp_data_str)

# load updated version of Payments for CPA codes
load("/Users/amantziou/Documents/R_directory_Turing_ONS/collab/data/public_updated18apr/IOT_flow_of_goods_Pluto_monthly_CPA_public.RData")

# remove counts of transactions
newd0 <- d0[,-which(grepl("_cnt_", colnames(d0), fixed = TRUE))]

# public data: remove NAs industry codes in From/To columns
newd0 <- newd0[-union(which(is.na(newd0$from)), which(is.na(newd0$to))),]

# check that the same CPA codes are shared in both dfs
all(unique(newd0$to)%in%unique(gdp_data_str_df$`CPA code`))
all(unique(newd0$from)%in%unique(gdp_data_str_df$`CPA code`))
all(unique(gdp_data_str_df$`CPA code`)%in%unique(newd0$to))
all(unique(gdp_data_str_df$`CPA code`)%in%unique(newd0$from))

# codes that are in GDP but not in payments
notin_codes_to <- setdiff(unique(gdp_data_str_df$`CPA code`),unique(newd0$to))
notin_codes_from <- setdiff(unique(gdp_data_str_df$`CPA code`),unique(newd0$from))
notin_codes_all <- unique(union(notin_codes_from,notin_codes_to))
gdp_data_str_df2 <- gdp_data_str_df[-which(gdp_data_str_df$`CPA code` %in% notin_codes_all),] 

# keep the columns/yearmonths of GDP corresponding to yearmonth columns in Payments data
starttime <- which(names(gdp_data_str_df2)=="2016-01") # start time of payments ts ("2015-08")
endtime <- which(names(gdp_data_str_df2)=="2023-09") # end time of payments ts is Oct 2023 but GDP data up to including Sep 2023
gdp_data_str_df2 <- gdp_data_str_df2[,c(1,starttime:endtime)]

# remove from payments data the time stamps for which no GDP (i.e. after Oct 2023)
newd0 <- newd0[,1:which(names(newd0)=="sep_amt_2023")]

# align names of timestamps between GDP ts and PAY ts
colnames(newd0)[-c(1,2)] <- colnames(gdp_data_str_df2[-c(1)])


########### NAs handling ########### 

# for public data: check which all NAs in rows
row_nas <- c()
for (i in 1:nrow(newd0)){
  if (all(is.na(newd0[i,3:ncol(newd0)]))){
    row_nas <- c(row_nas,i)
  }
}
# check which at least one NAs
#dim(newd0[!complete.cases(newd0[,3:ncol(newd0)]),])

# for public data: remove rows with all NAs
newd0 <- newd0[-row_nas,]
# assign 0 to NAs for the rest
newd0[is.na(newd0)] <- 0
# remove ts with more than 85 zeros
newd0 <- newd0[which(rowSums(newd0[,-c(1,2)]==0)<85),]

# for public data: remove all rows with at least one NA
newd0_nona <- newd0[complete.cases(newd0[,3:ncol(newd0)]),]

# remaining CPA codes after removal of transactions
CPA_codes_pub <- unique(union(unique(newd0_nona$to),unique(newd0_nona$from)))
# check which is not anymore and remove from GDP ts
gdp_data_str_df2 <- gdp_data_str_df2[-which(!(gdp_data_str_df2$`CPA code`%in%CPA_codes_pub)),]
# check that now CPA codes aligned between GDP ts and PAY ts
all(gdp_data_str_df2$`CPA code`%in%CPA_codes_pub,CPA_codes_pub%in%gdp_data_str_df2$`CPA code`)
# create numbering for CPA codes
CPA_node <- cbind(sort(CPA_codes_pub),seq(1:length(CPA_codes_pub)))
newd0_nona$to <- CPA_node[match(newd0_nona$to,CPA_node[,1]),2]
newd0_nona$from <- CPA_node[match(newd0_nona$from,CPA_node[,1]),2]
gdp_data_str_df2[,1] <- CPA_node[match(gdp_data_str_df2[,1],CPA_node[,1]),2]

########## create data frame for model input #########

ts_pay_gdp <- rbind(newd0_nona[,-c(1,2)],gdp_data_str_df2[,-c(1)])
ts_pay_gdp <- apply(ts_pay_gdp,2,as.numeric)
rownames(ts_pay_gdp) <- c(paste(newd0_nona[,1],newd0_nona[,2]),gdp_data_str_df2[,1])

######### create graph #########

g <- graph_from_edgelist(as.matrix(newd0_nona[,1:2]),directed = TRUE)
ecount(simplify(g,remove.loops = FALSE))==ecount(g) # check if non multiple edges
nnodes <- vcount(g)
nedges <- ecount(g)
n_edges_nodes <- nnodes+nedges
data_edges <- newd0_nona[,c(1,2)]
data_nodes <- gdp_data_str_df2[,1]



########## make data stationary #########

stationary_ts <- function(ts_data,diff=TRUE,normal=TRUE){
  if (diff & normal){
    new_ts <- apply(ts_data,1,diff)
    new_ts <- t(new_ts)
    sd_ts <- apply(new_ts,1,sd)
    new_ts <- apply(new_ts,1,function(x){x/sd(x)})
    return(list(new_ts,sd_ts))
  }
  if (normal & !diff){
    sd_ts <- apply(ts_data,1,sd)
    new_ts <- apply(ts_data,1,function(x){x/sd(x)})
    return(list(new_ts,sd_ts))
  }
  if (!normal & diff){
    new_ts <- apply(ts_data,1,diff)
    new_ts <- t(new_ts)
    return(new_ts)
  }
}

undo_diff1_sd <- function(ts_data,sd_data,pred){
  return((pred*sd_data)+ts_data[,ncol(ts_data)-1])
}


data <- data.matrix(ts_pay_gdp)
data_train <- data[,1:92]
globalalpha <- TRUE

# difference and normalise
dtn <- stationary_ts(data_train)
data_train_norm <- t(dtn[[1]])
data_train_sd <- dtn[[2]] # keep sd of data_train

#### STL #### 
### whole data ###
library(stats)
stl_list <- apply(data,1,function(x) stl(ts(as.vector(x),start = c(2016,1), end=c(2023,9),frequency = 12),s.window = "periodic"))
stl_mat <- t(sapply(stl_list, function(x) x$time.series[,3]))
seas <- t(sapply(stl_list, function(x) x$time.series[,1]))
trend <- t(sapply(stl_list, function(x) x$time.series[,2]))
colnames(stl_mat) <- colnames(data)[1:(ncol(data))]
rownames(stl_mat) <- rownames(data)
data_train_stl <- stl_mat[,1:92]
### only training ###
stl_list <- apply(data,1,function(x) stl(ts(as.vector(x),start = c(2016,1), end=c(2019,8),frequency = 12),s.window = "periodic"))
res <- t(sapply(stl_list, function(x) x$time.series[,3]))
seas <- t(sapply(stl_list, function(x) x$time.series[,1]))
trend <- t(sapply(stl_list, function(x) x$time.series[,2]))
data_train_stl <- res#+trend
colnames(data_train_stl) <- colnames(data)[1:44]
rownames(data_train_stl) <- rownames(data)
# colnames(res) <- colnames(data)[1:(ncol(data)-1)]
# rownames(res) <- rownames(data)
# colnames(trend) <- colnames(data)[1:(ncol(data)-1)]
# rownames(trend) <- rownames(data)
# last year same month Seasonal component (September)
seas_comp_sep <- t(sapply(stl_list, function(x) x$time.series[,1][9]))


rmse_mat_real <- matrix(nrow=9,ncol = 4)
rmse_mat_real_nodes <- matrix(nrow=9,ncol = 4)
dev_list <- lapply(1:9, function(x) list())
dev_list_stl <- lapply(1:9, function(x) list())
for (lagi in 1:9){
  ind_stag <- 1
  for (stag in 0:3){
    print(c("Lag iter",lagi," stag iter ",stag))
    alphaOrder <- lagi
    betaOrder <- rep(stag,lagi)
    gammaOrder <- lagi
    deltaOrder <- rep(stag,lagi)
    
    # fit GNAR-x
    fit_train <- gnar_x_fit(data_train_stl,nnodes,nedges,n_edges_nodes,data_edges,data_nodes,alphaOrder,betaOrder,
                            gammaOrder,deltaOrder,g,lead_lag_mat,globalalpha=TRUE,lead_lag_weights=FALSE,pay_now = FALSE,lag_0_sep = FALSE)
    print("hi")
    # Predict Sep 23
    sep_pred <- gnar_x_predict(fit_train,data_train_stl,alphaOrder,betaOrder,
                               gammaOrder,deltaOrder,n_edges_nodes,nnodes,nedges,1,fit_train$wei_mat,fit_train$data_loc_mat,
                               pay_now=FALSE,pay_data_now=NULL)
    
    # format data point (diff and norm) to compare to prediction
    #data_vec <- (data[,93]-data[,92])/data_train_sd
    #rmse_mat_real[lagi,ind_stag] <- sqrt(mean((sep_pred[,lagi+1]-data_vec)^2))
    #rmse_mat_real_nodes[lagi,ind_stag] <- sqrt(mean((tail(sep_pred[,lagi+1],n=nnodes)-tail(data_vec,n=nnodes))^2))
    dev_list_stl[[lagi]][[ind_stag]] <- sep_pred[,lagi+1]
    ind_stag <- ind_stag+1
  }
}


#saveRDS(dev_list,file = "/Users/amantziou/Documents/R_directory_Turing_ONS/collab/results/pred_maxlag9_maxstag3_ar.rds",version = 2)
#saveRDS(dev_list_stl,file = "/Users/amantziou/Documents/R_directory_Turing_ONS/collab/results/pred_maxlag9_maxstag3_stl.rds",version = 2)
#saveRDS(dev_list_stl,file = "/Users/amantziou/Documents/R_directory_Turing_ONS/collab/results/pred_maxlag9_maxstag3_stl_wholedata.rds",version = 2)
#saveRDS(dev_list_stl,file = "/Users/amantziou/Documents/R_directory_Turing_ONS/collab/results/pred_maxlag9_maxstag3_ar_stl_predSep19.rds",version = 2)

rmse_ar <- rmse_ar_nodes <- c()
rmse_ar_stl <- rmse_ar_stl_nodes <- c()
simtrainvar <- t(data_train_norm)
simtrainvar <- t(data_train_stl)
for (i in 1:9){
  simple_ar <- apply(simtrainvar, 2, function(x){forecast(auto.arima(x,d=0,D=0,max.p = i,max.q = 0,max.P = 0,max.Q = 0,stationary = TRUE,seasonal = FALSE,
                                                                     ic="bic",allowmean = FALSE,allowdrift = FALSE,trace = FALSE),h=1)$mean})
  #rmse_ar_nodes[i] <- sqrt(mean((simple_ar[(nedges+1):n_edges_nodes]-data_vec[(nedges+1):n_edges_nodes])^2))
  #rmse_ar[i] <- sqrt(mean((simple_ar-data_vec)^2))
  dev_list_stl[[i]][[5]] <- simple_ar
}
# get pred with ar for trend
x_trend <- seq(1,ncol(trend))
fitlist <- apply(trend, 1, function(x) lm(x~x_trend))
pred_trend <- sapply(fitlist, function(x) predict(x,data.frame(x_trend=45),type="response"))

### FOR DIFF & STAND DATA ###
rmse_mat <- matrix(nrow=9,ncol = 5)
rmse_mat_nodes <- matrix(nrow=9,ncol = 5)
for (i in 1:9){
  for (j in 1:5){
    data_vec <- (data[,93]-data[,92])/data_train_sd
    pred <- dev_list[[i]][[j]]
    rmse_mat[i,j] <- sqrt(mean((pred-data_vec)^2))
    rmse_mat_nodes[i,j] <- sqrt(mean((tail(pred,n=nnodes)-tail(data_vec,n=nnodes))^2))
  }
}
### FOR STL ###
# create dataframe for rmse from  devlist_stl
rmse_mat_stl <- matrix(nrow=9,ncol = 5)
rmse_mat_stl_nodes <- matrix(nrow=9,ncol = 5)
for (i in 1:9){
  for (j in 1:5){
    pred <- dev_list_stl[[i]][[j]]+seas[,93]+trend[,93]
    rmse_mat_stl[i,j] <- sqrt(mean((pred-data[,93])^2))
    rmse_mat_stl_nodes[i,j] <- sqrt(mean((tail(pred,n=nnodes)-tail(data[,93],n=nnodes))^2))
  }
}
# alternative: model fitted to seas adjusted data= T+R and add only S component from last year same month of pred
rmse_mat_stl <- matrix(nrow=9,ncol = 5)
rmse_mat_stl_nodes <- matrix(nrow=9,ncol = 5)
for (i in 1:9){
  for (j in 1:5){
    pred <- as.vector(dev_list_stl[[i]][[j]]+seas_comp_sep)
    rmse_mat_stl[i,j] <- sqrt(mean((pred-data[,93])^2))
    rmse_mat_stl_nodes[i,j] <- sqrt(mean((tail(pred,n=nnodes)-tail(data[,93],n=nnodes))^2))
  }
}
# ALTERNATIVE ADD pred_trend and seas_comp
rmse_mat_stl <- matrix(nrow=9,ncol = 5)
rmse_mat_stl_nodes <- matrix(nrow=9,ncol = 5)
for (i in 1:9){
  for (j in 1:5){
    pred <- as.vector(dev_list_stl[[i]][[j]]+pred_trend+seas_comp_sep)
    rmse_mat_stl[i,j] <- sqrt(mean((pred-data[,45])^2))
    rmse_mat_stl_nodes[i,j] <- sqrt(mean((tail(pred,n=nnodes)-tail(data[,45],n=nnodes))^2))
  }
}


rmse_df <- data.frame(rmse_mat_stl_nodes)
#rmse_df <- cbind(rmse_df,rmse_ar)
colnames(rmse_df) <- c("GNARx-0","GNARx-1","GNARx-2","GNARx-3","AR")

rmse_df_melt <- melt(rmse_df)
rmse_df_melt["lag"] <- rep(1:9,5)
colnames(rmse_df_melt) <- c("model","RMSE","lag")
ggplot(rmse_df_melt, aes(x = lag, y = RMSE, colour =model, group = model)) +  geom_line() + geom_point(size=1)+
  ylab("RMSE")+
  theme(axis.text.x = element_text(size=12),axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=12),axis.title.y = element_text(size=12),legend.text = element_text(size=12),legend.title = element_text(size=12))+
  scale_x_continuous("lag", labels = as.character(1:9), breaks = 1:9)



# Plot for errors per industry (node level) for lag=1, stage=1

load("/Users/amantziou/Documents/R_directory_Turing_ONS/collab/data/public_updated20apr/sector_names_public.RData")
CPA_ind_name <- conc_public[,3:4]
CPA_ind_name <- unique(CPA_ind_name)
codes_notin <- setdiff(unique(conc_public[,3]),CPA_node[,1])
CPA_ind_name <- CPA_ind_name[-which(CPA_ind_name[,1] %in% codes_notin),]
all(CPA_ind_name[,1]==CPA_node[,1])

# Create data frame 
errordf <- data.frame(
  x=CPA_ind_name[,2],
  y=abs(tail(dev_list[[2]][[1]],n=nnodes)- tail(data_vec,n=nnodes))
)
### FOR STL, lag=7, stage=3 ###
errordf <- data.frame(
  x=CPA_ind_name[,2],
  y=abs(tail(dev_list_stl[[7]][[4]]+seas[,93]+trend[,93],n=nnodes)- tail(data[,93],n=nnodes))
)
### FOR STL up to Aug 2023, sum GDP lag=2, stage=0 ###
errordf <- data.frame(
  x=CPA_ind_name[,2],
  y=abs((tail(as.vector(dev_list_stl[[2]][[1]]+pred_trend+seas_comp_sep),n=nnodes)- tail(data[,93],n=nnodes))/tail(data[,93],n=nnodes))
)
### FOR STL up to Aug 2019, sum GDP lag=3, stage=1 ###
errordf <- data.frame(
  x=CPA_ind_name[,2],
  y=abs((tail(as.vector(dev_list_stl[[2]][[1]]+pred_trend+seas_comp_sep),n=nnodes)- tail(data[,45],n=nnodes))/tail(data[,45],n=nnodes))
)

#errordf <- errordf[order(errordf$y),]

errordf %>% 
  arrange(y) %>%
  mutate(x=factor(x, levels=x)) %>%  
  ggplot(aes(x=x, y=y))+
  geom_segment( aes(x=x, xend=x, y=0, yend=y), color="skyblue") +
  geom_point( color="blue", size=4, alpha=0.6) +
  theme_light() +
  coord_flip() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
  )+ylab("relative error")+xlab("industry")

# Horizontal version
ggplot(errordf, aes(x=x, y=y)) +
  geom_segment( aes(x=x, xend=x, y=0, yend=y), color="skyblue") +
  geom_point( color="blue", size=4, alpha=0.6) +
  theme_light() +
  coord_flip() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  )



###### plot fot GDP sum for all lags all stages ###### 

tail(dev_list[[1]][[1]],n=nnodes)
gdp_sum <- sum(tail(data_vec,n=nnodes))
gdp_sum_mat <- matrix(nrow = 9,ncol = 5)
for (i in 1:9){
  for (j in 1:5){
    gdp_sum_mat[i,j] <- sum(tail(dev_list[[i]][[j]],n=nnodes))
  }
}

### FOR STL ###
gdp_sum <- sum(tail(data[,93],n=nnodes))
gdp_sum_mat <- matrix(nrow = 9,ncol = 5)
for (i in 1:9){
  for (j in 1:5){
    gdp_sum_mat[i,j] <- sum(tail(dev_list_stl[[i]][[j]],n=nnodes))
  }
}
# for pred trend
gdp_sum <- sum(tail(data[,45],n=nnodes))
gdp_sum_mat <- matrix(nrow = 9,ncol = 5)
for (i in 1:9){
  for (j in 1:5){
    gdp_sum_mat[i,j] <- sum(tail(as.vector(dev_list_stl[[i]][[j]]+pred_trend+seas_comp_sep),n=nnodes))
  }
}

gdp_sum_err <- abs(gdp_sum_mat-gdp_sum)
gdp_sum_rel_err <- abs((gdp_sum_mat-gdp_sum)/gdp_sum)
gdp_sum_err <- sqrt(mean((gdp_sum_mat-gdp_sum)^2))


rmse_df <- data.frame(gdp_sum_rel_err)
#rmse_df <- cbind(rmse_df,rmse_ar)
colnames(rmse_df) <- c("GNARx-0","GNARx-1","GNARx-2","GNARx-3","AR")

rmse_df_melt <- melt(rmse_df)
rmse_df_melt["lag"] <- rep(1:9,5)
colnames(rmse_df_melt) <- c("model","error","lag")
ggplot(rmse_df_melt, aes(x = lag, y = error, colour =model, group = model)) +  geom_line() + geom_point(size=1)+
  ylab("error")+
  theme(axis.text.x = element_text(size=12),axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=12),axis.title.y = element_text(size=12),legend.text = element_text(size=12),legend.title = element_text(size=12))+
  scale_x_continuous("lag", labels = as.character(1:9), breaks = 1:9)+ylab("relative error")

######## ######## ######## ######## ######## 
######## UNDO 1st DIFF and STAND ######## 
######## ######## ######## ######## ######## 

dev_list_real <- lapply(dev_list,function(x) lapply(x, function(y) undo_diff1_sd(data,data_train_sd,y)))

rmse_mat_realval <- matrix(nrow=9,ncol = 5)
rmse_mat_realval_nodes <- matrix(nrow=9,ncol = 5)
for (lagi in 1:9){
  for (stag in 1:5){
    rmse_mat_realval[lagi,stag] <- sqrt(mean((dev_list_real[[lagi]][[stag]]-data[,93])^2))
    rmse_mat_realval_nodes[lagi,stag] <- sqrt(mean((tail(dev_list_real[[lagi]][[stag]],n=nnodes)-tail(data[,93],n=nnodes))^2))
  }
}

rmse_df <- data.frame(rmse_mat_realval_nodes)
colnames(rmse_df) <- c("GNARx-0","GNARx-1","GNARx-2","GNARx-3","AR")

rmse_df_melt <- melt(rmse_df)
rmse_df_melt["lag"] <- rep(1:9,5)
colnames(rmse_df_melt) <- c("model","RMSE","lag")
ggplot(rmse_df_melt, aes(x = lag, y = RMSE, colour =model, group = model)) +  geom_line() + geom_point(size=1)+
  ylab("RMSE")+
  theme(axis.text.x = element_text(size=12),axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=12),axis.title.y = element_text(size=12),legend.text = element_text(size=12),legend.title = element_text(size=12))+
  scale_x_continuous("lag", labels = as.character(1:9), breaks = 1:9)


###### plot fot GDP sum for all lags all stages ###### 

gdp_sum <- sum(tail(data[,93],n=nnodes))
gdp_sum_mat <- matrix(nrow = 9,ncol = 5)
for (i in 1:9){
  for (j in 1:5){
    gdp_sum_mat[i,j] <- sum(tail(dev_list_real[[i]][[j]],n=nnodes))
  }
}

gdp_sum_err <- abs(gdp_sum_mat-gdp_sum)

rmse_df <- data.frame(gdp_sum_err)
#rmse_df <- cbind(rmse_df,rmse_ar)
colnames(rmse_df) <- c("GNARx-0","GNARx-1","GNARx-2","GNARx-3","AR")

rmse_df_melt <- melt(rmse_df)
rmse_df_melt["lag"] <- rep(1:9,5)
colnames(rmse_df_melt) <- c("model","error","lag")
ggplot(rmse_df_melt, aes(x = lag, y = error, colour =model, group = model)) +  geom_line() + geom_point(size=1)+
  ylab("error")+
  theme(axis.text.x = element_text(size=12),axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=12),axis.title.y = element_text(size=12),legend.text = element_text(size=12),legend.title = element_text(size=12))+
  scale_x_continuous("lag", labels = as.character(1:9), breaks = 1:9)
