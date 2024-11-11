######## ######## ######## ######## ######## ######## ######## ######## ######## 
######## CORRELATION ANALYSIS FOR NODE REMOVAL AND NETWORK SPARSIFICATION ######## 
######## ######## ######## ######## ######## ######## ######## ######## ######## 
dir_path <- "/Users/u2371456/Documents/R_directory_Turing_ONS/collab"

######## FIRST STEP: NODE REMOVAL ######## 

noderem_step1 <- function(prep_dat,node_vec){
  # prep_dat: data set with for version i of GDP release, with network structure aligned among releases
  # node_vec: vector of CPA codes corresponding to nodes 
  # return: data set with removed nodes and corresponding edges, according to node_vec vector
  
  # node labels to be removed
  node_rem <- prep_dat$CPA_node[which(prep_dat$CPA_node[,1] %in% node_vec),2]
  
  # check that df_gva uses labels from CPA_node 
  if (!all(prep_dat$CPA_node[,2]==prep_dat$df_gva[,1])){
    warning("labels of nodes in CPA_node not aligned with labels of nodes in df_gva, for prep_dat")
  }
  
  # check if rownames of ts_pay_gdp in alignment with df_gva and df_pay (REMEMBER: in Kerstin's data "to" means customer, and "from" means supplier thus reverse edge meaning)
  if(!all(rownames(prep_dat$ts_pay_gdp)==c(paste(prep_dat$df_pay$to,prep_dat$df_pay$from),prep_dat$df_gva$`CPA code`))){
    warning("labels of edges and nodes in ts_pay_gdp not aligned with df_pay and df_gva, for prep_dat")
  }
  
  new_dfpay <- prep_dat$df_pay[-which((prep_dat$df_pay[,1] %in% node_rem) | (prep_dat$df_pay[,2] %in% node_rem) ),]
  new_dfgva <- prep_dat$df_gva[-which(prep_dat$df_gva[,1] %in% node_rem),]
  
  new_dataedges <- prep_dat$data_edges[-which((prep_dat$data_edges[,1] %in% node_rem) | (prep_dat$data_edges[,2] %in% node_rem) ),]
  new_datanodes <- prep_dat$data_nodes[-which(prep_dat$data_nodes %in% node_rem)]
  new_cpanodes <- prep_dat$CPA_node[-which(prep_dat$CPA_node[,2] %in% node_rem),]
  
  new_graph <-  graph_from_data_frame(as.data.frame(new_dataedges),directed = TRUE,vertices = new_datanodes)
  new_nedges <-  ecount(new_graph)
  new_nnodes <- vcount(new_graph)
  new_nedgesnodes <-  ecount(new_graph)+vcount(new_graph)
  
  new_ts_pay_gdp <- rbind(new_dfpay[,-c(1,2)],new_dfgva[,-c(1)])
  new_ts_pay_gdp <- apply(new_ts_pay_gdp,2,as.numeric)
  rownames(new_ts_pay_gdp) <- c(paste(new_dfpay[,1],new_dfpay[,2]),new_dfgva[,1])
  
  return(list(df_pay_2=new_dfpay,df_gva_2=new_dfgva,ts_pay_gdp_2=new_ts_pay_gdp,CPA_node_2=new_cpanodes,
              graph_2=new_graph,nnodes_2=new_nnodes,nedges_2=new_nedges,
              n_edges_nodes_2=new_nedgesnodes,data_edges_2=new_dataedges,data_nodes_2=new_datanodes))
}


# Load data for all versions with aligned network structure
load("/Users/amantziou/Documents/R_directory_Turing_ONS/collab/data/data_preprocessing/data_allver_alignednet.RData")

# Node removal from economic interpretation: 
# O84 Public administration, G45 Wholesale vehicles,G46 Wholesale except vehicles, G47 Retail except vehicles, Q86 Human health,Q87-88 Residential care
node_rem_vec <- c("O84","G45","G46","G47","Q86","Q87 & Q88")

for (rev in 1:10){
  assign(paste("sparse_",rev,"_noderem",sep = ""),noderem_step1(get(paste("prep_data_rev_",rev,"uni",sep = "")),node_rem_vec))
}

# save(sparse_1_noderem,sparse_2_noderem,sparse_3_noderem,sparse_4_noderem,
#      sparse_5_noderem,sparse_6_noderem,sparse_7_noderem,sparse_8_noderem,
#      sparse_9_noderem,sparse_10_noderem,file = "/Users/amantziou/Documents/R_directory_Turing_ONS/collab/data/data_preprocessing/sparse_data_allver_noderem.RData")

# AFTER REVISED PREPROCESSING OF PREP_DATA_REV_UNI ALIGNED GRAPH VERSION OF DATA
# save(sparse_1_noderem,sparse_2_noderem,sparse_3_noderem,sparse_4_noderem,
#      sparse_5_noderem,sparse_6_noderem,sparse_7_noderem,sparse_8_noderem,
#      sparse_9_noderem,sparse_10_noderem,file = "/Users/amantziou/Documents/R_directory_Turing_ONS/collab/data/data_preprocessing/revised_preprocessing/sparse_data_allver_noderem_rev.RData")

load(paste(dir_path,"/data/data_preprocessing/revised_preprocessing/sparse_data_allver_noderem_rev.RData",sep=""))

######### STEP 2: EDGE REMOVAL ######## 
######## ACCORDING TO PEARSON'S CROSS INDUSTRY CORRELATIONS BETWEEN PAYMENTS AND GDP, ACROSS ALL RELEASES

#thres_p <- 0.4
thres_p_vec <- c(.3,.5,.6,.7)
ind_th <- 1
removepaym_pos <- list()
for (th in thres_p_vec){
  cormatrev <- list()
  paymremrev <- list()
  for (rev in 1:10){
    prep_data <- get(paste("sparse_",rev,"_noderem",sep = ""))$ts_pay_gdp_2
    
    edgerem <- pay_gdp_corr_crossind_egerem(prep_data,get(paste("sparse_",rev,"_noderem",sep = ""))$nnodes_2,get(paste("sparse_",rev,"_noderem",sep = ""))$nedges_2,th)
    
    paymremrev[[rev]] <- edgerem$paymrem
    cormatrev[[rev]] <- edgerem$cormat
    
    # hist(edgerem$cormat)
    # lowerb <- quantile(edgerem$cormat,.025)
    # upperb <- quantile(edgerem$cormat,.975)
    # abline(v=upperb,col="red")
    # abline(v=lowerb,col="red")
    #thres <- quantile(cormatrev[[rev]],.95)
  }
  
  removepaym_pos[[ind_th]] <- unique(unlist(paymremrev))
  ind_th <- ind_th+1
}


################################################
# graph related changes common for all versions
################################################

edgerem_step2 <- function(prep_dat,payrem_pos){
  # prep_dat: data set 
  # payrem_pos: positions of time series to be removed
  # return: updated data set after removal of edges (payments) and (potentially) corresponding nodes
  data_edges3 <-  prep_dat$data_edges_2[-payrem_pos,]
  # not in common nodes (after removing edges, removing also some nodes)
  removegdp_lab <- setdiff(prep_dat$data_nodes_2,unique(c(unique(data_edges3[,2]),unique(data_edges3[,1]))))
  removegdp_pos <- which(prep_dat$CPA_node_2[,2] %in% removegdp_lab)
  if (length(removegdp_pos)!=0){
    CPA_node3 <- prep_dat$CPA_node_2[-removegdp_pos,]
    data_nodes3 <- prep_dat$data_nodes_2[-removegdp_pos]
    df_gva3 <- prep_dat$df_gva_2[-removegdp_pos,]  
    ts_pay_gdp3 <-  prep_dat$ts_pay_gdp_2[-c(payrem_pos,removegdp_pos+prep_dat$nedges_2),]
  }else{
    CPA_node3 <- prep_dat$CPA_node_2
    data_nodes3 <- prep_dat$data_nodes_2
    df_gva3 <- prep_dat$df_gva_2
    ts_pay_gdp3 <-  prep_dat$ts_pay_gdp_2[-payrem_pos,]
  }
  graph3 <-  graph_from_data_frame(as.data.frame(data_edges3),directed = TRUE,vertices = data_nodes3)
  nedges3 <-  ecount(graph3)
  nnodes3 <- vcount(graph3)
  n_edges_nodes3 <-  ecount(graph3)+vcount(graph3)
  
  df_pay3 <- prep_dat$df_pay_2[-payrem_pos,]
  
  return(list(df_pay_2=df_pay3,df_gva_2=df_gva3,ts_pay_gdp_2=ts_pay_gdp3,CPA_node_2=CPA_node3,
              graph_2=graph3,nnodes_2=nnodes3,nedges_2=nedges3,
              n_edges_nodes_2=n_edges_nodes3,data_edges_2=data_edges3,data_nodes_2=data_nodes3))
}

for (th in 1:4){
  for (rev in 1:10){
    assign(paste("sparse_",rev,"_nodedgerem_thres",sep = ""),edgerem_step2(get(paste("sparse_",rev,"_noderem",sep = "")),removepaym_pos[[th]]))
  }
  save(sparse_1_nodedgerem_thres,sparse_2_nodedgerem_thres,sparse_3_nodedgerem_thres,sparse_4_nodedgerem_thres,
       sparse_5_nodedgerem_thres,sparse_6_nodedgerem_thres,sparse_7_nodedgerem_thres,sparse_8_nodedgerem_thres,
       sparse_9_nodedgerem_thres,sparse_10_nodedgerem_thres,
       file = paste("/Users/amantziou/Documents/R_directory_Turing_ONS/collab/data/data_preprocessing/revised_preprocessing/sparse_data_allver_nodedgerem_rev_pears",thres_p_vec[th],".RData"))
}


# DATA: AFTER REVISED PREPROCESSING OF PREP_DATA_REV_UNI ALIGNED GRAPH VERSION OF DATA
# save(sparse_1_nodedgerem_thres0.4,sparse_2_nodedgerem_thres0.4,sparse_3_nodedgerem_thres0.4,sparse_4_nodedgerem_thres0.4,
#      sparse_5_nodedgerem_thres0.4,sparse_6_nodedgerem_thres0.4,sparse_7_nodedgerem_thres0.4,sparse_8_nodedgerem_thres0.4,
#      sparse_9_nodedgerem_thres0.4,sparse_10_nodedgerem_thres0.4,file = "/Users/amantziou/Documents/R_directory_Turing_ONS/collab/data/data_preprocessing/revised_preprocessing/sparse_data_allver_nodedgerem_rev_pears0.4.RData")
load(paste(dir_path,"/data/data_preprocessing/revised_preprocessing/sparse_data_allver_nodedgerem_rev_pears0.4.RData",sep = ""))

#### RESULTS: GNAR and AR ####
# below  results for SPARSE network (pearson 0.4), grwoth rates and stl and evaluating predictions on original scale, after node and edge removal
# save(rmse_mat_real,
#      rmse_mat_real_nodes ,
#      bicres4,
#      pred_ver,
#      predsd_ver,
#      covresmat_list,
#      resmat_list,
#      pred_trend_list,
#      seas_comp_list,
#      data_vec_list,file = "/Users/amantziou/Documents/R_directory_Turing_ONS/collab/results/revised_preprocessing/res_public_updated13may_allversions_growthstl_sparse_0.4pearson_nodedgerem.RData")
# 
load("/Users/amantziou/Documents/R_directory_Turing_ONS/collab/results/revised_preprocessing/res_public_updated13may_allversions_growthstl_sparse_0.4pearson_nodedgerem.RData")

# below  results for SPARSE network (pearson 0.4), grwoth rates and stl and evaluating predictions on original scale, after node and edge removal, MAX LAG=12
# save(rmse_mat_real,
#      rmse_mat_real_nodes ,
#      bicres4,
#      pred_ver,
#      predsd_ver,
#      covresmat_list,
#      resmat_list,
#      pred_trend_list,
#      seas_comp_list,
#      data_vec_list,file = "/Users/amantziou/Documents/R_directory_Turing_ONS/collab/results/revised_preprocessing/res_public_updated13may_allversions_growthstl_sparse_0.4pearson_nodedgerem_maxlag12.RData")

#below  results for SPARSE network (pearson 0.4), grwoth rates and stl and evaluating predictions on original scale, after node and edge removal, INCLUDING PREDICTION INTERVALS
# save(rmse_mat_real,
#      rmse_mat_real_nodes ,
#      bicres4,
#      pred_ver,
#      predsd_ver,
#      covresmat_list,
#      resmat_list,
#      pred_trend_list,
#      seas_comp_list,
#      data_vec_list,
#      upper_or,
#      lower_or,
#      pred_incl_or,
#      upper_pre,
#      lower_pre,
#      pred_incl_pre,file = "/Users/amantziou/Documents/R_directory_Turing_ONS/collab/results/revised_preprocessing/res_public_updated13may_allversions_growthstl_sparse_0.4pearson_nodedgerem_predint.RData")
# 
load(paste(dir_path,"/results/revised_preprocessing/res_public_updated13may_allversions_growthstl_sparse_0.4pearson_nodedgerem_predint.RData",sep = ""))

# #below  results for SPARSE network (pearson 0.4), grwoth rates and stl and evaluating predictions on original scale, after node and edge removal, INCLUDING PREDICTION INTERVALS AND PRE-COVID PERIOD (TRAINING ON DATA UP TO AND INCL DEC 2019 ONLY FOR LATEST RELEASE DEC 2023)
# save(rmse_mat_real,
#      rmse_mat_real_nodes ,
#      bicres4,
#      pred_ver,
#      predsd_ver,
#      covresmat_list,
#      resmat_list,
#      pred_trend_list,
#      seas_comp_list,
#      data_vec_list,
#      upper_or,
#      lower_or,
#      pred_incl_or,
#      upper_pre,
#      lower_pre,
#      pred_incl_pre,file = paste(dir_path,"/results/revised_preprocessing/res_public_updated13may_versiondec2023_growthstl_sparse_0.4pearson_nodedgerem_predint_precovid.RData",sep = ""))


### RESULTS: AUTO.ARIMA on RAW DATA
# # below results from auto.arima on raw data (no stl, no grwoth rates) FOR SPARSE network (pearson 0.4) node and edge removal
# save(fit_arima,ord,pred_arima,file = "/Users/amantziou/Documents/R_directory_Turing_ONS/collab/results/revised_preprocessing/res_public_updated13may_allversions_arima_rawdata_sparse_0.4pearson_nodedgerem.RData")
load(paste(dir_path,"/results/revised_preprocessing/res_public_updated13may_allversions_arima_rawdata_sparse_0.4pearson_nodedgerem.RData",sep = ""))

# # below  results for SPARSE network (pearson 0.3), grwoth rates and stl and evaluating predictions on original scale, after node and edge removal
# save(rmse_mat_real,
#      rmse_mat_real_nodes ,
#      bicres4,
#      pred_ver,
#      predsd_ver,
#      covresmat_list,
#      resmat_list,
#      pred_trend_list,
#      seas_comp_list,
#      data_vec_list,file = "/Users/amantziou/Documents/R_directory_Turing_ONS/collab/results/revised_preprocessing/res_public_updated13may_allversions_growthstl_sparse_0.3pearson_nodedgerem.RData")
