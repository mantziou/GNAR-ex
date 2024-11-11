as_edgelist(graph_allv)
identical_graphs(intersection(prep_data_rev_1$graph_2,prep_data_rev_2$graph_2,prep_data_rev_3$graph_2,prep_data_rev_4$graph_2,prep_data_rev_5$graph_2,prep_data_rev_6$graph_2,prep_data_rev_7$graph_2,prep_data_rev_8$graph_2,prep_data_rev_9$graph_2,prep_data_rev_10$graph_2),graph_allv)
identical_graphs(intersection(graph_rev_1,graph_rev_2,graph_rev_3,graph_rev_4,graph_rev_5,graph_rev_7,graph_rev_6,graph_rev_8,graph_rev_9,graph_rev_10),graph_allv)


for (i in 1:10){
  prep_data_local <- get(paste("prep_data_rev_",i,sep = ""))
  edgelist_diff <- as_edgelist(igraph::union(difference(prep_data_local$graph_2,graph_allv),difference(graph_allv,prep_data_local$graph_2)))
  if(dim(edgelist_diff)[1]!=0){
    edge_paste <- paste(edgelist_diff[,1],edgelist_diff[,2])
    edge_pos <- which(rownames(prep_data_local$ts_pay_gdp_2) %in% edge_paste)
    ifelse((all(edge_paste %in% paste(prep_data_local$df_pay_2[edge_pos,1],prep_data_local$df_pay_2[edge_pos,2]) ) &
            all(paste(prep_data_local$df_pay_2[edge_pos,1],prep_data_local$df_pay_2[edge_pos,2]) %in% edge_paste) &
            all(edge_paste %in% rownames(prep_data_local$ts_pay_gdp_2)[edge_pos]) &
            all( rownames(prep_data_local$ts_pay_gdp_2)[edge_pos] %in% edge_paste) &
              all(edge_paste %in% paste(prep_data_local$data_edges_2[edge_pos,1],prep_data_local$data_edges_2[edge_pos,2])) &
              all(paste(prep_data_local$data_edges_2[edge_pos,1],prep_data_local$data_edges_2[edge_pos,2]) %in% edge_paste)),print("right indexing of edges"),warning("issue with edge indexing"))
    df_pay2 <- prep_data_local$df_pay[-edge_pos,]
    data_edges2 <-  prep_data_local$data_edges[-edge_pos,]
    ts_pay_gdp2 <-  prep_data_local$ts_pay_gdp[-edge_pos,]
    #graph2 <-  graph_from_edgelist(as.matrix(data_edges2),directed = TRUE)
    graph2 <-  graph_from_data_frame(as.data.frame(data_edges2),directed = TRUE,vertices = prep_data_local$data_nodes)
    nedges2 <-  ecount(graph2)
    n_edges_nodes2 <-  ecount(graph2)+vcount(graph2)
    assign(paste("prep_data_rev_",i,"uni",sep = ""),list(df_pay=df_pay2,df_gva=prep_data_local$df_gva,ts_pay_gdp=ts_pay_gdp2,CPA_node=prep_data_local$CPA_node,graph=graph2,nnodes=prep_data_local$nnodes,nedges=nedges2,
                                                         n_edges_nodes=n_edges_nodes2,data_edges=data_edges2,data_nodes=prep_data_local$data_nodes))
  }else{
    assign(paste("prep_data_rev_",i,"uni",sep = ""),get(paste("prep_data_rev_",i,sep = "")))
    print(c("Version ",i," does not require changes"))
  }
}

# save(prep_data_rev_1uni,prep_data_rev_2uni,prep_data_rev_3uni,prep_data_rev_4uni,prep_data_rev_5uni,prep_data_rev_6uni,
#      prep_data_rev_7uni,prep_data_rev_8uni,prep_data_rev_9uni,prep_data_rev_10uni,file="/Users/amantziou/Documents/R_directory_Turing_ONS/collab/data/data_preprocessing/revised_preprocessing/data_allver_alignednet_revised.RData")
