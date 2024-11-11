######### add to already existing lists results from lag=10,11,12
for (rev in 1:9){
  rmse_mat_real[[rev]] <- rbind(rmse_mat_real[[rev]],rmse_mat_real2[[rev]][10:12,])
  rmse_mat_real_nodes[[rev]] <- rbind(rmse_mat_real_nodes[[rev]],rmse_mat_real_nodes2[[rev]][10:12,])
  bicres4[[rev]] <- rbind(bicres4[[rev]],bicres42[[rev]][10:12,])
}
pred_ver3 <- lapply(1:length(1:9), function(x) vector("list",12))
predsd_ver3 <- lapply(1:length(1:9), function(x) vector("list",12))
covresmat_list3 <- resmat_list3 <- lapply(1:length(1:9), function(x) vector("list",12))

for (rev in 1:9){
  for (lagi in 1:12){
    for (stag in 1:5){
      if (lagi<=9){
        pred_ver3[[rev]][[lagi]][[stag]] <- pred_ver[[rev]][[lagi]][[stag]]
        
        if (stag<=4){
          predsd_ver3[[rev]][[lagi]][[stag]] <- predsd_ver[[rev]][[lagi]][[stag]]
          covresmat_list3[[rev]][[lagi]][[stag]] <- covresmat_list[[rev]][[lagi]][[stag]]
          resmat_list3[[rev]][[lagi]][[stag]] <- resmat_list[[rev]][[lagi]][[stag]]
        }
      }else{
        pred_ver3[[rev]][[lagi]][[stag]] <- pred_ver2[[rev]][[lagi]][[stag]]
        
        if (stag<=4){
          predsd_ver3[[rev]][[lagi]][[stag]] <- predsd_ver2[[rev]][[lagi]][[stag]]
          covresmat_list3[[rev]][[lagi]][[stag]] <- covresmat_list2[[rev]][[lagi]][[stag]]
          resmat_list3[[rev]][[lagi]][[stag]] <- resmat_list2[[rev]][[lagi]][[stag]]
        }
      }
    }
  }
}

pred_ver <- pred_ver3
predsd_ver <- predsd_ver3
covresmat_list <- covresmat_list3
resmat_list <- resmat_list3

for (i in 1:9){
  print(all(data_vec_list[[i]]==data_vec_list2[[i]]))
  print(all(seas_comp_list[[i]]==seas_comp_list2[[i]]))
  print(all(pred_trend_list[[i]]==pred_trend_list2[[i]]))
}