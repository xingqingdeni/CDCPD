df_PELT = matrix(,nrow=1000,ncol=20);df_RMosumMerged = matrix(,nrow=1000,ncol=20);df_RMosumPvalue = matrix(,nrow=1000,ncol=20);
df_WBS = matrix(,nrow=1000,ncol=200);df_cumseg = matrix(,nrow=1000,ncol=200);df_RANK = matrix(,nrow=1000,ncol=20);
df_ECP = matrix(,nrow=1000,ncol=20);df_NMCD = matrix(,nrow=1000,ncol=20);df_MosumMerged = matrix(,nrow=1000,ncol=20);df_MosumPvalue = matrix(,nrow=1000,ncol=20);
options (warn = - 1)
for(i in 1:1000){
  td1 <- testData(lengths = rep(10,14),means = c(0,1,0,1,0,1,0,1,0,1,0,1,0,1),sds = rep(0,14))
  set.seed(i+3000)
  x = td1$x + 0.2*rt(140,3)
  y1=length(cpts(cpt.mean(x/sd(x),penalty = 'BIC',method = 'PELT')))
  df_PELT[i,ncol(df_PELT)]=y1
  if (y1<1)
  {df_PELT[i,]=0}
  else{df_PELT[i,c(1:y1)]=cpts(cpt.mean(x/sd(x),penalty = 'BIC',method = 'PELT'))}#PELT
  write.csv(df_PELT,"D:/shicuo/PELT.csv")
  mub = multiscale.bottomUp(rank(x),G = c(8,10,20),eta = 2/3,var.est.method = 'mosum')$cpts
  y2=length(mub)
  df_RMosumMerged[i,ncol(df_RMosumMerged)]=y2
  if (y2<1)
  {df_RMosumMerged[i,]=0}
  else{df_RMosumMerged[i,c(1:y2)] = mub}  
  write.csv(df_RMosumMerged,"D:/shicuo/Merged RMOSUM (bandwidth).csv")
  mlp = multiscale.localPrune(rank(x), G = c(8,10,20),eta = 2/3,var.est.method = 'mosum')$cpts
  y3=length(mlp)
  df_RMosumPvalue[i,ncol(df_RMosumPvalue)]=y3
  if (y3<1)
  {df_RMosumPvalue[i,]=0}
  else{df_RMosumPvalue[i,c(1:y3)]=mlp}
  write.csv(df_RMosumPvalue,"D:/shicuo/Merged RMOSUM (p-value).csv")
  w <- wbs(x)
  w.cpt = as.numeric(unlist(changepoints(w,penalty="bic.penalty")$cpt.th))
  y4=length(w.cpt)
  df_WBS[i,ncol(df_WBS)]=y4
  if (y4<1)
  {df_WBS[i,]=0}
  else
  {df_WBS[i,c(1:y4)]=w.cpt}
  write.csv(df_WBS,"D:/shicuo/WBS.csv")
  cumseg <- jumpoints(x,output="2")
  cumseg_cpts = as.numeric(unlist(cumseg[5]))
  y5=length(cumseg_cpts)
  df_cumseg[i,ncol(df_cumseg)]=y5
  if (y5<1)
  {df_cumseg[i,]=0}
  else
  {df_cumseg[i,c(1:y5)]=cumseg_cpts}
  write.csv(df_cumseg,"D:/shicuo/Cumseg.csv")
  y6=length(cpts(cpt.mean(rank(x),method = 'PELT',penalty = 'Manual',pen.value = 0.1*500^2*log(500))))
  df_RANK[i,ncol(df_RANK)]=y6
  if (y6<1)
  {df_RANK[i,]=0}
  else
  {df_RANK[i,c(1:y6)]=cpts(cpt.mean(rank(x),method = 'PELT',penalty = 'Manual',pen.value = 0.1*500^2*log(500)))}  
  write.csv(df_RANK,"D:/shicuo/RANK.csv")
  y8=length(cpts(cpt.np(x,penalty = 'SIC',minseglen = 10)))
  df_NMCD[i,ncol(df_NMCD)]=y8
  if (y8<1)
  {df_NMCD[i,]=0}
  else{df_NMCD[i,c(1:y8)]=cpts(cpt.np(x,penalty = 'SIC',minseglen = 10))}
  write.csv(df_NMCD,"D:/shicuo/NMCD.csv")
  mub = multiscale.bottomUp(x,G = c(8,10,20),eta = 2/3,var.est.method = 'mosum')$cpts
  y9=length(mub)
  df_MosumMerged[i,ncol(df_MosumMerged)]=y9
  if (y9<1)
  {df_MosumMerged[i,]=0}
  else{df_MosumMerged[i,c(1:y9)] = mub}  
  write.csv(df_MosumMerged,"D:/shicuo/Merged MOSUM (bandwidth).csv")
  mlp = multiscale.localPrune(x,G = c(8,10,20),eta = 2/3,var.est.method = 'mosum')$cpts
  y10=length(mlp)
  df_MosumPvalue[i,ncol(df_MosumPvalue)]=y10
  if (y10<1)
  {df_MosumPvalue[i,]=0}
  else{df_MosumPvalue[i,c(1:y10)]=mlp}
  write.csv(df_MosumPvalue,"D:/shicuo/Merged MOSUM (p-value).csv")
}
for(i in 1:1000){
  td1 <- testData(lengths = rep(10,14),means = c(0,1,0,1,0,1,0,1,0,1,0,1,0,1),sds = rep(0,14))
  set.seed(i+3000)
  x = td1$x + 0.2*rt(140,3)
  y7=length(e.cp3o_delta(as.matrix(x),K =15,delta=8)$estimates)
  df_ECP[i,ncol(df_ECP)]=y7
  if (y7<1)
  {df_ECP[i,]=0}
  else{df_ECP[i,c(1:y7)]=e.cp3o_delta(as.matrix(x),K =15,delta=8)$estimates}
  write.csv(df_ECP,"D:/shicuo/ECP.csv")
}
table(df_PELT[,ncol(df_PELT)] - 13)/1000
table(df_RMosumMerged[,ncol(df_RMosumMerged)] - 13)/1000
table(df_RMosumPvalue[,ncol(df_RMosumPvalue)] - 13)/1000
table(df_MosumMerged[,ncol(df_MosumMerged)] - 13)/1000
table(df_MosumPvalue[,ncol(df_MosumPvalue)] - 13)/1000
table(df_WBS[,ncol(df_WBS)] - 13)/1000
table(df_cumseg[,ncol(df_cumseg)] - 13)/1000
table(df_RANK[,ncol(df_RANK)] - 13)/1000
table(df_ECP[,ncol(df_ECP)] - 13)/1000
table(df_NMCD[,ncol(df_NMCD)] - 13)/1000



#df_PELT
df_PELT1 = t(df_PELT)
True_changpoints = c(11,21,31,41,51,61,71,81,91,101,111,121,131)
True_changpointsnum = length(True_changpoints)
estimate_changpointsnum = df_PELT1[nrow(df_PELT1),]
estimate_changpointstimes = length(estimate_changpointsnum)
Use_data_nrow <- max(estimate_changpointsnum,True_changpointsnum,na.rm = TRUE)
df_PELT1 = df_PELT1[-c(Use_data_nrow+1:ncol(df_PELT1)),]
for(i in 1:estimate_changpointstimes){
  df_PELT1[,i] = sort(df_PELT1[,i],na.last = TRUE)
}
Use_changpointsnum = c()
for(i in 1:estimate_changpointstimes){
  Use_changpointsnum[i] <- min(sum(!is.na(df_PELT1[,i])),True_changpointsnum)
}
per_l1 = c()
mean_l1 = c()
per_hausdorff = c()
for(i in 1:estimate_changpointstimes){
  per_error = c()
  for(j in 1:True_changpointsnum){
    a = abs(df_PELT1[c(1:Use_changpointsnum[i]),i]-True_changpoints[j])
    b = min(abs(df_PELT1[c(1:Use_changpointsnum[i]),i]-True_changpoints[j]))
    min(which(a==b))
    per_error[j] = abs(df_PELT1[min(which(a==b)),i] - True_changpoints[j])
  }
  per_l1[i] = sum(per_error)
  mean_l1[i] = per_l1[i]/(Use_changpointsnum[i])
}
per_hausdorff = c()
for(i in 1:estimate_changpointstimes){
  h1 = c()
  h2 = c()
  for(j in 1:True_changpointsnum){
    h1[j] = min(abs(df_PELT1[c(1:Use_changpointsnum[i]),i]-True_changpoints[j]))
  }
  for(k in 1:Use_changpointsnum[i]){
    h2[k] = min(abs(df_PELT1[c(1:Use_changpointsnum[i]),i][k]-True_changpoints))
  }
  per_hausdorff[i] = max(max(h1),max(h2))
}
median(mean_l1)
mean(mean_l1)
mad(mean_l1)
sd(mean_l1)

median(per_hausdorff)
mean(per_hausdorff)
mad(per_hausdorff)
sd(per_hausdorff)

#df_RMosumMerged
df_RMosumMerged1= t(df_RMosumMerged)
True_changpoints = c(11,21,31,41,51,61,71,81,91,101,111,121,131)
True_changpointsnum = length(True_changpoints)
estimate_changpointsnum = df_RMosumMerged1[nrow(df_RMosumMerged1),]
estimate_changpointstimes = length(estimate_changpointsnum)
Use_data_nrow <- max(estimate_changpointsnum,True_changpointsnum,na.rm = TRUE)
df_RMosumMerged1= df_RMosumMerged1[-c(Use_data_nrow+1:ncol(df_RMosumMerged1)),]
for(i in 1:estimate_changpointstimes){
  df_RMosumMerged1[,i] = sort(df_RMosumMerged1[,i],na.last = TRUE)
}
Use_changpointsnum = c()
for(i in 1:estimate_changpointstimes){
  Use_changpointsnum[i] <- min(sum(!is.na(df_RMosumMerged1[,i])),True_changpointsnum)
}
per_l1 = c()
mean_l1 = c()
per_hausdorff = c()
for(i in 1:estimate_changpointstimes){
  per_error = c()
  for(j in 1:True_changpointsnum){
    a = abs(df_RMosumMerged1[c(1:Use_changpointsnum[i]),i]-True_changpoints[j])
    b = min(abs(df_RMosumMerged1[c(1:Use_changpointsnum[i]),i]-True_changpoints[j]))
    min(which(a==b))
    per_error[j] = abs(df_RMosumMerged1[min(which(a==b)),i] - True_changpoints[j])
  }
  per_l1[i] = sum(per_error)
  mean_l1[i] = per_l1[i]/(Use_changpointsnum[i])
}
per_hausdorff = c()
for(i in 1:estimate_changpointstimes){
  h1 = c()
  h2 = c()
  for(j in 1:True_changpointsnum){
    h1[j] = min(abs(df_RMosumMerged1[c(1:Use_changpointsnum[i]),i]-True_changpoints[j]))
  }
  for(k in 1:Use_changpointsnum[i]){
    h2[k] = min(abs(df_RMosumMerged1[c(1:Use_changpointsnum[i]),i][k]-True_changpoints))
  }
  per_hausdorff[i] = max(max(h1),max(h2))
}
median(mean_l1)
mean(mean_l1)
mad(mean_l1)
sd(mean_l1)

median(per_hausdorff)
mean(per_hausdorff)
mad(per_hausdorff)
sd(per_hausdorff)
#df_RMosumPvalue
df_RMosumPvalue1= t(df_RMosumPvalue)
True_changpoints = c(11,21,31,41,51,61,71,81,91,101,111,121,131)
True_changpointsnum = length(True_changpoints)
estimate_changpointsnum = df_RMosumPvalue1[nrow(df_RMosumPvalue1),]
estimate_changpointstimes = length(estimate_changpointsnum)
Use_data_nrow <- max(estimate_changpointsnum,True_changpointsnum,na.rm = TRUE)
df_RMosumPvalue1= df_RMosumPvalue1[-c(Use_data_nrow+1:ncol(df_RMosumPvalue1)),]
for(i in 1:estimate_changpointstimes){
  df_RMosumPvalue1[,i] = sort(df_RMosumPvalue1[,i],na.last = TRUE)
}
Use_changpointsnum = c()
for(i in 1:estimate_changpointstimes){
  Use_changpointsnum[i] <- min(sum(!is.na(df_RMosumPvalue1[,i])),True_changpointsnum)
}
per_l1 = c()
mean_l1 = c()
per_hausdorff = c()
for(i in 1:estimate_changpointstimes){
  per_error = c()
  for(j in 1:True_changpointsnum){
    a = abs(df_RMosumPvalue1[c(1:Use_changpointsnum[i]),i]-True_changpoints[j])
    b = min(abs(df_RMosumPvalue1[c(1:Use_changpointsnum[i]),i]-True_changpoints[j]))
    min(which(a==b))
    per_error[j] = abs(df_RMosumPvalue1[min(which(a==b)),i] - True_changpoints[j])
  }
  per_l1[i] = sum(per_error)
  mean_l1[i] = per_l1[i]/(Use_changpointsnum[i])
}
per_hausdorff = c()
for(i in 1:estimate_changpointstimes){
  h1 = c()
  h2 = c()
  for(j in 1:True_changpointsnum){
    h1[j] = min(abs(df_RMosumPvalue1[c(1:Use_changpointsnum[i]),i]-True_changpoints[j]))
  }
  for(k in 1:Use_changpointsnum[i]){
    h2[k] = min(abs(df_RMosumPvalue1[c(1:Use_changpointsnum[i]),i][k]-True_changpoints))
  }
  per_hausdorff[i] = max(max(h1),max(h2))
}
median(mean_l1)
mean(mean_l1)
mad(mean_l1)
sd(mean_l1)

median(per_hausdorff)
mean(per_hausdorff)
mad(per_hausdorff)
sd(per_hausdorff)

#df_MosumMerged
df_MosumMerged1 = t(df_MosumMerged)
True_changpoints = c(11,21,31,41,51,61,71,81,91,101,111,121,131)
True_changpointsnum = length(True_changpoints)
estimate_changpointsnum = df_MosumMerged1[nrow(df_MosumMerged1),]
estimate_changpointstimes = length(estimate_changpointsnum)
Use_data_nrow <- max(estimate_changpointsnum,True_changpointsnum,na.rm = TRUE)
df_MosumMerged1 = df_MosumMerged1[-c(Use_data_nrow+1:ncol(df_MosumMerged1)),]
for(i in 1:estimate_changpointstimes){
  df_MosumMerged1[,i] = sort(df_MosumMerged1[,i],na.last = TRUE)
}
Use_changpointsnum = c()
for(i in 1:estimate_changpointstimes){
  Use_changpointsnum[i] <- min(sum(!is.na(df_MosumMerged1[,i])),True_changpointsnum)
}
per_l1 = c()
mean_l1 = c()
per_hausdorff = c()
for(i in 1:estimate_changpointstimes){
  per_error = c()
  for(j in 1:True_changpointsnum){
    a = abs(df_MosumMerged1[c(1:Use_changpointsnum[i]),i]-True_changpoints[j])
    b = min(abs(df_MosumMerged1[c(1:Use_changpointsnum[i]),i]-True_changpoints[j]))
    min(which(a==b))
    per_error[j] = abs(df_MosumMerged1[min(which(a==b)),i] - True_changpoints[j])
  }
  per_l1[i] = sum(per_error)
  mean_l1[i] = per_l1[i]/(Use_changpointsnum[i])
}
per_hausdorff = c()
for(i in 1:estimate_changpointstimes){
  h1 = c()
  h2 = c()
  for(j in 1:True_changpointsnum){
    h1[j] = min(abs(df_MosumMerged1[c(1:Use_changpointsnum[i]),i]-True_changpoints[j]))
  }
  for(k in 1:Use_changpointsnum[i]){
    h2[k] = min(abs(df_MosumMerged1[c(1:Use_changpointsnum[i]),i][k]-True_changpoints))
  }
  per_hausdorff[i] = max(max(h1),max(h2))
}
median(mean_l1)
mean(mean_l1)
mad(mean_l1)
sd(mean_l1)

median(per_hausdorff)
mean(per_hausdorff)
mad(per_hausdorff)
sd(per_hausdorff)

#df_MosumPvalue
df_MosumPvalue1 = t(df_MosumPvalue)
True_changpoints = c(11,21,31,41,51,61,71,81,91,101,111,121,131)
True_changpointsnum = length(True_changpoints)
estimate_changpointsnum = df_MosumPvalue1[nrow(df_MosumPvalue1),]
estimate_changpointstimes = length(estimate_changpointsnum)
Use_data_nrow <- max(estimate_changpointsnum,True_changpointsnum,na.rm = TRUE)
df_MosumPvalue1 = df_MosumPvalue1[-c(Use_data_nrow+1:ncol(df_MosumPvalue1)),]
for(i in 1:estimate_changpointstimes){
  df_MosumPvalue1[,i] = sort(df_MosumPvalue1[,i],na.last = TRUE)
}
Use_changpointsnum = c()
for(i in 1:estimate_changpointstimes){
  Use_changpointsnum[i] <- min(sum(!is.na(df_MosumPvalue1[,i])),True_changpointsnum)
}
per_l1 = c()
mean_l1 = c()
per_hausdorff = c()
for(i in 1:estimate_changpointstimes){
  per_error = c()
  for(j in 1:True_changpointsnum){
    a = abs(df_MosumPvalue1[c(1:Use_changpointsnum[i]),i]-True_changpoints[j])
    b = min(abs(df_MosumPvalue1[c(1:Use_changpointsnum[i]),i]-True_changpoints[j]))
    min(which(a==b))
    per_error[j] = abs(df_MosumPvalue1[min(which(a==b)),i] - True_changpoints[j])
  }
  per_l1[i] = sum(per_error)
  mean_l1[i] = per_l1[i]/(Use_changpointsnum[i])
}
per_hausdorff = c()
for(i in 1:estimate_changpointstimes){
  h1 = c()
  h2 = c()
  for(j in 1:True_changpointsnum){
    h1[j] = min(abs(df_MosumPvalue1[c(1:Use_changpointsnum[i]),i]-True_changpoints[j]))
  }
  for(k in 1:Use_changpointsnum[i]){
    h2[k] = min(abs(df_MosumPvalue1[c(1:Use_changpointsnum[i]),i][k]-True_changpoints))
  }
  per_hausdorff[i] = max(max(h1),max(h2))
}
median(mean_l1)
mean(mean_l1)
mad(mean_l1)
sd(mean_l1)

median(per_hausdorff)
mean(per_hausdorff)
mad(per_hausdorff)
sd(per_hausdorff)

#df_WBS
df_WBS1= t(df_WBS)
True_changpoints = c(11,21,31,41,51,61,71,81,91,101,111,121,131)
True_changpointsnum = length(True_changpoints)
estimate_changpointsnum = df_WBS1[nrow(df_WBS1),]
estimate_changpointstimes = length(estimate_changpointsnum)
Use_data_nrow <- max(estimate_changpointsnum,True_changpointsnum,na.rm = TRUE)
df_WBS1= df_WBS1[-c(Use_data_nrow+1:ncol(df_WBS1)),]
for(i in 1:estimate_changpointstimes){
  df_WBS1[,i] = sort(df_WBS1[,i],na.last = TRUE)
}
Use_changpointsnum = c()
for(i in 1:estimate_changpointstimes){
  Use_changpointsnum[i] <- min(sum(!is.na(df_WBS1[,i])),True_changpointsnum)
}
per_l1 = c()
mean_l1 = c()
per_hausdorff = c()
for(i in 1:estimate_changpointstimes){
  per_error = c()
  for(j in 1:True_changpointsnum){
    a = abs(df_WBS1[c(1:Use_changpointsnum[i]),i]-True_changpoints[j])
    b = min(abs(df_WBS1[c(1:Use_changpointsnum[i]),i]-True_changpoints[j]))
    min(which(a==b))
    per_error[j] = abs(df_WBS1[min(which(a==b)),i] - True_changpoints[j])
  }
  per_l1[i] = sum(per_error)
  mean_l1[i] = per_l1[i]/(Use_changpointsnum[i])
}
per_hausdorff = c()
for(i in 1:estimate_changpointstimes){
  h1 = c()
  h2 = c()
  for(j in 1:True_changpointsnum){
    h1[j] = min(abs(df_WBS1[c(1:Use_changpointsnum[i]),i]-True_changpoints[j]))
  }
  for(k in 1:Use_changpointsnum[i]){
    h2[k] = min(abs(df_WBS1[c(1:Use_changpointsnum[i]),i][k]-True_changpoints))
  }
  per_hausdorff[i] = max(max(h1),max(h2))
}
median(mean_l1)
mean(mean_l1)
mad(mean_l1)
sd(mean_l1)

median(per_hausdorff)
mean(per_hausdorff)
mad(per_hausdorff)
sd(per_hausdorff)
#df_cumseg
df_cumseg1= t(df_cumseg)
True_changpoints = c(11,21,31,41,51,61,71,81,91,101,111,121,131)
True_changpointsnum = length(True_changpoints)
estimate_changpointsnum = df_cumseg1[nrow(df_cumseg1),]
estimate_changpointstimes = length(estimate_changpointsnum)
Use_data_nrow <- max(estimate_changpointsnum,True_changpointsnum,na.rm = TRUE)
df_cumseg1= df_cumseg1[-c(Use_data_nrow+1:ncol(df_cumseg1)),]
for(i in 1:estimate_changpointstimes){
  df_cumseg1[,i] = sort(df_cumseg1[,i],na.last = TRUE)
}
Use_changpointsnum = c()
for(i in 1:estimate_changpointstimes){
  Use_changpointsnum[i] <- min(sum(!is.na(df_cumseg1[,i])),True_changpointsnum)
}
per_l1 = c()
mean_l1 = c()
per_hausdorff = c()
for(i in 1:estimate_changpointstimes){
  per_error = c()
  for(j in 1:True_changpointsnum){
    a = abs(df_cumseg1[c(1:Use_changpointsnum[i]),i]-True_changpoints[j])
    b = min(abs(df_cumseg1[c(1:Use_changpointsnum[i]),i]-True_changpoints[j]))
    min(which(a==b))
    per_error[j] = abs(df_cumseg1[min(which(a==b)),i] - True_changpoints[j])
  }
  per_l1[i] = sum(per_error)
  mean_l1[i] = per_l1[i]/(Use_changpointsnum[i])
}
per_hausdorff = c()
for(i in 1:estimate_changpointstimes){
  h1 = c()
  h2 = c()
  for(j in 1:True_changpointsnum){
    h1[j] = min(abs(df_cumseg1[c(1:Use_changpointsnum[i]),i]-True_changpoints[j]))
  }
  for(k in 1:Use_changpointsnum[i]){
    h2[k] = min(abs(df_cumseg1[c(1:Use_changpointsnum[i]),i][k]-True_changpoints))
  }
  per_hausdorff[i] = max(max(h1),max(h2))
}
median(mean_l1)
mean(mean_l1)
mad(mean_l1)
sd(mean_l1)

median(per_hausdorff)
mean(per_hausdorff)
mad(per_hausdorff)
sd(per_hausdorff)

#df_RANK
df_RANK1= t(df_RANK)
True_changpoints = c(11,21,31,41,51,61,71,81,91,101,111,121,131)
True_changpointsnum = length(True_changpoints)
estimate_changpointsnum = df_RANK1[nrow(df_RANK1),]
estimate_changpointstimes = length(estimate_changpointsnum)
Use_data_nrow <- max(estimate_changpointsnum,True_changpointsnum,na.rm = TRUE)
df_RANK1= df_RANK1[-c(Use_data_nrow+1:ncol(df_RANK1)),]
for(i in 1:estimate_changpointstimes){
  df_RANK1[,i] = sort(df_RANK1[,i],na.last = TRUE)
}
Use_changpointsnum = c()
for(i in 1:estimate_changpointstimes){
  Use_changpointsnum[i] <- min(sum(!is.na(df_RANK1[,i])),True_changpointsnum)
}
per_l1 = c()
mean_l1 = c()
per_hausdorff = c()
for(i in 1:estimate_changpointstimes){
  per_error = c()
  for(j in 1:True_changpointsnum){
    a = abs(df_RANK1[c(1:Use_changpointsnum[i]),i]-True_changpoints[j])
    b = min(abs(df_RANK1[c(1:Use_changpointsnum[i]),i]-True_changpoints[j]))
    min(which(a==b))
    per_error[j] = abs(df_RANK1[min(which(a==b)),i] - True_changpoints[j])
  }
  per_l1[i] = sum(per_error)
  mean_l1[i] = per_l1[i]/(Use_changpointsnum[i])
}
per_hausdorff = c()
for(i in 1:estimate_changpointstimes){
  h1 = c()
  h2 = c()
  for(j in 1:True_changpointsnum){
    h1[j] = min(abs(df_RANK1[c(1:Use_changpointsnum[i]),i]-True_changpoints[j]))
  }
  for(k in 1:Use_changpointsnum[i]){
    h2[k] = min(abs(df_RANK1[c(1:Use_changpointsnum[i]),i][k]-True_changpoints))
  }
  per_hausdorff[i] = max(max(h1),max(h2))
}
median(mean_l1)
mean(mean_l1)
mad(mean_l1)
sd(mean_l1)

median(per_hausdorff)
mean(per_hausdorff)
mad(per_hausdorff)
sd(per_hausdorff)


#df_ECP
df_ECP1 = t(df_ECP)
True_changpoints = c(11,21,31,41,51,61,71,81,91,101,111,121,131)
True_changpointsnum = length(True_changpoints)
estimate_changpointsnum = df_ECP1[nrow(df_ECP1),]
estimate_changpointstimes = length(estimate_changpointsnum)
Use_data_nrow <- max(estimate_changpointsnum,True_changpointsnum,na.rm = TRUE)
df_ECP1 = df_ECP1[-c(Use_data_nrow+1:ncol(df_ECP1)),]
for(i in 1:estimate_changpointstimes){
  df_ECP1[,i] = sort(df_ECP1[,i],na.last = TRUE)
}
Use_changpointsnum = c()
for(i in 1:estimate_changpointstimes){
  Use_changpointsnum[i] <- min(sum(!is.na(df_ECP1[,i])),True_changpointsnum)
}
per_l1 = c()
mean_l1 = c()
per_hausdorff = c()
for(i in 1:estimate_changpointstimes){
  per_error = c()
  for(j in 1:True_changpointsnum){
    a = abs(df_ECP1[c(1:Use_changpointsnum[i]),i]-True_changpoints[j])
    b = min(abs(df_ECP1[c(1:Use_changpointsnum[i]),i]-True_changpoints[j]))
    min(which(a==b))
    per_error[j] = abs(df_ECP1[min(which(a==b)),i] - True_changpoints[j])
  }
  per_l1[i] = sum(per_error)
  mean_l1[i] = per_l1[i]/(Use_changpointsnum[i])
}
per_hausdorff = c()
for(i in 1:estimate_changpointstimes){
  h1 = c()
  h2 = c()
  for(j in 1:True_changpointsnum){
    h1[j] = min(abs(df_ECP1[c(1:Use_changpointsnum[i]),i]-True_changpoints[j]))
  }
  for(k in 1:Use_changpointsnum[i]){
    h2[k] = min(abs(df_ECP1[c(1:Use_changpointsnum[i]),i][k]-True_changpoints))
  }
  per_hausdorff[i] = max(max(h1),max(h2))
}
median(mean_l1)
mean(mean_l1)
mad(mean_l1)
sd(mean_l1)

median(per_hausdorff)
mean(per_hausdorff)
mad(per_hausdorff)
sd(per_hausdorff)

#df_NMCD
df_NMCD1 = t(df_NMCD)
True_changpoints = c(11,21,31,41,51,61,71,81,91,101,111,121,131)
True_changpointsnum = length(True_changpoints)
estimate_changpointsnum = df_NMCD1[nrow(df_NMCD1),]
estimate_changpointstimes = length(estimate_changpointsnum)
Use_data_nrow <- max(estimate_changpointsnum,True_changpointsnum,na.rm = TRUE)
df_NMCD1 = df_NMCD1[-c(Use_data_nrow+1:ncol(df_NMCD1)),]
for(i in 1:estimate_changpointstimes){
  df_NMCD1[,i] = sort(df_NMCD1[,i],na.last = TRUE)
}
Use_changpointsnum = c()
for(i in 1:estimate_changpointstimes){
  Use_changpointsnum[i] <- min(sum(!is.na(df_NMCD1[,i])),True_changpointsnum)
}
per_l1 = c()
mean_l1 = c()
per_hausdorff = c()
for(i in 1:estimate_changpointstimes){
  per_error = c()
  for(j in 1:True_changpointsnum){
    a = abs(df_NMCD1[c(1:Use_changpointsnum[i]),i]-True_changpoints[j])
    b = min(abs(df_NMCD1[c(1:Use_changpointsnum[i]),i]-True_changpoints[j]))
    min(which(a==b))
    per_error[j] = abs(df_NMCD1[min(which(a==b)),i] - True_changpoints[j])
  }
  per_l1[i] = sum(per_error)
  mean_l1[i] = per_l1[i]/(Use_changpointsnum[i])
}
per_hausdorff = c()
for(i in 1:estimate_changpointstimes){
  h1 = c()
  h2 = c()
  for(j in 1:True_changpointsnum){
    h1[j] = min(abs(df_NMCD1[c(1:Use_changpointsnum[i]),i]-True_changpoints[j]))
  }
  for(k in 1:Use_changpointsnum[i]){
    h2[k] = min(abs(df_NMCD1[c(1:Use_changpointsnum[i]),i][k]-True_changpoints))
  }
  per_hausdorff[i] = max(max(h1),max(h2))
}
median(mean_l1)
mean(mean_l1)
mad(mean_l1)
sd(mean_l1)

median(per_hausdorff)
mean(per_hausdorff)
mad(per_hausdorff)
sd(per_hausdorff)

