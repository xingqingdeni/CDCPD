#df_PELT
df_PELT1 = t(df_PELT)
True_changpoints = c(11,21,31,41,51,61,71,81,91,101,111,121,131,141)
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
per_error = c()
per_l1 = c()
mean_l1 = c()
for(i in 1:estimate_changpointstimes){
  for(j in 1:True_changpointsnum){
    a = abs(df_PELT1[c(1:Use_changpointsnum[i]),i]-True_changpoints[j])
    b = min(abs(df_PELT1[c(1:Use_changpointsnum[i]),i]-True_changpoints[j]))
    min(which(a==b))
    per_error[j] = abs(df_PELT1[min(which(a==b)),i] - True_changpoints[j])
  }
  per_l1[i] = sum(per_error)
  mean_l1[i] = per_l1[i]/(Use_changpointsnum[i])
}
mean(mean_l1)
median(mean_l1)
mad(mean_l1)
sd(mean_l1)
#df_MergedMosum
df_MergedMosum1 = t(df_MergedMosum)
True_changpoints = c(11,21,31,41,51,61,71,81,91,101,111,121,131,141)
True_changpointsnum = length(True_changpoints)
estimate_changpointsnum = df_MergedMosum1[nrow(df_MergedMosum1),]
estimate_changpointstimes = length(estimate_changpointsnum)
Use_data_nrow <- max(estimate_changpointsnum,True_changpointsnum,na.rm = TRUE)
df_MergedMosum1 = df_MergedMosum1[-c(Use_data_nrow+1:ncol(df_MergedMosum1)),]
for(i in 1:estimate_changpointstimes){
  df_MergedMosum1[,i] = sort(df_MergedMosum1[,i],na.last = TRUE)
}
Use_changpointsnum = c()
for(i in 1:estimate_changpointstimes){
  Use_changpointsnum[i] <- min(sum(!is.na(df_MergedMosum1[,i])),True_changpointsnum)
}
per_error = c()
per_l1 = c()
mean_l1 = c()
for(i in 1:estimate_changpointstimes){
  for(j in 1:True_changpointsnum){
    a = abs(df_MergedMosum1[c(1:Use_changpointsnum[i]),i]-True_changpoints[j])
    b = min(abs(df_MergedMosum1[c(1:Use_changpointsnum[i]),i]-True_changpoints[j]))
    min(which(a==b))
    per_error[j] = abs(df_MergedMosum1[min(which(a==b)),i] - True_changpoints[j])
  }
  per_l1[i] = sum(per_error)
  mean_l1[i] = per_l1[i]/(Use_changpointsnum[i])
}
mean(mean_l1)
median(mean_l1)
mad(mean_l1)
sd(mean_l1)






#df_MosumPvalue
df_MosumPvalue1 = t(df_MosumPvalue)
True_changpoints = c(11,21,31,41,51,61,71,81,91,101,111,121,131,141)
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
per_error = c()
per_l1 = c()
mean_l1 = c()
for(i in 1:estimate_changpointstimes){
  for(j in 1:True_changpointsnum){
    a = abs(df_MosumPvalue1[c(1:Use_changpointsnum[i]),i]-True_changpoints[j])
    b = min(abs(df_MosumPvalue1[c(1:Use_changpointsnum[i]),i]-True_changpoints[j]))
    min(which(a==b))
    per_error[j] = abs(df_MosumPvalue1[min(which(a==b)),i] - True_changpoints[j])
  }
  per_l1[i] = sum(per_error)
  mean_l1[i] = per_l1[i]/(Use_changpointsnum[i])
}
mean(mean_l1)
median(mean_l1)
mad(mean_l1)
sd(mean_l1)


#df_WBS
df_WBS1 = t(df_WBS)
True_changpoints = c(11,21,31,41,51,61,71,81,91,101,111,121,131,141)
True_changpointsnum = length(True_changpoints)
estimate_changpointsnum = df_WBS1[nrow(df_WBS1),]
estimate_changpointstimes = length(estimate_changpointsnum)
Use_data_nrow <- max(estimate_changpointsnum,True_changpointsnum,na.rm = TRUE)
df_WBS1 = df_WBS1[-c(Use_data_nrow+1:ncol(df_WBS1)),]
for(i in 1:estimate_changpointstimes){
  df_WBS1[,i] = sort(df_WBS1[,i],na.last = TRUE)
}
Use_changpointsnum = c()
for(i in 1:estimate_changpointstimes){
  Use_changpointsnum[i] <- min(sum(!is.na(df_WBS1[,i])),True_changpointsnum)
}
per_error = c()
per_l1 = c()
mean_l1 = c()
for(i in 1:estimate_changpointstimes){
  for(j in 1:True_changpointsnum){
    a = abs(df_WBS1[c(1:Use_changpointsnum[i]),i]-True_changpoints[j])
    b = min(abs(df_WBS1[c(1:Use_changpointsnum[i]),i]-True_changpoints[j]))
    min(which(a==b))
    per_error[j] = abs(df_WBS1[min(which(a==b)),i] - True_changpoints[j])
  }
  per_l1[i] = sum(per_error)
  mean_l1[i] = per_l1[i]/(Use_changpointsnum[i])
}
mean(mean_l1)
median(mean_l1)
mad(mean_l1)
sd(mean_l1)


#df_cumseg
df_cumseg1 = t(df_cumseg)
True_changpoints = c(11,21,31,41,51,61,71,81,91,101,111,121,131,141)
True_changpointsnum = length(True_changpoints)
estimate_changpointsnum = df_cumseg1[nrow(df_cumseg1),]
estimate_changpointstimes = length(estimate_changpointsnum)
Use_data_nrow <- max(estimate_changpointsnum,True_changpointsnum,na.rm = TRUE)
df_cumseg1 = df_cumseg1[-c(Use_data_nrow+1:ncol(df_cumseg1)),]
for(i in 1:estimate_changpointstimes){
  df_cumseg1[,i] = sort(df_cumseg1[,i],na.last = TRUE)
}
Use_changpointsnum = c()
for(i in 1:estimate_changpointstimes){
  Use_changpointsnum[i] <- min(sum(!is.na(df_cumseg1[,i])),True_changpointsnum)
}
per_error = c()
per_l1 = c()
mean_l1 = c()
for(i in 1:estimate_changpointstimes){
  for(j in 1:True_changpointsnum){
    a = abs(df_cumseg1[c(1:Use_changpointsnum[i]),i]-True_changpoints[j])
    b = min(abs(df_cumseg1[c(1:Use_changpointsnum[i]),i]-True_changpoints[j]))
    min(which(a==b))
    per_error[j] = abs(df_cumseg1[min(which(a==b)),i] - True_changpoints[j])
  }
  per_l1[i] = sum(per_error)
  mean_l1[i] = per_l1[i]/(Use_changpointsnum[i])
}
mean(mean_l1)
median(mean_l1)
mad(mean_l1)
sd(mean_l1)

#df_RANK
df_RANK1 = t(df_RANK)
True_changpoints = c(11,21,31,41,51,61,71,81,91,101,111,121,131,141)
True_changpointsnum = length(True_changpoints)
estimate_changpointsnum = df_RANK1[nrow(df_RANK1),]
estimate_changpointstimes = length(estimate_changpointsnum)
Use_data_nrow <- max(estimate_changpointsnum,True_changpointsnum,na.rm = TRUE)
df_RANK1 = df_RANK1[-c(Use_data_nrow+1:ncol(df_RANK1)),]
for(i in 1:estimate_changpointstimes){
  df_RANK1[,i] = sort(df_RANK1[,i],na.last = TRUE)
}
Use_changpointsnum = c()
for(i in 1:estimate_changpointstimes){
  Use_changpointsnum[i] <- min(sum(!is.na(df_RANK1[,i])),True_changpointsnum)
}
per_error = c()
per_l1 = c()
mean_l1 = c()
for(i in 1:estimate_changpointstimes){
  for(j in 1:True_changpointsnum){
    a = abs(df_RANK1[c(1:Use_changpointsnum[i]),i]-True_changpoints[j])
    b = min(abs(df_RANK1[c(1:Use_changpointsnum[i]),i]-True_changpoints[j]))
    min(which(a==b))
    per_error[j] = abs(df_RANK1[min(which(a==b)),i] - True_changpoints[j])
  }
  per_l1[i] = sum(per_error)
  mean_l1[i] = per_l1[i]/(Use_changpointsnum[i])
}
mean(mean_l1)
median(mean_l1)
mad(mean_l1)
sd(mean_l1)

#df_ECP
df_ECP1 = t(df_ECP)
True_changpoints = c(11,21,31,41,51,61,71,81,91,101,111,121,131,141)
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
per_error = c()
per_l1 = c()
mean_l1 = c()
for(i in 1:estimate_changpointstimes){
  for(j in 1:True_changpointsnum){
    a = abs(df_ECP1[c(1:Use_changpointsnum[i]),i]-True_changpoints[j])
    b = min(abs(df_ECP1[c(1:Use_changpointsnum[i]),i]-True_changpoints[j]))
    min(which(a==b))
    per_error[j] = abs(df_ECP1[min(which(a==b)),i] - True_changpoints[j])
  }
  per_l1[i] = sum(per_error)
  mean_l1[i] = per_l1[i]/(Use_changpointsnum[i])
}
mean(mean_l1)
median(mean_l1)
mad(mean_l1)
sd(mean_l1)

#df_NMCD
df_NMCD1 = t(df_NMCD)
True_changpoints = c(11,21,31,41,51,61,71,81,91,101,111,121,131,141)
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
per_error = c()
per_l1 = c()
mean_l1 = c()
for(i in 1:estimate_changpointstimes){
  for(j in 1:True_changpointsnum){
    a = abs(df_NMCD1[c(1:Use_changpointsnum[i]),i]-True_changpoints[j])
    b = min(abs(df_NMCD1[c(1:Use_changpointsnum[i]),i]-True_changpoints[j]))
    min(which(a==b))
    per_error[j] = abs(df_NMCD1[min(which(a==b)),i] - True_changpoints[j])
  }
  per_l1[i] = sum(per_error)
  mean_l1[i] = per_l1[i]/(Use_changpointsnum[i])
}
mean(mean_l1)
median(mean_l1)
mad(mean_l1)
sd(mean_l1)










