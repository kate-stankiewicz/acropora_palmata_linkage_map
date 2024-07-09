args<-commandArgs(TRUE)
chrom<-args[1]
mg <- args[2]

#read in files for all replicates
ldf <- list()
all_dosage <- list.files(pattern = paste0(chrom,"_",mg,"_",".*\\.ps21.txt"))
for (k in 1:length(all_dosage)){
  ldf[[k]] <- read.table(all_dosage[k], header = F)
}

#get the average of the replicates
avg <- Reduce(`+`, ldf) / length(ldf)

#write out the the result as a table
write.table(avg, paste0("avg_hic_scaffold_",chrom,"_mg",mg,".txt"), quote = F)
