#load libraries
library(ggplot2)
library(reshape)
library(dplyr)
library(tidyr)


args<-commandArgs(TRUE)
chrom<-args[1]

#read in files

#sample order for hybrids
hybr_sampIDs <- read.csv(paste0("samps_",chrom,".txt"), sep=",", header = F)

#global ancestry files for each mg
glanc_mg1 <- read.table(paste0("mg1_",chrom,"_glanc.txt"), header=F)

glanc_mg2 <- read.table(paste0("mg2_",chrom,"_glanc.txt"), header=F)

glanc_mg5 <- read.table(paste0("mg5_",chrom,"_glanc.txt"), header=F)

#dosage files for each mg
dosage_mg1 <- read.table(paste0("k_v_mg1_",chrom,".ps21.txt"), header = F )

dosage_mg2 <- read.table(paste0("k_v_mg2_",chrom,".ps21.txt"), header = F )

dosage_mg5 <- read.table(paste0("k_v_mg5_",chrom,".ps21.txt"), header = F )


#get average of each dataset

#Global Ancestry
global_anc_avg <- Reduce(`+`, mget(paste0("glanc_mg", c(1:2,5))))/3

#Dosage
dosage_avg <- Reduce(`+`, mget(paste0("dosage_mg", c(1:2,5))))/3




#format plot for global ancestry plotting
colnames(global_anc_avg) <-c("palm", "cerv")
hyb_ganc_chr1 <- global_anc_avg[51:72,]
rownames(hyb_ganc_chr1) <- hybr_sampIDs[1,]
hyb_ganc_chr1$sampID <- rownames(hyb_ganc_chr1)
hyb_ganc_chr1_m <- melt(hyb_ganc_chr1)




#format dfs for local ancestry plotting
t_dosage_avg <- as.data.frame(t(dosage_avg))
colnames(t_dosage_avg) <- hybr_sampIDs[1,]
palm_dosage_chr1 <- t_dosage_avg[c(TRUE, FALSE),]
cerv_dosage_chr1 <- t_dosage_avg[c(FALSE, TRUE),]
avg_palm_dos_chr1 <- as.data.frame(rowMeans(palm_dosage_chr1))
avg_cerv_dos_chr1 <- as.data.frame(rowMeans(cerv_dosage_chr1))
colnames(avg_palm_dos_chr1) <- c("palm")
colnames(avg_cerv_dos_chr1) <- c("cerv")
avg_palm_dos_chr1$pos <- seq.int(nrow(avg_palm_dos_chr1))
avg_cerv_dos_chr1$pos <- seq.int(nrow(avg_cerv_dos_chr1))
all_chr1_avg_dose<- merge(avg_palm_dos_chr1,avg_cerv_dos_chr1)
m_all_chr1_avg_dose <- melt(all_chr1_avg_dose, id=c("pos")) 

#sanity check on if the global averages in the log files match what you calculate from full data
global_chr1 <- as.data.frame(colMeans(palm_dosage_chr1))
colnames(global_chr1) <- c("palm_dos")
global_chr1$cerv_dos <- (colMeans(cerv_dosage_chr1))
global_chr1$palm <- global_chr1$palm_dos / 2
global_chr1$cerv <- global_chr1$cerv_dos / 2
global_chr1$match <- hyb_ganc_chr1[match(global_chr1$palm, hyb_ganc_chr1$palm), 1]

all_avg_dos_chr1 <- as.data.frame(rowMeans(t_dosage_avg))
colnames(all_avg_dos_chr1) <- c("avg_dose")
all_avg_dos_chr1$pos <- seq.int(nrow(all_avg_dos_chr1))

all_avg_dos_chr1$species <- ifelse(all_avg_dos_chr1$pos %% 2 == 0, "cerv", "palm")



#format df for local ancestry plots
palm_dosage_chr1$pos <- seq.int(nrow(palm_dosage_chr1))
cerv_dosage_chr1$pos <- seq.int(nrow(cerv_dosage_chr1))
palm_dosage_chr1$species <- "palm"
cerv_dosage_chr1$species <- "cerv"
chr1_dosage_all <- rbind(palm_dosage_chr1,cerv_dosage_chr1)
chr1_dosage_all$u_id <- seq.int(nrow(chr1_dosage_all))
m_chr1_dosage_all <- melt(chr1_dosage_all, id=c("pos","species", "u_id"), variable_name =  "ID")
um_chr1 <- unite(m_chr1_dosage_all, col = "id_spec", c("species","ID"), sep = "_", remove = F)


#set means for local ancestry plot
gd <- um_chr1 %>% 
  group_by(species, pos) %>% 
  summarise(value = mean(value))


#print local ancestry plot
pdf(paste0("dosage_avg_",chrom,".pdf"), paper = "a4r", width = 0, height = 0)
g2 <- ggplot(um_chr1, aes(x = pos, y = value, color = species)) +
  geom_line(aes(group = id_spec), alpha =0.3) +
  geom_line(data = gd, alpha = 0.8, size = 2) +
  theme_bw() +
  labs(
    title = "Local Ancestry",
    x = "Position",
    y = "Dosage",
    color = NULL
  ) + scale_colour_manual(values = c("cerv" = "#FBB116", "palm" = "#B33B77"))
print(g2)
dev.off


#print global ancestry plot
pdf(paste0("glanc_avg_",chrom,".pdf"), paper = "a4r", width = 0, height = 0)
p1 <- ggplot(hyb_ganc_chr1_m, aes(fill=variable, y=value, x=sampID)) + 
  geom_bar(position="fill", stat="identity") + theme_bw() + xlab("Sample ID") + ylab("Proportion Ancestry") + ggtitle("Global Ancestry") + scale_fill_manual(values = c("cerv" = "#FBB116", "palm" = "#B33B77"))
print(p1)
dev.off









