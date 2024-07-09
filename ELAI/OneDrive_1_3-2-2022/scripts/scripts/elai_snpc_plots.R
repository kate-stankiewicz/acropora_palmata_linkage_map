#elai plots for snp chip samples 
library(ggplot2)
library(reshape)
library(dplyr)
library(tidyr)

args<-commandArgs(TRUE)
chrom<-args[1]
mg<-args[2]

#read in file for global ancestry
global_anc_chr1 <- read.table(paste0("glanc_nou_mg",mg,"_",chrom,".log.txt"), header=F)
colnames(global_anc_chr1) <-c("palm", "cerv")
hyb_ganc_chr1 <- global_anc_chr1[204:263,]

#read in file for hybrid sample IDs
hybr_sampIDs <- read.csv(paste0("samps_snp_",chrom,".txt"), sep=",", header = F, check.names = F)
hyb_ganc_chr1$sampID <- t(hybr_sampIDs)
hyb_ganc_chr1_m <- melt(hyb_ganc_chr1, id=c("sampID"))

#read in the dosage file
dosage_chr1 <- read.table(paste0("nou_mg",mg,"_",chrom,".ps21.txt"), header = F )

#process it for plotting
t_dosage_chr1 <- as.data.frame(t(dosage_chr1))
colnames(t_dosage_chr1) <- hybr_sampIDs[1,]
palm_dosage_chr1 <- t_dosage_chr1[c(TRUE, FALSE),]
cerv_dosage_chr1 <- t_dosage_chr1[c(FALSE, TRUE),]
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

#don't need this part anymore
all_avg_dos_chr1 <- as.data.frame(rowMeans(t_dosage_chr1))
colnames(all_avg_dos_chr1) <- c("avg_dose")
all_avg_dos_chr1$pos <- seq.int(nrow(all_avg_dos_chr1))

all_avg_dos_chr1$species <- ifelse(all_avg_dos_chr1$pos %% 2 == 0, "cerv", "palm")


#set up processing for avg plot
palm_dosage_chr1$pos <- seq.int(nrow(palm_dosage_chr1))
cerv_dosage_chr1$pos <- seq.int(nrow(cerv_dosage_chr1))
palm_dosage_chr1$species <- "palm"
cerv_dosage_chr1$species <- "cerv"
chr1_dosage_all <- rbind(palm_dosage_chr1,cerv_dosage_chr1)
chr1_dosage_all$u_id <- seq.int(nrow(chr1_dosage_all))

m_chr1_dosage_all <- melt(chr1_dosage_all, id=c("pos","species", "u_id"), variable_name =  "ID")
um_chr1 <- unite(m_chr1_dosage_all, col = "id_spec", c("species","ID"), sep = "_", remove = F)

#now make plots

#global ancestry
pdf(paste0("glanc_snp_mg",mg,"_",chrom,".pdf"), paper = "a4r", width = 0, height = 0)
p1 <- ggplot(hyb_ganc_chr1_m, aes(fill=variable, y=value, x=sampID)) + 
  geom_bar(position="fill", stat="identity") + theme_bw() + xlab("Sample ID") + ylab("Proportion Ancestry") + 
ggtitle("Global proportion ancestry") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
scale_fill_manual(values = c("cerv" = "#FBB116", "palm" = "#B33B77"))
print(p1)
dev.off()

#set up avg for avg plot
gd <- um_chr1 %>% 
  group_by(species, pos) %>% 
  summarise(value = mean(value))


#plot loacl ancestry with averages bolded
pdf(paste0("dosage_snp_mg",mg,"_",chrom,".pdf"), paper = "a4r", width = 0, height = 0)
g1 <- ggplot(um_chr1, aes(x = pos, y = value, color = species)) +
  geom_line(aes(group = id_spec), alpha =0.3) +
  geom_line(data = gd, alpha = 0.8, size = 2) +
  theme_bw() +
  labs(
    title = "Dosage",
    x = "Position",
    y = "Dosage",
    color = NULL
  ) + scale_colour_manual(values = c("cerv" = "#FBB116", "palm" = "#B33B77"))
print(g1)
dev.off()
