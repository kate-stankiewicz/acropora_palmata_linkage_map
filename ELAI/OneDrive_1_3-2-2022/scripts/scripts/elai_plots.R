library(ggplot2)
library(reshape)
library(dplyr)
library(tidyr)

args<-commandArgs(TRUE)
chrom<-args[1]
mg<-args[2]

#read in files for global ancestry
global_anc_chr1 <- read.table(paste0(mg,"_",chrom,"_glanc.txt"), header=F)
colnames(global_anc_chr1) <-c("palm", "cerv")

#add hybrid IDs
hyb_ganc_chr1 <- global_anc_chr1[51:72,]
hybr_sampIDs <- read.csv(paste0("samps_",chrom,".txt"), sep=",", header = F)
rownames(hyb_ganc_chr1) <- hybr_sampIDs[1,]
hyb_ganc_chr1$sampID <- rownames(hyb_ganc_chr1)
hyb_ganc_chr1_m <- melt(hyb_ganc_chr1)

#read in the dosage file
dosage_chr1 <- read.table(paste0("k_v_",mg,"_",chrom,".ps21.txt"), header = F )

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

#process per species, per individual
palm_dosage_chr1$pos <- seq.int(nrow(palm_dosage_chr1))
cerv_dosage_chr1$pos <- seq.int(nrow(cerv_dosage_chr1))
palm_dosage_chr1$species <- "palm"
cerv_dosage_chr1$species <- "cerv"
chr1_dosage_all <- rbind(palm_dosage_chr1,cerv_dosage_chr1)
chr1_dosage_all$u_id <- seq.int(nrow(chr1_dosage_all))

m_chr1_dosage_all <- melt(chr1_dosage_all, id=c("pos","species", "u_id"), variable_name =  "ID")
um_chr1 <- unite(m_chr1_dosage_all, col = "id_spec", c("species","ID"), sep = "_", remove = F)

#now plot

#set the averaged
gd <- um_chr1 %>% 
  group_by(species, pos) %>% 
  summarise(value = mean(value))

#write out the plot
pdf(paste0("dosage_all_",mg,"_",chrom,".pdf"), paper = "a4r", width = 0, height = 0)
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

#now plot the global ancestry
pdf(paste0("glanc_",mg,"_",chrom,".pdf"), paper = "a4r", width = 0, height = 0)
g2<-ggplot(hyb_ganc_chr1_m, aes(fill=variable, y=value, x=sampID)) + 
  geom_bar(position="fill", stat="identity") + theme_bw() + xlab("Sample ID") + ylab("Proportion Ancestry") + ggtitle("Global ancestry")
ggplot(data=m_all_chr1_avg_dose, aes(x=pos, y=value, color=variable)) + geom_line() + ggtitle("Average allele dosage") + scale_fill_manual(values = c("cerv" = "#FBB116", "palm" = "#B33B77"))
print(g2)
dev.off()
