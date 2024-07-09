library(ggplot2)
library(reshape)
library(dplyr)
library(tidyr)
library(stringr)
library(ggpubr)

#chrom sizes
chrS<-read.table("E:/PSU/Hybridization Project/Genome_paper/Analysis/Circos_plot/AP_chrom.txt",
                 check.names=FALSE, header=T, 
                 na.strings=c("", "NA"), stringsAsFactors=FALSE, sep="\t",quote="")
chrom<-chrS$chr
chrom<-c("hic_scaffold_15","hic_scaffold_17", "hic_scaffold_20", "hic_scaffold_21", "hic_scaffold_30", "hic_scaffold_31", "hic_scaffold_35")

mg<-c(2)

#SNPchip samples
sc<-c("E:/PSU/Hybridization Project/Genome_paper/Analysis/HybridCalling/ELAI/apalm_snpchip_plots/apalm_snpchip/")
#genome samples
gen<-c("E:/PSU/Hybridization Project/Genome_paper/Analysis/HybridCalling/ELAI/apalm_genome_plots/apalm_genome/")

repeats<-c("13704","13706","13736","13748","13762","1834","13764","1837","1839",
           "13782","13821","13823","13833","13843","13807","13778","6779","6791",
           "8939","13913","1303","1667")

repeats<-c("1667","1839")

for(i in chrom){
  print(i)
  for(j in mg){
    print(j)
    for (k in repeats){
      print(k)
    #read in snpchip data
    hybr_sampIDs<-read.csv(paste0(sc,"/samps_snp_",i,".txt"), sep=",", header = F)
    thy<-as.data.frame(t(hybr_sampIDs[1,]))
    colnames(thy)<-c("ID")
    thy$ID<-as.character(thy$ID)
    dosage_chr<-read.table(paste0(sc,"rr_",j,"_",i,".ps21.txt"), header = F)
    pos_chr <- read.table(paste0(sc,"rr_",j,"_",i,".snpinfo.txt"), header = T)
    thy$ID<-str_replace(thy$ID, "x", ".")
    t_dosage_chr <- as.data.frame(t(dosage_chr))
    colnames(t_dosage_chr) <- thy$ID
    
    #separate out the bc and hybrids
    rep_hybrid_sc <- as.data.frame(t_dosage_chr[, k])
    
    #process the bc palms for plotting, get average dosage
    rep_hybrid_sc_palm_dosage_chr <- as.data.frame(rep_hybrid_sc[c(TRUE, FALSE),])
    colnames(rep_hybrid_sc_palm_dosage_chr)<-k
    rep_hybrid_sc_cerv_dosage_chr <- as.data.frame(rep_hybrid_sc[c(FALSE, TRUE),])
    colnames(rep_hybrid_sc_cerv_dosage_chr)<-k
  
    #process bc palms per species, per individual
    rep_hybrid_sc_palm_dosage_chr$pos <- pos_chr$pos
    rep_hybrid_sc_cerv_dosage_chr$pos <- pos_chr$pos
    rep_hybrid_sc_palm_dosage_chr$species <- "palm"
    rep_hybrid_sc_cerv_dosage_chr$species <- "cerv"
    rep_hybrid_sc_chr_dosage_all <- rbind(rep_hybrid_sc_palm_dosage_chr,rep_hybrid_sc_cerv_dosage_chr)
    rep_hybrid_sc_chr_dosage_all$u_id <- seq.int(nrow(rep_hybrid_sc_chr_dosage_all))
    
    rep_hybrid_sc_m_chr_dosage_all <- melt(rep_hybrid_sc_chr_dosage_all, id=c("pos","species", "u_id"), variable_name =  "ID")
    rep_hybrid_sc_um_chr <- unite(rep_hybrid_sc_m_chr_dosage_all, col = "id_spec", c("species","ID"), sep = "_", remove = F)
    rep_hybrid_sc_um_chr$method <- "snpchip"
    
    #read in genome data
    hybr_sampIDs_gen<-read.csv(paste0(gen,"samps_",i,".txt"), sep=",", header = F)
    thy_gen<-as.data.frame(t(hybr_sampIDs_gen[1,]))
    colnames(thy_gen)<-c("ID")
    thy_gen$ID<-as.character(thy_gen$ID)
    dosage_chr_gen<-read.table(paste0(gen,"/rr__",j,"_",i,".ps21.txt"), header = F)
    pos_chr_gen <- read.table(paste0(gen,"/rr__",j,"_",i,".snpinfo.txt"), header = T)
    t_dosage_chr_gen <- as.data.frame(t(dosage_chr_gen))
    colnames(t_dosage_chr_gen) <- thy_gen$ID
    
    #separate out the bc and hybrids
    rep_hybrid_gen <- as.data.frame(t_dosage_chr_gen[, k])
    
    #process the bc palms for plotting, get average dosage
    rep_hybrid_gen_palm_dosage_chr <- as.data.frame(rep_hybrid_gen[c(TRUE, FALSE),])
    colnames(rep_hybrid_gen_palm_dosage_chr)<-k
    rep_hybrid_gen_cerv_dosage_chr <- as.data.frame(rep_hybrid_gen[c(FALSE, TRUE),])
    colnames(rep_hybrid_gen_cerv_dosage_chr)<-k
    
    #process bc palms per species, per individual
    rep_hybrid_gen_palm_dosage_chr$pos <- pos_chr_gen$pos
    rep_hybrid_gen_cerv_dosage_chr$pos <- pos_chr_gen$pos
    rep_hybrid_gen_palm_dosage_chr$species <- "palm"
    rep_hybrid_gen_cerv_dosage_chr$species <- "cerv"
    rep_hybrid_gen_chr_dosage_all <- rbind(rep_hybrid_gen_palm_dosage_chr,rep_hybrid_gen_cerv_dosage_chr)
    rep_hybrid_gen_chr_dosage_all$u_id <- seq.int(nrow(rep_hybrid_gen_chr_dosage_all))
    
    rep_hybrid_gen_m_chr_dosage_all <- melt(rep_hybrid_gen_chr_dosage_all, id=c("pos","species", "u_id"), variable_name =  "ID")
    rep_hybrid_gen_um_chr <- unite(rep_hybrid_gen_m_chr_dosage_all, col = "id_spec", c("species","ID"), sep = "_", remove = F)
    rep_hybrid_gen_um_chr$method <- "genome"
    
    #read in the remarkable interval file for plotting
    RI <- read.table(paste0("E:/PSU/Hybridization Project/Genome_paper/Analysis/Circos_plot/Apalm_RI.BED"), header = T )
    RI_chr1<- RI %>% filter(scaffold %in% i)
    
    #set the average for bc palms
    sc_gd <- rep_hybrid_sc_um_chr %>% 
      group_by(species, pos) %>% 
      summarise(value = mean(value))
    
    gd<-bind_rows(rep_hybrid_sc_um_chr,rep_hybrid_gen_um_chr)
    
    #plot bc
    pdf(paste0(sc,"/genome_snpchip_sameID/",k,"_",j,"_",i,".pdf"), paper = "a4r", width = 8, height = 8)
    g1 <- ggplot(rep_hybrid_sc_um_chr, aes(x = (pos/1000000), y = value, color = species)) +
      geom_line(aes(group = id_spec), alpha =0.5, lwd=1.2, linetype="dashed") +
      geom_line(data = rep_hybrid_gen_um_chr, alpha = 0.5, lwd=0.9, linetype="solid") +
      theme_bw() +
      labs(
        title=k,
        x = paste0(i," Position (Mb)"),
        y = "Dosage",
        color = NULL
      ) + scale_colour_manual(name="species",values = c("cerv" = "#FBB116", "palm" = "#B33B77"))+
      scale_linetype_manual(name="linetype", values = c("genome" = "solid", "snpchip" = "dashed")) +
      scale_y_continuous(limits=c(0, 2))
    g2<-g1 +
      geom_rect(data=RI_chr1,inherit.aes=FALSE, mapping=aes(xmin=(start/1000000), xmax=(end/1000000), ymin=0,ymax=0.05), color="darkgrey")
    
    print(g2)
    dev.off()
    }
  }
}


