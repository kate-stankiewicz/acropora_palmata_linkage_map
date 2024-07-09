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

#SNPchip samples average 10 runs
sc<-c("E:/PSU/Hybridization Project/Genome_paper/Analysis/HybridCalling/ELAI/apalm_snpchip_avg10/avg_dosage_elai_runs/")
bc<-c("1839","1667","15727")
cross<-c("P199.C1821","P1017.C433","P199.C1822","P1017.C434","P199.C1823","P1017.C435","P1017.C436","P199.C1825","P199.C1826")
pop_sc<-c("Belize","Belize","Curacao","Belize","Curacao","Belize","Curacao","Curacao","Belize","Belize","Belize",
          "Belize","Belize","Cuba","Cuba","Bahamas","USVI","Belize","Belize","Antigua","Antigua","Antigua","Antigua",
          "Antigua","Antigua","Antigua","Antigua","USVI","USVI","USVI","USVI","USVI","USVI","Belize","Belize","Belize",
          "Belize","USVI","Florida","USVI","Bahamas","Belize","Curacao","Belize","Belize","Belize","USVI","USVI","Belize",
          "Belize","Belize","Belize","Belize","Belize","Belize","Belize","Belize","Belize","Belize","Bahamas")

mg<-c(1,2,3)

for(i in chrom){
  print(i)
  for(j in mg){
    print(j)
    hybr_sampIDs<-read.csv(paste0(sc,"samps_snp_",i,".txt"), sep=",", header = F)
    thy<-as.data.frame(t(hybr_sampIDs[1,]))
    thy$pop <-pop_sc
    colnames(thy)<-c("ID", "pop")
    thy$ID<-as.character(thy$ID)
    dosage_chr<-read.table(paste0(sc,"avg_hic_scaffold_",i,"_mg",j,".txt"), header = T,fill = TRUE,stringsAsFactors=FALSE, na.strings="NA")
    pos_chr <- read.table(paste0(sc,"rr_",j,"_",i,".snpinfo.txt"), header = T)
    thy$ID<-str_replace(thy$ID, "x", ".")
    t_dosage_chr <- as.data.frame(t(dosage_chr))
    colnames(t_dosage_chr) <- thy$ID
    
    #separate out the bc and hybrids
    bc_palms <- t_dosage_chr[, c("1839", "1667","15727")]
    hybrids_only <- t_dosage_chr[, -which(names(t_dosage_chr) %in% c("1839","1667","15727"))]
    
    #process the bc palms for plotting, get average dosage
    bc_palm_dosage_chr <- bc_palms[c(TRUE, FALSE),]
    bc_cerv_dosage_chr <- bc_palms[c(FALSE, TRUE),]
    bc_avg_palm_dos_chr <- as.data.frame(rowMeans(bc_palm_dosage_chr))
    bc_avg_cerv_dos_chr <- as.data.frame(rowMeans(bc_cerv_dosage_chr))
    colnames(bc_avg_palm_dos_chr) <- c("palm")
    colnames(bc_avg_cerv_dos_chr) <- c("cerv")
    bc_avg_palm_dos_chr$pos <- pos_chr$pos
    bc_avg_cerv_dos_chr$pos <- pos_chr$pos
    bc_all_chr_avg_dose<- merge(bc_avg_palm_dos_chr,bc_avg_cerv_dos_chr)
    bc_m_all_chr_avg_dose <- melt(bc_all_chr_avg_dose, id=c("pos")) 
    
    #process bc palms per species, per individual
    bc_palm_dosage_chr$pos <- pos_chr$pos
    bc_cerv_dosage_chr$pos <- pos_chr$pos
    bc_palm_dosage_chr$species <- "palm"
    bc_cerv_dosage_chr$species <- "cerv"
    bc_chr_dosage_all <- rbind(bc_palm_dosage_chr,bc_cerv_dosage_chr)
    bc_chr_dosage_all$u_id <- seq.int(nrow(bc_chr_dosage_all))
    
    bc_m_chr_dosage_all <- melt(bc_chr_dosage_all, id=c("pos","species", "u_id"), variable_name =  "ID")
    bc_um_chr <- unite(bc_m_chr_dosage_all, col = "id_spec", c("species","ID"), sep = "_", remove = F)
    
    #now repeat the process for hybrids
    #process the hybrids for plotting, get average dosage
    hy_palm_dosage_chr <- hybrids_only[c(TRUE, FALSE),]
    hy_cerv_dosage_chr <- hybrids_only[c(FALSE, TRUE),]
    hy_avg_palm_dos_chr <- as.data.frame(rowMeans(hy_palm_dosage_chr))
    hy_avg_cerv_dos_chr <- as.data.frame(rowMeans(hy_cerv_dosage_chr))
    colnames(hy_avg_palm_dos_chr) <- c("palm")
    colnames(hy_avg_cerv_dos_chr) <- c("cerv")
    hy_avg_palm_dos_chr$pos <- pos_chr$pos
    hy_avg_cerv_dos_chr$pos <- pos_chr$pos
    hy_all_chr_avg_dose<- merge(hy_avg_palm_dos_chr,hy_avg_cerv_dos_chr)
    hy_m_all_chr_avg_dose <- melt(hy_all_chr_avg_dose, id=c("pos")) 
    
    #process hybrids per species, per individual
    hy_palm_dosage_chr$pos <- pos_chr$pos
    hy_cerv_dosage_chr$pos <- pos_chr$pos
    hy_palm_dosage_chr$species <- "palm"
    hy_cerv_dosage_chr$species <- "cerv"
    hy_chr_dosage_all <- rbind(hy_palm_dosage_chr,hy_cerv_dosage_chr)
    hy_chr_dosage_all$u_id <- seq.int(nrow(hy_chr_dosage_all))
    
    hy_m_chr_dosage_all <- melt(hy_chr_dosage_all, id=c("pos","species", "u_id"), variable_name =  "ID")
    hy_um_chr <- unite(hy_m_chr_dosage_all, col = "id_spec", c("species","ID"), sep = "_", remove = F)
    
    hy_um_chr<-hy_um_chr %>%
      full_join(thy, by="ID")
    
    #read in the remarkable interval file for plotting
    RI <- read.table(paste0("E:/PSU/Hybridization Project/Genome_paper/Analysis/Circos_plot/Apalm_RI.BED"), header = T )
    RI_chr1<- RI %>% filter(scaffold %in% i)
    
    #set the average for bc palms
    bc_gd <- bc_um_chr %>% 
      group_by(species, pos) %>% 
      summarise(value = mean(value))
    
    #plot bc
    g1 <- ggplot(bc_um_chr, aes(x = (pos/1000000), y = value, color = species)) +
      geom_line(aes(group = id_spec,linetype=ID), alpha =0.5) +
      geom_line(data = bc_gd, alpha = 0.8, size = 2) +
      theme_bw() +
      labs(
        title = "Dosage Backcross A. palmata",
        x = paste0(i," Position (Mb)"),
        y = "Dosage",
        color = NULL
      ) + scale_colour_manual(values = c("cerv" = "#FBB116", "palm" = "#B33B77"))+
      scale_y_continuous(limits=c(0, 2))
    g2<-g1 +
      geom_rect(data=RI_chr1,inherit.aes=FALSE, mapping=aes(xmin=(start/1000000), xmax=(end/1000000), ymin=0,ymax=0.05), color="darkgrey")
    
    #set the average for F1-cross
    x_gd <- hy_um_chr %>% 
      filter(ID %in% cross) %>%
      group_by(species, pos) %>% 
      summarise(value = mean(value))
    
    #write out the plot
    g3 <- ggplot(hy_um_chr %>% filter(ID %in% cross), aes(x = (pos/1000000), y = value, color = species)) +
      geom_line(aes(group = id_spec,linetype=ID), alpha =0.5) +
      geom_line(data = x_gd, alpha = 0.8, size = 2) +
      theme_bw() +
      labs(
        title = "Dosage F1 outcross",
        x = paste0(i," Position (Mb)"),
        y = "Dosage",
        color = NULL
      ) + scale_colour_manual(values = c("cerv" = "#FBB116", "palm" = "#B33B77")) +
          scale_y_continuous(limits=c(0, 2))
    g4<-g3 +
      geom_rect(data=RI_chr1,inherit.aes=FALSE, mapping=aes(xmin=(start/1000000), xmax=(end/1000000), ymin=0,ymax=0.05), color="darkgrey")
    
    #set the average for F1-cross
    hy_gd <- hy_um_chr %>% 
      filter(!ID %in% cross) %>%
      group_by(species, pos) %>% 
      summarise(value = mean(value))
    
    #write out the plot
    g5 <- ggplot(hy_um_chr %>% filter(!ID %in% cross), aes(x = (pos/1000000), y = value, color = species)) +
      geom_line(aes(group = id_spec,linetype=pop), alpha =0.5) +
      geom_line(data = hy_gd, alpha = 0.8, size = 2) +
      theme_bw() +
      labs(
        title = "Dosage Wild F1 hybrids",
        x = paste0(i," Position (Mb)"),
        y = "Dosage",
        color = NULL
      ) + scale_colour_manual(values = c("cerv" = "#FBB116", "palm" = "#B33B77"))+
      scale_y_continuous(limits=c(0, 2))
    g6<-g5 +
      geom_rect(data=RI_chr1,inherit.aes=FALSE, mapping=aes(xmin=(start/1000000), xmax=(end/1000000), ymin=0,ymax=0.05), color="darkgrey")
  
    pdf(paste0(sc,"../combined_plots/all_",j,"_",i,".pdf"), paper = "a4r", width = 0, height = 0)
    print(ggarrange(g2+ rremove("x.title"),g4+ rremove("x.title"),g6,nrow = 3, align="v"))
    dev.off()
  }
}


