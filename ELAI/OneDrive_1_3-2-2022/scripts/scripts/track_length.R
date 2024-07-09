library(ggplot2)
library(reshape)
library(dplyr)
library(tidyr)
library(stringr)

#chrom sizes
chrS<-read.table("E:/PSU/Hybridization Project/Genome_paper/Analysis/Circos_plot/AP_chrom.txt",
                 check.names=FALSE, header=T, 
                 na.strings=c("", "NA"), stringsAsFactors=FALSE, sep="\t",quote="")
chrom<-chrS$chr

bc<-c("1839","1667")

mg<-c(1,2,3)

#genome samples
setwd("E:/PSU/Hybridization Project/Genome_paper/Analysis/HybridCalling/ELAI/apalm_genome_plots/apalm_genome/")

## Read in all data frames using a loop
datalist_palm_cnt = list()
datalist_palm_len = list()
datalist_cerv_cnt = list()
datalist_cerv_len = list()

for(i in chrom){
  print(i)
  for(j in mg){
    print(j)
    hybr_sampIDs<-read.csv(paste0("samps_",i,".txt"), sep=",", header = F)
    thy<-as.data.frame(t(hybr_sampIDs[1,]))
    colnames(thy)<-c("ID")
    thy$ID<-as.character(thy$ID)
    dosage_chr<-read.table(paste0("rr__",j,"_",i,".ps21.txt"), header = F)
    pos_chr <- read.table(paste0("rr__",j,"_",i,".snpinfo.txt"), header = T)
    t_dosage_chr <- as.data.frame(t(dosage_chr))
    colnames(t_dosage_chr) <- hybr_sampIDs[1,]

    #separate out the bc and hybrids
    bc_palms <- t_dosage_chr[, c("1839", "1667")]
    hybrids_only <- t_dosage_chr[, -which(names(t_dosage_chr) %in% c("1839","1667"))]
  
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
  
    #count/length of Ap SNPs
    bc_palm_len <- bc_um_chr %>% 
      group_by(ID) %>% 
      filter(species=="palm") %>% 
      mutate(grp=ifelse(value > 1.5, cumsum(c(value[1]>1.5, diff(value>1.5)) == 1), 0)) %>% 
      group_by(ID, grp) %>%
      summarize(len=max(pos)- min(pos), min=min(pos),max=max(pos)) %>%
      filter(!grp == 0) %>%
      mutate(hyb="Bc", method="genome",cnt=length(len), chr=i, mg=j)
  
    hy_palm_len <- hy_um_chr %>% 
      group_by(ID) %>% 
      filter(species=="palm") %>% 
      mutate(grp=ifelse(value > 1.5, cumsum(c(value[1]>1.5, diff(value>1.5)) == 1), 0)) %>% 
      group_by(ID, grp) %>%
      summarize(len=max(pos)- min(pos), min=min(pos),max=max(pos)) %>%
      filter(!grp == 0) %>%
      mutate( hyb="F1",method="genome",cnt=length(len),chr=i, mg=j)
 
    palm_ct<-bc_palm_len %>% 
      bind_rows(hy_palm_len) %>%
      full_join(thy, by="ID") %>%
      mutate(cnt=ifelse(is.na(cnt), 0, cnt),chr=ifelse(is.na(chr), i, chr), 
           mg=ifelse(is.na(mg), j, mg),
           method=ifelse(is.na(method),"genome", method),
           hyb=ifelse(is.na(hyb), ifelse(ID %in% bc, "Bc", "F1"), hyb)) %>%
      select(ID, cnt, hyb, method,mg, chr)
  
    palm_len<-bc_palm_len %>% 
      bind_rows(hy_palm_len) 
    
    name <- paste(i,"-mg",j, sep='')

    datalist_palm_cnt[[name]] <- palm_ct
    
    datalist_palm_len[[name]] <- palm_len
    
    #count/length of Ap SNPs
    bc_cerv_len <- bc_um_chr %>% 
      group_by(ID) %>% 
      filter(species=="cerv") %>% 
      mutate(grp=ifelse(value > 1.5, cumsum(c(value[1]>1.5, diff(value>1.5)) == 1), 0)) %>% 
      group_by(ID, grp) %>%
      summarize(len=max(pos)- min(pos), min=min(pos),max=max(pos)) %>%
      filter(!grp == 0) %>%
      mutate(hyb="Bc", method="genome",cnt=length(len), chr=i, mg=j)
    
    hy_cerv_len <- hy_um_chr %>% 
      group_by(ID) %>% 
      filter(species=="cerv") %>% 
      mutate(grp=ifelse(value > 1.5, cumsum(c(value[1]>1.5, diff(value>1.5)) == 1), 0)) %>% 
      group_by(ID, grp) %>%
      summarize(len=max(pos)- min(pos), min=min(pos),max=max(pos)) %>%
      filter(!grp == 0) %>%
      mutate( hyb="F1",method="genome",cnt=length(len),chr=i, mg=j)
    
    cerv_ct<-bc_cerv_len %>% 
      bind_rows(hy_cerv_len) %>%
      full_join(thy, by="ID") %>%
      mutate(cnt=ifelse(is.na(cnt), 0, cnt),chr=ifelse(is.na(chr), i, chr), 
             mg=ifelse(is.na(mg), j, mg),
             method=ifelse(is.na(method),"genome", method),
             hyb=ifelse(is.na(hyb), ifelse(ID %in% bc, "Bc", "F1"), hyb)) %>%
      select(ID, cnt, hyb, method,mg, chr)
    
    cerv_len<-bc_cerv_len %>% 
      bind_rows(hy_cerv_len) 
    
    datalist_cerv_cnt[[name]] <- cerv_ct
    
    datalist_cerv_len[[name]] <- cerv_len
  }
}



big_data_palm_cnt <- dplyr::bind_rows(datalist_palm_cnt) %>% distinct() %>% mutate(species="palm")
big_data_palm_len <- dplyr::bind_rows(datalist_palm_len) %>% distinct()%>% mutate(species="palm")

write.table(big_data_palm_cnt, "palmTracks_cnt_genome.txt", sep="\t")
write.table(big_data_palm_len, "palmTracks_len_genome.txt", sep="\t")

big_data_cerv_cnt <- dplyr::bind_rows(datalist_cerv_cnt) %>% distinct()%>% mutate(species="cerv")
big_data_cerv_len <- dplyr::bind_rows(datalist_cerv_len) %>% distinct()%>% mutate(species="cerv")

write.table(big_data_cerv_cnt, "cervTracks_cnt_genome.txt", sep="\t")
write.table(big_data_cerv_len, "cervTracks_len_genome.txt", sep="\t")

#SNPchip samples
setwd("E:/PSU/Hybridization Project/Genome_paper/Analysis/HybridCalling/ELAI/apalm_snpchip_plots/apalm_snpchip/")
bc<-c("1839","1667","15727")
cross<-c("P199.C1821","P1017.C433","P199.C1822","P1017.C434","P199.C1823","P1017.C435","P1017.C436","P199.C1825","P199.C1826")

## Read in all data frames using a loop
datalist_palm_sc_cnt = list()
datalist_palm_sc_len = list()
datalist_cerv_sc_cnt = list()
datalist_cerv_sc_len = list()

for(i in chrom){
  print(i)
  for(j in mg){
    print(j)
    hybr_sampIDs<-read.csv(paste0("samps_snp_",i,".txt"), sep=",", header = F)
    thy<-as.data.frame(t(hybr_sampIDs[1,]))
    colnames(thy)<-c("ID")
    thy$ID<-as.character(thy$ID)
    dosage_chr<-read.table(paste0("rr_",j,"_",i,".ps21.txt"), header = F)
    pos_chr <- read.table(paste0("rr_",j,"_",i,".snpinfo.txt"), header = T)
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
    
    #count/length of Ap SNPs
    bc_palm_len <- bc_um_chr %>% 
      group_by(ID) %>% 
      filter(species=="palm") %>% 
      mutate(grp=ifelse(value > 1.5, cumsum(c(value[1]>1.5, diff(value>1.5)) == 1), 0)) %>% 
      group_by(ID, grp) %>%
      summarize(len=max(pos)- min(pos), min=min(pos),max=max(pos)) %>%
      filter(!grp == 0) %>%
      mutate(hyb="Bc", method="snpchip",cnt=length(len), chr=i, mg=j)
    
    hy_palm_len <- hy_um_chr %>% 
      group_by(ID) %>% 
      filter(species=="palm") %>% 
      mutate(grp=ifelse(value > 1.5, cumsum(c(value[1]>1.5, diff(value>1.5)) == 1), 0)) %>% 
      group_by(ID, grp) %>%
      summarize(len=max(pos)- min(pos), min=min(pos),max=max(pos)) %>%
      filter(!grp == 0) %>%
      mutate(hyb=ifelse(ID %in% cross, "F1-x", "F1"),method="snpchip",cnt=length(len),chr=i, mg=j)
    
    palm_ct<-bc_palm_len %>% 
      bind_rows(hy_palm_len) %>%
      full_join(thy, by="ID") %>%
      mutate(cnt=ifelse(is.na(cnt), 0, cnt),chr=ifelse(is.na(chr), i, chr), 
             mg=ifelse(is.na(mg), j, mg),
             method=ifelse(is.na(method),"snpchip", method),
             hyb=ifelse(is.na(hyb), ifelse(ID %in% bc, "Bc", ifelse(ID %in% cross, "F1-x","F1")), hyb)) %>%
      select(ID, cnt, hyb, method,mg, chr)
    
    palm_len<-bc_palm_len %>% 
      bind_rows(hy_palm_len) 
    
    name <- paste(i,"-mg",j, sep='')
    
    datalist_palm_sc_cnt[[name]] <- palm_ct
    
    datalist_palm_sc_len[[name]] <- palm_len
    
    #count/length of Ap SNPs
    bc_cerv_len <- bc_um_chr %>% 
      group_by(ID) %>% 
      filter(species=="cerv") %>% 
      mutate(grp=ifelse(value > 1.5, cumsum(c(value[1]>1.5, diff(value>1.5)) == 1), 0)) %>% 
      group_by(ID, grp) %>%
      summarize(len=max(pos)- min(pos), min=min(pos),max=max(pos)) %>%
      filter(!grp == 0) %>%
      mutate(hyb="Bc", method="snpchip",cnt=length(len), chr=i, mg=j)
    
    hy_cerv_len <- hy_um_chr %>% 
      group_by(ID) %>% 
      filter(species=="cerv") %>% 
      mutate(grp=ifelse(value > 1.5, cumsum(c(value[1]>1.5, diff(value>1.5)) == 1), 0)) %>% 
      group_by(ID, grp) %>%
      summarize(len=max(pos)- min(pos), min=min(pos),max=max(pos)) %>%
      filter(!grp == 0) %>%
      mutate(hyb=ifelse(ID %in% cross, "F1-x", "F1"),method="snpchip",cnt=length(len),chr=i, mg=j)
    
    cerv_ct<-bc_cerv_len %>% 
      bind_rows(hy_cerv_len) %>%
      full_join(thy, by="ID") %>%
      mutate(cnt=ifelse(is.na(cnt), 0, cnt),chr=ifelse(is.na(chr), i, chr), 
             mg=ifelse(is.na(mg), j, mg),
             method=ifelse(is.na(method),"snpchip", method),
             hyb=ifelse(is.na(hyb), ifelse(ID %in% bc, "Bc", ifelse(ID %in% cross, "F1-x","F1")), hyb)) %>%
      select(ID, cnt, hyb, method,mg, chr)
    
    cerv_len<-bc_cerv_len %>% 
      bind_rows(hy_cerv_len) 
    
    datalist_cerv_sc_cnt[[name]] <- cerv_ct
    
    datalist_cerv_sc_len[[name]] <- cerv_len
  }
}


big_data_sc_palm_cnt <- dplyr::bind_rows(datalist_palm_sc_cnt) %>% distinct() %>% mutate(species="palm")
big_data_sc_palm_len <- dplyr::bind_rows(datalist_palm_sc_len) %>% distinct()%>% mutate(species="palm")

write.table(big_data_sc_palm_cnt, "palmTracks_cnt_snpchip.txt", sep="\t")
write.table(big_data_sc_palm_len, "palmTracks_len_snpchip.txt", sep="\t")

big_data_sc_cerv_cnt <- dplyr::bind_rows(datalist_cerv_sc_cnt) %>% distinct()%>% mutate(species="cerv")
big_data_sc_cerv_len <- dplyr::bind_rows(datalist_cerv_sc_len) %>% distinct()%>% mutate(species="cerv")

write.table(big_data_sc_cerv_cnt, "cervTracks_cnt_snpchip.txt", sep="\t")
write.table(big_data_sc_cerv_len, "cervTracks_len_snpchip.txt", sep="\t")


#combine genome and snpchip data
big_data_cnt<- big_data_palm_cnt %>% bind_rows(big_data_cerv_cnt,big_data_sc_palm_cnt,big_data_sc_cerv_cnt)
big_data_len<- big_data_palm_len %>% bind_rows(big_data_cerv_len,big_data_sc_palm_len,big_data_sc_cerv_len)


data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

#now plot
# count
g5 <-ggplot(big_data_cnt %>% filter(mg==1), aes(x=species, y=cnt, color=factor(hyb), fill=factor(hyb)), group=factor(method)) +
  geom_point(aes(shape=method, color=hyb),position=position_jitterdodge(jitter.height=0.1,dodge.width=0.8),alpha=0.1)+ 
  theme_bw()+
  theme(panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
        strip.background = element_rect(colour="grey30", fill="grey95",
                                        size=1.1, linetype="solid"),
        strip.text.x = element_text(size=12, color="black",
                                    face="bold")) +
  labs(y="Count of Parental Ancestry Tracks",
       x="")+
  scale_color_manual(values = c("F1" = "#6F81BC","F1-x" = "#485b99", "Bc" = "#6D2B75"))+
  scale_fill_manual(values = c("F1" = "#ACB6D8", "F1-x" = "#485b99","Bc" = "#6D2B75"))+
  scale_shape_manual(values=c(21,24))+
  facet_wrap(~ method)
g6<-g5 + stat_summary(fun.data=data_summary, geom="pointrange",shape=c(18),lwd=1,position=position_dodge(width=0.8))
print(g6)


# length
g7 <- ggplot(big_data_len %>% filter(mg==1), aes(x=species, y=(len/1000000), color=factor(hyb), fill=factor(hyb)), group=factor(method)) +
  geom_point(aes(shape=method, color=hyb),position=position_jitterdodge(jitter.height=0.2,dodge.width=0.8),alpha=0.1) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(colour="grey30", fill="grey95", 
                                        size=1.1, linetype="solid"),
        strip.text.x = element_text(size=12, color="black",
                                    face="bold")) +
  labs(y="Length of Parental Ancestry Tracks (Mb)",
       x="") +
  scale_color_manual(values = c("F1" = "#6F81BC","F1-x" = "#485b99", "Bc" = "#6D2B75"))+
  scale_fill_manual(values = c("F1" = "#ACB6D8","F1-x" = "#485b99", "Bc" = "#6D2B75"))+
  facet_wrap(~ method)
g8<-g7 + stat_summary(fun.data=data_summary, geom="pointrange",shape=18, lwd=1,position=position_dodge(width=0.8))
print(g8)

library(ggpubr)

ggarrange(g6, g8, labels = c("A", "B"),nrow = 2)

