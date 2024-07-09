library(ggplot2)
library(reshape)
library(dplyr)
library(tidyr)
library(stringr)
library(GenomicRanges)

#chrom sizes
chrS<-read.table("E:/PSU/Hybridization Project/Genome_paper/Analysis/Circos_plot/AP_chrom.txt",
                 check.names=FALSE, header=T, 
                 na.strings=c("", "NA"), stringsAsFactors=FALSE, sep="\t",quote="")
chrom<-chrS$chr

bc<-c("1839","1667")

mg<-c(1,2,3)
pop_gen<-c("USVI","USVI","USVI","USVI","USVI","USVI","USVI","USVI","USVI","USVI","Belize","Belize","Belize",
           "Belize","Belize","Belize","Belize","Belize","Belize","Curacao","Florida","Bahamas")

#genome samples
#average 10 runs
gen<-c("E:/PSU/Hybridization Project/Genome_paper/Analysis/HybridCalling/ELAI/apalm_genome_avg10/avg_dosage_elai_runs/")

## Read in all data frames using a loop
datalist_palm_cnt = list()
datalist_palm_len = list()
datalist_cerv_cnt = list()
datalist_cerv_len = list()
datalist_hy_gen_avg = list()
datalist_bc_gen_avg = list()

for(i in chrom){
  print(i)
  for(j in mg){
    print(j)
    hybr_sampIDs<-read.csv(paste0(gen,"samps_",i,".txt"), sep=",", header = F)
    thy<-as.data.frame(t(hybr_sampIDs[1,]))
    thy$pop <-pop_gen
    colnames(thy)<-c("ID", "pop")
    thy$ID<-as.character(thy$ID)
    dosage_chr<-read.table(paste0(gen,"avg_hic_scaffold_",i,"_mg",j,".txt"), header = T,fill = TRUE,stringsAsFactors=FALSE, na.strings="NA")
    pos_chr <- read.table(paste0(gen,"rr__",j,"_",i,".snpinfo.txt"), header = T)
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
    
    #bc_avg_dose <- bc_m_all_chr_avg_dose %>%
    #  mutate(chr=i, mg=j) %>%
    #  group_by(variable) %>% 
    # mutate(grp=ifelse(value > 1.5, cumsum(c(value[1]>1.5, diff(value>1.5)) == 1), 0)) %>% 
    #  group_by(grp, variable) %>%
    #  summarize(len=max(pos)- min(pos), min=min(pos),max=max(pos)) %>%
    #  filter(!grp == 0) %>%
    #  mutate(method="genome",cnt=length(len), chr=i, mg=j)
  
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
    
    #hy_avg_dose <- hy_m_all_chr_avg_dose %>%
    #  mutate(chr=i, mg=j) %>%
    # group_by(variable) %>% 
    #  mutate(grp=ifelse(value > 1.5, cumsum(c(value[1]>1.5, diff(value>1.5)) == 1), 0)) %>% 
    #  group_by(grp, variable) %>%
    # summarize(len=max(pos)- min(pos), min=min(pos),max=max(pos)) %>%
    #  filter(!grp == 0) %>%
    #  mutate(method="genome",cnt=length(len), chr=i, mg=j)
  
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
      select(ID, pop,cnt, hyb, method,mg, chr)
  
    palm_len<-bc_palm_len %>% 
      bind_rows(hy_palm_len)  %>%
      full_join(thy, by="ID") %>%
      mutate(len=ifelse(is.na(len), 0, len),chr=ifelse(is.na(chr), i, chr), 
             mg=ifelse(is.na(mg), j, mg),
             method=ifelse(is.na(method),"genome", method),
             hyb=ifelse(is.na(hyb), ifelse(ID %in% bc, "Bc", "F1"), hyb)) %>%
      select(ID, pop,len, hyb, method,mg, chr) 
    
    hy_avg_dose <- hy_um_chr %>%
      group_by(species,ID) %>%
      mutate(i = value > 1.5 & lead(value,default = FALSE) > 1.5) %>%
      filter(i =="TRUE") %>% 
      ungroup() %>%
      group_by(pos,species)%>%
      mutate(pos_count = n()) %>%
      filter(!pos_count == 1) %>%
      ungroup() %>%
      group_by(species,pos_count,ID) %>%
      summarize(len=max(pos)- min(pos), Min=min(pos),Max=max(pos), samples=str_c(unique(ID), collapse = ", ")) %>%
      mutate(method="genome",chr=i, mg=j, rng=paste0(Min,"-",Max))%>%
      ungroup()
    
    hy_avg_dose_filt <-as.data.frame(disjoin(GRanges(hy_avg_dose$species, IRanges(hy_avg_dose$Min, hy_avg_dose$Max))))
    colnames(hy_avg_dose_filt)<-c("species","Min","Max","length","strand")
    
    hy_avg_dose_filt2 <- hy_avg_dose_filt %>%
      mutate(method="genome",chr=i, mg=j,rng=paste0(Min,"-",Max)) %>%
      full_join(hy_avg_dose %>%
                  select(rng, ID,Min, Max), by= c("rng"))%>%
      filter(!is.na(ID)) %>%
      group_by(Min.y)%>%
      mutate(samples=str_c(unique(ID), collapse = ", ")) %>%
      ungroup ()%>%
      group_by(Max.y)%>%
      mutate(samples2=str_c(unique(ID), collapse = ", ")) %>%
      unite("samps",samples:samples2, sep=", ") %>%
      mutate(samps = as.character(sapply(samps, function(x) paste(unique(unlist(str_split(x,", "))), collapse = ", ")))) %>%
      drop_na() %>%
      ungroup() %>%
      select(-Min.y,-Max.y, -ID, -strand) %>%
      distinct()

    bc_avg_dose <- bc_um_chr %>%
      group_by(species,ID) %>%
      mutate(i = value > 1.5 & lead(value,default = FALSE) > 1.5) %>%
      filter(i =="TRUE") %>% 
      ungroup() %>%
      group_by(pos,species)%>%
      mutate(pos_count = n()) %>%
      filter(!pos_count == 1) %>%
      ungroup() %>%
      group_by(species,pos_count,ID) %>%
      summarize(len=max(pos)- min(pos), Min=min(pos),Max=max(pos), samples=str_c(unique(ID), collapse = ", ")) %>%
      mutate(method="genome",chr=i, mg=j, rng=paste0(Min,"-",Max))%>%
      ungroup()
    
    bc_avg_dose_filt <-as.data.frame(disjoin(GRanges(bc_avg_dose$species, IRanges(bc_avg_dose$Min, bc_avg_dose$Max))))
    colnames(bc_avg_dose_filt)<-c("species","Min","Max","length","strand")
    
    bc_avg_dose_filt2 <- bc_avg_dose_filt %>%
      mutate(method="genome",chr=i, mg=j,rng=paste0(Min,"-",Max)) %>%
      full_join(bc_avg_dose %>%
                  select(rng, ID,Min, Max), by= c("rng"))%>%
      filter(!is.na(ID)) %>%
      group_by(Min.y)%>%
      mutate(samples=str_c(unique(ID), collapse = ", ")) %>%
      ungroup ()%>%
      group_by(Max.y)%>%
      mutate(samples2=str_c(unique(ID), collapse = ", ")) %>%
      unite("samps",samples:samples2, sep=", ") %>%
      mutate(samps = as.character(sapply(samps, function(x) paste(unique(unlist(str_split(x,", "))), collapse = ", ")))) %>%
      drop_na() %>%
      ungroup() %>%
      select(-Min.y,-Max.y, -ID, -strand) %>%
      distinct()

    name <- paste(i,"-mg",j, sep='')
    
    datalist_hy_gen_avg[[name]] <-hy_avg_dose_filt2
    
    datalist_bc_gen_avg[[name]] <-bc_avg_dose_filt2

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
      select(ID, pop,cnt, hyb, method,mg, chr)
  
    cerv_len<-bc_cerv_len %>% 
      bind_rows(hy_cerv_len)  %>%
      full_join(thy, by="ID") %>%
      mutate(len=ifelse(is.na(len), 0, len),chr=ifelse(is.na(chr), i, chr), 
             mg=ifelse(is.na(mg), j, mg),
             method=ifelse(is.na(method),"genome", method),
             hyb=ifelse(is.na(hyb), ifelse(ID %in% bc, "Bc", "F1"), hyb)) %>%
      select(ID, pop,len, hyb, method,mg, chr) 
    
    datalist_cerv_cnt[[name]] <- cerv_ct
    
    datalist_cerv_len[[name]] <- cerv_len
  }
}

big_data_palm_cnt <- dplyr::bind_rows(datalist_palm_cnt) %>% distinct() %>% mutate(species="palm")
#big_data_palm_cnt<-read.table("palmTracks_cnt_genome.txt",sep="\t")
big_data_palm_cnt$ID<-as.character(big_data_palm_cnt$ID)

big_data_palm_len <- dplyr::bind_rows(datalist_palm_len) %>% distinct()%>% mutate(species="palm")
#big_data_palm_len <-read.table("palmTracks_len_genome.txt",sep="\t")
big_data_palm_len$ID<-as.character(big_data_palm_len$ID)

palm_len<-big_data_palm_len %>% 
  full_join(thy, by="ID") %>%
  mutate(len=ifelse(is.na(len), 0, len)) %>%
  select(ID, pop,len, hyb, method,mg, chr, species) 

write.table(big_data_palm_cnt, "E:/PSU/Hybridization Project/Genome_paper/Analysis/HybridCalling/ELAI/apalm_genome_avg10/avg_dosage_elai_runs/palmTracks_cnt_genome.txt", sep="\t",row.names=FALSE)
write.table(big_data_palm_len, "E:/PSU/Hybridization Project/Genome_paper/Analysis/HybridCalling/ELAI/apalm_genome_avg10/avg_dosage_elai_runs/palmTracks_len_genome.txt", sep="\t",row.names=FALSE)

big_data_cerv_cnt <- dplyr::bind_rows(datalist_cerv_cnt) %>% distinct()%>% mutate(species="cerv")
#big_data_cerv_cnt<-read.table("cervTracks_cnt_genome.txt",sep="\t")
big_data_cerv_cnt$ID<-as.character(big_data_cerv_cnt$ID)

big_data_cerv_len <- dplyr::bind_rows(datalist_cerv_len) %>% distinct()%>% mutate(species="cerv")
#big_data_cerv_len <-read.table("cervTracks_len_genome.txt",sep="\t")
big_data_cerv_len$ID<-as.character(big_data_cerv_len$ID)

cerv_len<-big_data_cerv_len %>% 
  full_join(thy, by="ID") %>%
  mutate(len=ifelse(is.na(len), 0, len)) %>%
  select(ID, pop,len, hyb, method,mg, chr,species) 

write.table(big_data_cerv_cnt, "E:/PSU/Hybridization Project/Genome_paper/Analysis/HybridCalling/ELAI/apalm_genome_avg10/avg_dosage_elai_runs/cervTracks_cnt_genome.txt", sep="\t",row.names=FALSE)
write.table(big_data_cerv_len, "E:/PSU/Hybridization Project/Genome_paper/Analysis/HybridCalling/ELAI/apalm_genome_avg10/avg_dosage_elai_runs/cervTracks_len_genome.txt", sep="\t",row.names=FALSE)

big_data_hy_gen_avg <- dplyr::bind_rows(datalist_hy_gen_avg) %>% distinct() 
write.table(big_data_hy_gen_avg, "E:/PSU/Hybridization Project/Genome_paper/Analysis/HybridCalling/ELAI/apalm_genome_avg10/avg_dosage_elai_runs/shared_F1_Tracks_genome.txt", sep="\t",row.names=FALSE)

big_data_bc_gen_avg <- dplyr::bind_rows(datalist_bc_gen_avg) %>% distinct() 
write.table(big_data_bc_gen_avg, "E:/PSU/Hybridization Project/Genome_paper/Analysis/HybridCalling/ELAI/apalm_genome_avg10/avg_dosage_elai_runs/shared_BC_Tracks_genome.txt", sep="\t",row.names=FALSE)

#SNPchip samples
#average of 10 runs
sc<-c("E:/PSU/Hybridization Project/Genome_paper/Analysis/HybridCalling/ELAI/apalm_snpchip_avg10/avg_dosage_elai_runs/")

bc<-c("1839","1667","15727")
cross<-c("P199.C1821","P1017.C433","P199.C1822","P1017.C434","P199.C1823","P1017.C435","P1017.C436","P199.C1825","P199.C1826")
pop_sc<-c("Belize","Belize","Curacao","Belize","Curacao","Belize","Curacao","Curacao","Belize","Belize","Belize",
       "Belize","Belize","Cuba","Cuba","Bahamas","USVI","Belize","Belize","Antigua","Antigua","Antigua","Antigua",
       "Antigua","Antigua","Antigua","Antigua","USVI","USVI","USVI","USVI","USVI","USVI","Belize","Belize","Belize",
       "Belize","USVI","Florida","USVI","Bahamas","Belize","Curacao","Belize","Belize","Belize","USVI","USVI","Belize",
       "Belize","Belize","Belize","Belize","Belize","Belize","Belize","Belize","Belize","Belize","Bahamas")
#pop_sc<-as.data.frame(pop_sc)

## Read in all data frames using a loop
datalist_palm_sc_cnt = list()
datalist_palm_sc_len = list()
datalist_cerv_sc_cnt = list()
datalist_cerv_sc_len = list()
datalist_hy_sc_avg = list()
datalist_bc_sc_avg = list()

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
    bc_m_all_chr_avg_dose <- reshape2::melt(bc_all_chr_avg_dose, id=c("pos")) 
    
    #bc_avg_dose <- bc_m_all_chr_avg_dose %>%
    #  mutate(chr=i, mg=j) %>%
    #  group_by(variable) %>% 
    #  mutate(grp=ifelse(value > 1.5, cumsum(c(value[1]>1.5, diff(value>1.5)) == 1), 0)) %>% 
    #  group_by(grp, variable) %>%
    #  summarize(len=max(pos)- min(pos), min=min(pos),max=max(pos)) %>%
    #  filter(!grp == 0) %>%
    #  mutate(method="snpchip",cnt=length(len), chr=i, mg=j)
    
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
    hy_m_all_chr_avg_dose <- reshape2::melt(hy_all_chr_avg_dose, id=c("pos")) 
    
    #hy_avg_dose <- hy_m_all_chr_avg_dose %>%
    #  mutate(chr=i, mg=j) %>%
    #  group_by(variable) %>% 
    #  mutate(grp=ifelse(value > 1.5, cumsum(c(value[1]>1.5, diff(value>1.5)) == 1), 0)) %>% 
    #  group_by(grp,variable) %>%
    # summarize(len=max(pos)- min(pos), min=min(pos),max=max(pos)) %>%
    #  filter(!grp == 0) %>%
    # mutate(method="snpchip",cnt=length(len), chr=i, mg=j)
    
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
      mutate(grp=ifelse(value > 1.5, cumsum(c(value[1]>=1.5, diff(value>=1.5)) == 1), 0)) %>% 
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
      mutate(hyb=ifelse(ID %in% cross, "F1-x", "F1"),method="snpchip",cnt=length(len),chr=i, mg=j)%>%
      mutate_if(is.logical, as.character)
    
    palm_ct<-bc_palm_len %>% 
      bind_rows(hy_palm_len) %>%
      full_join(thy, by="ID") %>%
      mutate(cnt=ifelse(is.na(cnt), 0, cnt),chr=ifelse(is.na(chr), i, chr), 
             mg=ifelse(is.na(mg), j, mg),
             method=ifelse(is.na(method),"snpchip", method),
             hyb=ifelse(is.na(hyb), ifelse(ID %in% bc, "Bc", ifelse(ID %in% cross, "F1-x","F1")), hyb)) %>%
      select(ID, pop,cnt, hyb, method,mg, chr)
    
    palm_len <-bc_palm_len %>% 
      bind_rows(hy_palm_len)  %>%
      full_join(thy, by="ID") %>%
      mutate(len=ifelse(is.na(len), 0, len),chr=ifelse(is.na(chr), i, chr), 
           mg=ifelse(is.na(mg), j, mg),
           method=ifelse(is.na(method),"snpchip", method),
           hyb=ifelse(is.na(hyb), ifelse(ID %in% bc, "Bc", ifelse(ID %in% cross, "F1-x","F1")), hyb)) %>%
      select(ID, pop,len, hyb, method,mg, chr)
    
    hy_avg_dose <- hy_um_chr %>%
      group_by(species,ID) %>%
      mutate(i = value > 1.5 & lead(value,default = FALSE) > 1.5) %>%
      filter(i =="TRUE") %>% 
      ungroup() %>%
      group_by(pos,species)%>%
      mutate(pos_count = n()) %>%
      filter(!pos_count == 1) %>%
      ungroup() %>%
      group_by(species,pos_count,ID) %>%
      summarize(len=max(pos)- min(pos), Min=min(pos),Max=max(pos), samples=str_c(unique(ID), collapse = ", ")) %>%
      mutate(method="snpchip",chr=i, mg=j, rng=paste0(Min,"-",Max))%>%
      ungroup()
    
    hy_avg_dose_filt <-as.data.frame(disjoin(GRanges(hy_avg_dose$species, IRanges(hy_avg_dose$Min, hy_avg_dose$Max))))
    colnames(hy_avg_dose_filt)<-c("species","Min","Max","length","strand")
    
    hy_avg_dose_filt2 <- hy_avg_dose_filt %>%
      mutate(method="snpchip",chr=i, mg=j,rng=paste0(Min,"-",Max)) %>%
      full_join(hy_avg_dose %>%
                  select(rng, ID,Min, Max), by= c("rng"))%>%
      filter(!is.na(ID)) %>%
      group_by(Min.y)%>%
      mutate(samples=str_c(unique(ID), collapse = ", ")) %>%
      ungroup ()%>%
      group_by(Max.y)%>%
      mutate(samples2=str_c(unique(ID), collapse = ", ")) %>%
      unite("samps",samples:samples2, sep=", ") %>%
      mutate(samps = as.character(sapply(samps, function(x) paste(unique(unlist(str_split(x,", "))), collapse = ", ")))) %>%
      drop_na() %>%
      ungroup() %>%
      select(-Min.y,-Max.y, -ID, -strand) %>%
      distinct()
    
    bc_avg_dose <- bc_um_chr %>%
      group_by(species,ID) %>%
      mutate(i = value > 1.5 & lead(value,default = FALSE) > 1.5) %>%
      filter(i =="TRUE") %>% 
      ungroup() %>%
      group_by(pos,species)%>%
      mutate(pos_count = n()) %>%
      filter(!pos_count == 1) %>%
      ungroup() %>%
      group_by(species,pos_count,ID) %>%
      summarize(len=max(pos)- min(pos), Min=min(pos),Max=max(pos), samples=str_c(unique(ID), collapse = ", ")) %>%
      mutate(method="snpchip",chr=i, mg=j, rng=paste0(Min,"-",Max))%>%
      ungroup()
    
    bc_avg_dose_filt <-as.data.frame(disjoin(GRanges(bc_avg_dose$species, IRanges(bc_avg_dose$Min, bc_avg_dose$Max))))
    colnames(bc_avg_dose_filt)<-c("species","Min","Max","length","strand")
    
    bc_avg_dose_filt2 <- bc_avg_dose_filt %>%
      mutate(method="snpchip",chr=i, mg=j,rng=paste0(Min,"-",Max)) %>%
      full_join(bc_avg_dose %>%
                  select(rng, ID,Min, Max), by= c("rng"))%>%
      filter(!is.na(ID)) %>%
      group_by(Min.y)%>%
      mutate(samples=str_c(unique(ID), collapse = ", ")) %>%
      ungroup ()%>%
      group_by(Max.y)%>%
      mutate(samples2=str_c(unique(ID), collapse = ", ")) %>%
      unite("samps",samples:samples2, sep=", ") %>%
      mutate(samps = as.character(sapply(samps, function(x) paste(unique(unlist(str_split(x,", "))), collapse = ", ")))) %>%
      drop_na() %>%
      ungroup() %>%
      select(-Min.y,-Max.y, -ID, -strand) %>%
      distinct()
    
    name <- paste(i,"-mg",j, sep='')
    
    datalist_hy_sc_avg[[name]] <-hy_avg_dose_filt2
    
    datalist_bc_sc_avg[[name]] <-bc_avg_dose_filt2
    
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
      mutate(hyb=ifelse(ID %in% cross, "F1-x", "F1"),method="snpchip",cnt=length(len),chr=i, mg=j) %>%
      mutate_if(is.logical, as.character)
    
    cerv_ct<-bc_cerv_len %>% 
      bind_rows(hy_cerv_len) %>%
      full_join(thy, by="ID") %>%
      mutate(cnt=ifelse(is.na(cnt), 0, cnt),chr=ifelse(is.na(chr), i, chr), 
             mg=ifelse(is.na(mg), j, mg),
             method=ifelse(is.na(method),"snpchip", method),
             hyb=ifelse(is.na(hyb), ifelse(ID %in% bc, "Bc", ifelse(ID %in% cross, "F1-x","F1")), hyb)) %>%
      select(ID, pop, cnt, hyb, method,mg, chr)
    
    cerv_len<-bc_cerv_len %>% 
      bind_rows(hy_cerv_len) %>% 
      full_join(thy, by="ID") %>%
      mutate(len=ifelse(is.na(len), 0, len),chr=ifelse(is.na(chr), i, chr), 
             mg=ifelse(is.na(mg), j, mg),
             method=ifelse(is.na(method),"snpchip", method),
             hyb=ifelse(is.na(hyb), ifelse(ID %in% bc, "Bc", ifelse(ID %in% cross, "F1-x","F1")), hyb)) %>%
      select(ID, pop,len, hyb, method,mg, chr)
    
    datalist_cerv_sc_cnt[[name]] <- cerv_ct
    
    datalist_cerv_sc_len[[name]] <- cerv_len
  }
}


big_data_sc_palm_cnt <- dplyr::bind_rows(datalist_palm_sc_cnt) %>% distinct() %>% mutate(species="palm")
big_data_sc_palm_len <- dplyr::bind_rows(datalist_palm_sc_len) %>% distinct()%>% mutate(species="palm")

#write.table(big_data_sc_palm_cnt, "palmTracks_cnt_snpchip.txt", sep="\t",row.names=FALSE)
#write.table(big_data_sc_palm_len, "palmTracks_len_snpchip.txt", sep="\t",row.names=FALSE)

big_data_sc_cerv_cnt <- dplyr::bind_rows(datalist_cerv_sc_cnt) %>% distinct()%>% mutate(species="cerv")
big_data_sc_cerv_len <- dplyr::bind_rows(datalist_cerv_sc_len) %>% distinct()%>% mutate(species="cerv")

#write.table(big_data_sc_cerv_cnt, "cervTracks_cnt_snpchip.txt", sep="\t",row.names=FALSE)
#write.table(big_data_sc_cerv_len, "cervTracks_len_snpchip.txt", sep="\t",row.names=FALSE)

big_data_hy_sc_avg <- dplyr::bind_rows(datalist_hy_sc_avg) %>% distinct() 
write.table(big_data_hy_sc_avg, "shared_F1_Tracks_snpchip.txt", sep="\t",row.names=FALSE)

big_data_bc_sc_avg <- dplyr::bind_rows(datalist_bc_sc_avg) %>% distinct() 
write.table(big_data_bc_sc_avg, "shared_BC_Tracks_snpchip.txt", sep="\t",row.names=FALSE)

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

pchr<-big_data_cnt %>% filter(mg==1)%>% group_by(chr,species,hyb,method) %>% mutate(cnt=sum(cnt)) %>%
  select(cnt,hyb,method,chr,species)%>% distinct()

g5 <-ggplot(big_data_cnt %>% filter(mg==1)%>% group_by(chr,species,hyb,method) %>% mutate(cnt=sum(cnt)) %>%
              select(cnt,hyb,method,chr,species)%>% distinct(),
            aes(x=species, y=cnt, color=factor(hyb), fill=factor(hyb)), group=factor(method)) +
  geom_point(aes(shape=method, color=hyb),position=position_jitterdodge(jitter.height=0.1,dodge.width=0.8),alpha=0.4)+ 
  theme_bw()+
  theme(panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
        strip.background = element_rect(colour="grey30", fill="grey95",
                                        size=1.1, linetype="solid"),
        strip.text.x = element_text(size=12, color="black",
                                    face="bold")) +
  labs(y="Count of Parental Ancestry Tracks / Chromosome",
       x="")+
  scale_color_manual(values = c("F1" = "#6F81BC","F1-x" = "#485b99", "Bc" = "#6D2B75"))+
  scale_fill_manual(values = c("F1" = "#ACB6D8", "F1-x" = "#485b99","Bc" = "#6D2B75"))+
  scale_shape_manual(values=c(21,24))+
  facet_wrap(~ method)
g6<-g5 + stat_summary(fun.data=data_summary, geom="pointrange",shape=c(18),lwd=1,position=position_dodge(width=0.8))
g6

# length
g7 <- ggplot(big_data_len %>% filter(mg==1), aes(x=species, y=(len/1000000), color=factor(hyb), fill=factor(hyb)), group=factor(method)) +
  geom_point(aes(shape=method, color=hyb),position=position_jitterdodge(jitter.height=0.2,dodge.width=0.8),alpha=0.4) +
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
g8

library(ggpubr)
pdf(paste0(gen,"../../track_lengths_avg10/tracks_genome_snpchip.pdf"), width = 5, height = 5)
ggarrange(g6, g8, labels = c("b", "c"),nrow = 2)
dev.off()

# strip-plot by population, snpchip samples
library(forcats)

big_data_sc_palm_cnt_sum <- dplyr::bind_rows(datalist_palm_sc_cnt) %>% distinct() %>% mutate(species="palm") %>%
  group_by(ID,species,mg,hyb,pop) %>%
  summarize(cnt=sum(cnt)) %>%
  ungroup()%>%
  mutate(ID=fct_reorder(ID,pop))

big_data_sc_palm_len_sum <- dplyr::bind_rows(datalist_palm_sc_len) %>% distinct()%>% mutate(species="palm") %>%
  group_by(ID,species,mg,hyb,pop) %>%
  summarize(len=sum(len)) %>%
  ungroup() %>%
  mutate(ID=fct_reorder(ID,pop))

big_data_sc_cerv_cnt_sum <- dplyr::bind_rows(datalist_cerv_sc_cnt) %>% distinct()%>% mutate(species="cerv")%>%
  group_by(ID,species,mg,hyb,pop) %>%
  summarize(cnt=sum(cnt)) %>%
  ungroup()%>%
  mutate(ID=fct_reorder(ID,pop))

big_data_sc_cerv_len_sum <- dplyr::bind_rows(datalist_cerv_sc_len) %>% distinct()%>% mutate(species="cerv")%>%
  group_by(ID,species,mg,hyb,pop) %>%
  summarize(len=sum(len)) %>%
  ungroup()%>%
  mutate(ID=fct_reorder(ID,pop)) 

big_data_cnt_sum <- big_data_sc_palm_cnt_sum %>% bind_rows(big_data_sc_cerv_cnt_sum)

big_data_len_sum <- big_data_sc_palm_len_sum %>% bind_rows(big_data_sc_cerv_len_sum) 

g9 <-ggplot(big_data_cnt_sum %>% filter(mg==1, hyb=="F1-x"), aes(x=ID, y=cnt, color=factor(species))) +
  geom_point(stat='identity',position=position_dodge(width=0.3),alpha=0.6, size=2.2) + 
  theme_bw()+coord_flip()+ylim(0,15)+
  theme(plot.margin=margin(0, 5, 0, 2, unit = "pt"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  labs(y="Count of Parental Ancestry Tracks",
       x="")+
  scale_color_manual(values = c("palm" = "#B33B77","cerv" = "#FBB116"))+
  scale_fill_manual(values = c("palm" = "#B33B77","cerv" = "#FBB116"))

g9b <-ggplot(big_data_len_sum %>% filter(mg==1, hyb=="F1-x",!len == 0), aes(x=ID, y=len/1000000, color=factor(species))) +
  geom_point(stat='identity',position=position_dodge(width=0.3),alpha=0.6, size=2.2) + 
  theme_bw()+coord_flip()+ylim(0,15)+
  theme(plot.margin=margin(0, 5, 1, 2, unit = "pt"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  labs(y="Length of Parental Ancestry Tracks",
       x="")+
  scale_color_manual(values = c("palm" = "#B33B77","cerv" = "#FBB116"))+
  scale_fill_manual(values = c("palm" = "#B33B77","cerv" = "#FBB116"))

g10 <-ggplot(big_data_cnt_sum %>% filter(mg==1,hyb=="F1"), aes(x=ID, y=cnt, color=factor(species))) +
  geom_point(stat='identity',position=position_dodge(width=0.3), alpha=0.6, size=2.2) + 
  theme_bw()+coord_flip()+ylim(0,15)+
  theme(plot.margin=margin(0, 5, 0, 2, unit = "pt"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  labs(y="Count of Parental Ancestry Tracks",
       x="")+
  scale_color_manual(values = c("palm" = "#B33B77","cerv" = "#FBB116"))+
  scale_fill_manual(values = c("palm" = "#B33B77","cerv" = "#FBB116"))

g10b <-ggplot(big_data_len_sum %>% filter(mg==1, hyb=="F1",!len == 0) , aes(x=ID, y=len/1000000, color=factor(species))) +
  geom_point(stat='identity',position=position_dodge(width=0.3), alpha=0.6, size=2.2) + 
  theme_bw()+coord_flip()+ylim(0,15)+
  theme(plot.margin=margin(0, 5, 1, 2, unit = "pt"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())+
  labs(y="Length of Parental Ancestry Tracks",
       x="")+
  scale_color_manual(values = c("palm" = "#B33B77","cerv" = "#FBB116"))+
  scale_fill_manual(values = c("palm" = "#B33B77","cerv" = "#FBB116"))

g11 <-ggplot(big_data_cnt_sum %>% filter(mg==2,hyb=="Bc"), aes(x=ID, y=cnt, color=factor(species))) +
  geom_point(stat='identity',position=position_dodge(width=0.3),alpha=0.6, size=2.2) + 
  theme_bw()+coord_flip()+ ylim(0,15)+
  theme(plot.margin=margin(0, 5, 2, 2, unit = "pt"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  labs(y="Count of Parental Ancestry Tracks",
       x="")+
  scale_color_manual(values = c("palm" = "#B33B77","cerv" = "#FBB116"))+
  scale_fill_manual(values = c("palm" = "#B33B77","cerv" = "#FBB116"))

g11b <-ggplot(big_data_len_sum %>% filter(mg==2,hyb=="Bc",!len == 0), aes(x=ID, y=len/1000000, color=factor(species))) +
  geom_point(stat='identity',position=position_dodge(width=0.3),alpha=0.6, size=2.2) + 
  geom_hline(yintercept =15, size=0.8)+
  theme_bw()+coord_flip()+
  theme(plot.margin=margin(0, 5, 2, 2, unit = "pt"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  labs(y="Total Length of Parental Ancestry Tracks (Mb)",
       x="")+
  scale_color_manual(values = c("palm" = "#B33B77","cerv" = "#FBB116"))+
  scale_fill_manual(values = c("palm" = "#B33B77","cerv" = "#FBB116"))

ggarrange(g9 + rremove("x.title") + rremove("x.text")+ rremove("x.ticks"), 
          g9b + rremove("x.title"),
          g10 +  rremove("x.title") +rremove("x.text")+ rremove("x.ticks"),
          g10b +  rremove("x.title") ,
          g11, g11b,
          nrow = 3,ncol=2, align="v", common.legend = TRUE, legend="top", 
          heights = c(0.3,1.7,0.25))


### genome samples
big_data_gen_palm_cnt_sum <- big_data_palm_cnt %>%
  group_by(ID,species,mg,hyb,pop) %>%
  summarize(cnt=sum(cnt)) %>%
  ungroup()%>%
  mutate(ID=fct_reorder(ID,pop))

big_data_gen_palm_len_sum <- big_data_palm_len %>%
  group_by(ID,species,mg,hyb,pop) %>%
  summarize(len=sum(len)) %>%
  ungroup() %>%
  mutate(ID=fct_reorder(ID,pop))

big_data_gen_cerv_cnt_sum <- big_data_cerv_cnt %>%
  group_by(ID,species,mg,hyb,pop) %>%
  summarize(cnt=sum(cnt)) %>%
  ungroup()%>%
  mutate(ID=fct_reorder(ID,pop))

big_data_gen_cerv_len_sum <- big_data_cerv_len %>% 
  group_by(ID,species,mg,hyb,pop) %>%
  summarize(len=sum(len)) %>%
  ungroup()%>%
  mutate(ID=fct_reorder(ID,pop)) 

big_data_cnt_sum <- big_data_gen_palm_cnt_sum %>% bind_rows(big_data_gen_cerv_cnt_sum)

big_data_len_sum <- big_data_gen_palm_len_sum %>% bind_rows(big_data_gen_cerv_len_sum)

g12 <-ggplot(big_data_cnt_sum %>% filter(mg==1, hyb=="F1"), aes(x=ID, y=cnt, color=factor(species))) +
  geom_point(stat='identity',position=position_dodge(width=0.3),alpha=0.6, size=2.2) + 
  theme_bw()+coord_flip()+ylim(0,30)+
  theme(plot.margin=margin(0, 5, 0, 2, unit = "pt"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  labs(y="Count of Parental Ancestry Tracks",
       x="")+
  scale_color_manual(values = c("palm" = "#B33B77","cerv" = "#FBB116"))+
  scale_fill_manual(values = c("palm" = "#B33B77","cerv" = "#FBB116"))

g12b <-ggplot(big_data_len_sum %>% filter(mg==1, hyb=="F1",!len == 0), aes(x=ID, y=len/1000000, color=factor(species))) +
  geom_point(stat='identity',position=position_dodge(width=0.3),alpha=0.6, size=2.2) + 
  theme_bw()+coord_flip()+ylim(0,15)+
  theme(plot.margin=margin(0, 5, 1, 2, unit = "pt"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  labs(y="Length of Parental Ancestry Tracks",
       x="")+
  scale_color_manual(values = c("palm" = "#B33B77","cerv" = "#FBB116"))+
  scale_fill_manual(values = c("palm" = "#B33B77","cerv" = "#FBB116"))

g13 <-ggplot(big_data_cnt_sum %>% filter(mg==2,hyb=="Bc"), aes(x=ID, y=cnt, color=factor(species))) +
  geom_point(stat='identity',position=position_dodge(width=0.3),alpha=0.6, size=2.2) + 
  theme_bw()+coord_flip()+ ylim(0,30)+
  theme(plot.margin=margin(0, 5, 2, 2, unit = "pt"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  labs(y="Count of Parental Ancestry Tracks",
       x="")+
  scale_color_manual(values = c("palm" = "#B33B77","cerv" = "#FBB116"))+
  scale_fill_manual(values = c("palm" = "#B33B77","cerv" = "#FBB116"))

g13b <-ggplot(big_data_len_sum %>% filter(mg==2,hyb=="Bc",!len == 0), aes(x=ID, y=len/1000000, color=factor(species))) +
  geom_point(stat='identity',position=position_dodge(width=0.3),alpha=0.6, size=2.2) + 
  geom_hline(yintercept =15, size=0.8)+
  theme_bw()+coord_flip()+
  theme(plot.margin=margin(0, 5, 2, 2, unit = "pt"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  labs(y="Total Length of Parental Ancestry Tracks (Mb)",
       x="")+
  scale_color_manual(values = c("palm" = "#B33B77","cerv" = "#FBB116"))+
  scale_fill_manual(values = c("palm" = "#B33B77","cerv" = "#FBB116"))

ggarrange(g12 + rremove("x.title") , 
          g12b + rremove("x.title"),
          g13, g13b,
          nrow = 2,ncol=2, align="v", common.legend = TRUE, legend="top", 
          heights = c(1.75,0.25))
