#read in command line args
args<-commandArgs(TRUE)
phen<-args[1]
spec<-args[2]

#Summary Statistics for Hyperparameters and construction of Manhattan Plots.
#Based on code available at: http://romainvilloutreix.alwaysdata.net/romainvilloutreix/wp-content/uploads/2017/01/gwas_gemma-2017-01-17.pdf

library("data.table")
#read in the hyp file
hyperParameters <- read.table(paste0("resids_allchr_",spec,"_",phen,".hyp.txt"),header=T)

#get the 95% equal tailed probability interval (credible interval)
#h denotes the approximation to the proportion of phenotypic variance explained by the PVE
h <- c("h",mean(hyperParameters$h),quantile(hyperParameters$h,probs = c(0.5,0.025,0.975)))

#pve denotes the proportion of phenotypic variance explained by sparse and random effects
pve <- c("pve",mean(hyperParameters$pve),quantile(hyperParameters$pve,probs=c(0.5,0.025,0.975)))

#pge denotes the proportion of PVE explained by the sparse effects only
pge <- c("pge",mean(hyperParameters$pge),quantile(hyperParameters$pge,probs=c(0.5,0.025,0.975)))

#rho denotes the approximation to the proportion of genetic variance explained by the variants with a major effect (PGE)
#rho = 0 has a highly polygenic basis
#rho = 1 means pure BVSR, few major effect loci
rho <- c("rho",mean(hyperParameters$rho),quantile(hyperParameters$rho,probs=c(0.5,0.025,0.975)))

#pi denotes the proportion of variants with nonzero effects.
pi <- c("pi",mean(hyperParameters$pi),quantile(hyperParameters$pi,probs=c(0.5,0.025,0.975)))

#nGamma denotes the number of variants with a major effect.
nGamma <- c("nGamma",mean(hyperParameters$n_gamma),quantile(hyperParameters$n_gamma,probs=c(0.5,0.025,0.975)))

#Put together a data table containing summary statistics for the hyperparamters
summaryTable <- as.data.frame(rbind(h,pve,pge,rho,pi,nGamma),row.names=F)
colnames(summaryTable) <- c("hyperparameter","mean","median","2.5%","97.5%")

#Let's see which SNPs have important/nontrivial effects.

#Now we read the parameters of interest. 
params <- as.data.frame(fread(paste0("resids_allchr_", spec, "_", phen, ".param.txt"), header=T, sep="\t"))

#total effect size
params["eff"] <- abs(params$alpha + params$beta*params$gamma)

#get all with effect size > 0
paramsEffects <- params[params$eff>0,]

#number of important effects
nrow(paramsEffects)

#sort the important effects by decreasing size
sortEffects <- paramsEffects[order(-paramsEffects$eff),]
head(sortEffects,10)

#get top 0.01% of variants
top001<-sortEffects[sortEffects$eff>quantile(sortEffects$eff,0.9999),]

#get variants with high Posterior Inclusion Probability
pipSort <- params[order(-params$gamma),]

#plot variant PIPs across chromosomes - can't do this if GEMMA is run chromosome by chromosome
chrs <- gsub("lg|_.+","",params$rs)
params["chr"]<-chrs
linkGroupSort <- params[order(as.numeric(params$chr),params$rs),]
#get list of linkage groups
chromosomeList <- sort(as.numeric(unique(gsub("^.*?p","",chrs)))) #get snps in terms of numbers

par(mfrow=c(1,1))

library(tidyverse)
don <- pipSort %>% 
  
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(ps)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(pipSort, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, ps) %>%
  mutate(POScum=ps+tot)

axisdf = don %>% group_by(chr) %>% summarize(center=( max(POScum) + min(POScum) ) / 2 )

pdf(paste0("man_resids_anc_", spec, "_phen", phen, ".pdf"), paper = "a4r", width = 0, height = 0)
gg <- ggplot(don, aes(x=POScum, y=eff)) +
  
  # Show all points
  geom_point( aes(color=as.factor(chr)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  
  # custom X axis:
  scale_x_continuous(label = axisdf$chr, breaks= axisdf$center) +
  scale_y_continuous(expand = c(0, 0)) +     # remove space between plot area and x axis
  coord_cartesian(ylim = c(0, signif(max(don$eff)+(max(don$eff)*0.1), digits = 6))) +
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +

  # Custom axis labels
  labs(title=paste0("Manhattan Plot for ancestry ", spec, " phenotype ", phen),
      x ="Chromosome", y = "Effect Size") +
  
  # Label the top results
  geom_text(aes(label=ifelse(eff>min(top001$eff),as.character(rs),'')),hjust=0,vjust=0,size=3)
print(gg)
dev.off

write.table(top001, paste0("top001_resids_",spec, "_phen", phen, ".txt"), quote = F, sep = "\t", row.names = F)
write.table(summaryTable, paste0("summaryTable_resids_",spec,"_phen",phen, ".txt"), quote = F, sep = "\t", row.names = F)
