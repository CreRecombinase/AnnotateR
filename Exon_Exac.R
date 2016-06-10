#New calculation of Exac Exon-level constraint.
#This script is written so that it may call code from other files
#But will not require any preprocessing of data (as much as possible)
#3/19/2016
#NWK
library(dplyr)
library(lazyeval)
source("~/Dropbox/AnnotateR/SetOps.R")
source("~/Dropbox/BayesianDA/Exac/ExacFuncs.R")

hdir <- "~/Desktop/Exac/temp/"

dbdir = "/media/nwknoblauch/Data/annovar/humandb/"
ExacDataDir <- "~/Desktop/Exac/temp/"
ResultDir <- "~/Dropbox/BayesianDA/Exac/"
covdir <- "/media/nwknoblauch/Data/ExacCoverage/"

genef <- file.path(hdir,"Kaitlin_esp6500_ac10_Zdata_fordist.csv")
synonf <- file.path(hdir,"MutRate_exomeWide.synonymous.txt")
nsynonf <- file.path(hdir,"MutRate_exomeWide.nonsynonymous.txt")
Exac_All_VCF <- file.path(hdir,"ExAC.r0.3.sites.vep.vcf.gz")
Exac_All_of <- file.path(hdir,"ExAC_All")
Exac_All_info <- paste0(Exac_All_of,".INFO")



Exac_VCF <- Exac_All_VCF
Exac_of <- Exac_All_of
Exac_info <- Exac_All_info
#Exac_VCF <- file.path(ExacDataDir,paste0("ExAC",VariantClass,".vcf.gz"))
# Exac_of <- file.path(ExacDataDir,paste0("ExAC_",VariantClass))
#Exac_info <- paste0(Exac_of,".INFO")

#VCFFunc(VCFF,gzipped=T,fields = c("AC_Adj","AN_Adj","VQSLOD"),outfile = Exac_All_of)
VCFFunc(Exac_All_VCF,gzipped = T,fields = c("AC_Adj","AN_Adj","VQSLOD"),outfile = Exac_All_of)
#VCFFunc(Exac_NP_VCF,gzipped = T,fields = c("AC_Adj","AN_Adj","VQSLOD"),outfile = Exac_NP_of)

if(RegionClass=="Exon"){
  synondf <- read.table(synonf,header=F,stringsAsFactors = F)
  nsynondf <- read.table(nsynonf,header=F,stringsAsFactors = F)
  colnames(synondf) <- c("chr","start","end","gene","syn_rate")
  colnames(nsynondf) <- c("chr","start","end","gene","nsyn_rate")
  region_df <- full_join(synondf,nsynondf)
}else{
  genedf <- read.table(genef,header=T,stringsAsFactors = F,sep="\t")
  sgenedf <- select(genedf,chr,start,end=stop,gene,syn_rate=p_syn,nsyn_rate=p_mis) %>% mutate(chr=paste0("chr",chr))
  region_df  <- sgenedf
}

Exac_bed_full <- info.bed(Exac_info) %>% mutate(af_adj=as.numeric(ac_adj)/as.numeric(an_adj))
Exac_bedd <-filter(Exac_bed_full,af_adj<0.001)
Exac_rdsf <- file.path(ExacDataDir,paste0("Reference_Exac_full_df.RDS"))
saveRDS(Exac_bed_full,file = Exac_rdsf)
ref_f <- file.path(ExacDataDir,paste0("Reference_Exac_",RegionClass,".txt"))
covmf = file.path(covdir,paste0(RegionClass,"_covmap.txt"))

region_df <-  semi_join(region_df,rename(FindMax(select(bdf,chrom=region_df,start,end,name=gene)),gene=name)) %>%
  mutate(start=format(start,scientific=F,trim=T),end=format(end,scientific=F,trim=T))

write.table(region_df,ref_f,col.names=T,sep="\t",row.names=F,quote=F)

cov_df <- ExacCov(covdir = covdir,exonf = ref_f,covmf = covmf)
Exac_varf <- annovar_bed(bdf = Exac_bed,dbdir=dbdir,tmpdir = ExacDataDir ,name = VariantClass)

Exac_synvarf <- file.path(ExacDataDir,paste0(VariantClass,"_synon.txt"))
Exac_nsynvarf <- file.path(ExacDataDir,paste0(VariantClass,"_nonsynon.txt"))

splitdf(Exac_varf,Exac_synvarf,Exac_nsynvarf)


Exac_cvf <- varcount(Exac_synvarf,Exac_nsynvarf,ref_f,cov_exdf = cov_df)
Exac_ref_ctf <- file.path(ResultDir,paste0(VariantClass,"_Reference_Exac_",RegionClass,".txt"))


write.table(Exac_cvf,Exac_ref_ctf,col.names=T,row.names = F,quote=F)



exclude.exons <- data_frame(start=c(42976313,
                                    105404399,
                                    100633911,
                                    103381800,
                                    103381800,
                                    144990344,
                                    9056172),
                            gene=c("STARD9",
                                   "AHNAK2",
                                   "MUC12",
                                   "CCDC168",
                                   "CCDC168",
                                   "PLEC",
                                   "MUC16"))

exclude.gene=c("TTN")


compute_exac_z <- function(df,include.sex=FALSE,exclude.gene=c("TTN"),exlude.exons=NULL){

  if(!include.sex){
    ndf <- filter(df,chrom %in% paste0("chr",1:22))
  }else{
    ndf <- df
  }
  ndf <- filter(ndf,!gene %in% exclude.gene)
  ndf <- filter(ndf,n_syn!=0,n_syn!=0)
  ndf <- mutate(ndf,bp=end-start)
  if(!is.null(exclude.exons)){
    ndf <- anti_join(ndf,exclude.exons)
  }
  ndf <- filter(ndf,!is.na(syn_rate),!is.na(nsyn_rate))


  ndf <- group_by(ndf,gene) %>% mutate(rn=paste0(gene,":",1:n())) %>% ungroup()
  rownames(ndf) <- ndf$rn
  tndf <- filter(ndf,median>50)
  syn.pred.lm.bp <- lm(n_syn~syn_rate+bp,data=tndf)
  ndf$pred_syn <- predict(syn.pred.lm.bp,newdata=ndf)
  ndf <- mutate(ndf,pred_syn)
  ndf$pred_nsyn <- predict(syn.pred.lm.bp,newdata = data_frame(syn_rate=ndf$nsyn_rate,bp=ndf$bp))
  ndf <- mutate(ndf,med_cut_bin=cut(median,breaks=c(0,seq(to=max(median),by=2.0),max(median)+1),labels = F,include.lowest = T))
  ndf <- mutate(ndf,med_cut_bin=(med_cut_bin/max(med_cut_bin,na.rm = T))*max(median))
  #From Samocha et al method (slightly different values from the above fit)
  #Also see page 69 of supplement
  ndf <- mutate(ndf,adj_p_syn=ifelse(median>=50,pred_syn,pred_syn*(0.09189+0.27086*log(pmax(1,median)))))
  ndf <- mutate(ndf,adj_p_nsyn=ifelse(median>=50,pred_nsyn,pred_nsyn*(0.09189+0.27086*log(pmax(1,median)))))
  ndf <- mutate(ndf,adj_p_nsyn=pmax(1,adj_p_nsyn))
  ndf <- mutate(ndf,mis_cs=(n_nsyn-adj_p_nsyn)^2/adj_p_nsyn,
                z_sign_mis=sqrt(mis_cs)*ifelse(n_nsyn>adj_p_nsyn,-1,1))
  ndf <- select(ndf,-med_cut_bin,-mis_cs)
  return(ndf)

}

Exac_X <- filter(NP_nexondf,chrom=="chrX")

Exac_final_auto <- compute_exac_z(df=Exac_cvf,   include.sex = F,exclude.gene=c("TTN"),exlude.exons = exclude.exons)
Exac_final_X <-    compute_exac_z(df=Exac_X,include.sex = T,exclude.gene=c("TTN"),exlude.exons = NULL)

Exac_final <- bind_rows(Exac_final_auto,NP_final_X)

write.table(Exac_final,file.path(ResultDir,paste0(VariantClass,"Exac_",RegionClass,"_Cons.txt")),
            col.names=T,row.names=F,sep="\t",quote=F)

saveRDS(final_exon,file.path(ResultDir,paste0(VariantClass,"Exac_",RegionClass,"_Cons.RDS")))


