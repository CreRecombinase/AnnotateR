#New Schizophrenia dataset
source("~/Dropbox/AnnotateR/SetOps.R")
varf <- "~/Downloads/swedexgf.counts.gz"
exonf <- "~/Dropbox/BayesianDA/Exac/Exon_Convservation.RDS"
tblf <- "~/Dropbox/BayesianDA/VariantAnnotation/swedschiz/swed_schiz_anno.db"


totdf <- read.table(varf,header=T,sep="\t",stringsAsFactors = F)
totdf <- mutate(totdf,chrom= sapply(strsplit(VAR,split = ":",fixed = T),"[",1),POS=sapply(strsplit(VAR,split=":",fixed=T),"[",2))
totdf <- mutate(totdf,start=as.numeric(sapply(strsplit(POS,"..",fixed = T),"[",1)),
                ref=sapply(strsplit(REF.ALT,"/",fixed=T),"[",1),
                end=start+nchar(ref)-1)
totdf <- mutate(totdf,alt=sapply(strsplit(REF.ALT,"/",fixed=T),"[",2)) %>% mutate(alt=sapply(strsplit(alt,",",fixed=T),"[",1))%>%
  select(chrom,start,end,ref,alt,cnta=CNTA,tota=TOTA,cntu=CNTU,totu=TOTU)
totdf <- mutate(totdf,evid=paste0(chrom,":",start,"-",end,"_",ref,"_",alt))
casedf <- select(totdf,evid,cnta,tota)
ctrldf <- select(totdf,evid,cntu,totu)
posdf <- select(totdf,chrom,start,end,ref,alt)

exondf <- readRDS(exonf)
exondf <- group_by(exondf,gene) %>% mutate(name=paste0(gene,":",1:n()))%>% ungroup
zexondf <- select(exondf,name,z_sign_mis,syn_rate,nsyn_rate,n_nsyn_exac=n_nsyn,n_syn_exac=n_syn)
texondf <- select(exondf,chrom,start,end,name)
#texondf <- group_by(texondf,gene) %>% mutate(name=paste0(gene,":",1:n())) %>% ungroup

posdf <- findintervar(texondf,posdf,tempdir="/tmp")
posdf <- mutate(posdf,map=ifelse(map=="","inter",map),evid=paste0(chrom,":",start,"-",end,"_",ref,"_",alt))
consdf <- annovar_bed(bdf=posdf,name="swed")
consdf <- mutate(consdf,evid=paste0(chrom,":",start,"-",end,"_",ref,"_",alt))
consdf <- select(consdf,-one_of(c("chrom","start","end","ref","alt")))
saveRDS(posdf,annomatf)
saveRDS(consdf,fpredf)
gcdf <- GCstarch("~/Desktop/VariantPrioritization/GChg19.starch",posdf)
gcdf <- mutate(gcdf,evid=paste0(chrom,":",start,"-",end,"_",ref,"_",alt))
saveRDS(gcdf,gcf)
gcdf <- select(gcdf,evid,gc=map)


casedf <- mutate(casedf,
                 alts=cnta,
                 refs=tota-alts,Pheno="Schizophrenia") %>% select(evid,alts,refs,Pheno)
ctrldf <- mutate(ctrldf,
                 alts=cntu,
                 refs=totu-alts,Pheno="Control") %>% select(evid,alts,refs,Pheno)


swedtbl <- src_sqlite(tblf,create = T)
                 
copy_to(swedtbl,casedf,"CaseDat",temporary = F,indexes = list(c("evid")))
copy_to(swedtbl,ctrldf,"CtrlDat",temporary=F,indexes=list(c("evid")))
copy_to(swedtbl,gcdf,"GCdat",temporary=F,indexes=list(c("evid")))
copy_to(swedtbl,consdf,"MutationConsequence",temporary=F,indexes=list(c("evid")))
copy_to(swedtbl,posdf,"VariantPosition",temporary=F,indexes=list("evid","map"))
copy_to(swedtbl,zexondf,"ExacExonCons",temporary=F,indexes=list("name"))



