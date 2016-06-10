#Set operations using BEDOPS
library(dplyr)
library(ggplot2)


FindMax <- function(df){
  stopifnot(!is.null(df$chrom),!is.null(df$start),!is.null(df$end),!is.null(df$name))
  nndf <- CollapseNested(df,concat.name = F)

  nndf <- arrange(nndf,chrom,start,end) %>% group_by(chrom) %>% mutate(prevend=lag(end,default = 0),
                                                                       prevstart=lag(start,default = 0),
                                                                       nstart=lead(start,default=0),
                                                                       nend=lead(end,default=0),
                                                                       prevsize=prevend-prevstart,
                                                                       nsize=nend-nstart,
                                                                       size=end-start) %>% ungroup()
  i <- 0
  nndf <- nndf[!duplicated(select(nndf,chrom,start,end)),]
  while(sum((nndf$prevend>nndf$start)|((nndf$nstart<nndf$end)&(nndf$nstart>0)))>0){
    # print(i)
    # print(sum((nndf$prevend>nndf$start)|((nndf$nstart<nndf$end)&(nndf$nstart>0))))
    nndf <-  group_by(nndf,chrom) %>% mutate(kill=ifelse(((prevend>start)&(prevsize>=size))|((nstart<end)&(nsize>=size)&(nstart>0)),TRUE,FALSE)) %>% ungroup()
    nndf <- filter(nndf,kill==FALSE)
    nndf <- arrange(nndf,chrom,start,end) %>% group_by(chrom) %>% mutate(prevend=lag(end,default = 0),
                                                                         prevstart=lag(start,default = 0),
                                                                         nstart=lead(start,default=0),
                                                                         nend=lead(end,default=0),
                                                                         prevsize=prevend-prevstart,
                                                                         nsize=nend-nstart) %>% ungroup()
    i <- i+1

  }
  nndf <- select(nndf,chrom,start,end,name)
  return(nndf)
}



CollapseNested <- function(df,concat.name=T){
  stopifnot(!is.null(df$chrom),!is.null(df$start),!is.null(df$end),!is.null(df$name))
  df <- arrange(df,chrom,start,end) %>% group_by(chrom) %>% mutate(prevend=lag(end,default = 0),prevstart=lag(start,default = 0)) %>% ungroup()
  sum(df$prevend>df$end)
  i <- 0
  while(sum(df$prevend>df$end)>0){
    #    print(i)
    #    print(sum(df$prevend>df$end))
    df <- group_by(df,chrom) %>% mutate(start=ifelse(prevend>end,prevstart,start),
                                        end=ifelse(prevend>end,prevend,end),
                                        kill=ifelse(prevend>end,TRUE,FALSE)) %>% ungroup()

    if(concat.name){
      df <- group_by(df,chrom,start,end) %>% summarise(name=paste(name,collapse=",")) %>% arrange(chrom,start,end)
    }else{
      df <- filter(df,kill==FALSE)
    }
    df <- group_by(df,chrom) %>% mutate(prevend=lag(end,default=0),prevstart=lag(start,default=0)) %>% ungroup()
    i <- i+1
    if(i>50){
      print("Quit!")
      break
    }
  }

  df <- mutate(df,name=sapply(lapply(strsplit(name,split = ",",fixed = T),unique),paste,collapse=",")) %>% arrange(chrom,start,end)
  df <- select(df,chrom,start,end,name)
  return(df)
}

MergeReg <- function(df){
  stopifnot(!is.null(df$chrom),!is.null(df$start),!is.null(df$end),!is.null(df$name))
  nndf <- CollapseNested(df)

  nndf <- arrange(nndf,chrom,start,end) %>% group_by(chrom) %>% mutate(prevend=lag(end,default = 0),prevstart=lag(start,default = 0),nstart=lead(start,default=0),nend=lead(end,default=0)) %>% ungroup()
  i <- 0
  sum(nndf$prevend>nndf$start)
  while(sum(nndf$prevend>nndf$start)>0){
    #  print(i)
    #   print(sum(nndf$prevend>nndf$start))
    nndf <-  group_by(nndf,chrom) %>% mutate(start=ifelse(prevend>start,prevstart,start)) %>% ungroup()
    nndf <- arrange(nndf,chrom,start,end) %>% group_by(chrom) %>% mutate(prevend=lag(end,default = 0),prevstart=lag(start,default = 0),nstart=lead(start,default=0),nend=lead(end,default=0)) %>% ungroup()
    nndf <- group_by(nndf,chrom) %>% mutate(end=ifelse(nstart==start&nend>end,nend,end)) %>% ungroup()
    nndf <- group_by(nndf,chrom,start,end) %>% summarise(name=paste(name,collapse=",")) %>% arrange(chrom,start,end)
    nndf <- arrange(nndf,chrom,start,end) %>% group_by(chrom) %>% mutate(prevend=lag(end,default = 0),prevstart=lag(start,default = 0),nstart=lead(start,default=0),nend=lead(end,default=0)) %>% ungroup()
    i <- i+1
  }
  nndf <- mutate(nndf,name=sapply(lapply(strsplit(name,split = ",",fixed = T),unique),paste,collapse=",")) %>% arrange(chrom,start,end)
  nndf <- select(nndf,chrom,start,end,name)
  return(nndf)
}



findintervar <- function(regdf,vardf,tempdi=tempdir()){
  varf <- file.path(tempdi,"VARDAT.bed")
  regf <- file.path(tempdi,"REGDAT.bed")
  svarf <- file.path(tempdi,"SVARDAT.bed")
  sregf <- file.path(tempdi,"SREGDAT.bed")
  vardat <- select(vardf,chrom,start,end,ref,alt) %>% mutate(start=format(start,scientific=F,trim=T),
                                                             end=format(end+1,scientific = F,trim=T))
  regdat <- select(regdf,chrom,start,end,name) %>% mutate(start=format(start,scientific=F,trim=T),
                                                          end=format(end,scientific=F,trim=T))
  write.table(vardat,
              file=varf,
              sep="\t",
              col.names=F,
              row.names=F,
              quote=F)
  write.table(regdat,
              file=regf,
              sep="\t",
              col.names=F,
              row.names=F,
              quote=F)
  system(paste0("sort-bed ",varf, " > ",svarf))
  system(paste0("sort-bed ",regf, " > ",sregf))
  varmap <- system(paste0("bedmap --ec --faster --echo-map-id ",svarf," ",sregf),intern = T)
  nvd <- read.table(svarf,header=F,sep="\t",stringsAsFactors = F)
  nvd <- rename(nvd,chrom=V1,start=V2,end=V3,ref=V4,alt=V5) %>% mutate(end=end-1,map=varmap)
  system(paste("rm ",varf,regf,sregf))
  return(nvd)
}


findinter <- function(genedf,cnvdat,tmpdir=tempdir(),bedmap=NULL){
  stopifnot(!is.null(genedf$chrom),
            !is.null(genedf$start),
            !is.null(genedf$end),
            !is.null(genedf$name),
            !is.null(cnvdat$sampleid),
            !is.null(cnvdat$cnvstart),
            !is.null(cnvdat$cnvstop),
            !is.null(cnvdat$chrom),
            !is.null(cnvdat$cnvid))
  if(is.null(bedmap)){
    bedmap <- system("which bedmap",intern = T)
    bedsort <- system("which sort-bed",intern=T)
  }
  cnvf <- file.path(tmpdir,"CNVDAT.bed")
  scnvf <- file.path(tmpdir,"SCNVDAT.bed")
  genef <- file.path(tmpdir,"GENEDAT.bed")
  sgenef <- file.path(tmpdir,"SGENEDAT.bed")
  cnvdat <- select(cnvdat,chrom,cnvstart,cnvstop,cnvid,sampleid)
  genedf <- select(genedf,chrom,start,end,name)
  write.table(cnvdat,
              file=cnvf,
              sep="\t",
              col.names=F,
              row.names=F,
              quote=F)
  write.table(genedf,
              file=genef,
              sep="\t",
              col.names=F,
              row.names=F,
              quote=F)
  system(paste0(bedsort," ",cnvf, " > ",scnvf))
  system(paste0(bedsort," ",genef, " > ",sgenef))
  cnvm <- system(paste0(bedmap," --faster --echo-map-id ",scnvf," ",sgenef),intern = T)
  cnvdat <- read.table(scnvf,header=F,sep="\t",stringsAsFactors = F)
  cnvdat <- rename(cnvdat,chrom=V1,cnvstart=V2,cnvstop=V3,cnvid=V4,sampleid=V5)

  genedf <- read.table(sgenef,header=F,sep="\t",stringsAsFactors = F)
  genedf <- rename(genedf,chrom=V1,start=V2,end=V3,name=V4)
  file.remove(sgenef,scnvf)
  gene.id <- bind_rows(mapply(FUN = function(x,y,z){
    genes <- strsplit(y,split = ";",fixed = T)
    if(length(genes[[1]])==0){
      genes[[1]] <- "NA"
    }
    return(data.frame(cnvid=x,sampleid=z,gene=genes[[1]],stringsAsFactors = F))
  },cnvdat[["cnvid"]],cnvm,cnvdat[["sampleid"]],SIMPLIFY = F))
  return(gene.id)
}


annovar_bed <- function(bdf,dbdir="/media/nwknoblauch/Data/annovar/humandb/",tmpdir,name){

  outfile <- file.path(tmpdir,paste0(name,"_Results_EV"))
  outfilet <- file.path(tmpdir,paste0(name,"_Results_T"))
  aexonf<- paste0(outfile,".exonic_variant_function")
  aconsf <- paste0(outfilet,".hg19_multianno.txt")

  if(!any(file.exists(c(aexonf,aconsf)))){

    infile <- file.path(tmpdir,"Input.bed")
    write.table(bdf,infile,sep="\t",col.names=F,row.names=F,quote=F)
    cmd <- paste("/media/nwknoblauch/Data/annovar/annotate_variation.pl -outfile",outfile,"-build hg19",infile,dbdir,"--geneanno","-dbtype refGene")
    system(cmd)
    cmdt <- paste("/media/nwknoblauch/Data/annovar/table_annovar.pl",
                  infile,
                  dbdir,
                  "--outfile ",outfilet,
                  " -protocol ljb26_all -operation f -build hg19 -otherinfo -nastring .")
    system(cmdt)

  }
  ccol <- scan(aconsf,what=character(),nlines=1,sep="\t")
  allcol <- c(ccol[-length(ccol)],colnames(bdf)[-(1:5)])
  consdf <- read.table(aconsf,header=F,skip = 1,sep="\t",stringsAsFactors = F)
  colnames(consdf) <- allcol
  consdf <- mutate_each(consdf,funs(ifelse(.==".",NA,.)))
  consdf <- rename(consdf,chrom=Chr,start=Start,end=End,ref=Ref,alt=Alt)
  vardf <- read.table(aexonf,header=F,sep="\t",stringsAsFactors=F)
  vardf <- select(vardf,chrom=V4,start=V5,end=V6,ref=V7,alt=V8,consequence=V2)


  bdf <- left_join(bdf,vardf)
  bdf <- mutate(bdf,consequence=ifelse(is.na(consequence),"inter",consequence))
  bdf <- mutate(bdf,ac_adj=as.numeric(ac_adj))

  bdf <- arrange(bdf,chrom,start,end,ref,alt)
  consdf <- arrange(consdf,chrom,start,end,ref,alt)
  rbdf <- bind_cols(bdf,select(consdf,-chrom,-start,-end,-ref,-alt,-ac_adj,-an_adj,-af_adj))
  rbdf <- select(rbdf,-one_of("Otherinfo"))
  return(rbdf)
}

splitdf <- function(df,synonf,nsynonf){
  nsyndf <- filter(df,consequence=="nonsynonymous SNV") %>% mutate(end=format(end+1,trim = T,scientific = F),
                                                                   start=format(start,trim=T,scientific=F)) %>%
    select(chrom,start,end,ref,alt)
  syndf <- filter(df,consequence=="synonymous SNV") %>%  mutate(end=format(end+1,trim = T,scientific = F),
                                                                start=format(start,trim=T,scientific=F)) %>%
    select(chrom,start,end,ref,alt)
  write.table(nsyndf,nsynonf,col.names = F,sep="\t",row.names=F,quote=F)
  write.table(syndf,synonf,col.names = F,sep="\t",row.names=F,quote=F)
}







BedC <- function(varfile,regionfile,regionheader=T,outdir=tempdir(),fieldname){
  bedsort <- system("which sort-bed",intern=T)
  bedmap <- system("which bedmap",intern=T)
  csf <- file.path(outdir,"tvarfile.txt")
  crf <- file.path(outdir,"trfile.txt")
  system(paste0(bedsort," ",varfile, " > ",csf))
  if(regionheader){
    system(paste0("tail -n +2 ",regionfile," | ",bedsort," - > ",crf))
  }else{
    system(paste0("cat ",regionfile," | ",bedsort," - > ",crf))
  }
  cnvm <- system(paste0(bedmap," --faster --count ",crf," ",csf),intern = T)
  regdf <- read.table(crf,header=F,sep="\t",stringsAsFactors = F)
  regdf <- rename(regdf,chrom=V1,start=V2,end=V3)
  regdf[[fieldname]] <- cnvm
  return(regdf)
}

varcount <- function(synvarf,nsynvarf,regf,cov_exdf){
  synncount <- BedC(synvarf,regf,regionheader = T,fieldname="n_syn")
  synncount <- rename(synncount,gene=V4,syn_rate=V5,nsyn_rate=V6) %>% mutate(n_syn= as.numeric(n_syn))
  nsynncount <- BedC(nsynvarf,regf,regionheader = T,fieldname="n_nsyn")
  nsynncount <- rename(nsynncount,gene=V4,syn_rate=V5,nsyn_rate=V6) %>% mutate(n_nsyn=as.numeric(n_nsyn))
  nexondf <- inner_join(nsynncount,synncount)
  nexondf <- left_join(nexondf,cov_exdf)
  return(nexondf)
}




GCstarch <- function(starchfile,vardf,tmpdir=tempdir()){

  varf <- file.path(tmpdir,"VARDAT.bed")
  svarf <- file.path(tmpdir,"SVARDAT.bed")
  vardat <- select(vardf,chrom,start,end,ref,alt) %>% mutate(start=format(start,scientific=F,trim=T),
                                                             end=format(end+1,scientific = F,trim=T))
  write.table(vardat,
              file=varf,
              sep="\t",
              col.names=F,
              row.names=F,
              quote=F)
  system(paste0("sort-bed ",varf, " > ",svarf))
  varmap <- system(paste0("bedmap --faster --mean ",svarf," ",starchfile),intern = T)
  nvd <- read.table(svarf,header=F,sep="\t",stringsAsFactors = F)
  nvd <- rename(nvd,chrom=V1,start=V2,end=V3,ref=V4,alt=V5) %>% mutate(end=end-1,map=varmap)
#  system(paste("rm ",svarf,varf,regf,sregf))
  return(nvd)
}


ExacCov <- function(covdir,exonf,covmf){
  tempd <- tempdir()
  oexonf <- file.path(tempd,"sexf.txt")
  system(paste("tail -n +2 ",exonf,"| sort-bed - >",oexonf))
  starch.file <- file.path(covdir,"all_cov.starch")
  bedfiles <- file.path(covdir,paste0("file_",c(1:22,"X","Y"),".bed"))
  gzfiles <- dir(covdir,pattern = ".gz",full.names = T)
  if(!file.exists(starch.file)){
    if(length(gzfiles)>0){
      sapply(gzfiles,function(x){
        cat(x)
        cat("\n")
        filechrom <- gsub("^.+chr(.+).coverage.txt.gz","\\1",x)
        filedir <- gsub("(.+)Panel.chr.+.txt.gz","\\1",x)
        outfile <- file.path(filedir,paste0("file_",filechrom,".bed"))

        tf <- read.table(x,header=F,sep="\t",col.names=c(c("chrom","start","mean"),4:13),
                                      colClasses = c("character","integer","numeric",rep("NULL",10)))
        nr <- nrow(tf)
        tf <- mutate(tf,end=format(start+1,scientific = F,trim=T),start=format(start,scientific=F,trim=T),chrom=paste0("chr",chrom),id=paste0(chrom,1:nr)) %>% select(chrom,start,end,id,mean)
        write.table(tf,outfile,col.names=F,row.names=F,quote=F,sep="\t")
      })
    }
    system(paste("sort-bed ",paste0(bedfiles,collapse=" "),"| starch - >",starch.file))
  }
  if(!file.exists(covmf)){
    system(paste("bedmap --faster --mean --median --min --variance ",oexonf,starch.file,"> ",covmf))
  }
  exdf <- read.table(oexonf,header=F,sep="\t",stringsAsFactors = F)
  exdf <- rename(exdf,chrom=V1,start=V2,end=V3,gene=V4,syn_rate=V5,nsyn_rate=V6)
  covm <- read.table(covmf,sep="|",header=F,stringsAsFactors=F,colClasses = rep("numeric",4),na.strings = "NAN")
  covm <- rename(covm,mean=V1,median=V2,min=V3,variance=V4)
  exdf <- bind_cols(exdf,covm)
  return(exdf)
}



