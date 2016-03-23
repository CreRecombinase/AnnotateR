#Set operations using BEDOPS
library(dplyr)
library(ggplot2)


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
  

findintervar <- function(regdf,vardf,tempdir="/tmp"){
  varf <- file.path(tempdir,"VARDAT.bed")
  regf <- file.path(tempdir,"REGDAT.bed")
  svarf <- file.path(tempdir,"SVARDAT.bed")
  sregf <- file.path(tempdir,"SREGDAT.bed")
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
  varmap <- system(paste0("bedmap --faster --echo-map-id ",svarf," ",sregf),intern = T)
  nvd <- read.table(svarf,header=F,sep="\t",stringsAsFactors = F)
  nvd <- rename(nvd,chrom=V1,start=V2,end=V3,ref=V4,alt=V5) %>% mutate(end=end-1,map=varmap)
  system(paste("rm ",varf,regf,sregf))
  return(nvd)
}


findinter <- function(genedf,cnvdat,tempdir="/tmp/",bedmap=NULL){
  stopifnot(!is.null(genedf$chrom),
            !is.null(genedf$start),
            !is.null(genedf$end),
            !is.null(genedf$name),
            !is.null(cnvdat$sid),
            !is.null(cnvdat$cnvstart),
            !is.null(cnvdat$cnvstop),
            !is.null(cnvdat$chrom),
            !is.null(cnvdat$cnvid))
  if(is.null(bedmap)){
    bedmap <- system("which bedmap",intern = T)
    bedsort <- system("which sort-bed",intern=T)
  }
  cnvf <- file.path(tempdir,"CNVDAT.bed")
  scnvf <- file.path(tempdir,"SCNVDAT.bed")
  genef <- file.path(tempdir,"GENEDAT.bed")
  sgenef <- file.path(tempdir,"SGENEDAT.bed")
  cnvdat <- select(cnvdat,chrom,cnvstart,cnvstop,cnvid,sid)
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
  cnvdat <- rename(cnvdat,chrom=V1,cnvstart=V2,cnvstop=V3,cnvid=V4,sid=V5)
  
  genedf <- read.table(sgenef,header=F,sep="\t",stringsAsFactors = F)
  genedf <- rename(genedf,chrom=V1,start=V2,end=V3,name=V4)
  gene.id <- bind_rows(mapply(FUN = function(x,y,z){
    genes <- strsplit(y,split = ";",fixed = T)
    if(length(genes[[1]])==0){
      genes[[1]] <- "NA"
    }
    return(data.frame(cnvid=x,sid=z,gene=genes[[1]],stringsAsFactors = F))
  },cnvdat[["cnvid"]],cnvm,cnvdat[["sid"]],SIMPLIFY = F))
  return(gene.id)
}


annovar_bed <- function(bdf,dbdir="/media/nwknoblauch/Data/annovar/humandb/",tempdir="/tmp/",name){
  
  outfile <- file.path(tempdir,paste0(name,"_Results_EV"))
  outfilet <- file.path(tempdir,paste0(name,"_Results_T"))
  aexonf<- paste0(outfile,".exonic_variant_function")
  aconsf <- paste0(outfilet,".hg19_multianno.txt")
  
  if(!any(file.exists(c(aexonf,aconsf)))){

    infile <- file.path(tempdir,"Input.bed")
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
  consdf <- read.table(aconsf,header=T,sep="\t",stringsAsFactors = F,row.names=NULL)
  cn <- colnames(consdf)
  consdf <- consdf[,-ncol(consdf)]
  colnames(consdf) <- cn[-1]
  consdf <- mutate_each(consdf,funs(ifelse(.==".",NA,.)))
  consdf <- rename(consdf,chrom=Chr,start=Start,end=End,ref=Ref,alt=Alt)
  vardf <- read.table(aexonf,header=F,sep="\t",stringsAsFactors=F)
  vardf <- select(vardf,chrom=V4,start=V5,end=V6,ref=V7,alt=V8,consequence=V2)
  
  bdf <- left_join(bdf,vardf)
  bdf <- mutate(bdf,consequence=ifelse(is.na(consequence),"inter",consequence))
  bdf <- left_join(bdf,consdf)
  bdf <- select(bdf,-one_of("Otherinfo"))
  return(bdf)
}






BedC <- function(varfile,regionfile,regionheader=T,outdir="/tmp/",fieldname){
  
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

GCstarch <- function(starchfile,vardf,tempdir="/tmp/"){
  
  varf <- file.path(tempdir,"VARDAT.bed")
  svarf <- file.path(tempdir,"SVARDAT.bed")
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


ExacCov <- function(covdir,exonfile){
  starch.file <- file.path(covdir,"all_cov.starch")
  bedfiles <- file.path(covdir,paste0("file_",1:22,".bed"))
  gzfiles <- dir(covdir,pattern = ".gz",full.names = T)
  if(length(gzfiles)>0){
  sapply(gzfiles,function(x){
    cat(x)
    cat("\n")
    filechrom <- gsub(".+chr(.+).coverage.txt","\1",x)
    filedir <- gsub("(.+)Panel.chr.+.txt","\\1")
    outfile <- file.path(filedir,paste0("file_",filechrom,".bed"))
    tf <- read.table(x,header=F,sep="\t",stringsAsFactors = F)
    tf <- select(tf,chrom=V1,start=V2,mean=V3) %>% mutate(end=start+1)
    write.table(tf,outfile,col.names=F,row.names=F,quote=F,sep="\t")
  })
  }
  if(!file.exists(starch.file)){
  system(paste("sort-bed ",paste0(bedfiles,collapse=" "),"| starch - >",starch.file))
  }
  
  covmf <- file.path(covdir,"covmap.txt")
  if(!file.exists(covmf)){
  system(paste("bedmap --faster --mean --median --min --variance ",exonfile,starch.file,"> ",covmf))
  }
  exdf <- read.table(exonfile,header=F,sep="\t",stringsAsFactors = F)
  exdf <- rename(exdf,chrom=V1,start=V2,end=V3,gene=V4,syn_rate=V5,nsyn_rate=V6)
  covm <- read.table(covmf,sep="|",header=F,stringsAsFactors=F,colClasses = rep("numeric",4),na.strings = "NAN")
  covm <- rename(covm,mean=V1,median=V2,min=V3,variance=V4)
  exdf <- bind_cols(exdf,covm)
  return(exdf)
}


# tbf <- read.table(bedfiles[1],header=F,sep=" ",stringsAsFactors = F)
# ttbf <- filter(tbf,V2>=69090,V2<=70008)




  # bedsort <- system("which sort-bed",intern=T)
  # bedmap <- system("which bedmap",intern=T)
  # csf <- file.path(outdir,"tvarfile.txt")
  # crf <- file.path(outdir,"trfile.txt")
  # system(paste0(bedsort," ",varfile, " > ",csf))
  # if(regionheader){
  #   system(paste0("tail -n +2 ",regionfile," | ",bedsort," - > ",crf))    
  # }else{
  #   system(paste0("cat ",regionfile," | ",bedsort," - > ",crf))    
  # }
  # cnvm <- system(paste0(bedmap," --faster --count ",crf," ",csf),intern = T)
  # regdf <- read.table(crf,header=F,sep="\t",stringsAsFactors = F)
  # regdf <- rename(regdf,chrom=V1,start=V2,end=V3)
  # regdf[[fieldname]] <- cnvm
  # return(regdf)
# }


# compCovereage <- function(genedf,cnvdat,windowsize=10000,overlap=1,bedmap=NULL){
#   stopifnot(!is.null(genedf$chrom),
#             !is.null(genedf$start),
#             !is.null(genedf$end),
#             !is.null(genedf$name),
#             !is.null(cnvdat$sid),
#             !is.null(cnvdat$cnvstart),
#             !is.null(cnvdat$cnvstop),
#             !is.null(cnvdat$chrom),
#             !is.null(cnvdat$cnvid))
#   if(is.null(bedmap)){
#     bedmap <- system("which bedmap",intern = T)
#     bedsort <- system("which sort-bed",intern=T)
#     bedops <- system("which bedops",intern=T)
#     
#   }
#   cnvf <- file.path(tempdir,"CNVDAT.bed")
#   scnvf <- file.path(tempdir,"SCNVDAT.bed")
#   genef <- file.path(tempdir,"GENEDAT.bed")
#   sgenef <- file.path(tempdir,"SGENEDAT.bed")
#   cgenef <- file.path(tempdir,"CSGENEDAT.bed")
#   #ccnvf <- file.path(tempdir,"CCNVDAT.bed")
#   cnvdat <- select(cnvdat,chrom,cnvstart,cnvstop,cnvid,sid)
# 
#   genedf <- select(genedf,chrom,start,end,name)
#   write.table(cnvdat,
#               file=cnvf,
#               sep="\t",
#               col.names=F,
#               row.names=F,
#               quote=F)
#   write.table(genedf,
#               file=genef,
#               sep="\t",
#               col.names=F,
#               row.names=F,
#               quote=F)
#   system(paste0(bedsort," ",cnvf, " > ",scnvf))
#   system(paste0(bedsort," ",genef, " > ",sgenef))
#   system(paste0(bedops," -w ",windowsize," ",sgenef," > ",cgenef))
#   cnvm <- system(paste0(bedmap," --faster --echo-map-id ",cgenef," ",sgenef),intern = T)
#   
#   cnvc <- as.numeric(system(paste0(bedmap," --count --bp-ovr ",overlap," ",cgenef," ",scnvf),intern=T))
#   ctdf <- data_frame(name=cnvm)
#   posd <- read.table(cgenef,header=F,stringsAsFactors = F)
#   posd <- mutate(posd,ct=cnvc)
#   
#   ggplot(posd)+geom_segment(aes(x=V2,y=ct,xend=V3,yend=ct),size=9)  
# }





