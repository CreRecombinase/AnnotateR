res.double <- function(df){
  comcols <- apply(df,2,function(x)any(grepl(pattern = ",",x)))
  ccl <- apply(df[,comcols],2,function(x)unlist(strsplit(x,",")))
  dupc <- lengths(strsplit(df$ALT,","))
  interp("x=y",x=colnames(df)[1])
  rdf <- data.frame(unlist(mapply(rep,df[[colnames(df)[1]]],dupc)))
  colnames(rdf) <- colnames(df)[1]
  for(i in 2:ncol(df)){
    if(comcols[i]){
      rdf[[colnames(df)[i]]] <-ccl[,colnames(df)[i]]
    }else {
      rdf[[colnames(df)[i]]] <-unlist(mapply(rep,df[[colnames(df)[i]]],dupc))
    }
  }
  rdf <- as_data_frame(rdf)
  return(rdf)
}

VCFFunc <- function(vcffile,gzipped=T,fields,outfile){
  vcfpath <- system("which vcftools",intern=T)
  if(gzipped){
    fpath <- paste0("--gzvcf ",vcffile)
  }else{
    fpath <- paste0("--vcf ",vcffile)
  }
  infos <- paste0("--get-INFO ",fields,collapse = " ")
  outs <- paste0("--out ",outfile)
  system(paste(vcfpath,fpath,infos,outs))
}

info.bed <- function(infofile,VQLcut=-2.632){
  adat <- read.table(infofile,header=T,stringsAsFactors = F,sep="\t")
  ddf <- filter(adat,grepl(",",ALT))
  ndf <- res.double(ddf)
  nadat <- bind_rows(filter(adat,!grepl(",",ALT)),ndf)
  nadat <- filter(nadat,VQSLOD>VQLcut)
  nadat <- select(nadat,chrom=CHROM,start=POS,ref=REF,alt=ALT,ac_adj=AC_Adj,an_adj=AN_Adj)
  nadat <- mutate(nadat,end=as.integer(start+nchar(ref)-1),chrom=paste0("chr",chrom)) %>% select(chrom,
                                                                                                  start,end,ref,alt,ac_adj,an_adj)
  return(nadat)
}






