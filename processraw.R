agfun <- function(x) {
  if (length(x)==1) return(x)
  else {
    if (length(unique(x))==1) return(unique(x))
    else return(as.integer(NA))
  }
}
source('~/Dropbox/monkeybris/rscript/SNPpreproc.R')
source('~/Dropbox/monkeybris/rscript/setup.R')
d <- "~/Dropbox/Manuscripts_Cayo/GeneralSNPpaper/SNP_Data/SNPplateData/"
fname <- grep("RESULTS",dir(d),value = T)

dat <- list()
for (i in 1:length(fname))
{
  dat[[i]] <- fread(paste0(d,fname[i]),skip = 9)
  dat[[i]] <- dat[[i]][,.(SNP=`SNP Name`,ID=`Sample ID`,geno=(`Allele1 - AB`=="B") + (`Allele2 - AB`=="B"))]
}
dat <- as.data.table(plyr::ldply(dat))
# dat[ID=="Zl4 (L)",ID:="Z14"]
# dat[ID=="9.00E+05",ID:="9E5"]
# dat[ID=="9.00E+03",ID:="9E3"]
# dat[ID=="8.00E+02",ID:="8E2"]
# dat[ID=="7.00E+00",ID:="7E0"]
# dat[ID=="7.00E+01",ID:="7E1"]
# dat[ID=="7.00E+03",ID:="7E3"]
# dat[ID=="5.00E+03",ID:="5E3"]
# dat[ID=="2.00E+04",ID:="2E4"]
# dat[ID=="2.00E+08",ID:="2E8"]
# dat[ID=="1.00E+03",ID:="1E3"]
# dat[ID=="1.00E+00",ID:="1E0"]
# dat[ID=="1.00E+07",ID:="1E7"]
# dat[ID=="1.00E+06",ID:="1E6"]
# dat[ID=="O1K",ID:="01K"]
# dat[ID=="OD3",ID:="0D3"]
# # dat[grepl("0\\.00E\\+00",ID),ID:="0E0"]
# dat[grepl("1D5",ID),ID:="1D5"]
# dat[grepl("6C0",ID),ID:="6C0"]
# dat[grepl("O3D",ID),ID:="03D"]
# dat[,ID:=toupper(ID)]

dat[,geno := if (mean(geno,na.rm = T)>1) as.integer(geno*-1+2) else geno,by=SNP]
#dat[SNP=="HTR2A1425",SNP:="HTR2A-1425"]
#dat[SNP=="SLC6A22495",SNP:="SLC6A2-2495"]
dat <- dcast(melt(dat),ID ~SNP,fun.aggregate = agfun)


readyparallel()
m <- foreach(i = 1:418,.combine=cbind) %dopar% {
  x <- unlist(dat[i,-1,with=F])[1:301]
  m <- vector()
  for (j in 1:397)
    m[j] <- cor(unlist(xSNP[j,]),x,use="pairwise") %>% abs()
  return(m)
}

dat <- dat[-((apply(m,2,max) < 0.85) %>% which())]
m <- m[,-((apply(m,2,max) < 0.85) %>% which())]
dat$ID <- rownames(xSNP)[apply(m,2,which.max)]

extra <- (table(dat$ID) %>% as.data.table())[N>1,V1]
for (id in extra)
{
  i <- which(dat$ID %in% id)[-1]
  dat <- dat[-i]
}
setkey(dat,"ID")

nfSNP <- rbind(fread("~/Dropbox/Manuscripts_Cayo/GeneralSNPpaper/SNP_Data/GGWG_SNPfiles032811/110328 final_Illumina_AIMS.csv",select=c(1,4,5,3)),
               fread("~/Dropbox/Manuscripts_Cayo/GeneralSNPpaper/SNP_Data/GGWG_SNPfiles032811/110328 final_Illumina_Paternity.csv",select=c(1,4,5,3)))
nfSNP$Functional <- F
nfSNP$Locus_Name <- str_replace_all(nfSNP$Locus_Name,"\\.","-")
#fSNP <- fread("~/Dropbox/Manuscripts_Cayo/GeneralSNPpaper/SNP_Data/IlluminaSubmittedFiles040111/HorvathCandidateGeneSNPs040111.csv")
fSNP <- fread("~/Dropbox/monkeybris/oldata/SNPs_with locations_full.csv",select=c(1,3,4,2))
fSNP$Functional <- T
SNPdat <- rbind(fSNP,nfSNP)
SNPdat <- SNPdat[Locus_Name %in% names(dat)]
SNPdat[,Coordinate:= str_replace_all(Coordinate,",","") %>% as.numeric()]
SNPdat[,Sequence:=toupper(Sequence)]
#setnames(SNPdat,1:3,names(fv$SNP)[1:3])
setcolorder(SNPdat,c(1,3,4,2,5))
SNPdat[,Chromosome:= str_replace(Chromosome,"chr","")]
setkey(SNPdat,"Locus_Name")
dupe <- SNPdat[,as.data.table(table(Locus_Name))][N>1,Locus_Name]
for (l in dupe) SNPdat <- SNPdat[-which(Locus_Name==l)[2]]
dat <- dat[,-which(!(names(dat) %in% SNPdat$Locus_Name))[-1],with=F]