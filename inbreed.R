library(rstan)
getpars <- c("maf","f","f_mu","f_s")
model <- stan_model("~/code/inbreed/flipperking_logit.stan")

cdat <- cullSNP(as.data.frame(dat[,-1,with=F]),SNPdat,mthresh=0.25,lthresh = 0.1)
fX <- cdat$X[cdat$SNP$Functional]
fSNP <- cdat$SNP[functional==T,]
fX <- fX[-SNPdelink(fSNP)]

inds <- which(!is.na(fX),arr.ind = T)
standat <- list(N=nrow(inds),n=nrow(fX),m=ncol(fX),
                ii=as.vector(inds[,1]),jj=as.vector(inds[,2]),X=fX[inds],
                mu0=10,s0=10)

fit <- sampling(model,data=standat,pars = getpars,
             iter=300,chains = 3,cores = 3)

nX <- cdat$X[!cdat$SNP$functional]
nSNP <- cdat$SNP[functional==F]
#nX <- nX[-SNPdelink(nSNP)]

inds <- which(!is.na(nX),arr.ind = T)
standat <- list(N=nrow(inds),n=nrow(nX),m=ncol(nX),
                ii=as.vector(inds[,1]),jj=as.vector(inds[,2]),X=nX[inds],
                mu0=10,s0=10)
nfit <- sampling(model,data=standat,pars = getpars,
                iter=150,chains = 3,cores = 3)

fXc <- fX[sapply(fX,mean,na.rm=T)>0.4]
inds <- which(!is.na(fXc),arr.ind = T)
standat <- list(N=nrow(inds),n=nrow(fXc),m=ncol(fXc),
                ii=as.vector(inds[,1]),jj=as.vector(inds[,2]),X=fXc[inds],
                mu0=10,s0=10)

fcfit <- sampling(model,data=standat,pars = getpars,
                 iter=300,chains = 3,cores = 3)


hwe <- function(x) {
  p <- mean(x,na.rm=T)/2
  2*p*(1-p)
}
qplot(sapply(fX,hwe),sapply(fX,function(x) mean(x==1,na.rm=T))) + geom_abline()
qplot(sapply(nX,hwe),sapply(nX,function(x) mean(x==1,na.rm=T))) + geom_abline()

fXc2 <- fX[abs(sapply(fX,hwe)-sapply(fX,function(x) mean(x==1,na.rm=T))) < .04]
inds <- which(!is.na(fXc2),arr.ind = T)
standat <- list(N=nrow(inds),n=nrow(fXc2),m=ncol(fXc2),
                ii=as.vector(inds[,1]),jj=as.vector(inds[,2]),X=fXc2[inds],
                mu0=10,s0=10)
fc2fit2 <- sampling(model,data=standat,pars = getpars,
                  iter=150,chains = 3,cores = 3)



##eigenanalysis
guh <- as.matrix(fXc2)
pcout <- pca(guh,"ppca",scale = "uv",center = T,nPcs = 50)

Xp <- guh
p <- colMeans(guh,na.rm=T)/2
r2 <- foreach (j = 1:100,.combine = cbind) %dopar% {
  library(pcaMethods)
  for (i in 1:ncol(guh)) Xp[,i] <- rbinom(392,2,p[i])#sample(guh[,i],replace = F)
  return(sDev(pca(as.matrix(Xp),"ppca",scale = "uv",center = T,nPcs = 50)))
}
plt <- rbind(data.table(r2=sDev(pcout),lb=sDev(pcout),ub=sDev(pcout),model="data"),
             data.table(r2=rowMeans(r2),lb=apply(r2,1,quantile,probs=0.025),ub=apply(r2,1,quantile,probs=0.975),model="null"))

ggplot(plt,aes(x=rep(1:50,2),y=r2,ymin=lb,ymax=ub,color=model,fill=model)) + geom_line() + geom_ribbon(alpha=0.5)

