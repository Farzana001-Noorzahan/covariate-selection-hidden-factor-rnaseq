library(csrnaseq)
FixCov <- csrnaseq::FixCov
VarCov <- csrnaseq::VarCov
counts <- csrnaseq::counts


filtered_id <- apply(counts, 1, function(x) length(x[x>=5]) >= 2)
SimCnt_filtered <-counts[filtered_id,]

lib.size <- apply(SimCnt_filtered, 2, quantile, .75)
y <- t(log2(t(SimCnt_filtered + 0.5)/(lib.size + 1) * 1e+06))

SimCntOut<-FSRAnalysisBS(
counts,
FixCov,
VarCov,
option = "OWN",
B = 100,
m = 3,
lambdamax = 10,
alpha0 = 0.05,
ncores =1,
print.progress= FALSE,
saveall = TRUE
)

BestCovOut <- VarCov[names(SimCntOut$BestCovOut$BestRE)]


# FSR and then sva--------
if(ncol(BestCovOut)>0){
  mod0 <- model.matrix(~., data = BestCovOut)
  mod <- model.matrix(~., data = cbind(FixCov, BestCovOut))
}else{
  mod0 <- model.matrix(~1, data = BestCovOut)
  mod <- model.matrix(~., data = FixCov)
}
FSRsvaout <- sva::sva(dat = y, mod = mod, mod0 = mod0)
FSRsvacov <- data.frame(FSRsvaout$sv)

if(ncol(FSRsvacov) >=1){names(FSRsvacov) <- paste0("sva", 1:ncol(FSRsvacov))

svacov_on_FSR <- cbind(BestCovOut, FSRsvacov)} else{
  svacov_on_FSR<-BestCovOut
  
}


# sva on nothing-----
mod1 <- model.matrix(~., data = FixCov)
svaout0 <- sva::sva(dat = y, mod = mod1)
svacov0 <- data.frame(svaout0$sv)
if(ncol(svacov0) >=1){
  names(svacov0) <- paste0("sva", 1:ncol(svacov0))
  FixCovOut_sva_FSR <- cbind(FixCov,svacov0)
  
}

# sva on everything then FSR
mod0 <- model.matrix(~., data =  VarCov)
mod <- model.matrix(~., data = cbind(FixCov, VarCov))
# svaoutall <- sva::svaseq(dat = SimCnt_filtered, mod = mod, mod0 = mod0)
svaoutall <- sva::sva(dat = y, mod = mod, mod0 = mod0)
svacovall <- data.frame(svaoutall$sv)
if(ncol(svacovall) >=1){
  names(svacovall) <- paste0("sva", 1:ncol(svacovall))
  AllCovOut_sva_FSR <- cbind(FixCov, VarCov, svacovall)
  
}
cbind(VarCov,svacovall)

SVAall_FSRvar<-csrnaseq::FSRAnalysisBS(
  counts,
  FixCov,
  VarCov=cbind(VarCov,svacovall),
  option = "OWN",
  B = 100,
  m = 3,
  lambdamax = 10,
  alpha0 = 0.05,
  ncores = 1,
  print.progress = FALSE,
  saveall = TRUE
)

Allsvaallcov<-cbind(VarCov,svacovall)

BestCovOutFSR_all <- Allsvaallcov[names(SVAall_FSRvar$BestCovOut$BestRE)]


#### DE including only Line
lineonly <- csrnaseq:::VoomPv(counts = csrnaseq::counts, AllCov = csrnaseq::FixCov)
deline <- sum(csrnaseq:::jabes.q(lineonly$pvs[,1])<= 0.05)


##### DE Genes including sva and line
linesva <- csrnaseq:::VoomPv(counts = csrnaseq::counts, AllCov = cbind(csrnaseq::FixCov,svacov0))
delinesva <- sum(csrnaseq:::jabes.q(linesva$pvs[,1])<= 0.05)


#### DE including only the bestcov
Bestcovonly <- csrnaseq:::VoomPv(counts = csrnaseq::counts, AllCov = cbind(csrnaseq::FixCov ,BestCovOut))
debestcov <- sum(csrnaseq:::jabes.q(Bestcovonly$pvs[,1])<= 0.05)


#### DE including only the bestcov and sva
Bestcovsav <- csrnaseq:::VoomPv(counts = csrnaseq::counts, AllCov = cbind(csrnaseq::FixCov ,BestCovOut,FSRsvacov))
debestcovsva <- sum(csrnaseq:::jabes.q(Bestcovsav$pvs[,1])<= 0.05)


### DE genes including all variables
Allcov <- csrnaseq:::VoomPv(counts = csrnaseq::counts, AllCov = cbind(csrnaseq::FixCov,VarCov))
deallcov <- sum(csrnaseq:::jabes.q(Allcov$pvs[,1])<= 0.05)

### DE genes including all variables and estimated surrogate variables
Allcovsva <- csrnaseq:::VoomPv(counts = csrnaseq::counts, AllCov = cbind(csrnaseq::FixCov,BestCovOutFSR_all ))
deallsva <- sum(csrnaseq:::jabes.q(Allcovsva$pvs[,1])<= 0.05)

