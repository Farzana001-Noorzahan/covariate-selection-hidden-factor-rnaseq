## This function run all simulation for total 12 scenarios of
## ModelSize (6: 0, 1, 2, 6, 7, 8)
## alpha0 (1:.05)
## m (1: 3)
## ideal (2: TRUE, FALSE)


#rm(list = ls())
#pm1 <- proc.time()
#library("csrnaseq")
#env <-  as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
#print(env)
#nGene <- 2000
#nSim <- 1:500
#B <- 100
#lambdamax <- 10
#ncores <- 20
#ModelSizevec <- c( 1, 2, 3, 4, 5, 6, 7, 8)
#alpha0vec <- 0.05
#mvec <- 3
#idealvec <- c(TRUE, FALSE)
#Scenarios <- expand.grid(alpha0 = alpha0vec, m = mvec, ModelSize = ModelSizevec, ideal = idealvec)
#print.progress <- F
#saveall <- FALSE
#savesim <- TRUE

#RealDataOut <- readRDS(file = paste0("/home/fnoor001/cs-hidden-factor-rnaseq/cs-hidden-factor-rnaseq-main/analysis/RealDataOutBS/ModelSize_",Scenarios$ModelSize[env], ".rds" ))
#alpha0 <- Scenarios$alpha0[env]; m <- Scenarios$m[env]; ideal <- Scenarios$ideal[env]
#FSRSimOut_sva_ruv <- csrnaseq:::FSRnSimBS_sva_ruv(RealDataOut, nGene, nSim, B, m,
#                       lambdamax, alpha0, ncores, print.progress,
#                       saveall, savesim, ideal)
#proc.time() -pm1


#####
# Get the SLURM array task ID
env <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (is.na(env)) {
  stop("This script must be run as a SLURM array job")
}

# Start timing
pm1 <- proc.time()

# Load required library
library("csrnaseq")

# Print the scenario index
print(paste("Running scenario", env))

# Simulation parameters
nGene <- 2000
nSim <- 1:100
B <- 100
lambdamax <- 10
ncores <- 20 

# Experimental design
ModelSizevec <- c(1, 2, 3, 4, 5, 6, 7, 8)
alpha0vec <- 0.05
mvec <- 3
idealvec <- c(TRUE, FALSE)
Scenarios <- expand.grid(alpha0 = alpha0vec, m = mvec, ModelSize = ModelSizevec, ideal = idealvec)

# Additional parameters
print.progress <- FALSE
saveall <- FALSE
savesim <- TRUE

# Verify scenario exists
if (env > nrow(Scenarios)) {
  stop(paste("Invalid scenario number", env, ". Max scenario is", nrow(Scenarios)))
}

# Read real data
RealDataOut <- readRDS(file = paste0("RealDataOutBS/ModelSize_", Scenarios$ModelSize[env], ".rds"))

# Set parameters for this scenario
alpha0 <- Scenarios$alpha0[env]
m <- Scenarios$m[env]
ideal <- Scenarios$ideal[env]

# Create output directory if it doesn't exist
dir.create("results", showWarnings = FALSE)

#FSRnSimBS_sva
FSR1SimBS_sva <- function(RealDataOut, nGene, nrep, B, m, lambdamax, alpha0, ncores, print.progress, saveall = FALSE, savesim = TRUE, ideal = TRUE){
  FixCov <- RealDataOut$FixCov # the main factor of interest
  VarCov0 <- RealDataOut$VarCov # all available known covariates
  if(ideal){
    VarCov <- VarCov0
  }else{
    VarCov <- VarCov0[, !names(VarCov0)%in%c("neut", "lymp", "mono", "eosi",  "baso")]
  }
  TrueSelCov <- names(RealDataOut$BestCovOut$BestER)
  SimCnt <- csrnaseq:::SimCounts(RealDataOut, nGene,nrep)
  saveRDS(SimCnt, "SimCnterr.rds")
  filtered_id <- apply(SimCnt$SimCnt, 1, function(x) length(x[x>=5]) >= 2)
  SimCnt_filtered <-SimCnt$SimCnt[filtered_id,]
  
  lib.size <- apply(SimCnt_filtered, 2, quantile, .75)
  y <- t(log2(t(SimCnt_filtered + 0.5)/(lib.size + 1) * 1e+06))
  
  option <- "OWN"
  SimCntOut <- csrnaseq::FSRAnalysisBS(#counts = SimCnt$SimCnt,
    counts = SimCnt_filtered,
    FixCov, VarCov,
    option, B, m,lambdamax,
    alpha0, ncores,print.progress, saveall)
  
  # Calculate the empirical FSR, S, U, R for FSRBS
  
  SimSelCov.RE <- names(SimCntOut$BestCovOut$BestRE)
  S.RE <-  length(intersect(TrueSelCov, SimSelCov.RE))
  R.RE <-  length(SimSelCov.RE)
  U.RE <- R.RE - S.RE
  FSP.RE <- U.RE/(R.RE+1) # i
  FSROut <- c(S.RE, R.RE, U.RE, FSP.RE)
  names(FSROut) <- paste0("FSR.", c("S", "R", "U", "FSP"), ".OWN.RE")
  
  # DEA result for 1Sim
  lab <-as.factor(as.numeric(!(SimCnt$IndSample[filtered_id]%in%SimCnt$EEGene$Line)))
  # lab <-as.factor(as.numeric(!(SimCnt$IndSample%in%SimCnt$EEGene$Line)))
  # OnlyLine
  OnlyLine <- csrnaseq:::DEAeval(p = csrnaseq:::VoomPv(counts = SimCnt_filtered, AllCov = FixCov)$pvs[,"Line"], lab = lab,  method = "OnlyLine")
  # AllCov
  All <- csrnaseq:::DEAeval(p = csrnaseq:::VoomPv(counts = SimCnt_filtered, AllCov = cbind(FixCov, VarCov))$pvs[,"Line"], lab = lab,  method = "All")
  #Oracle
  TrueCov <- VarCov0[TrueSelCov]
  Oracle <- csrnaseq:::DEAeval(p = csrnaseq:::VoomPv(counts = SimCnt_filtered, AllCov = cbind(FixCov,TrueCov))$pvs[,"Line"], lab = lab,  method = "Oracle")
  
  # FSR only-------
  BestCovOut <- VarCov[names(SimCntOut$BestCovOut$BestRE)]
  OnlyFSR<- csrnaseq:::DEAeval(p = csrnaseq:::VoomPv(counts = SimCnt_filtered,
                                                     AllCov = cbind(FixCov,BestCovOut))$pvs[,"Line"],
                               lab = lab,  method = "OnlyFSR")
  # FSRDEA <- csrnaseq:::DEAeval(p = csrnaseq:::VoomPv(counts = SimCnt$SimCnt,
  #                                                    AllCov = cbind(FixCov,BestCovOut))$pvs[,"Line"],
  #                              lab = lab,  method = "OWN.RE")
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
  if(ncol(FSRsvacov) >=1)names(FSRsvacov) <- paste0("sva", 1:ncol(FSRsvacov))
  svacov_on_FSR <- cbind(BestCovOut, FSRsvacov)
  FSRsva <- csrnaseq:::DEAeval(p = csrnaseq:::VoomPv(counts = SimCnt_filtered,
                                                     AllCov = cbind(FixCov,svacov_on_FSR))$pvs[,"Line"],
                               lab = lab,  method = "FSRsva")
  
  
  FSRsvacov
  
  
  # sva on nothing-----
  
  mod <- model.matrix(~., data = FixCov)
  svaout0 <- sva::sva(dat = y, mod = mod)
  svacov0 <- data.frame(svaout0$sv)
  if(ncol(svacov0) >=1)names(svacov0) <- paste0("sva", 1:ncol(svacov0))
  SVA0 <- csrnaseq:::DEAeval(p = csrnaseq:::VoomPv(counts = SimCnt_filtered,
                                                   AllCov = cbind(FixCov,svacov0))$pvs[,"Line"],
                             lab = lab,  method = "SVA0")
 svacov0 

  # sva on everything then FSR-----
  
  mod0 <- model.matrix(~., data =  VarCov)
  mod <- model.matrix(~., data = cbind(FixCov, VarCov))
  # svaoutall <- sva::svaseq(dat = SimCnt_filtered, mod = mod, mod0 = mod0)
  svaoutall <- sva::sva(dat = y, mod = mod, mod0 = mod0)
  svacovall <- data.frame(svaoutall$sv)
  if(ncol(svacovall) >=1){
    names(svacovall) <- paste0("sva", 1:ncol(svacovall))
    if(m + ncol(svacovall) + ncol(mod) >= nrow(mod)) m <- max(1, nrow(mod) - ncol(mod) - ncol(svacovall) - 1)
    SimCntOut_sva_FSR <- csrnaseq::FSRAnalysisBS(#counts = SimCnt$SimCnt,
      counts = SimCnt_filtered,
      FixCov, cbind(VarCov, svacovall),
      option, B, m,lambdamax,
      alpha0, ncores,print.progress, saveall)
    
    BestCovOut_sva_FSR <- cbind(VarCov, svacovall)[names(SimCntOut_sva_FSR$BestCovOut$BestRE)]
    SVAall_FSR <- csrnaseq:::DEAeval(p = csrnaseq:::VoomPv(counts = SimCnt_filtered,
                                                           AllCov = cbind(FixCov,BestCovOut_sva_FSR))$pvs[,"Line"],
                                     lab = lab,  method = "SVAall_FSR")
                                     
    
  }else{ # SVAall_FSR <-OnlyFSR
    BestCovOut_sva_FSR <- VarCov[names(SimCntOut$BestCovOut$BestRE)]
    SVAall_FSR<- csrnaseq:::DEAeval(p = csrnaseq:::VoomPv(counts = SimCnt_filtered,
                                                          AllCov = cbind(FixCov,BestCovOut_sva_FSR))$pvs[,"Line"],
                                    lab = lab,  method = "SVAall_FSR")
  }
   svacovall 
   

# Calculate the empirical FSR, S, U, R for FSRBS

  SimSelCov.SVAallFSR <- names(BestCovOut_sva_FSR)
  SimSelCov.SVAallFSR
  S.SVAallFSR <-  length(intersect(TrueSelCov, SimSelCov.SVAallFSR))
  R.SVAallFSR <-  length(intersect(names(VarCov), SimSelCov.SVAallFSR))  #max(length(SimSelCov.SVAallFSR) - ncol(svacovall), 0)
  U.SVAallFSR <- R.SVAallFSR - S.SVAallFSR
  FSP.SVAallFSR <- U.SVAallFSR/(R.SVAallFSR+1) # i
  FSROut.SVAallFSR <- c(S.SVAallFSR, R.SVAallFSR, U.SVAallFSR, FSP.SVAallFSR)
  names(FSROut.SVAallFSR) <- paste0("FSR.", c("S", "R", "U", "FSP"), ".OWN.RE")
   
   
 
  DEAevalOut <- c(OnlyFSR, FSRsva, SVA0, SVAall_FSR,
                  OnlyLine, All, Oracle)
  
  SimOut <- list(SimCnt = SimCnt,
                 # SimCntOut = SimCntOut,
                 TrueSelCov = TrueSelCov,
                 SimSelCov.RE = SimSelCov.RE,
                 FSROut = FSROut,
                 FSROut.SVAallFSR = FSROut.SVAallFSR,
                 DEAevalOut = DEAevalOut,
                 FixCov = FixCov,
                 VarCov0 = VarCov0,
                 VarCov = VarCov,
                 FSRsvacov = FSRsvacov,
                 svacov0 = svacov0,
                 svacovall = svacovall
                 )
  
  ModelSize <- length(RealDataOut$BestCovOut$BestER)
  if(savesim){
    SimPath <- paste0("SimulationOutsva/ModelSize_", ModelSize, "_nGene_", nGene, "_B_", B,  "_alpha0_", alpha0, "_ideal_", ideal)
    dir.create(path = SimPath, showWarnings = F, recursive = T)
    saveRDS(SimOut, file = paste0(SimPath, "/nrep_", nrep, ".rds"))
  }
  
  
  
  
  SimSelCov <- list(TrueSelCov = TrueSelCov,
                    SimSelCov.RE = SimSelCov.RE,
                    FSROut = FSROut,
                    FSROut.SVAallFSR = FSROut.SVAallFSR,
                    DEAevalOut = DEAevalOut)
  SimSelCov
}

FSRnSimBS_sva<- function(RealDataOut, nGene, nSim, B, m, lambdamax, alpha0, ncores, print.progress,saveall = FALSE, savesim = TRUE, ideal = TRUE){
  nSimSelCov <- plyr::ldply(nSim, function(nrep){
    FSR1SimBSOut <- FSR1SimBS_sva(RealDataOut, nGene, nrep, B, m, lambdamax, alpha0, ncores, print.progress, saveall, savesim, ideal)
    c(FSR1SimBSOut$FSROut, FSR1SimBSOut$DEAevalOut)
  })
  res <- list(TrueSelCov = names(RealDataOut$BestCovOut$BestER),
              nSimSelCov = nSimSelCov)
  ModelSize <- length(RealDataOut$BestCovOut$BestER)
  SimPath <- paste0("SimulationOutsva/ModelSize_", ModelSize, "_nGene_", nGene, "_B_", B, "_m_", m,  "_alpha0_", alpha0, "_ideal_", ideal)
  dir.create(path = SimPath, showWarnings = F, recursive = T)
  saveRDS(res, file = paste0(SimPath, "/nSim_", min(nSim), "_", max(nSim), ".rds"))
  saveRDS(res, file = paste0(SimPath, "_nSim_", min(nSim), "_", max(nSim), ".rds"))
  res
}


# Run simulation
FSRSimOut_sva<- FSRnSimBS_sva(RealDataOut, nGene, nSim, B, m, lambdamax, alpha0, ncores, print.progress, saveall, savesim, ideal)

FSRSimOut_sva

means <- colMeans(FSRSimOut_sva$nSimSelCov)
sds <- apply(FSRSimOut_sva$nSimSelCov, 2, sd)

# Save results with scenario number in filename
saveRDS(FSRSimOut_sva, file = paste0("results/FSRSimOut_sva_scenario_", env, ".rds"))

# Calculate and print execution time
runtime <- proc.time() - pm1
print(runtime)

# Save runtime information
write.csv(data.frame(scenario = env, 
                    user_time = runtime[1], 
                    system_time = runtime[2], 
                    elapsed_time = runtime[3]),
          file = paste0("results/runtime_scenario_", env, ".csv"))