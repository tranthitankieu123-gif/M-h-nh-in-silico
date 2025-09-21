# Script to reproduce P. aeruginose biofilm simulation
# File S1 is needed (contains metabolic model and diet)

library(BacArena)
library(parallel)

medium = read.csv('./poa_minimal_medium.csv')
load("./poa2.RData")

ex <- findExchReact(poa2) # change exchange names (only for better readability)
dict <- poa2@met_name[ex@met_pos]; names(dict) <- ex@react_id
poa2@react_id[ex@react_pos] <- paste0("EX_", poa2@met_name[ex@met_pos])
min <- poa2@react_id[ex@react_pos[which(ex@react_id!=0)]]

replicates <- 10; timesteps=192 # Attention: will use 10 cores!, timesteps in 0.5h
simlist <- list()

cores <- ifelse(replicates>detectCores(), detectCores(), replicates)
cl <- makeCluster(cores, type="PSOCK")
clusterExport(cl, c("poa2", "medium", "timesteps"))

print(system.time(simlist <- parLapply(cl, 1:replicates, function(i){
  print(sybil::SYBIL_SETTINGS()$SOLVER)
  bac = BacArena::Bac(model=poa2, setExInf=TRUE)
  bac = BacArena::setKinetics(bac, exchangeR="EX_D-Glucose[e]", Km=0.01, vmax=7.56) # Gosset2005
  arena = BacArena::Arena(tstep=0.5)
  start_bacs = 900; start_bacs_tmp = sqrt(start_bacs)
  middle = arena@n/2 - expand.grid(1:start_bacs_tmp,1:start_bacs_tmp) + start_bacs_tmp/2
  arena = BacArena::addOrg(arena, bac, amount=start_bacs, x=middle[1], y=middle[2])
  arena <- BacArena::addSubs(arena, smax=medium$Concentration, difspeed=medium$diffusion.constant, unit='mM', 
          mediac=paste0("EX_",medium$Substance,"[e]"))
  arena <- BacArena::addSubs(arena, smax=0, mediac="EX_CO2[e]", add=FALSE)
  BacArena::simEnv(arena, time=timesteps, sec_obj = "mtf")
}) ))
stopCluster(cl)

