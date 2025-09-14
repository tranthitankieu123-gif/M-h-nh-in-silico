#load library
library(BacArena)
library(sybil)
library(glpkAPI)
library(data.table)
library(ggplot2)
sybil::SYBIL_SETTINGS("SOLVER","cplexAPI")

pa <- readRDS("./../../gem-models/pa.RDS")
sa <- readRDS("./../../gem-models/sa.RDS")

#from data to bacterial
bacpa <- Bac(pa, type="Pseudomonas aeruginosa", chem = "EX_cpd00007_e0")
bacsa <- Bac(sa, type="MRSA")

#setting enviroment culture
arena = Arena(n=100, m=100, Lx=0.025, Ly=0.025)

#add bacteria
arena <- addOrg(arena, bacpa, amount=50)
arena <- addOrg(arena, bacsa,  amount=50)

arena_subs <- fread("TSBmed.csv") #dỉrectory
arena_subs[, ex.rxn := paste0("EX_", compounds, "_e0")]
arena <- addSubs(arena, smax = arena_subs$maxFlux, mediac = arena_subs$ex.rxn, unit = "mM", addAnyway = T)
arena <- createGradient(arena,smax=200,mediac=c("EX_cpd00007_e0"), position='top',steep=0.5, add=TRUE)
chemotaxis(bacpa, arena, 1, "EX_cpd00007_e0", arena@occupyM)

#start simulation
sim <- simEnv(arena,time=10, sec_obj='mtf',)
save(sim, file="sim_oxy.RData")

#show result
plotCurves2(sim, legendpos = "topleft")
evalArena(sim)
plotSpecActivity(sim)
#end