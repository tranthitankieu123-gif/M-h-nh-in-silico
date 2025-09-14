#install package BacArena, sybil, glpkAPI, data.table

#load library
library(BacArena)
library(sybil)
library(glpkAPI)
library(data.table)
library(ggplot2)
sybil::SYBIL_SETTINGS("SOLVER","cplexAPI")

#input genome metabolic model
pa <- readRDS("./../../gem-models/pa.RDS") #directory file in your computer
sa <- readRDS("./../../gem-models/sa.RDS") #directory file in your computer

#from data to bacterial
bacpa <- Bac(pa, type="Pseudomonas aeruginosa", chem = "EX_cpd00007_e0")
bacsa <- Bac(sa, type="MRSA")

#setting enviroment culture
arena = Arena(n=100, m=100, Lx=0.025, Ly=0.025)

#add bacteria
arena <- addOrg(arena, bacsa, amount=2)

arena_subs <- fread("TSBmed.csv") #dỉrectory
arena_subs[, ex.rxn := paste0("EX_", compounds, "_e0")]
arena <- addSubs(arena, smax = arena_subs$maxFlux, 
                 mediac = arena_subs$ex.rxn, unit = "mM", addAnyway = T)
arena <- createGradient(arena,smax=50,mediac=c("EX_cpd00007_e0"),
                    position='top',steep=0.5, add=TRUE)
chemotaxis(bacpa,arena,1, "EX_cpd00007_e0", arena@occupyM)

#start simulation
sim <- simEnv(arena, time=24, sec_obj='mtf')

arena <- addOrg(arena, bacpa, amount=2)

sim2 <- simEnv(arena, time=24, (sec_obj='mtf')

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
save(sim, file=paste0("sim_after_24h_v2", timestamp,".RData"))
save(sim2, file=paste0("sim2_after_24h_v2", timestamp,".RData"))
#show result
pdf(paste0("MSARFirstSA24HourLater_first",timestamp,".pdf")
evalArena(sim, legend_pos = "bottom")
plotCurves2(sim, legendpos = "topleft")

pdf(paste0("MSARFirstSA24HourLater",timestamp,".pdf")
evalArena(sim2, legend_pos = "bottom")
plotCurves2(sim2, legendpos = "topleft")