#install package BacArena, sybil, glpkAPI, data.table

#load library
library(BacArena)
library(sybil)
library(glpkAPI)
library(data.table)
library(ggplot2)
sybil::SYBIL_SETTINGS("SOLVER","cplexAPI")

#input genome metabolic model
pa <- readRDS("pa.RDS") #directory file in your computer
sa <- readRDS("sa.RDS") #directory file in your computer

#from data to bacterial
bacpa <- Bac(pa, type="Pseudomonas aeruginosa", chem = "EX_cpd00007_e0")
bacsa <- Bac(sa, type="MRSA")

#setting enviroment culture
arena = Arena(n=100, m=100, Lx=0.025, Ly=0.025)

#add bacteria
arena <- addOrg(arena, bacsa, amount=60)

arena_subs <- fread("TSBmed.csv") #dỉrectory
arena_subs[, ex.rxn := paste0("EX_", compounds, "_e0")]
arena <- addSubs(arena, smax = arena_subs$maxFlux, mediac = arena_subs$ex.rxn, unit = "mM", addAnyway = T)
chemotaxis(bacpa,arena,1, "EX_cpd00007_e0", arena@occupyM)

#start simulation
sim <- simEnv(arena, time=12, sec_obj='mtf')
arena <- addOrg(arena, bacpa, amount=10)
sim2 <- simEnv(arena, time=12, sec_obj='mtf')
save(sim, file="biomass_check_6_1_sim.RData")
save(sim2, file="biomass_check_6_1_sim2.RData")

#show result
pdf("biomass_check_6_1_sim.pdf")
evalArena(sim, legend_pos = "bottom")
plotCurves2(sim, legendpos = "topleft")
plotSpecActivity(sim)

pdf("biomass_check_6_2_sim.pdf")
evalArena(sim2, legend_pos = "bottom")
plotCurves2(sim2, legendpos = "topleft")
plotSpecActivity(sim2)