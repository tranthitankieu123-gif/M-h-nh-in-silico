#install package BacArena, sybil, glpkAPI, data.table

#load library
library(BacArena)
library(sybil)
library(glpkAPI)
library(data.table)
library(ggplot2)
library(sybilSBML)
sybil::SYBIL_SETTINGS("SOLVER","cplexAPI")

#input genome metabolic model
pa <- readSBMLmod("GCF_000009645.1_ASM964v1_genomic.xml") #directory file in your computer
sa <- readSBMLmod("GCF_000006765.1_ASM676v1_genomic.xml")

#from data to bacterial
bacpa <- Bac(pa, type="Pseudomonas aeruginosa", chem = "EX_cpd00007_e0")
bacsa <- Bac(sa, type="MRSA")

#setting enviroment culture
arena = Arena(n=100, m=100, Lx=0.025, Ly=0.025)

#add bacteria
arena <- addOrg(arena, bacpa, amount=10, x = c(51, 45, 60, 38, 50, 53, 42, 59, 48, 63), y = c(48, 55, 62, 41, 50, 39, 60, 44, 52, 57))
arena <- addOrg(arena, bacsa,  amount=10, x = c(40, 56, 49, 61, 44, 57, 36, 54, 47, 64), y = c(40, 58, 63, 37, 46, 51, 59, 42, 49, 53))

arena_subs <- fread("TSBmed.csv") #dỉrectory
arena_subs[, ex.rxn := paste0("EX_", compounds, "_e0")]
arena <- addSubs(arena, smax = arena_subs$maxFlux, mediac = arena_subs$ex.rxn, unit = "mM", addAnyway = T)
arena <- addSubs(arena, smax = 200, mediac = "EX_cpd00007_e0", unit = "mM", addAnyway = T)
arena <- createGradient(arena,smax=200,mediac=c("EX_cpd00007_e0"), position='top',steep=0.5, add=TRUE)
chemotaxis(bacpa, arena, 10, "EX_cpd00007_e0", arena@occupyM)


sim <- simEnv(arena,time=20, sec_obj='mtf')
save(sim, file="sim_oxy_v3_chemotaxis.RData")

#show result
pdf("sim_oxy_v3_chemotaxis.pdf")
plotCurves2(sim, legendpos = "topleft")
evalArena(sim)
plotSpecActivity(sim)