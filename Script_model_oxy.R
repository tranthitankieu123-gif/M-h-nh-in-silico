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

arena_end <- getArena(sim, time=4)
df <- arena_end@orgdat

# Đổi tên type cho dễ nhìn nếu cần
df$type <- factor(df$type, 
                  levels = c("Pseudomonas aeruginosa", "MRSA"))

ggplot(df, aes(x = x, y = y, color = type, shape = type)) +
  geom_point(alpha = 0.7, size = 1.2) +
  scale_color_manual(values = c("Pseudomonas aeruginosa" = "#1f77b4", # xanh
                                "MRSA" = "#d62728")) +                # đỏ
  scale_shape_manual(values = c("Pseudomonas aeruginosa" = 4,
                                "MRSA" = 3)) +
  labs(x = "Không gian X", y = "Không gian Y", color = "Species", shape = "Species") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  ggtitle("Phân bố population sau 15 giờ mô phỏng")

#end