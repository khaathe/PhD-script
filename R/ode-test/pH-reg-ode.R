library(deSolve)
library(plotly)

rm(list=ls())

odes <- function(t, state, parameters){
  with(as.list(state, parameters),{
    V_mct_in <- a2ch * V_max_mct * ( He/(Ve*Km_mct+He) )
    V_mct_out <- c2ah * V_max_mct * ( Hi/(Vi*Km_mct+Hi) )
    dHi.dt <- V_mct_in - V_mct_out
    dHe.dt <- V_mct_out - V_mct_in
    return( list( c(dHi.dt, dHe.dt) ) )
  })
}

pHi <- 7.2
pHe <- 7.2
Hi <- 10^(-pHi)
He <- 10^(-pHe)
rc <- 4.5
Vi <- 4/3 * pi * (rc^3)
Ve <- 10^12

Vmax_AcL <- 9.58e-5 # pg/s/um^2
MW_h <- 1 # g/mol
MW_AcL <- 90.1 # g/mol
Si <- 1
V_max_mct <- Vmax_AcL * (MW_h/MW_AcL) * Si

Km_AcL <-0.405e-3 # pg/um^3 
Km_mct <- Km_AcL * (MW_h/MW_AcL)

a2ch_slope <- 1.5 # dimensionless
a2ch_thr <- 7 # dimensionless
a2ch <- 2 - tanh( a2ch_slope * pHi - a2ch_thr)

c2ah_slope <- 1.5 # dimensionless
c2ah_thr <- 7 # dimensionless
c2ah <- 2 - tanh( c2ah_slope * pHe - c2ah_thr)

parameters <- list(V_max_mct = V_max_mct, Km_mct=Km_mct, Vi=Vi, Ve=Ve)
state <- c(Hi=Hi, He=He )
times <- seq(from = 0.0, to = 10, by = 0.01)

out <- ode(y = state, times = times, func = odes, parms =  parameters, method = 'euler')
out <- as.data.frame(out)

plot_ly(data = out, x=~times, y=~Hi, name = "Protons", type="scatter", mode="lines") %>%
  add_trace(y=~He, name="He")
