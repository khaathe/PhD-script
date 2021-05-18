library(deSolve)
library(plotly)

rm(list=ls())

odes <- function(t, state, parameters){
  with(as.list(c(state, parameters, t=t)),{
    fo <- -Vo*(O/(O+Ko))
    fg <- -( (1/2)*(pg*A0) + (27/10)*fo ) * (G/(G+Kg))
    fa <- -(2*fg + ( (27/5)*fo ) )
    fh <- Kh*( (29/5)*(pg*Vo+fo) )
    return( list( c(fo, fg, fa, fh) ) )
  })
}

t.start <- 0
t.end <- 60 # s
dt <- 0.01
times <- seq(t.start, t.end, by = dt)
O <- 0.05*760*1.39e-3 # 0.05*760*1.39e-3 mmol/L
G <- 5 # mmol/L
ATP <- 0 # 
H <- 7.4 # pH
Vo <- 0.012 # 0.012 mmol/L/s
Ko <- 0.005 # mmol/L
A0 <- (29/5)*Vo # 
pg <- 1 # dimensionless
Kg <- 0.04 # mmol/L
Kh <- 0.00025 # dimensionless

parameters <- list(Vo=Vo, Ko = Ko, pg=pg, A0 = A0, Kg = Kg, Kh=Kh)
state <- c(O=O, G=G, ATP=ATP, H=H)

out <- ode(y = state, times = times, func = odes, parms =  parameters, method = 'rk4')

plot_ly(data = as.data.frame(out), x=~time, y=~O, name = "Oxygen (mmol/L)", type="scatter", mode='lines') %>%
  add_trace(y=~G, name = "Glucose (mmol/L)") %>%
  add_trace(y=~ATP, name = "ATP (mmol/L)") %>%
  add_trace(y=~H, name = "pH" ) %>%
  layout(title=paste0(O, " Oxygen (mmol/L), ", G, " Glucose (mmol/L)"),
         xaxis = list(title="Time (s)"),
         yaxis = list(title="Concentration")
  )
