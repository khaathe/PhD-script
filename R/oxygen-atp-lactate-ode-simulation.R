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

times <- seq(0, 100, by = 0.01)
O <- 0.052 # 0.05*760*1.39e-3
G <- 5 # 4e-11
ATP <- 1
H <- 7.4
parameters <- list(Vo=0.012, Ko = 0.005, pg=1, A0 = (29/5)*0.012, Kg = 0.04, Kh=0.00025)
state <- c(O=O, G=G, ATP=ATP, H=H)

out <- ode(y = state, times = times, func = odes, parms =  parameters, method = 'rk4')

plot_ly(data = as.data.frame(out), x=~time, y=~O, name = "Oxygen (mmol/L)", type="scatter", mode='lines') %>%
  add_trace(y=~G, name = "Glucose (mmol/L)") %>%
  add_trace(y=~ATP, name = "ATP (mmol/L)") %>%
  add_trace(y=~H, name = "pH" ) %>%
  layout(title=paste0(O, " Oxygen (mmol/L), ", G, " Glucose (mmol/L)"),
         xaxis = list(title="Time (min)"),
         yaxis = list(title="Concentration")
  )
