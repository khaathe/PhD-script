library(deSolve)
library(plotly)

rm(list=ls())
source('~/PhD/Code/R/function/model-gene-level.R')

t.end <- 7200
t <- seq(0, t.end, by = 6)
O2.cos <- 80 * cos( 0.001 * t) + 80
names(O2.cos) <- t

hypoxia.start.time <- 1000
hypoxia.end.time <- 5200
O2.loss <- 160/(hypoxia.end.time - hypoxia.start.time)
O2.decrease <- sapply(t, function (t) {
  O2 <- 160
  if(t>hypoxia.start.time){
    O2 <- O2 - (t-hypoxia.start.time)  * O2.loss
  }
  if (t>hypoxia.end.time){
    O2 <- 0
  }
  O2
})
names(O2.decrease) <- t

# Default Parameters
# c.h <- 2.0 * 3.0 * 0.4 * 0.36 * 1.5 * 10.94 * 10
# c.ldh <-  2.6
# c.pdk <- 0.8 * 29.6
# c.pdh <- 0.09

# Custom Parameters
c.h <- 2.0 * 3.0 * 0.4 * 0.36 * 1.5 * 10.94 * 10
c.ldh <-  1.0
c.pdk <- 1.0
c.pdh <- 1.0

parameters <- list(A=0.005, D=0.005, s=0.02*760, n=4, gamma.h=40, pO2=O2.cos, c.h = c.h, c.ldh = c.ldh, c.pdk = c.pdk, c.pdh = c.pdh)
y0 <- c(H = 0.0, LDH = 0.0, PDK = 0.0, PDH = 0.0)

dt <- data.frame(T=t)
gamma.values <- c(50, 40, 30, 20, 10, 5,1)
nb.col <- length(gamma.values) + 2
dt[1:length(t), 2:nb.col] <- 0
colnames(dt)[2:nb.col] <- c("O2", gamma.values)
dt$O2 <- parameters[["pO2"]]

for ( i in gamma.values){
  parameters[["gamma.h"]] <- i
  out <- ode(y = y0, times = t, func = hif.model.odes, parms = parameters, method = 'euler')
  dt[1:length(t), as.character(i)] <- out[1:length(t), 2]
}

plot.hif.reg <- plot_ly(dt, x= ~T, y=~ O2, name="O2 pressure (mmHg)", type='scatter', mode='lines')
for (i in gamma.values){
  plot.hif.reg <- plot.hif.reg %>% add_trace(y= dt[1:length(t), as.character(i)], name = paste0("γ = ", i) )
}
plot.hif.reg <- plot.hif.reg %>% layout(title="Evolution of HIF1 level according to the O2 pressure (initail HIF1 level = 0.0 ) ",
                        xaxis = list(title="Time (min)"),
                        yaxis = list(title="Level (dimensionless)")
                        )
plot.hif.reg

parameters[["gamma.h"]] <- 40
ode.res <- ode(y = y0, times = t, func = hif.model.odes, parms = parameters, method = 'euler')
ode.res <- as.data.frame(ode.res)
ode.res$O2 <- parameters[["pO2"]]
plot.reg <- plot_ly(data = ode.res, x= ~time, y= ~H, name = "HIF1", type="scatter", mode='lines') %>% 
  add_trace(y= ~LDH, name = "LDH") %>% 
  add_trace(y= ~PDK, name = "PDK") %>% 
  add_trace(y= ~PDH, name = "PDH") %>% 
  add_trace(y= ~O2, name = "O2") %>% 
  layout(title=paste0("Gamma O2 = ", parameters[["gamma.h"]]),
         xaxis = list(title="Time (min)"),
         yaxis = list(title="Level (dimensionless)")
  )
plot.reg
