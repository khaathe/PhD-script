library(deSolve)
library(plotly)
hif.reg.model <- function(t, y0, parameters){
  with(as.list(c(y0, parameters)),{
    dO2 <- 0.0
    if (t>hypoxia.start.time && t < hypoxia.end.time) {
      dO2 <- -1* O2.time.loss 
    }
    h.O2 <- (s^n/(s^n+O2^n))+gamma*(O2^n/(s^n+O2^n))
    dH <- A * R - D * H * h.O2 
    list(c(dH, dO2))
  })
}

t.end <- 7200
hypoxia.start.time <- 2200
hypoxia.end.time <- 5200
O2.time.loss <- (0.21*760)/(hypoxia.end.time-hypoxia.start.time)
R <- 2.0 * 3.0 * 0.4 * 0.36 * 1.5 * 10.94 * 10
parameters <- c(A=0.005, D=0.005, s=0.02*760, n=4, gamma=40, hypoxia.start.time=hypoxia.start.time, O2.time.loss=O2.time.loss,
                hypoxia.end.time = hypoxia.end.time, R)
y0 <- c(H = 0.0, O2=160)
t <- seq(0, t.end, by = 6)

dt <- data.frame(T=t)
gamma.values <- c(100, 50, 40, 30, 20, 10, 5,1)
nb.col <- length(gamma.values) + 2
dt[1:length(t), 2:nb.col] <- 0
colnames(dt)[2:nb.col] <- c("O2", gamma.values)

for ( i in gamma.values){
  parameters[["gamma"]] <- i
  out <- ode(y = y0, times = t, func = hif.reg.model, parms = parameters)
  dt[1:length(t), as.character(i)] <- out[1:length(t), 2]
}
dt$O2 <- out[1:length(t), 3]

plot <- plot_ly(dt, x= ~T, y=~O2, name="O2 pressure (mmHg)", type='scatter', mode='lines')
for (i in gamma.values){
  plot <- plot %>% add_trace(y= dt[1:length(t), as.character(i)], name = paste0("γ = ", i) )
}
plot <- plot %>% layout(title="Evolution of HIF1 level according to the O2 pressure (initail HIF1 level = 0.0 ) ",
                        xaxis = list(title="Time (min)"),
                        yaxis = list(title="Level (dimensionless)")
                        )
plot

