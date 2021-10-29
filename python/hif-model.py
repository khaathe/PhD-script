#!/usr/bin

import numpy as np
from numpy.typing import _128Bit
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def H(y,s,n,gamma):
    return s**n / ( s**n + y**n ) + gamma * y**n / (s**n + y**n)

def model(t, x, A, D, N, hif_coeff_prod, S, gamma, Vo, Ko, pg, A0,  Kg, Kh):
    dxdt = [None] * 8
    hif, ldh, pdk, pdh, o, g, atp, h = x 
    Co = Vo * ( o/(o+Ko) )
    Cg = ( (pg*A0/2.0) - (27.0*Co/10.0) ) * ( g/(g+Kg) )
    Pa = ( (2.0 * Cg) + ( (27.0/5.0)*Co ) )
    Ph = Kh*( (29.0*(pg*Vo-Co))/5.0 )
    dxdt[0] = A * hif_coeff_prod - D * H(o, S[0][1], N, gamma[0][1]) * hif
    dxdt[1] = A * H(hif, S[3][4], N, gamma[3][4]) - D * ldh
    dxdt[2] = A * H(hif, S[3][4], N, gamma[3][4]) - D * pdk
    dxdt[3] = A * H(pdk, S[3][4], N, gamma[3][4]) - D * pdh
    dxdt[4] = -Co
    dxdt[5] = -Cg
    dxdt[6] = Pa
    dxdt[7] = Ph
    return dxdt

A = 0.005
D = 0.005
N = 4.0
hif_coeff_prod = 140.0
Vo = 2.1e-11
Ko = 0.005
pg = 1.0
A0 = 29.0/5.0 * Vo
Kg = 0.04
Kh = 2.5e-4
S = [ 
    [0.0,20.85e-3,0.0,0.0,0.0],
    [0.0,0.0,4.64,4.0,0.0],
    [0.0,0.0,0.0,0.0,0.0],
    [0.0,0.0,0.0,0.0,2.2],
    [0.0,0.0,0.0,0.0,0.0]
]
gamma = [
    [0.0,40.0,0.0,0.0,0.0],
    [0.0,0.0,3.81,6.97,0.0],
    [0.0,0.0,0.0,0.0,0.0],
    [0.0,0.0,0.0,0.0,0.14],
    [0.0,0.0,0.0,0.0,0.0]
]
hypoxia = 0.010564
normoxia = 0.056
gluNormal = 5.0
gluLow = 1.0

x0 = [ hif_coeff_prod / gamma[0][1], gamma[1][2], gamma[1][3], gamma[3][4], normoxia, gluNormal, 0.0, 10**(-7.4)/1000]
tspan = [0.0, 480.0]
dt = 0.1
p = [A, D, N, hif_coeff_prod, S, gamma, Vo, Ko, pg, A0,  Kg, Kh]

sol = solve_ivp(model, tspan, x0, t_eval= np.arange(0.0, 480.0, 0.1), method='RK45', args=p, rtol = 1e-6, atol=1e-10)

labels = ["HIF", "LDH", "PDK", "PDH", "Oxygen", "Glucose", "ATP", "H+"]

for i in range(7):
    plt.plot(sol.t, sol.y[i], label=labels[i])
plt.legend()
plt.show()