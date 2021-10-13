#!/usr/bin/python
import math
import matplotlib
import numpy
from scipy.optimize import least_squares
from scipy.integrate import solve_ivp
import pandas as pd
import matplotlib.pyplot as plt

MW_H = 1.0                    
MW_CO2 = 44.0                 
MW_O2 = 32.0                  
MW_HCO3 = 61.0                
MW_AcL = 90.1                 
PM_CO2 = 3.2 * pow(10, 4)     
gAcL = 3.8 * pow(10,-4)       
q_O2 = 3.5 * pow(10,-5)       
k1 = 0.144                    
k2 = 1.9 * pow(10, 5)         
VMAXAcL = 9.58 * pow(10,-5)   
K_mAcL = 0.405 * pow(10,-3)   
a2cH_slope = 1.5             
a2cH_thr = 7.0               
c2aH_slope = 1.5              
c2aH_thr = 7.0                
VMAXNHE = 5.15 * pow(10,-7)   
K_mNHE = 0.196 * pow(10,-6)   
a = 2.67                      
l_NHE = 0.076                 
pH0_NHE = 7.1                 
VMAXTHCO3 = 2.02 * pow(10,-5) 
K_mTHCO3 = 7.38 * pow(10,-3)  
l_THCO3 = 1.63                
pHe0_THCO3 = 6.85             
g_THCO3 = 4.2                 
pHi0_THCO3 = 6.90             
VMAXCA9 = 9.47 * pow(10,-2)   
K_mCA9 = 7.2 * pow(10,-3)     
d_CA9 = 7.3
SensO2 =  1.0
SensATP =  1.0
r_C = 6.55       
pKa = -math.log10(k1 / k2)          
V_c = 2e11

"""
    VMAXAcL = 0.0057480000000000005
    VMAXNHE = 3.0900000000000006e-05
    VMAXTHCO3 = 0.0012120000000000002
    VMAXCA9 = 5.6820000000000004
"""

def ph_ode (t, x):
    dxdt = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    # intracellular carbondioxide dynamics
    dxdt[0] = (
        SensO2 * q_O2 * MW_CO2 / MW_O2 # internal rate
            -k1 * x[0] + k2 * x[1] * x[2] * 1000 * MW_CO2 / ((4.0 / 3.0 * math.pi * pow(r_C, 3.0)) * MW_H * MW_HCO3) # chemical equilibrium
            + PM_CO2 * (x[3] / V_c-x[0] / (4.0 / 3.0 * math.pi * pow(r_C, 3.0))) * (4.0 * math.pi * pow(r_C, 2.0))
    )
            

    # intracellular hydrogen dynamics
    dxdt[1] = (
        SensATP * gAcL * MW_H / MW_AcL # internal rate
        + k1 * x[0] * MW_H / MW_CO2 - k2 *x[1] * x[2] * 1000 / ((4.0 / 3.0 * math.pi * pow(r_C, 3.0)) * MW_HCO3) # chemical equilibrium
        - (2.0 - math.tanh(c2aH_slope * (-math.log10(1000 * x[1] / (4.0 / 3.0 * math.pi * pow(r_C, 3.0) * MW_H))) - c2aH_thr)) * VMAXAcL * MW_H / MW_AcL 
        * (4.0 * math.pi * pow(r_C, 2.0)) * x[1] / ((4.0 / 3.0 * math.pi * pow(r_C, 3.0)) * K_mAcL * MW_H / MW_AcL + x[1])
        + (2.0 - math.tanh(a2cH_slope * (-math.log10(1000 * x[4] / (V_c * MW_H))) - a2cH_thr)) * VMAXAcL * MW_H / MW_AcL 
        * (4.0 * math.pi * pow(r_C, 2.0)) * x[4] / (V_c * K_mAcL * MW_H / MW_AcL + x[4])
        - SensO2 * SensATP * 0.5 * (1.0 + ((-math.log10(1000 * x[4] / (V_c * MW_H))) - pH0_NHE) / (l_NHE + abs((-math.log10(1000 * x[4] / (V_c * MW_H))) - pH0_NHE))) 
        * VMAXNHE * (4.0 * math.pi * pow(r_C, 2.0)) * pow(x[1], a) / (pow(4.0 / 3.0 * math.pi * pow(r_C, 3.0) * MW_H * K_mNHE / 1000, a) + pow(x[1], a))
    )

    # intracellular bicarbonate ions dynamics
    dxdt[2] = (
        k1 * MW_HCO3 / MW_CO2 * x[0] - k2 * x[1] * x[2] * 1000 / ((4.0 / 3.0 * math.pi * pow(r_C, 3.0)) * MW_H) # chemical equilibrium
        + SensATP * (0.5) * (1.0 + math.tanh(l_THCO3 * (-math.log10(1000 * x[4] / (V_c * MW_H)) - pHe0_THCO3))) * (0.5) * (1.0 + math.tanh(g_THCO3 * (pHi0_THCO3 - (-math.log10(1000 * x[1] / (4.0 / 3.0 * math.pi * pow(r_C, 3.0) * MW_H)))))) 
        * VMAXTHCO3 * (4.0 * math.pi * pow(r_C, 2.0)) * x[5] / (V_c * K_mTHCO3 * MW_HCO3 / 1000 + x[5])
    )

    # extracellular carbondioxide dynamics
    dxdt[3] = 0.0

    # extracellular hydrogen dynamics
    dxdt[4] = (
        k1 * x[3] * MW_H / MW_CO2 - k2 * x[4] * x[5] * 1000 / (V_c * MW_HCO3) # chemical equilibrium
	    + (2.0 - math.tanh(c2aH_slope * (-math.log10(1000 * x[1] / (4.0 / 3.0 * math.pi * pow(r_C, 3.0) * MW_H))) - c2aH_thr)) * VMAXAcL * MW_H / MW_AcL 
        * (4.0 * math.pi * pow(r_C, 2.0)) * x[1] / ((4.0 / 3.0 * math.pi * pow(r_C, 3.0)) * K_mAcL * MW_H / MW_AcL + x[1])
        - (2.0 - math.tanh(a2cH_slope * (-math.log10(1000 * x[4] / (V_c * MW_H))) - a2cH_thr)) * VMAXAcL * MW_H / MW_AcL 
        * (4.0 * math.pi * pow(r_C, 2.0)) * x[4] / (V_c * K_mAcL * MW_H / MW_AcL + x[4]) 
        + SensATP * SensO2 * 0.5 * (1.0 + ((-math.log10(1000 * x[4] / (V_c * MW_H))) - pH0_NHE) / (l_NHE + abs((-math.log10(1000 * x[4] / (V_c * MW_H))) - pH0_NHE))) 
        * VMAXNHE * (4.0 * math.pi * pow(r_C, 2.0)) * pow(x[1], a) / (pow(4.0 / 3.0 * math.pi * pow(r_C, 3.0) * MW_H * K_mNHE / 1000, a) + pow(x[1], a)) 
        + (3.0 + 2.0 * math.tanh(-d_CA9 * SensO2)) * VMAXCA9 * 4.0 * math.pi * pow(r_C, 2.0) * x[3] / (V_c * K_mCA9 * MW_CO2 / 1000 + x[3]) 
        * MW_H / MW_CO2
    )
    
    # extracellular bicarbonate ions dynamics
    dxdt[5] = (
        k1 * x[3] * MW_HCO3 / MW_CO2 - k2 * x[4] * x[5] * 1000 / (V_c * MW_H) # chemical equilibrium
		- SensATP * (0.5) * (1.0 + math.tanh(l_THCO3 * (-math.log10(1000 * x[4] / (V_c * MW_H)) - pHe0_THCO3))) 
        * (0.5) * (1.0 + math.tanh(g_THCO3 * (pHi0_THCO3 - (-math.log10(1000 * x[1] / (4.0 / 3.0 * math.pi * pow(r_C, 3.0) * MW_H)))))) 
        * VMAXTHCO3 * (4.0 * math.pi * pow(r_C, 2.0)) * x[5] / (V_c * K_mTHCO3 * MW_HCO3 / 1000 + x[5])
        + (3.0 + 2.0 * math.tanh(-d_CA9 * SensO2)) * VMAXCA9 * 4.0 * math.pi * pow(r_C, 2.0) * x[3] / (V_c * K_mCA9 * MW_CO2 / 1000 + x[3]) * MW_HCO3 / MW_CO2
    )
    return(dxdt)
    
# Starting conditions
pH_cell = 7.4
pH_extra = 7.4

CO2_intra_initial = 5.39 * pow(10.0, -5) * (4.0 / 3.0 * math.pi * pow(r_C, 3.0))
H_intra_initial = pow(10.0, - pH_cell - 3.0) * (4.0 / 3.0 * math.pi * pow(r_C, 3.0))
HCO3_intra_initial = MW_HCO3 / MW_CO2 * CO2_intra_initial * pow(10.0, pH_cell - pKa)
CO2_extra_initial = 5.39 * pow(10.0,-5) * V_c
H_extra_initial = pow(10.0, -pH_extra - 3.0) * V_c
HCO3_extra_initial = MW_HCO3 / MW_CO2 * CO2_extra_initial * pow(10.0, pH_extra - pKa)
x_initial = [CO2_intra_initial, H_intra_initial, HCO3_intra_initial, CO2_extra_initial, H_extra_initial, HCO3_extra_initial]
t = numpy.arange(0, 28800, step=1e-2, dtype=float)
sol = solve_ivp(ph_ode, [0, 28800], x_initial, t_eval=t, method='Radau', rtol=1e-10, atol=1e-12)

pHi_pred = []
for mHi in sol.y[1]: pHi_pred.append(-math.log10(1000 * mHi / ((4.0 / 3.0 * math.pi * pow(r_C, 3.0)) * MW_H) ))
pHe_pred = []
for mHe in sol.y[4]: pHe_pred.append(-math.log10(1000 * mHe / (V_c * MW_H) ))
plt.plot(sol.t, pHi_pred, label='pHi')
plt.plot(sol.t, pHe_pred, label='pHe')
plt.show()