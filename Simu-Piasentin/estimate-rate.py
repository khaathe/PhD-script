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
PM_CO2 = 1920000
gAcL = 0.022800000000000001
q_O2 = 0.0021000000000000003
k1 = 8.64
k2 = 11400000
K_mAcL = 0.405e-3
a2cH_slope = 1.5
a2cH_thr = 7.0
c2aH_slope = 1.5
c2aH_thr = 7.0
K_mNHE = 0.196e-6
a = 2.67
l_NHE = 0.076
pH0_NHE = 7.1
K_mTHCO3 = 7.38e-3
l_THCO3 = 1.63
pHe0_THCO3 = 6.85
g_THCO3 = 4.2
pHi0_THCO3 = 6.90
K_mCA9 = 7.2e-3
d_CA9 = 7.3
pKa = 6.12039110885758
SensO2 =  1.0
SensATP =  1.0
R_cell = 8.41
V_cell = 4/3 * math.pi * math.pow(R_cell, 3)
S_cell = 4 * math.pi * math.pow(R_cell, 2)
V_extracell = 2e11

"""
    VMAXAcL = 0.0057480000000000005
    VMAXNHE = 3.0900000000000006e-05
    VMAXTHCO3 = 0.0012120000000000002
    VMAXCA9 = 5.6820000000000004
"""
# starting conditions
pH_cell = 7.40
pH_extra = 7.4
CO2_intra_initial = 5.39 * pow(10.0, -5) * (4.0 / 3.0 * math.pi * pow(R_cell, 3.0))
H_intra_initial = pow(10.0, -pH_cell - 3.0) * (4.0 / 3.0 * math.pi * pow(R_cell, 3.0))
HCO3_intra_initial = MW_HCO3 / MW_CO2 * CO2_intra_initial * pow(10.0, pH_cell - pKa)

CO2_extra_initial = 5.39 * pow(10.0,-5) * V_extracell
H_extra_initial = pow(10.0, -pH_extra - 3.0) * V_extracell
HCO3_extra_initial = MW_HCO3 / MW_CO2 * CO2_extra_initial * pow(10.0, pH_extra - pKa)

x = [CO2_intra_initial, H_intra_initial, 
    HCO3_intra_initial, CO2_extra_initial, H_extra_initial, HCO3_extra_initial]

print(x)

def ph_ode (t, x, VMAXAcL, VMAXNHE, VMAXTHCO3, VMAXCA9):
    dxdt = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    diff_CO2_in_out = PM_CO2 * ( (x[3] / V_extracell) - (x[0] / V_cell) ) * S_cell
	
    nu_MCT_in_out = (
        (2.0 - math.tanh(c2aH_slope * (-math.log10(1000 * x[1] / (V_cell * MW_H))) - c2aH_thr)) * VMAXAcL * MW_H / MW_AcL 
        * (S_cell) * x[1] / ((V_cell) * K_mAcL * MW_H / MW_AcL + x[1])
    )
        
    nu_MCT_out_in = (
        (2.0 - math.tanh(a2cH_slope * (-math.log10(1000 * x[4] / (V_extracell * MW_H))) - a2cH_thr)) * VMAXAcL * MW_H / MW_AcL 
        * (S_cell) * x[4] / (V_extracell * K_mAcL * MW_H / MW_AcL + x[4])
    )
        
    nu_NHE_in_out = (
        SensO2 * SensATP * 0.5 * (1.0 + ((-math.log10(1000 * x[4] / (V_extracell * MW_H))) - pH0_NHE) / (l_NHE + abs((-math.log10(1000 * x[4] / (V_extracell * MW_H))) - pH0_NHE))) 
        * VMAXNHE * (S_cell) * pow(x[1], a) / (pow(V_cell * MW_H * K_mNHE / 1000, a) + pow(x[1], a))
    )
        
    nu_THCO3_out_in = (
        SensATP * (0.5) * 
        (1.0 + math.tanh(l_THCO3 * (-math.log10(1000 * x[4] / (V_extracell * MW_H)) - pHe0_THCO3))) * (0.5) 
        * (1.0 + math.tanh(g_THCO3 * (pHi0_THCO3 - (-math.log10(1000 * x[2] / (V_cell * MW_H)))))) 
        * VMAXTHCO3 * (S_cell) * x[5] / (V_extracell * K_mTHCO3 * MW_HCO3 / 1000 + x[5])
    )
    
    nu_CA9 = (3.0 + 2.0 * math.tanh(-d_CA9 * SensO2)) * VMAXCA9 * S_cell * x[3] / (V_extracell * K_mCA9 * MW_CO2 / 1000 + x[3])

    # intracellular carbondioxide dynamics
    dxdt[0] = (
        SensO2 * q_O2 * MW_CO2 / MW_O2 # internal rate
            -k1 * x[0] + k2 * x[1] * x[2] * 1000 * MW_CO2 / (V_cell * MW_H * MW_HCO3) # chemical equilibrium
            + diff_CO2_in_out
    )
            

    # intracellular hydrogen dynamics
    dxdt[1] = (
        SensATP * gAcL * MW_H / MW_AcL * MW_H # internal rate
        + k1 * x[0] * MW_H / MW_CO2 - k2 *x[1] * x[2] * 1000 / (V_cell * MW_HCO3) # chemical equilibrium
        - nu_MCT_in_out + nu_MCT_out_in - nu_NHE_in_out
    )

    # intracellular bicarbonate ions dynamics
    dxdt[2] = (
        k1 * MW_HCO3 / MW_CO2 * x[0] - k2 * x[1] * x[2] * 1000 / (V_cell * MW_H) # chemical equilibrium
        + nu_THCO3_out_in
    )

    # extracellular carbondioxide dynamics
    dxdt[3] = 0.0

    # extracellular hydrogen dynamics
    dxdt[4] = (
        k1 * x[3] * MW_H / MW_CO2 - k2 * x[4] * x[5] * 1000 / (V_extracell * MW_HCO3) # chemical equilibrium
	    + nu_MCT_in_out - nu_MCT_out_in + nu_NHE_in_out + nu_CA9  * MW_H / MW_CO2 # extracellular hydrogen dynamics
    )
    
    # extracellular bicarbonate ions dynamics
    dxdt[5] = (
        k1 * x[3] * MW_HCO3 / MW_CO2 - k2 * x[4] * x[5] * 1000 / (V_extracell * MW_H) # chemical equilibrium
		- nu_THCO3_out_in + nu_CA9 * MW_HCO3 / MW_CO2
    )
    return(dxdt)
    
def model (x, y, t, dt, VMAXAcL, VMAXNHE, VMAXTHCO3, VMAXCA9):
    print("Hello World")

def funfit (df, VMAXAcL, VMAXNHE, VMAXTHCO3, VMAXCA9):
    print("Hello World")

u87_file_path = "/home/spinicck/PhD/Data/alaa-experiment/U87_pH_time_regulation.txt"
f98_file_path = "/home/spinicck/PhD/Data/alaa-experiment/F98_pH_time_regulation.txt"
u87_df = pd.read_table(u87_file_path, sep="\t", index_col=0)

t_start = 0.0
t_end = 36000 # 8h -> 28800 s
dt = 1e-2
t = numpy.linspace(t_start, t_end, num=int( (t_end-t_start)/dt ) )
constant_rates = (9.58 * math.pow(10,-5), 5.15 * math.pow(10,-7),  2.02 * math.pow(10,-5), 9.47 * math.pow(10,-2) )
sol = solve_ivp(ph_ode, [t_start, t_end], x, t_eval=t, args=constant_rates, method='Radau' )

pHi = []
pHe = []
for mHi in sol.y[1]: pHi.append( -math.log10(1000 * mHi / (V_cell * MW_H) ) )
for mHe in sol.y[4]: pHe.append( -math.log10(1000 * mHe / (V_extracell * MW_H) ) )

plt.plot(t, pHi)
plt.plot(t, pHe)
plt.legend(['pHi', 'pHe'])
plt.show()