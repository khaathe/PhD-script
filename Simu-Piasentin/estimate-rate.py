#!/usr/bin/python
import math
import numpy
from scipy.optimize import least_squares
from scipy.integrate import solve_ivp
import pandas as pd
import matplotlib.pyplot as plt
import pickle

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
pKa = -math.log10(k1 / k2)         
SensO2 =  1.0
SensATP =  1.0
R_cell = 6.7       
V_cell = 4/3 * math.pi * math.pow(R_cell, 3)
S_cell = 4 * math.pi * math.pow(R_cell, 2)
V_extracell = 2e11

"""
    VMAXAcL = 0.0057480000000000005
    VMAXNHE = 3.0900000000000006e-05
    VMAXTHCO3 = 0.0012120000000000002
    VMAXCA9 = 5.6820000000000004
"""

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
        SensATP * (0.5) * (1.0 + math.tanh(l_THCO3 * (-math.log10(1000 * x[4] / (V_extracell * MW_H)) - pHe0_THCO3))) * (0.5) * (1.0 + math.tanh(g_THCO3 * (pHi0_THCO3 - (-math.log10(1000 * x[1] / V_cell * MW_H))))) 
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
        SensATP * gAcL * MW_H / MW_AcL # internal rate
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
    
def model (t_start, t_end, dt, constant_rates, pH_cell, pH_extra):
    t = numpy.arange(t_start, t_end+dt, step=dt,  dtype=float)
    
    # Starting conditions
    CO2_intra_initial = 5.39 * pow(10.0, -5) * (4.0 / 3.0 * math.pi * pow(R_cell, 3.0))
    H_intra_initial = pow(10.0, -pH_cell - 3.0) * (4.0 / 3.0 * math.pi * pow(R_cell, 3.0))
    HCO3_intra_initial = MW_HCO3 / MW_CO2 * CO2_intra_initial * pow(10.0, pH_cell - pKa)
    CO2_extra_initial = 5.39 * pow(10.0,-5) * V_extracell
    H_extra_initial = pow(10.0, -pH_extra - 3.0) * V_extracell
    HCO3_extra_initial = MW_HCO3 / MW_CO2 * CO2_extra_initial * pow(10.0, pH_extra - pKa)
    x_initial = [CO2_intra_initial, H_intra_initial, HCO3_intra_initial, CO2_extra_initial, H_extra_initial, HCO3_extra_initial]

    sol = solve_ivp(ph_ode, [t_start, t_end+dt], x_initial, t_eval=t, args=(constant_rates), method='LSODA', rtol=2.5e-14, atol=1e-16 )
    return(sol)

def funfit (constant_rates_init, pHi_data, pH_cell, pH_extra, t_start, t_end, dt):
    sol = model(t_start, t_end, dt, constant_rates_init, pH_cell, pH_extra)
    pHi = []
    for mHi in sol.y[1, ::int(time_step/dt)]: pHi.append(-math.log10(1000 * mHi / (V_cell * MW_H) ))
    return(pHi-pHi_data)

def make_fit(pHi_data, pHi0, pHe0, t_start, t_end, dt, constant_rates, constant_bounds ):
    arguments = (pHi_data, pHi0, pHe0, t_start, t_end, dt)
    fit = least_squares(funfit, constant_rates, args=arguments, bounds=constant_bounds)
    return(fit)

def to_pH (mH, V):
    pH = []
    for m in mH: pH.append(-math.log10(1000 * m / (V * MW_H) ))
    return( pH )

u87_file_path = "/home/spinicck/PhD/Data/alaa-experiment/U87_pH_time_regulation.txt"
# f98_file_path = "/home/spinicck/PhD/Data/alaa-experiment/F98_pH_time_regulation.txt"
u87_df = pd.read_table(u87_file_path, sep="\t", index_col=0)

t = numpy.arange(0,9)*60

constant_rates = [ 0.005748, 3.09e-05, 0.001212, 5.682 ]
constant_bounds = (
    [0.003, 1e-5, 1e-5, 3], 
    [1e-1, 1e-1, 1e-1, 8]
) 
time_step = 60
arguments = (u87_df["7.4"].to_numpy(), 7.484051505, 7.4, 0, 8*60, 1e-2/60)

k = 0

with open('fit.txt', 'w') as out:
    header = "pHe\tVMAXAcL\tVMAXNHE\tVMAXTHCO3\tVMAXCA9"
    print(header)
    out.write(header + "\n")
    for pHe in u87_df.columns:
        fit = make_fit(u87_df[pHe], 7.484051505, float(pHe), 0, 8*60, 1e-2/60, constant_rates, constant_bounds)
        pred = model(0, 8*60, 1e-2/60, fit.x, 7.484051505, float(pHe) )
        pHi_final = -math.log10(1000 * pred.y[1, -1] / (V_cell * MW_H) )
        pHe_final = -math.log10(1000 * pred.y[4, -1] / (V_extracell * MW_H) )
        res = f'{pHe}\t{fit.x[0]}\t{fit.x[1]}\t{fit.x[2]}\t{fit.x[3]}'
        print(res)
        out.write(res + "\n")
        k=k+1
        plt.subplot(6,3,k)
        plt.plot(pred.t, to_pH(pred.y[1], V_cell), label='pHi model')
        plt.plot(pred.t, to_pH(pred.y[4], V_extracell), label='pHe model')
        plt.plot(t, u87_df[pHe], 'o', label='pHi Data')
        plt.title(f'pHe : {pHe}', loc='left')

plt.show()