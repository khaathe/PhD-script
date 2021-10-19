#!/usr/bin/julia

using Plots
using DifferentialEquations
using LsqFit
using LSODA




function pow(x,y)
    return x^y
end

function ph_ode(dxdt, x, p, t)
    VMAXAcL = p[1]
    VMAXNHE = p[2]
    VMAXTHCO3 = p[3]
    VMAXCA9 = p[4]

    diff_CO2_in_out = PM_CO2 * ( (x[4] / V_extracell) - (x[1] / V_cell) ) * S_cell
	
    nu_MCT_in_out = (
        (2.0 - tanh(c2aH_slope * (-log10(1000 * x[2] / (V_cell * MW_H))) - c2aH_thr)) * VMAXAcL * MW_H / MW_AcL 
        * (S_cell) * x[2] / ((V_cell) * K_mAcL * MW_H / MW_AcL + x[2])
    )
        
    nu_MCT_out_in = (
        (2.0 - tanh(a2cH_slope * (-log10(1000 * x[5] / (V_extracell * MW_H))) - a2cH_thr)) * VMAXAcL * MW_H / MW_AcL 
        * (S_cell) * x[5] / (V_extracell * K_mAcL * MW_H / MW_AcL + x[5])
    )
        
    nu_NHE_in_out = (
        SensO2 * SensATP * 0.5 * (1.0 + ((-log10(1000 * x[5] / (V_extracell * MW_H))) - pH0_NHE) / (l_NHE + abs((-log10(1000 * x[5] / (V_extracell * MW_H))) - pH0_NHE))) 
        * VMAXNHE * (S_cell) * pow(x[2], a) / (pow(V_cell * MW_H * K_mNHE / 1000, a) + pow(x[2], a))
    )
        
    nu_THCO3_out_in = (
        SensATP * (0.5) * (1.0 + tanh(l_THCO3 * (-log10(1000 * x[5] / (V_extracell * MW_H)) - pHe0_THCO3))) * (0.5) * (1.0 + tanh(g_THCO3 * (pHi0_THCO3 - (-log10(1000 * x[2] / V_cell * MW_H))))) 
        * VMAXTHCO3 * (S_cell) * x[6] / (V_extracell * K_mTHCO3 * MW_HCO3 / 1000 + x[6])
    )

    nu_CA9 = (3.0 + 2.0 * tanh(-d_CA9 * SensO2)) * VMAXCA9 * S_cell * x[4] / (V_extracell * K_mCA9 * MW_CO2 / 1000 + x[4])

    # intracellular carbondioxide dynamics
    dxdt[1] = (
        SensO2 * q_O2 * MW_CO2 / MW_O2 # internal rate
            -k1 * x[1] + k2 * x[2] * x[3] * 1000 * MW_CO2 / (V_cell * MW_H * MW_HCO3) # chemical equilibrium
            + diff_CO2_in_out
    )
            

    # intracellular hydrogen dynamics
    dxdt[2] = (
        SensATP * gAcL * MW_H / MW_AcL # internal rate
        + k1 * x[1] * MW_H / MW_CO2 - k2 *x[2] * x[3] * 1000 / (V_cell * MW_HCO3) # chemical equilibrium
        - nu_MCT_in_out + nu_MCT_out_in - nu_NHE_in_out
    )

    # intracellular bicarbonate ions dynamics
    dxdt[3] = (
        k1 * MW_HCO3 / MW_CO2 * x[1] - k2 * x[2] * x[3] * 1000 / (V_cell * MW_H) # chemical equilibrium
        + nu_THCO3_out_in
    )

    # extracellular carbondioxide dynamics
    dxdt[4] = 0.0

    # extracellular hydrogen dynamics
    dxdt[5] = (
        k1 * x[4] * MW_H / MW_CO2 - k2 * x[5] * x[6] * 1000 / (V_extracell * MW_HCO3) # chemical equilibrium
	    + nu_MCT_in_out - nu_MCT_out_in + nu_NHE_in_out + nu_CA9  * MW_H / MW_CO2 # extracellular hydrogen dynamics
    )
    
    # extracellular bicarbonate ions dynamics
    dxdt[6] = (
        k1 * x[4] * MW_HCO3 / MW_CO2 - k2 * x[5] * x[6] * 1000 / (V_extracell * MW_H) # chemical equilibrium
		- nu_THCO3_out_in + nu_CA9 * MW_HCO3 / MW_CO2
    )
    return dxdt
end
    
function model(x, constant_rates)
    # Starting conditions
    CO2_intra_initial = 5.39 * pow(10.0, -5) * (4.0 / 3.0 * pi * pow(R_cell, 3.0))
    H_intra_initial = pow(10.0, -pH_cell - 3.0) * (4.0 / 3.0 * pi * pow(R_cell, 3.0))
    HCO3_intra_initial = MW_HCO3 / MW_CO2 * CO2_intra_initial * pow(10.0, pH_cell - pKa)
    CO2_extra_initial = 5.39 * pow(10.0,-5) * V_extracell
    H_extra_initial = pow(10.0, -pH_extra - 3.0) * V_extracell
    HCO3_extra_initial = MW_HCO3 / MW_CO2 * CO2_extra_initial * pow(10.0, pH_extra - pKa)
    x_initial = [CO2_intra_initial, H_intra_initial, HCO3_intra_initial, CO2_extra_initial, H_extra_initial, HCO3_extra_initial]
    
    prob = ODEProblem(ph_ode, x_initial, (t_start, t_end), constant_rates, saveat=60.0 )
    # pred = solve(prob, RadauIIA5(), reltol=1e-14, abstol=1e-16)
    pred = solve(prob, lsoda(), reltol=1e-12, abstol=1e-14)

    return [-log10(1000 * u[2] / (V_cell * MW_H) ) for (u) in pred.u]
end

const MW_H = 1.0
const MW_CO2 = 44.0
const MW_O2 = 32.0
const MW_HCO3 = 61.0
const MW_AcL = 90.1
const PM_CO2 = 1920000
const gAcL = 0.022800000000000001
const q_O2 = 0.0021000000000000003
const k1 = 8.64
const k2 = 11400000
const K_mAcL = 0.405e-3
const a2cH_slope = 1.5
const a2cH_thr = 7.0
const c2aH_slope = 1.5
const c2aH_thr = 7.0
const K_mNHE = 0.196e-6
const a = 2.67
const l_NHE = 0.076
const pH0_NHE = 7.1
const K_mTHCO3 = 7.38e-3
const l_THCO3 = 1.63
const pHe0_THCO3 = 6.85
const g_THCO3 = 4.2
const pHi0_THCO3 = 6.90
const K_mCA9 = 7.2e-3
const d_CA9 = 7.3
const pKa = -log10(k1 / k2)         
const SensO2 =  1.0
const SensATP =  1.0
const R_cell = 6.7       
const V_cell = 4/3 * pi * pow(R_cell, 3)
const S_cell = 4 * pi * pow(R_cell, 2)
const V_extracell = 2e11
const time_step = 60
const t_start = 0.0
const t_end = 480.0
const  pH_cell = 7.4
const  pH_extra = 7.4
    
 
constant_rates = [ 0.005748, 3.09e-05, 0.001212, 5.682 ]
t = collect(0.0:8.0).*60
lb = [0.0, 0.0, 0.0, 0.0]
ub = [1.0, 1.0, 1.0, 20.0]
pHi_data = [7.484051505, 7.16635527, 7.182903, 7.221744159, 7.253844111, 7.167054266, 7.192465359, 7.2350498, 7.253665358]

fit = curve_fit(model, t, pHi_data, constant_rates, lower=lb, upper=ub)

# println("Running Model ...")
# @time pred = model(0.0, 480.0, constant_rates, 7.4, 7.4)
# println("Plotting Result")
# plot(pred.t, [-log10(1000 * u[2] / (V_cell * MW_H) ) for (u) in pred.u], label="pHi" )