#!/usr/bin/python
from sympy import Matrix, MatrixSymbol
from sympy.core.symbol import Symbol, symbols
from sympy import tanh
from sympy.codegen.cfunctions import log10, Pow
from sympy.printing.c import C99CodePrinter
import re

def replace_x(line):
    pat = "X\[(\d+),(\d+)\]"
    match = re.search(pat, line)
    while match != None:
        replacement = f"pow({match.group(1)},{match.group(2)})"
        replacement = f'X[{match.group(1)}]'
        line = re.sub(pat, replacement, line, 1, re.IGNORECASE)
        match = re.search(pat, line)
    return line

def replace_power(line):
    member_of_power = "\w+|\w+\[\d+,?\d*\]"
    pat = f"({member_of_power})\*\*({member_of_power})"
    match = re.search(pat, line)
    while match != None:
        replacement = f"pow({match.group(1)},{match.group(2)})"
        line = re.sub(pat, replacement, line, 1, re.IGNORECASE)
        match = re.search(pat, line)
    return line

O,A,D,N,HIF_PROD,MW_H,MW_CO2,MW_O2,MW_HCO3,MW_AcL,V_C,V_c,S_C,PM_CO2,gAcL,q_O2,k1,k2,VMAXAcL,K_mAcL,a2cH_slope,a2cH_thr,c2aH_slope,c2aH_thr,VMAXNHE,K_mNHE,a,l_NHE,pH0_NHE,VMAXTHCO3,K_mTHCO3,l_THCO3,pHe0_THCO3,g_THCO3,pHi0_THCO3,VMAXCA9,K_mCA9,d_CA9,pKa,pH_cell,SensO2,SensATP = symbols('O,A,D,N,HIF_PROD,MW_H,MW_CO2,MW_O2,MW_HCO3,MW_AcL,V_C,V_c,S_C,PM_CO2,gAcL,q_O2,k1,k2,VMAXAcL,K_mAcL,a2cH_slope,a2cH_thr,c2aH_slope,c2aH_thr,VMAXNHE,K_mNHE,a,l_NHE,pH0_NHE,VMAXTHCO3,K_mTHCO3,l_THCO3,pHe0_THCO3,g_THCO3,pHi0_THCO3,VMAXCA9,K_mCA9,d_CA9,pKa,pH_cell,SensO2,SensATP')

m_CO2_C,m_H_C,m_HCO3_C,m_CO2_c_old, m_H_c, m_HCO3_c = symbols('m_CO2_C,m_H_C,m_HCO3_C,m_CO2_c_old, m_H_c, m_HCO3_c')

X = MatrixSymbol('x', 8, 1)
S = MatrixSymbol('S', 5, 5)
Ph = Symbol('Ph')
gamma_array = MatrixSymbol('gamma_array', 5, 5)
M1 = Matrix([
    [A * HIF_PROD - D * X[0] * ((Pow(S[0, 1], N) / (Pow(S[0, 1], N) + Pow(O, N))) + gamma_array[0, 1] * (Pow(O, N) / (Pow(S[0, 1], N) + Pow(O, N)))) ],
    [A * ((Pow(S[1, 2], N) / (Pow(S[1, 2], N) + Pow(X[0], N))) + gamma_array[1, 2] * (Pow(X[0], N) / (Pow(S[1, 2], N) + Pow(X[0], N)))) - D * X[1] ],
    [A * ((Pow(S[1, 3], N) / (Pow(S[1, 3], N) + Pow(X[0], N))) + gamma_array[1, 3] * (Pow(X[0], N) / (Pow(S[1, 3], N) + Pow(X[0], N)))) - D * X[2] ],
    [A * ((Pow(S[3, 4], N) / (Pow(S[3, 4], N) + Pow(X[2], N))) + gamma_array[3, 4] * (Pow(X[2], N) / (Pow(S[3, 4], N) + Pow(X[2], N)))) - D * X[3] ],
    [0.0],
    [SensO2 * q_O2 * MW_CO2 / MW_O2 - k1 * m_CO2_C + k2 * m_H_C * m_HCO3_C * 1000 * MW_CO2 / (V_C * MW_H * MW_HCO3) + PM_CO2 * (m_CO2_c_old / V_c-m_CO2_C / (V_C)) * (S_C)],
    [SensATP * Ph * V_C * 1e-6 + k1 * m_CO2_C * MW_H / MW_CO2 - k2 * m_H_C * m_HCO3_C * 1000 / (V_C * MW_HCO3) - (2.0 - tanh(c2aH_slope * (-log10(1000 * m_H_C / (V_C * MW_H))) - c2aH_thr)) * VMAXAcL * MW_H / MW_AcL * (S_C) * m_H_C / ((V_C) * K_mAcL * MW_H / MW_AcL + m_H_C) + (2.0 - tanh(a2cH_slope * (-log10(1000 * m_H_c / (V_c * MW_H))) - a2cH_thr)) * VMAXAcL * MW_H / MW_AcL * (S_C) * m_H_c / (V_c * K_mAcL * MW_H / MW_AcL + m_H_c) - SensO2 * SensATP * 0.5 * (1.0 + ((-log10(1000 * m_H_c / (V_c * MW_H))) - pH0_NHE) / (l_NHE + abs((-log10(1000 * m_H_c / (V_c * MW_H))) - pH0_NHE))) * VMAXNHE * (S_C) * Pow(m_H_C, a) / (Pow(V_C * MW_H * K_mNHE / 1000, a) + Pow(m_H_C, a))],
    [k1 * MW_HCO3 / MW_CO2 * m_CO2_C - k2 * m_H_C * m_HCO3_C * 1000 / (V_C * MW_H) + SensATP * (0.5) * (1.0 + tanh(l_THCO3 * (-log10(1000 * m_H_c / (V_c * MW_H)) - pHe0_THCO3))) * (0.5) * (1.0 +  tanh(g_THCO3 * (pHi0_THCO3 - (-log10(1000 * m_H_C / (V_C * MW_H)))))) * VMAXTHCO3 * (S_C) * m_HCO3_c / (V_c * K_mTHCO3 * MW_HCO3 / 1000 + m_HCO3_c)]
])
M2 = Matrix([
    [X[0]],
    [X[1]],
    [X[2]],
    [X[3]],
    [X[5]],
    [m_CO2_C],
    [m_H_C],
    [m_HCO3_C]
])
Jacobian = M1.jacobian(M2)
class MyCodePrinter(C99CodePrinter):        
    def _print_MatrixElement(self, expr):
        matrice_element = '{0}[{1}][{2}]'.format(expr.parent, expr.i, expr.j)
        if expr.parent.name.upper() == 'X':
            matrice_element = '{0}[{1}]'.format(expr.parent.name, expr.i)
        elif expr.parent.name.upper() == 'J':
            matrice_element = '{0}({1})({2})'.format(expr.parent, expr.i, expr.j)
        return self._print(matrice_element)

my_printer = MyCodePrinter()

cpp_code = '''void jacobi ( const state_type &x , matrix_type &J , const double &t , state_type &dfdt, SimuKevinParameters *p, vector<double> nearest_gradient, Cell *pCell )
{
    // Binding ODE parameters to Variables
	auto& [
		A, D, N, hif_coeff_prod, Vo, Ko, pg, A0,  Kg, Kh, Vmax_f, Vmax_r, Km_Hi, Km_He,
		MW_H, MW_CO2, MW_O2, MW_HCO3, MW_AcL, r_C, PM_CO2, gAcL, q_O2, k1, k2,  VMAXAcL, K_mAcL, a2cH_slope, a2cH_thr, c2aH_slope,
		c2aH_thr,  VMAXNHE, K_mNHE, a, l_NHE, pH0_NHE,  VMAXTHCO3, K_mTHCO3, l_THCO3, pHe0_THCO3, g_THCO3, pHi0_THCO3, VMAXCA9, 
		K_mCA9, d_CA9, pKa, pH_cell, SensO2, SensATP
	] = p->ode_parameters;
	auto S = p->S;
	auto gamma_array = p->gamma_array;

	const double m_CO2_C = x[5];
    const double m_H_C = x[6];
    const double m_HCO3_C = x[7];
    const double m_H_c = nearest_gradient[microenvironment.find_density_index("proton")]; 
	const double m_HCO3_c = nearest_gradient[microenvironment.find_density_index("HCO3")];
	const double m_CO2_c_old = nearest_gradient[microenvironment.find_density_index("CO2")];
	const double O = nearest_gradient[microenvironment.find_density_index("oxygen")];
	const double V_c = microenvironment.mesh.dV;
	const double V_C = pCell->phenotype.volume.total;
	const double S_C = pCell->phenotype.geometry.surface_area;

'''
i=0
for i in range(Jacobian.shape[0]):
    for j in range(Jacobian.shape[1]):
        cpp_code = cpp_code + "\t" + my_printer.doprint(Jacobian[i,j], assign_to=f'J({i},{j})') + '\n'

for i in range(Jacobian.shape[0]):
    cpp_code = cpp_code + "\t" + f'dfdt[{i}]=0.0;' + '\n'

cpp_code = cpp_code + "}"

with open('jacobi.cpp', 'w') as out:
    out.write(cpp_code)


# Just In case
    # [
    #     k1 * m_CO2_c_old * MW_H / MW_CO2 - k2 * m_H_c * m_HCO3_c * 1000 / (V_c * MW_HCO3)
    #     + (2.0 - tanh(c2aH_slope * (-log10(1000 * m_H_C / (V_C * MW_H))) - c2aH_thr)) * VMAXAcL * MW_H / MW_AcL * (S_C) * m_H_C / ((V_C) * K_mAcL * MW_H / MW_AcL + m_H_C)
    #     - (2.0 - tanh(a2cH_slope * (-log10(1000 * m_H_c / (V_c * MW_H))) - a2cH_thr)) * VMAXAcL * MW_H / MW_AcL * (S_C) * m_H_c / (V_c * K_mAcL * MW_H / MW_AcL + m_H_c) 
    #     + SensATP * SensO2 * 0.5 * (1.0 + ((-log10(1000 * m_H_c / (V_c * MW_H))) - pH0_NHE) / (l_NHE + abs((-log10(1000 * m_H_c / (V_c * MW_H))) - pH0_NHE))) * VMAXNHE * (S_C) * pow(m_H_C, a) / (pow(V_C * MW_H * K_mNHE / 1000, a) + pow(m_H_C, a)) 
    #     + (3.0 + 2.0 * tanh(-d_CA9 * SensO2)) * VMAXCA9 * S_C * m_CO2_c_old / (V_c * K_mCA9 * MW_CO2 / 1000 + m_CO2_c_old) * MW_H / MW_CO2
    # ],
    # [
    #     k1 * m_CO2_c_old * MW_HCO3 / MW_CO2 - k2 * m_H_c * m_HCO3_c * 1000 / (V_c * MW_H)
    #     - SensATP * (0.5) * (1.0 + tanh(l_THCO3 * (-log10(1000 * m_H_c / (V_c * MW_H)) - pHe0_THCO3))) * (0.5) * (1.0 + tanh(g_THCO3 * (pHi0_THCO3 - (-log10(1000 * m_H_C / (V_C * MW_H))))))  * VMAXTHCO3 * (S_C) * m_HCO3_c / (V_c * K_mTHCO3 * MW_HCO3 / 1000 + m_HCO3_c)
    #     + (3.0 + 2.0 * tanh(-d_CA9 * SensO2)) * VMAXCA9 * S_C * m_CO2_c_old / (V_c * K_mCA9 * MW_CO2 / 1000 + m_CO2_c_old) * MW_HCO3 / MW_CO2
    # ]