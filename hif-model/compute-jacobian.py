#!/usr/bin/python
from sympy import Matrix, MatrixSymbol
from sympy.core.function import Function
from sympy.core.symbol import Symbol, symbols
from sympy import tanh
from sympy.codegen.cfunctions import log10, Pow
from sympy.printing.c import C99CodePrinter
from sympy.functions import exp

O,G,N,HIF_PROD,MW_H,MW_CO2,MW_O2,MW_HCO3,MW_AcL,V_cell,V_extracell,S_cell,PM_CO2,gAcL,q_O2,k1,k2,VMAXAcL,K_mAcL,a2cH_slope,a2cH_thr,c2aH_slope,c2aH_thr,VMAXNHE,K_mNHE,a,l_NHE,pH0_NHE,VMAXTHCO3,K_mTHCO3,l_THCO3,pHe0_THCO3,g_THCO3,pHi0_THCO3,VMAXCA9,K_mCA9,d_CA9,pKa,pH_cell,SensO2,SensATP = symbols('O,G,N,HIF_PROD,MW_H,MW_CO2,MW_O2,MW_HCO3,MW_AcL,V_cell,V_extracell,S_cell,PM_CO2,gAcL,q_O2,k1,k2,VMAXAcL,K_mAcL,a2cH_slope,a2cH_thr,c2aH_slope,c2aH_thr,VMAXNHE,K_mNHE,a,l_NHE,pH0_NHE,VMAXTHCO3,K_mTHCO3,l_THCO3,pHe0_THCO3,g_THCO3,pHi0_THCO3,VMAXCA9,K_mCA9,d_CA9,pKa,pH_cell,SensO2,SensATP', real=True)

CO2_cell,H_cell,HCO3_cell,CO2_extracell, H_extracell, HCO3_extracell = symbols('CO2_cell,H_cell,HCO3_cell,CO2_extracell, H_extracell, HCO3_extracell', real=True)

Oxygen, Glucose, ATP, Proton = symbols("Oxygen,Glucose,ATP,H", real = True)

X = MatrixSymbol('x', 8, 1)
A = MatrixSymbol('A', 4, 1)
D = MatrixSymbol('D', 4, 1)
S = MatrixSymbol('S', 5, 5)
gamma_array = MatrixSymbol('gamma', 5, 5)

Vo, Ko, A0, Kg, Kh, pg_max, pg_min, k, ldh0, po_max, po_min, l, pdh0 = symbols('Vo, Ko, A0,  Kg, Kh, pg_max, pg_min, k, ldh0, po_max, po_min, l, pdh0', real=True)
po = (po_max-po_min)/(1+exp(-l*(X[3]-pdh0))) + po_min
pg = (pg_max-pg_min)/(1+exp(-k*(X[1]-ldh0))) + pg_min
Co = po * Vo * ( O/(O+Ko) )
Cg = ( (pg*A0/2.0) - (27.0*Co/10.0) ) * ( G/(G+Kg) )
Pa = ( (2.0 * Cg) + ( (27.0/5.0)*Co ) )
Ph = Kh*( (29.0*(pg*Vo-Co))/5.0 )

ode = Matrix([
    # HIF
    [A[0] - D[0] * X[0] * ((Pow(S[0, 1], N) / (Pow(S[0, 1], N) + Pow(O, N))) + gamma_array[0, 1] * (Pow(O, N) / (Pow(S[0, 1], N) + Pow(O, N)))) ],
    # LDH
    [A[1] * ((Pow(S[1, 2], N) / (Pow(S[1, 2], N) + Pow(X[0], N))) + gamma_array[1, 2] * (Pow(X[0], N) / (Pow(S[1, 2], N) + Pow(X[0], N)))) - D[1] * X[1] ],
    # PDK
    [A[2] * ((Pow(S[1, 3], N) / (Pow(S[1, 3], N) + Pow(X[0], N))) + gamma_array[1, 3] * (Pow(X[0], N) / (Pow(S[1, 3], N) + Pow(X[0], N)))) - D[2] * X[2] ],
    # PDK
    [A[3] * ((Pow(S[3, 4], N) / (Pow(S[3, 4], N) + Pow(X[2], N))) + gamma_array[3, 4] * (Pow(X[2], N) / (Pow(S[3, 4], N) + Pow(X[2], N)))) - D[3] * X[3] ],
    # ATP
    [Pa],
    # Oxygen
    [-Co],
    # Glucose
    [-Cg],
    # Protons
    [Ph]
    # # CO2 cell
    # [
    #     SensO2 * q_O2 * MW_CO2 / MW_O2 
    #     - k1 * CO2_cell + k2 * H_cell * HCO3_cell * 1000 * MW_CO2 / (V_cell * MW_H * MW_HCO3) 
    #     + PM_CO2 * (CO2_extracell / V_extracell-CO2_cell / (V_cell)) * (S_cell)
    # ],
    # # H cell
    # [
    #     SensATP * Ph * V_cell * 1e-6 
    #     + k1 * CO2_cell * MW_H / MW_CO2 - k2 * H_cell * HCO3_cell * 1000 / (V_cell * MW_HCO3) 
    #     - (2.0 - tanh(c2aH_slope * (-log10(1000 * H_cell / (V_cell * MW_H))) - c2aH_thr)) * VMAXAcL * MW_H / MW_AcL * (S_cell) * H_cell / ((V_cell) * K_mAcL * MW_H / MW_AcL + H_cell) 
    #     + (2.0 - tanh(a2cH_slope * (-log10(1000 * H_extracell / (V_extracell * MW_H))) - a2cH_thr)) * VMAXAcL * MW_H / MW_AcL * (S_cell) * H_extracell / (V_extracell * K_mAcL * MW_H / MW_AcL + H_extracell) 
    #     - SensO2 * SensATP * 0.5 * (1.0 + ((-log10(1000 * H_extracell / (V_extracell * MW_H))) - pH0_NHE) / (l_NHE + abs((-log10(1000 * H_extracell / (V_extracell * MW_H))) - pH0_NHE))) * VMAXNHE * (S_cell) * Pow(H_cell, a) / (Pow(V_cell * MW_H * K_mNHE / 1000, a) + Pow(H_cell, a))
    # ],
    # # HCO3 cell
    # [
    #     k1 * MW_HCO3 / MW_CO2 * CO2_cell - k2 * H_cell * HCO3_cell * 1000 / (V_cell * MW_H) 
    #     + SensATP * (0.5) * (1.0 + tanh(l_THCO3 * (-log10(1000 * H_extracell / (V_extracell * MW_H)) - pHe0_THCO3))) * (0.5) * (1.0 
    #     +  tanh(g_THCO3 * (pHi0_THCO3 - (-log10(1000 * H_cell / (V_cell * MW_H)))))) * VMAXTHCO3 * (S_cell) * HCO3_extracell / (V_extracell * K_mTHCO3 * MW_HCO3 / 1000 + HCO3_extracell)
    # ],
    # # H extracell
    # [
    #     k1 * CO2_extracell * MW_H / MW_CO2 - k2 * H_extracell * HCO3_extracell * 1000 / (V_extracell * MW_HCO3)
    #     + (2.0 - tanh(c2aH_slope * (-log10(1000 * H_cell / (V_cell * MW_H))) - c2aH_thr)) * VMAXAcL * MW_H / MW_AcL * (S_cell) * H_cell / ((V_cell) * K_mAcL * MW_H / MW_AcL + H_cell)
    #     - (2.0 - tanh(a2cH_slope * (-log10(1000 * H_extracell / (V_extracell * MW_H))) - a2cH_thr)) * VMAXAcL * MW_H / MW_AcL * (S_cell) * H_extracell / (V_extracell * K_mAcL * MW_H / MW_AcL + H_extracell) 
    #     + SensATP * SensO2 * 0.5 * (1.0 + ((-log10(1000 * H_extracell / (V_extracell * MW_H))) - pH0_NHE) / (l_NHE + abs((-log10(1000 * H_extracell / (V_extracell * MW_H))) - pH0_NHE))) * VMAXNHE * (S_cell) * pow(H_cell, a) / (pow(V_cell * MW_H * K_mNHE / 1000, a) + pow(H_cell, a)) 
    #     + (3.0 + 2.0 * tanh(-d_CA9 * SensO2)) * VMAXCA9 * S_cell * CO2_extracell / (V_extracell * K_mCA9 * MW_CO2 / 1000 + CO2_extracell) * MW_H / MW_CO2
    # ],
    # # HCO3 extracell
    # [
    #     k1 * CO2_extracell * MW_HCO3 / MW_CO2 - k2 * H_extracell * HCO3_extracell * 1000 / (V_extracell * MW_H)
    #     - SensATP * (0.5) * (1.0 + tanh(l_THCO3 * (-log10(1000 * H_extracell / (V_extracell * MW_H)) - pHe0_THCO3))) * (0.5) * (1.0 + tanh(g_THCO3 * (pHi0_THCO3 - (-log10(1000 * H_cell / (V_cell * MW_H))))))  * VMAXTHCO3 * (S_cell) * HCO3_extracell / (V_extracell * K_mTHCO3 * MW_HCO3 / 1000 + HCO3_extracell)
    #     + (3.0 + 2.0 * tanh(-d_CA9 * SensO2)) * VMAXCA9 * S_cell * CO2_extracell / (V_extracell * K_mCA9 * MW_CO2 / 1000 + CO2_extracell) * MW_HCO3 / MW_CO2
    # ]
])
var = Matrix([
    [X[0]],
    [X[1]],
    [X[2]],
    [X[3]],
    [ATP],
    [Oxygen],
    [Glucose],
    [Proton]
    # [CO2_cell],
    # [H_cell],
    # [HCO3_cell],
    # [H_extracell],
    # [HCO3_extracell]
])
Jacobian = ode.jacobian(var)
class MyCodePrinter(C99CodePrinter):        
    def _print_MatrixElement(self, expr):
        matrice_element = '{0}[{1}][{2}]'.format(expr.parent, expr.i, expr.j)
        if expr.parent.name.upper() in ['X', 'A', 'D']:
            matrice_element = '{0}[{1}]'.format(expr.parent.name, expr.i)
        elif expr.parent.name.upper() == 'J':
            matrice_element = '{0}({1})({2})'.format(expr.parent, expr.i, expr.j)
        return self._print(matrice_element)

my_printer = MyCodePrinter()

cpp_code = '''void jacobi ( const state_type &x , matrix_type &J , const double &t , state_type &dfdt, SimuKevinParameters *p, vector<double> nearest_gradient, Cell *pCell )
{
    // Binding ODE parameters to Variables
	const auto& [
		A, D, N, hif_coeff_prod, Vo, Ko, pg, A0,  Kg, Kh, Vmax_f, Vmax_r, Km_Hi, Km_He,
		MW_H, MW_CO2, MW_O2, MW_HCO3, MW_AcL, r_C, PM_CO2, gAcL, q_O2, k1, k2,  VMAXAcL, K_mAcL, a2cH_slope, a2cH_thr, c2aH_slope,
		c2aH_thr,  VMAXNHE, K_mNHE, a, l_NHE, pH0_NHE,  VMAXTHCO3, K_mTHCO3, l_THCO3, pHe0_THCO3, g_THCO3, pHi0_THCO3, VMAXCA9, 
		K_mCA9, d_CA9, pKa, pH_cell, SensO2, SensATP
	] = p->ode_parameters;
	auto S = p->S;
	auto gamma_array = p->gamma_array;

	const double O = p->nearest_density[microenvironment.find_density_index("oxygen")] / (microenvironment.mesh.dV * 1e-15); // Convert mmol to mmol/L
	const double H_extracell = p->nearest_density[microenvironment.find_density_index("proton")]; 
	const double HCO3_extracell = p->nearest_density[microenvironment.find_density_index("HCO3")];
	const double V_extracell = p->simu_parameters->V_extracell;
	const double V_cell = p->cell->phenotype.volume.total;
	const double S_cell = p->cell->phenotype.geometry.surface_area;
    double CO2_cell = x[5];
    double H_cell = x[6];
    double HCO3_cell = x[7];

'''
i=0
for i in range(Jacobian.shape[0]):
    for j in range(Jacobian.shape[1]):
        cpp_code = cpp_code + "\t" + my_printer.doprint(Jacobian[i,j], assign_to=f'J({i},{j})') + '\n'

for i in range(Jacobian.shape[0]):
    cpp_code = cpp_code + "\t" + f'dfdt[{i}]=0.0;' + '\n'

cpp_code = cpp_code + "}"

with open('/home/spinicck/PhD/Code/output/jacobi.cpp', 'w') as out:
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