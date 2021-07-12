from sympy import Matrix, MatrixSymbol
from sympy.core.symbol import symbols
from sympy import tanh
from sympy.codegen.cfunctions import log10

# A, D, N, S, Gamma, X, Y = symbols('A,D,N,S,Gamma,X,Y')
# M1 = Matrix([
#     [A * 140.0 - D * X * ((pow(S, N) / (pow(S, N) + pow(X, N))) + Gamma * (pow(X, N) / (pow(S, N) + pow(X, N)))) ]
# ])

A,D,N,X,HIF_PROD = symbols('A,D,N,X,HIF_PROD')

MW_H,MW_CO2,MW_O2,MW_HCO3,MW_AcL,V_C,V_c,S_C,PM_CO2,gAcL,q_O2,k1,k2,VMAXAcL,K_mAcL,a2cH_slope,a2cH_thr,c2aH_slope,c2aH_thr,VMAXNHE,K_mNHE,a,l_NHE,pH0_NHE,VMAXTHCO3,K_mTHCO3,l_THCO3,pHe0_THCO3,g_THCO3,pHi0_THCO3,VMAXCA9,K_mCA9,d_CA9,pKa,pH_cell,SensO2,SensATP = symbols('MW_H,MW_CO2,MW_O2,MW_HCO3,MW_AcL,V_C,V_c,S_C,PM_CO2,gAcL,q_O2,k1,k2,VMAXAcL,K_mAcL,a2cH_slope,a2cH_thr,c2aH_slope,c2aH_thr,VMAXNHE,K_mNHE,a,l_NHE,pH0_NHE,VMAXTHCO3,K_mTHCO3,l_THCO3,pHe0_THCO3,g_THCO3,pHi0_THCO3,VMAXCA9,K_mCA9,d_CA9,pKa,pH_cell,SensO2,SensATP')

m_CO2_C,m_H_C,m_HCO3_C,m_CO2_c_old, m_H_c, m_HCO3_c = symbols('m_CO2_C,m_H_C,m_HCO3_C,m_CO2_c_old, m_H_c, m_HCO3_c')

S = MatrixSymbol('X', 5, 5)
Gamma = MatrixSymbol('Gamma', 5, 5)
M1 = Matrix([
    [A * HIF_PROD - D * X * ((pow(S[0, 1], N) / (pow(S[0, 1], N) + pow(X, N))) + Gamma[0, 1] * (pow(X, N) / (pow(S[0, 1], N) + pow(X, N)))) ],
    [A * ((pow(S[1, 2], N) / (pow(S[1, 2], N) + pow(X, N))) + Gamma[1, 2] * (pow(X, N) / (pow(S[1, 2], N) + pow(X, N)))) - D * X ],
    [A * ((pow(S[1, 3], N) / (pow(S[1, 3], N) + pow(X, N))) + Gamma[1, 3] * (pow(X, N) / (pow(S[1, 3], N) + pow(X, N)))) - D * X ],
    [A * ((pow(S[3, 4], N) / (pow(S[3, 4], N) + pow(X, N))) + Gamma[3, 4] * (pow(X, N) / (pow(S[3, 4], N) + pow(X, N)))) - D * X ],
    [
        SensO2 * q_O2 * MW_CO2 / MW_O2 
        - k1 * m_CO2_C + k2 * m_H_C * m_HCO3_C * 1000 * MW_CO2 / (V_C * MW_H * MW_HCO3) 
        + PM_CO2 * (m_CO2_c_old / V_c-m_CO2_C / (V_C)) * (S_C)
    ],
    [    
        SensATP * gAcL * MW_H / MW_AcL
        + k1 * m_CO2_C * MW_H / MW_CO2 - k2 * m_H_C * m_HCO3_C * 1000 / (V_C * MW_HCO3)
        - (2.0 - tanh(c2aH_slope * (-log10(1000 * m_H_C / (V_C * MW_H))) - c2aH_thr)) * VMAXAcL * MW_H / MW_AcL * (S_C) * m_H_C / ((V_C) * K_mAcL * MW_H / MW_AcL + m_H_C)
        + (2.0 - tanh(a2cH_slope * (-log10(1000 * m_H_c / (V_c * MW_H))) - a2cH_thr)) * VMAXAcL * MW_H / MW_AcL * (S_C) * m_H_c / (V_c * K_mAcL * MW_H / MW_AcL + m_H_c)
        - SensO2 * SensATP * 0.5 * (1.0 + ((-log10(1000 * m_H_c / (V_c * MW_H))) - pH0_NHE) / (l_NHE + abs((-log10(1000 * m_H_c / (V_c * MW_H))) - pH0_NHE))) * VMAXNHE * (S_C) * pow(m_H_C, a) / (pow(V_C * MW_H * K_mNHE / 1000, a) + pow(m_H_C, a))
    ],
    [
        k1 * MW_HCO3 / MW_CO2 * m_CO2_C - k2 * m_H_C * m_HCO3_C * 1000 / (V_C * MW_H)
        + SensATP * (0.5) * (1.0 + tanh(l_THCO3 * (-log10(1000 * m_H_c / (V_c * MW_H)) - pHe0_THCO3))) * (0.5) * (1.0 + tanh(g_THCO3 * (pHi0_THCO3 - (-log10(1000 * m_H_C / (V_C * MW_H)))))) * VMAXTHCO3 * (S_C) * m_HCO3_c / (V_c * K_mTHCO3 * MW_HCO3 / 1000 + m_HCO3_c)
    ]
])
M2 = Matrix([
    [X],
    [X],
    [X],
    [X],
    [m_CO2_C],
    [m_H_C],
    [m_HCO3_C]
])
J = M1.jacobian(M2)
print(f'jacobian={J}')


# ode_metabolism = Matrix([
#     [A * 140.0 - D * X * ((pow(S(0, 1), N) / (pow(S(0, 1), N) + pow(X, N))) + Gamma(0, 1) * (pow(X, N) / (pow(S(0, 1), N) + pow(X, N)))) ],
#     [A * ((pow(S(1, 3)(1, 2), N) / (pow(S(1, 3)(1, 2), N) + pow(X, N))) + Gamma(1, 3)(1, 2) * (pow(X, N) / (pow(S(1, 3)(1, 2), N) + pow(X, N)))) - D * X ],
#     [A * ((pow(S(3, 4)(1, 3), N) / (pow(S(3, 4)(1, 3), N) + pow(X, N))) + Gamma(3, 4)(1, 3) * (pow(X, N) / (pow(S(3, 4)(1, 3), N) + pow(X, N)))) - D * X ],
#     [A * ((pow(S(3, 4), N) / (pow(S(3, 4), N) + pow(X, N))) + Gamma(3, 4) * (pow(X, N) / (pow(S(3, 4), N) + pow(X, N)))) - D * X ]    
# ])
# jacobian = ode_metabolism.jacobian()

# print("Hello World !")
# print(f'A={A}, D={D}')
# print(f'S={S}')
# print(f'Gamma={Gamma}')
# print(f'Metabolism={ode_metabolism}')


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