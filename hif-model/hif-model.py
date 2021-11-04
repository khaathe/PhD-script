#!/usr/bin

import numpy as np
from plotly.subplots import make_subplots
from scipy.integrate import solve_ivp
import plotly.graph_objects as go
import math

def H(y,s,n,gamma):
    return s**n / ( s**n + y**n ) + gamma * y**n / (s**n + y**n)

def base_model(t, x, p):
    dxdt = [None] * 8
    hif, ldh, pdk, pdh, o, g, atp, h = x
    A, D, N, hif_coeff_prod, S, gamma, Vo, Ko, pg, A0,  Kg, Kh, L, beta, k, ldh0 = p.values()
    Co = Vo * ( o/(o+Ko) )
    Cg = ( (pg*A0/2.0) - (27.0*Co/10.0) ) * ( g/(g+Kg) )
    Pa = ( (2.0 * Cg) + ( (27.0/5.0)*Co ) )
    Ph = Kh*( (29.0*(pg*Vo-Co))/5.0 )
    dxdt[0] = A * hif_coeff_prod - D * H(o, S[0][1], N, gamma[0][1]) * hif
    dxdt[1] = A * H(hif, S[1][2], N, gamma[1][2]) - D * ldh
    dxdt[2] = A * H(hif, S[1][3], N, gamma[1][3]) - D * pdk
    dxdt[3] = A * H(pdk, S[3][4], N, gamma[3][4]) - D * pdh
    dxdt[4] = -Co
    dxdt[5] = -Cg
    dxdt[6] = Pa
    dxdt[7] = Ph
    return dxdt

def new_model(t, x, p):
    dxdt = [None] * 8
    hif, ldh, pdk, pdh, o, g, atp, h = x
    A, D, N, hif_coeff_prod, S, gamma, Vo, Ko, pg, A0,  Kg, Kh, L, beta, k, ldh0 = p.values()
    
    pg = (L-beta)/(1+math.exp(-k*(ldh-ldh0))) + beta

    Co = Vo * ( o/(o+Ko) )
    Cg = ( (pg*A0/2.0) - (27.0*Co/10.0) ) * ( g/(g+Kg) )
    Pa = ( (2.0 * Cg) + ( (27.0/5.0)*Co ) )
    Ph = Kh * ( (29.0*(pg*Vo-Co))/5.0 )
    dxdt[0] = A * hif_coeff_prod - D * H(o, S[0][1], N, gamma[0][1]) * hif
    dxdt[1] = A * H(hif, S[1][2], N, gamma[1][2]) - D * ldh
    dxdt[2] = A * H(hif, S[1][3], N, gamma[1][3]) - D * pdk
    dxdt[3] = A * H(pdk, S[3][4], N, gamma[3][4]) - D * pdh
    dxdt[4] = -Co
    dxdt[5] = -Cg
    dxdt[6] = Pa
    dxdt[7] = Ph
    return dxdt

def o2_decreasing(t, x, p):
    dxdt = [None] * len(x)
    hif, ldh, pdk, pdh, o, g, atp, h, oConsumption, gConsumption = x
    A, D, N, hif_coeff_prod, S, gamma, Vo, Ko, pg, A0,  Kg, Kh, L, beta, k, ldh0 = p.values()
    
    pg = (L-beta)/(1+math.exp(-k*(ldh-ldh0))) + beta
    Co = Vo * ( o/(o+Ko) )
    Cg = ( (pg*A0/2.0) - (27.0*Co/10.0) ) * ( g/(g+Kg) )
    Pa = ( (2.0 * Cg) + ( (27.0/5.0)*Co ) )
    Ph = Kh * ( (29.0*(pg*Vo-Co))/5.0 )
    dxdt[0] = A * hif_coeff_prod - D * H(o, S[0][1], N, gamma[0][1]) * hif
    dxdt[1] = A * H(hif, S[1][2], N, gamma[1][2]) - D * ldh
    dxdt[2] = A * H(hif, S[1][3], N, gamma[1][3]) - D * pdk
    dxdt[3] = A * H(pdk, S[3][4], N, gamma[3][4]) - D * pdh
    dxdt[4] = -0.0009333333 if (t>(480-30) and t<(480+30) ) else 0.0
    dxdt[5] = 0.0
    dxdt[6] = Pa
    dxdt[7] = Ph
    dxdt[8] = Co
    dxdt[9] = Cg
    return dxdt
class Simu:
    NORMOXIA = 0.056 # ~5%
    HYPOXIA = 0.010564 # 1%
    ANOXIA = 1.0564e-05 # 0.001%
    GLUCOSE_NORMAL = 5.0
    GLUCOSE_LOW = 1.0
    GLUCOSE_VERY_LOW = 0.1

    LABELS_PLOT = ["HIF", "LDH", "PDK", "PDH", "Oxygen", "Glucose", "ATP", "H+"]
    COLOR_PLOT = ["blue", "green", "goldenrod", "grey", "red", "purple", "cyan", "magenta"]
    CONFIG = {'displaylogo': False}

    def __init__(self):
        self.initialCondition = {
            "normoxia+normal_glucose" : [ 3.565, 1.726, 3.31, 0.28, self.NORMOXIA, self.GLUCOSE_NORMAL, 0.0, 10**(-7.4)/1000],
            "hypoxia+normal_glucose" : [ 3.565, 1.726, 3.31, 0.28, self.HYPOXIA, self.GLUCOSE_NORMAL, 0.0, 10**(-7.4)/1000],
            "normoxia+low_glucose" : [ 3.565, 1.726, 3.31, 0.28, self.NORMOXIA, self.GLUCOSE_LOW, 0.0, 10**(-7.4)/1000],
            "hypoxia+low_glucose" : [ 3.565, 1.726, 3.31, 0.28, self.HYPOXIA, self.GLUCOSE_LOW, 0.0, 10**(-7.4)/1000]
        }
        self.pOde = {
            "A" : 0.005,
            "D" : 0.005,
            "N" : 4.0,
            "hif_coeff_prod" : 140.0,
            "S" : [ 
                [0.0,20.85e-3,0.0,0.0,0.0],
                [0.0,0.0,4.64,4.0,0.0],
                [0.0,0.0,0.0,0.0,0.0],
                [0.0,0.0,0.0,0.0,2.2],
                [0.0,0.0,0.0,0.0,0.0]
            ],
            "gamma" : [
                [0.0,40.0,0.0,0.0,0.0],
                [0.0,0.0,3.81,6.97,0.0],
                [0.0,0.0,0.0,0.0,0.0],
                [0.0,0.0,0.0,0.0,0.14],
                [0.0,0.0,0.0,0.0,0.0]
            ],
            "Vo" : 2.1e-11,
            "Ko" : 0.005,
            "pg" : 1.0,
            "A0" : 29.0/5.0 * 2.1e-11,
            "Kg" : 0.04,
            "Kh" : 2.5e-4,
            "L" : 50,
            "beta" : 1,
            "k" : 10,
            "ldh0" : 2.176675
        }
        self.pSimu = {
            "tspan" : [0.0, 1440.0],
            "dt" : 0.1,
            "subRows" : 4,
            "subCols" : 1
        }
        self.solutions = {}
        self.model = base_model
        self.modelName = "Base Model"
        self.description = ""

    def run(self):
        self.solutions= dict()
        tspan = self.pSimu["tspan"]
        dt = self.pSimu["dt"]
        t = np.arange(tspan[0], tspan[1], dt)
        for condition in self.initialCondition:
            self.solutions[condition] = solve_ivp(self.model, tspan, self.initialCondition[condition], t_eval= t, method='DOP853', args=[self.pOde], rtol = 1e-6, atol=1e-10)

    def setVo(self, Vo):
        self.pOde["Vo"] = Vo
        self.pOde["A0"] = 29.0/5.0 * Vo

    def plot(self):
        nbRows = self.pSimu["subRows"]
        nbCols = self.pSimu["subCols"]
        solKeys = list(self.solutions.keys())
        sub = make_subplots(rows = nbRows, cols=nbCols, subplot_titles=solKeys)
        for i in range(nbRows):
            for j in range(nbCols):
                self.subplot_solution(sub, i+1, j+1, self.solutions[ solKeys[i*nbCols+j] ], plot_legend=(i*nbCols+j == 0))
        sub.update_layout(title_text= "{} - {}".format(self.modelName, self.description) )
        sub.show(config = self.CONFIG)
        pass
    
    def subplot_solution(self, sub, row_index, col_index, sol, plot_legend):
        for i in range(len(sol.y)):
            sub.append_trace( 
                go.Scatter(
                    x = sol.t, 
                    y = sol.y[i], 
                    mode = "lines", 
                    name = self.LABELS_PLOT[i], 
                    legendgroup='group_{}'.format(self.LABELS_PLOT[i]), 
                    showlegend=plot_legend,
                    line_color = self.COLOR_PLOT[i]
                ), 
                row=row_index,
                col=col_index
            )
        sub.append_trace(
            go.Scatter(
                x = sol.t,
                y = -np.log10(sol.y[7]*1000),
                mode = "lines",
                name = 'pH',
                legendgroup='group_pH',
                showlegend=plot_legend,
                line_color = self.COLOR_PLOT[i]
            ), 
            row=row_index,
            col=col_index
        )
        return sub

    def setModel(self, modelType):
        if ( modelType.lower() == "base_model" ):
            self.model = base_model
            self.modelName = "Base Model"
        elif ( modelType.lower() == "new_model" ):
            self.model = new_model
            self.modelName = "New Model"
        elif ( modelType.lower() == "o2_decreasing" ):
            self.model = o2_decreasing
            self.modelName = "Model with O2 decreasing"
        else:
            self.model = base_model
            self.modelName = "Base Model"
    
    def updateParamOde(self, newVal):
        self.pOde.update(newVal)

    def updateParamSimu(self, newVal):
        self.pSimu.update(newVal)

    def setInitialCondition(self, initialConditions):
        self.initialCondition = initialConditions
        self.pSimu["subRows"] = len(self.initialCondition)

if __name__ == "__main__":
    # simulation = Simu()
    # simulation.run()
    # simulation.description = "Vo={}".format( simulation.pOde["Vo"] )
    # simulation.plot()

    # simulation = Simu()
    # simulation.setVo(1.16e-4)
    # simulation.description = "Vo={}".format( simulation.pOde["Vo"] )
    # simulation.run()
    # simulation.plot()
    
    # simulation = Simu()
    # simulation.setVo(1.16e-4)
    # simulation.updateParamOde({"pg" : 50})
    # simulation.description = "Vo = {}, pg = {}".format(simulation.pOde["Vo"], simulation.pOde["pg"])
    # simulation.run()
    # simulation.plot()

    # simulation = Simu()
    # simulation.setModel("new_model")
    # simulation.description = "Vo={}".format( simulation.pOde["Vo"] )
    # simulation.run()
    # simulation.plot()

    # simulation = Simu()
    # simulation.setModel("new_model")
    # simulation.setVo(1.16e-4)
    # simulation.description = "Vo={}".format( simulation.pOde["Vo"] )
    # simulation.run()
    # simulation.plot()

    simulation = Simu()
    simulation.setInitialCondition({
        "" : [ 3.565, 1.726, 3.31, 0.28, 0.056, 5.0, 0.0, 10**(-7.4)/1000, 0.0, 0.0]
    })
    simulation.LABELS_PLOT = ["HIF", "LDH", "PDK", "PDH", "Oxygen", "Glucose", "ATP", "H+", "Oxygen Consumption", "Glucose Consumption"]
    simulation.COLOR_PLOT.append("#7FFFD4")
    simulation.COLOR_PLOT.append("#FFE4C4")
    simulation.description = "Vo={}".format( simulation.pOde["Vo"] )
    simulation.setModel("o2_decreasing")
    simulation.run()
    simulation.plot()