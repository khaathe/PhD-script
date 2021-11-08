#!/usr/bin

from abc import ABC, abstractclassmethod
import numpy as np
from numpy.linalg import cond
from plotly.subplots import make_subplots
from scipy.integrate import solve_ivp
import plotly.graph_objects as go
import math

class Simu(ABC):
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
            "Kh" : 2.5e-4
        }
        self.pSimu = {
            "tspan" : [0.0, 1440.0],
            "dt" : 0.1,
            "subRows" : 4,
            "subCols" : 1
        }
        self.solutions = {}
        self.model = self.model
        self.modelName = ""
        self.description = ""

    @abstractclassmethod
    def model(self, t, x, p):
        pass

    def H(self,y,s,n,gamma):
        return s**n / ( s**n + y**n ) + gamma * y**n / (s**n + y**n)

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

    def plot_vars(self, condition, vars_index, vars_label):
        nbRows = len(vars_index)
        nbCols = 1
        sub = make_subplots(rows = nbRows, cols=nbCols)
        for i in range(nbRows):
            for j in range(nbCols):
                sub.append_trace(
                    go.Scatter(
                        x = self.solutions[condition].t, 
                        y = self.solutions[condition].y[ vars_index[i*nbCols+j] ], 
                        mode = "lines", 
                        name = vars_label[i*nbCols+j],
                        showlegend=True
                    ), 
                    row=i+1,
                    col=j+1
                )
        sub.update_layout(title_text= "{} - {}".format(self.modelName, self.description) )
        sub.show()

    def plot_all_condition_all_vars(self):
        nbRows = len(self.initialCondition)
        nbCols = 1
        solKeys = list(self.solutions.keys())
        sub = make_subplots(rows = nbRows, cols=nbCols, subplot_titles=solKeys)
        for i in range(nbRows):
            for j in range(nbCols):
                self.subplot_all_var(sub, i+1, j+1, self.solutions[ solKeys[i*nbCols+j] ], plot_legend=(i*nbCols+j == 0))
        sub.update_layout(title_text= "{} - {}".format(self.modelName, self.description) )
        sub.show(config = self.CONFIG)
    
    def subplot_all_var(self, sub, row_index, col_index, sol, plot_legend):
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
    
    def updateParamOde(self, newVal):
        self.pOde.update(newVal)

    def updateParamSimu(self, newVal):
        self.pSimu.update(newVal)

    def setInitialCondition(self, initialConditions):
        self.initialCondition = initialConditions
        self.pSimu["subRows"] = len(self.initialCondition)

class BaseModel(Simu):
    def __init__(self):
        super().__init__()
        self.modelName = "Base Model"
        self.description = ""
    
    def model(self, t, x, p):
        dxdt = [None] * 8
        hif, ldh, pdk, pdh, o, g, atp, h = x
        A, D, N, hif_coeff_prod, S, gamma, Vo, Ko, pg, A0,  Kg, Kh = p.values()
        Co = Vo * ( o/(o+Ko) )
        Cg = ( (pg*A0/2.0) - (27.0*Co/10.0) ) * ( g/(g+Kg) )
        Pa = ( (2.0 * Cg) + ( (27.0/5.0)*Co ) )
        Ph = Kh*( (29.0*(pg*Vo-Co))/5.0 )
        dxdt[0] = A * hif_coeff_prod - D * self.H(o, S[0][1], N, gamma[0][1]) * hif
        dxdt[1] = A * self.H(hif, S[1][2], N, gamma[1][2]) - D * ldh
        dxdt[2] = A * self.H(hif, S[1][3], N, gamma[1][3]) - D * pdk
        dxdt[3] = A * self.H(pdk, S[3][4], N, gamma[3][4]) - D * pdh
        dxdt[4] = -Co
        dxdt[5] = -Cg
        dxdt[6] = Pa
        dxdt[7] = Ph
        return dxdt

class NewModel(Simu):
    def __init__(self):
        super().__init__()
        self.pOde.update({
            "L" : 50,
            "beta" : 1,
            "k" : 10,
            "ldh0" : 2.176675
        })
        self.modelName = "New Model"
        self.description = ""
    
    def model(self, t, x, p):
        dxdt = [None] * 8
        hif, ldh, pdk, pdh, o, g, atp, h = x
        A, D, N, hif_coeff_prod, S, gamma, Vo, Ko, pg, A0,  Kg, Kh, L, beta, k, ldh0 = p.values()
        pg = (L-beta)/(1+math.exp(-k*(ldh-ldh0))) + beta
        Co = Vo * ( o/(o+Ko) )
        Cg = ( (pg*A0/2.0) - (27.0*Co/10.0) ) * ( g/(g+Kg) )
        Pa = ( (2.0 * Cg) + ( (27.0/5.0)*Co ) )
        Ph = Kh * ( (29.0*(pg*Vo-Co))/5.0 )
        dxdt[0] = A * hif_coeff_prod - D * self.H(o, S[0][1], N, gamma[0][1]) * hif
        dxdt[1] = A * self.H(hif, S[1][2], N, gamma[1][2]) - D * ldh
        dxdt[2] = A * self.H(hif, S[1][3], N, gamma[1][3]) - D * pdk
        dxdt[3] = A * self.H(pdk, S[3][4], N, gamma[3][4]) - D * pdh
        dxdt[4] = -Co
        dxdt[5] = -Cg
        dxdt[6] = Pa
        dxdt[7] = Ph
        return dxdt
        

class O2Decreasing(NewModel):
    def __init__(self):
        super().__init__()
        self.initialCondition = {
             "" : [ 3.565, 1.726, 3.31, 0.28, 0.056, 5.0, 0.0, 10**(-7.4)/1000, 0.0, 0.0]
        }
        self.pSimu.update({
            "subRows" : 1,
            "subCols" : 1
        })
        self.LABELS_PLOT.extend(["Oxygen Consumption", "Glucose Consumption"])
        self.COLOR_PLOT.extend(["#7FFFD4", "#FFE4C4"])
        self.modelName = "Model with Decreasing O2 over time"
        self.description = ""
    
    def model(self, t, x, p):
        dxdt = [None] * len(x)
        hif, ldh, pdk, pdh, o, g, atp, h, oConsumption, gConsumption = x
        A, D, N, hif_coeff_prod, S, gamma, Vo, Ko, pg, A0,  Kg, Kh, L, beta, k, ldh0 = p.values()
        
        pg = (L-beta)/(1+math.exp(-k*(ldh-ldh0))) + beta
        Co = Vo * ( o/(o+Ko) )
        Cg = ( (pg*A0/2.0) - (27.0*Co/10.0) ) * ( g/(g+Kg) )
        Pa = ( (2.0 * Cg) + ( (27.0/5.0)*Co ) )
        Ph = Kh * ( (29.0*(pg*Vo-Co))/5.0 )
        dxdt[0] = A * hif_coeff_prod - D * self.H(o, S[0][1], N, gamma[0][1]) * hif
        dxdt[1] = A * self.H(hif, S[1][2], N, gamma[1][2]) - D * ldh
        dxdt[2] = A * self.H(hif, S[1][3], N, gamma[1][3]) - D * pdk
        dxdt[3] = A * self.H(pdk, S[3][4], N, gamma[3][4]) - D * pdh
        dxdt[4] = -0.0009333333 if (t>(480-30) and t<(480+30) ) else 0.0
        dxdt[5] = 0.0
        dxdt[6] = Pa
        dxdt[7] = Ph
        dxdt[8] = Co
        dxdt[9] = Cg
        return dxdt

if __name__ == "__main__":
    # simulation = BaseModel()
    # simulation.run()
    # simulation.description = "Vo={}".format( simulation.pOde["Vo"] )
    # simulation.plot_all_condition_all_vars()

    # simulation = BaseModel()
    # simulation.setVo(1.16e-4)
    # simulation.description = "Vo={}".format( simulation.pOde["Vo"] )
    # simulation.run()
    # simulation.plot()
    
    # simulation = BaseModel()
    # simulation.setVo(1.16e-4)
    # simulation.updateParamOde({"pg" : 50})
    # simulation.description = "Vo = {}, pg = {}".format(simulation.pOde["Vo"], simulation.pOde["pg"])
    # simulation.run()
    # simulation.plot()

    # simulation = NewModel()
    # simulation.description = "Vo={}".format( simulation.pOde["Vo"] )
    # simulation.run()
    # simulation.plot_all_condition_all_vars()

    # simulation = NewModel()
    # simulation.setVo(1.16e-4)
    # simulation.description = "Vo={}".format( simulation.pOde["Vo"] )
    # simulation.run()
    # simulation.plot()

    # simulation = O2Decreasing()
    # simulation.description = "Vo={}".format( simulation.pOde["Vo"] )
    # simulation.run()
    # simulation.plot()

    simulation = NewModel()
    simulation.setInitialCondition({
        "Normoxia - Glucose = 5 mmol/L" : [ 3.565, 1.726, 3.31, 0.28, 0.056, 5.0, 0.0, 10**(-7.4)/1000],
        "Normoxia - Glucose = 1 mmol/L" : [ 3.565, 1.726, 3.31, 0.28,0.056, 1.0, 0.0, 10**(-7.4)/1000]
    })
    simulation.setVo(2.1e-11 * 5e6)
    simulation.description = "Vo={}".format( simulation.pOde["Vo"] )
    simulation.run()
    simulation.plot_vars("Normoxia - Glucose = 5 mmol/L", [0, 4, 5, 7], ["HIF", "Oxygen", "Glucose", "H+"])
    simulation.plot_vars("Normoxia - Glucose = 1 mmol/L", [0, 4, 5, 7], ["HIF", "Oxygen", "Glucose", "H+"])
