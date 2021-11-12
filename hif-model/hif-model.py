#!/usr/bin

from abc import ABC, abstractclassmethod
from operator import mod
import numpy as np
from numpy.linalg import cond
from plotly.subplots import make_subplots
from scipy.integrate import solve_ivp
import plotly.graph_objects as go
import plotly.express as px
import pandas as pd
import math

class Simu(ABC):
    NORMOXIA = 0.056 # ~5%
    HYPOXIA = 0.010564 # 1%
    ANOXIA = 1.0564e-05 # 0.001%
    GLUCOSE_NORMAL = 5.0
    GLUCOSE_LOW = 1.0
    GLUCOSE_VERY_LOW = 0.1

    DF_COLUMN_NAMES = ["HIF", "LDH", "PDK", "PDH", "Oxygen", "Glucose", "ATP", "H+"]
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
            sol = solve_ivp(self.model, tspan, self.initialCondition[condition], t_eval= t, method='DOP853', args=[self.pOde], rtol = 1e-6, atol=1e-10)
            data = { "Time" : sol.t}
            for i in range(len(self.DF_COLUMN_NAMES)):
                data[ self.DF_COLUMN_NAMES[i] ] = sol.y[i]
            df = pd.DataFrame(data)
            self.solutions[condition] = df

    def setVo(self, Vo):
        self.pOde["Vo"] = Vo
        self.pOde["A0"] = 29.0/5.0 * Vo

    def plot(self, condition, title, x_title="Time (min)", y_title="Value over Time"):
        data = self.solutions[condition]
        fig = px.scatter(data, x="Time", y=data.columns[1:data.shape[1]]) 
        fig.update_layout(
            title = title,
            xaxis_title = x_title,
            yaxis={
                "title" : y_title,
                "showexponent" : 'all',
                "exponentformat" : 'e'
            }
        )
        fig.show()
        return fig

    def plot_var_across_conditions(self, var, title):
        pass

    def plot_vars(self, condition, vars_index, vars_label):
        nbRows = len(vars_index)
        nbCols = 1
        sub = make_subplots(rows = nbRows, cols=nbCols, vertical_spacing=0.1)
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
        sub.update_layout(title_text= "{} - {}".format(self.modelName, self.description))
        sub.update_xaxes({ "title" : "Time (min)" })
        sub.update_yaxes({
            "showexponent" : 'all',
            "exponentformat" : 'e'
        })
        sub.update_yaxes(title_text = "Concentration", row=1, col=1)
        sub.update_yaxes(title_text = "Gene Level", row=2, col=1)
        sub.update_yaxes(title_text = "Consumption Rate", row=3, col=1)
        sub.update_yaxes(title_text = "Consumption Rate", row=4, col=1)
        sub.update_yaxes(title_text = "Production Rate", row=5, col=1)
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
            "pg_max" : 50,
            "pg_min" : 1,
            "k" : 10,
            "ldh0" : 2.176675,
            "po_max" : 1,
            "po_min" : 0.0,
            "l" : 30,
            "pdh0" : 0.21245
        })
        self.modelName = "New Model"
        self.description = ""
    
    def model(self, t, x, p):
        dxdt = [None] * 8
        hif, ldh, pdk, pdh, o, g, atp, h = x
        A, D, N, hif_coeff_prod, S, gamma, Vo, Ko, pg, A0,  Kg, Kh, pg_max, pg_min, k, ldh0, po_max, po_min, l, pdh0 = p.values()
        po = (po_max-po_min)/(1+math.exp(-l*(pdh-pdh0))) + po_min
        pg = (pg_max-pg_min)/(1+math.exp(-k*(ldh-ldh0))) + pg_min
        Co = po * Vo * ( o/(o+Ko) )
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
        
class DecreasingO2(Simu):
    DF_COLUMN_NAMES = ["HIF", "LDH", "PDK", "PDH", "Oxygen", "Glucose", "ATP", "H+", "Oxygen Consumed", "Glucose Consumed"]
    
    def __init__(self):
        super().__init__()

    def calculate_rates(self, substrate, dt):
        rates = [0.0] * len(substrate)
        for i in np.arange(1, len(substrate)):
            rates[i] = ( substrate[i] - substrate[i-1] ) / dt
        return rates

    def dO2Dt(self, t, o):
        dxdt = 0.0 
        if (t>(480-30) and t<(480+30) and o > Simu.ANOXIA):
            dxdt = -0.0009333333 
        elif (t>(960-30) and t<(960+30) and o < self.NORMOXIA):
            dxdt = 0.0009333333
        return dxdt
    
    def run(self):
        self.solutions= dict()
        tspan = self.pSimu["tspan"]
        dt = self.pSimu["dt"]
        t = np.arange(tspan[0], tspan[1], dt)
        for condition in self.initialCondition:
            sol = solve_ivp(self.model, tspan, self.initialCondition[condition], t_eval= t, method='DOP853', args=[self.pOde], rtol = 1e-6, atol=1e-10)
            data = { "Time" : sol.t}
            for i in range(len(self.DF_COLUMN_NAMES)):
                data[ self.DF_COLUMN_NAMES[i] ] = sol.y[i]
            data["Oxygen Consumption Rate"] = self.calculate_rates(data["Oxygen Consumed"], self.pSimu["dt"])
            data["Glucose Consumption Rate"] = self.calculate_rates(data["Glucose Consumed"], self.pSimu["dt"])
            data["H+ Production Rate"] = self.calculate_rates(data["H+"], self.pSimu["dt"])
            data["ATP Production Rate"] = self.calculate_rates(data["ATP"], self.pSimu["dt"])
            data["pH"] = [ -math.log10(1000*x) for x in data["H+"]]
            df = pd.DataFrame(data)
            self.solutions[condition] = df
    

class NewModelWithO2Decreasing(DecreasingO2, NewModel):

    def __init__(self):
        DecreasingO2.__init__(self)
        NewModel.__init__(self)
        self.initialCondition = {
             "" : [ 3.565, 1.726, 3.31, 0.28, 0.056, 5.0, 0.0, 10**(-7.4)/1000, 0.0, 0.0]
        }
        self.pSimu.update({
            "subRows" : 1,
            "subCols" : 1
        })
        self.LABELS_PLOT.extend(["Oxygen Consumption Rate (mmol/L/min)", "Glucose Consumption Rate (mmol/L/min)", "Protons Productions Rate (mmol/L/min)"])
        self.COLOR_PLOT.extend(["#7FFFD4", "#FFE4C4", "#F25A07"])
        self.modelName = "New Model with Varying O2 over time"
        self.description = ""
    
    def model(self, t, x, p):
        dxdt = [None] * len(x)
        hif, ldh, pdk, pdh, o, g, atp, h, oConsumption, gConsumption = x
        A, D, N, hif_coeff_prod, S, gamma, Vo, Ko, pg, A0,  Kg, Kh, pg_max, pg_min, k, ldh0, po_max, po_min, l, pdh0 = p.values()
        po = (po_max-po_min)/(1+math.exp(-l*(pdh-pdh0))) + po_min
        pg = (pg_max-pg_min)/(1+math.exp(-k*(ldh-ldh0))) + pg_min
        Co = po * Vo * ( o/(o+Ko) )
        Cg = ( (pg*A0/2.0) - (27.0*Co/10.0) ) * ( g/(g+Kg) )
        Pa = ( (2.0 * Cg) + ( (27.0/5.0)*Co ) )
        Ph = Kh * ( (29.0*(pg*Vo-Co))/5.0 )
        dxdt[0] = A * hif_coeff_prod - D * self.H(o, S[0][1], N, gamma[0][1]) * hif
        dxdt[1] = A * self.H(hif, S[1][2], N, gamma[1][2]) - D * ldh
        dxdt[2] = A * self.H(hif, S[1][3], N, gamma[1][3]) - D * pdk
        dxdt[3] = A * self.H(pdk, S[3][4], N, gamma[3][4]) - D * pdh
        dxdt[4] = self.dO2Dt(t,o)
        dxdt[5] = 0.0
        dxdt[6] = Pa
        dxdt[7] = Ph
        dxdt[8] = Co
        dxdt[9] = Cg
        return dxdt

class BaseModelWithO2Decreasing(DecreasingO2,BaseModel):
    def __init__(self):
        DecreasingO2.__init__(self)
        BaseModel.__init__(self)
        self.initialCondition = {
             "" : [ 3.565, 1.726, 3.31, 0.28, 0.056, 5.0, 0.0, 10**(-7.4)/1000, 0.0, 0.0]
        }
        self.pSimu.update({
            "subRows" : 1,
            "subCols" : 1
        })
        self.LABELS_PLOT.extend(["Oxygen Consumption", "Glucose Consumption"])
        self.COLOR_PLOT.extend(["#7FFFD4", "#FFE4C4"])
        self.modelName = "Base Model with Varying O2 over time"
        self.description = ""
    
    def model(self, t, x, p):
        dxdt = [None] * len(x)
        hif, ldh, pdk, pdh, o, g, atp, h, oConsumption, gConsumption = x
        A, D, N, hif_coeff_prod, S, gamma, Vo, Ko, pg, A0,  Kg, Kh = p.values()
        Co = Vo * ( o/(o+Ko) )
        Cg = ( (pg*A0/2.0) - (27.0*Co/10.0) ) * ( g/(g+Kg) )
        Pa = ( (2.0 * Cg) + ( (27.0/5.0)*Co ) )
        Ph = Kh*( (29.0*(pg*Vo-Co))/5.0 )
        dxdt[0] = A * hif_coeff_prod - D * self.H(o, S[0][1], N, gamma[0][1]) * hif
        dxdt[1] = A * self.H(hif, S[1][2], N, gamma[1][2]) - D * ldh
        dxdt[2] = A * self.H(hif, S[1][3], N, gamma[1][3]) - D * pdk
        dxdt[3] = A * self.H(pdk, S[3][4], N, gamma[3][4]) - D * pdh
        dxdt[4] = self.dO2Dt(t, o)
        dxdt[5] = 0.0
        dxdt[6] = Pa
        dxdt[7] = Ph
        dxdt[8] = Co
        dxdt[9] = Cg
        return dxdt


if __name__ == "__main__":
    simulation = NewModelWithO2Decreasing()
    simulation.run()
    simulation.plot_vars_V2("", ["HIF", "Oxygen"], "test")
    
    # simulation = NewModelWithO2Decreasing()
    # simulation.setVo(0.012*60)
    # simulation.run()
    # simulation.modelName = "Preliminary Result"
    
    # simulation.description = "Impact of Oxygen on Gene Levels and Metabolism according to HIF levels"
    # simulation.plot_vars("", [4, 0, 1, 2, 3, 8, 9, 10], ["Oxygen", "HIF", "LDH", "PDK", "PDH", "Oxygen Consumption Rate (mmol/L/min)", "Glucose Consumption Rate (mmol/L/min)", "Protons Productions Rate (mmol/L/min)"])

    # simulation.description = "Impact of Oxygen on Gene Levels"
    # simulation.plot_vars("", [4, 0, 1, 2, 3], ["Oxygen", "HIF", "LDH", "PDK", "PDH"])


    # simulation.description = "Impact of HIF on Oxygen and Glucose consumption rate and H+ production rate"
    # simulation.plot_vars("", [4, 0, 8, 9, 10], ["Oxygen (mmol/L)", "HIF (Adimentional)", "Oxygen Consumption Rate (mmol/L/min)", "Glucose Consumption Rate (mmol/L/min)", "Protons Productions Rate (mmol/L/min)"])

    # simulation = NewModel()
    # simulation.run()
    # simulation.plot_all_condition_all_vars()
    # data = simulation.solutions[""]
    # sub = make_subplots(rows = 8, cols=1)
    # sub.append_trace(
    #     go.Scatter(
    #     x = data.t, 
    #     y = data.y[4], 
    #     mode = "lines", 
    #     name = "Oxygen",
    #     showlegend=True        
    # ), 
    #     row=1,
    #     col=1
    # )
    # sub.update_layout(title_text= "Impact of Oxygen on Gene Levels and Metabolism according to HIF levels")
    # sub.update_xaxes({
    #     "title" : 'Time (min)'
    # })
    # sub.update_yaxes({
    #     "showexponent" : 'all',
    #     "exponentformat" : 'e'
    # })
    # sub.show()