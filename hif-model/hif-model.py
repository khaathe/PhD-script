#!/usr/bin/python3

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
from scipy import optimize

# Example of initial Condition
# Order -> [HIF, LDH, PDK, PDH, Oxygen, Glucose, ATP, H+, Oxygen Consumed, Glucose Consumed]
# {
#             "normoxia+normal_glucose" : [ 3.565, 1.726, 3.31, 0.28, 0.056, 5.0, 0.0, 10**(-7.4)/1000, 0.0, 0.0],
#             "hypoxia+normal_glucose" : [ 3.565, 1.726, 3.31, 0.28, 0.010564, 5.0, 0.0, 10**(-7.4)/1000, 0.0, 0.0],
#             "normoxia+low_glucose" : [ 3.565, 1.726, 3.31, 0.28, 0.056, 1.0, 0.0, 10**(-7.4)/1000, 0.0, 0.0],
#             "hypoxia+low_glucose" : [ 3.565, 1.726, 3.31, 0.28, 0.010564, 1.0, 0.0, 10**(-7.4)/1000, 0.0, 0.0]
# }

class Simu(ABC):
    NORMOXIA = 0.056 # ~5%
    HYPOXIA = 0.010564 # 1%
    ANOXIA = 1.0564e-05 # 0.001%
    GLUCOSE_NORMAL = 5.0
    GLUCOSE_LOW = 1.0
    GLUCOSE_VERY_LOW = 0.1

    DF_COLUMN_NAMES = ["HIF", "LDH", "PDK", "PDH", "Oxygen Consumed", "Glucose Consumed", "ATP", "H+"]
    LABELS_PLOT = ["HIF", "LDH", "PDK", "PDH", "Oxygen", "Glucose", "ATP", "H+"]
    COLOR_PLOT = ["blue", "green", "goldenrod", "grey", "red", "purple", "cyan", "magenta"]
    CONFIG = {'displaylogo': False}

    def __init__(self):
        self.initialCondition = None
        self.pOde = {
            "A" : [0.7, 0.005, 0.005, 0.005, 0.005,], # We multiply HIF prod rate by 140 or we should change gamma_O2-> as well
            "D" : [0.005, 0.005, 0.005, 0.005, 0.005,],
            "N" : 4.0,
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
            "subCols" : 1,
            "oExtraInitial" : 0.056,
            "gExtraInitial" : 5.0,
            "protonExtraInitial" : 10**(-7.4 + 3)
        }
        self.solutions = {}
        self.model = self.model
        self.modelName = ""
        self.description = ""

    @abstractclassmethod
    def model(self, t, x, p):
        pass

    def getO2Extra(self, t):
        return 0.056

    def getGlucoseExtra(self, t):
        return 5.0

    def H(self,y,s,n,gamma):
        return s**n / ( s**n + y**n ) + gamma * y**n / (s**n + y**n)

    def calculate_rates(self, substrate, dt):
        rates = [0.0] * len(substrate)
        for i in np.arange(1, len(substrate)):
            rates[i] = abs( substrate[i] - substrate[i-1] ) / dt
        return rates

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
            data["Oxygen Extracellular"] = [self.getO2Extra(t) for t in sol.t ]
            data["Glucose Extracellular"] = [self.getGlucoseExtra(t) for t in sol.t ]
            data["H+ Extracellular"] = [ self.pSimu["protonExtraInitial"] + x for x in sol.y[7] ]
            data["pH"] = [ -math.log10(1e-3*x) for x in data["H+ Extracellular"]]
            df = pd.DataFrame(data)
            self.solutions[condition] = df

    def setVo(self, Vo):
        self.pOde["Vo"] = Vo
        self.pOde["A0"] = 29.0/5.0 * Vo

    def plot(self, condition, title, x_title="Time (min)", y_title="Value over Time"):
        data = self.solutions[condition]
        fig = px.line(data, x="Time", y=data.columns[1:data.shape[1]]) 
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

    def plot_var_across_conditions(self, var, title, x_title="Time (min)", y_title="Value over Time"):
        i = 1
        sub = make_subplots(rows = len(self.solutions), cols = 1, subplot_titles = list(self.solutions.keys()) )
        for k in self.solutions.keys():
            sub.append_trace( go.Scatter(x=self.solutions[k]["Time"], y=self.solutions[k][var], mode='lines', showlegend=False), row=i, col=1)
            i = i + 1
        sub.update_layout(title = title)
        sub.update_xaxes(title_text = x_title)
        sub.update_yaxes(title_text = y_title, showexponent = 'all', exponentformat = 'e')
        sub.show()
        return sub

    def plot_vars(self, condition, vars, title, x_title="Time (min)", y_title="Value over Time"):
        i = 1
        sub = make_subplots(rows = len(vars), cols = 1, subplot_titles = vars )
        for v in vars:
            sub.append_trace( go.Scatter(x=self.solutions[condition]["Time"], y=self.solutions[condition][v], mode='lines', showlegend=False), row=i, col=1)
            i = i + 1
        sub.update_layout(title = title)
        sub.update_xaxes(title_text = x_title)
        sub.update_yaxes(title_text = y_title, showexponent = 'all', exponentformat = 'e')
        sub.show()
        return sub

    def setInitialCondition(self, initialConditions):
        self.initialCondition = initialConditions
        self.pSimu["subRows"] = len(self.initialCondition)

    def stable_state(self, title, xaxis, yaxis, t=0.0, nbInterval = 30):
        wrapper = lambda x: self.model(t, x, self.pOde)
        x0 = np.linspace(0,1,nbInterval).reshape(nbInterval,1)
        max = np.array([140.0, 10.0, 10.0, 1.0, 1.0, 5.0, 10.0, 10.0])
        min = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        x0 = (max - min) * x0 + min
        sol = [optimize.root(wrapper, x, method='lm') for x in x0]
        xplot = [ s.x[1] for s in sol if s.success]
        yplot = [ s.x[3] for s in sol if s.success]
        fig = go.Figure(data=go.Scatter(
            x=xplot,
            y=yplot,
            mode='markers'
        ))
        for x in x0:
            fig.add_vline(x=x[1],line_dash="dash", line_width=1, opacity=0.5)
            fig.add_hline(y=x[3], line_dash="dash", line_width=1, opacity=0.5)
        fig.update_layout(
            title=title,
            xaxis_title = xaxis,
            yaxis_title = yaxis
        )
        fig.show()
        return sol

class BaseModel(Simu):
    def __init__(self):
        super().__init__()
        self.modelName = "Base Model"
        self.description = ""
    
    def model(self, t, x, p):
        dxdt = [None] * len(x)
        hif, ldh, pdk, pdh, oConsumption, gConsumption, atp, h  = x
        A, D, N, S, gamma, Vo, Ko, pg, A0,  Kg, Kh = p.values()
        oExtra = self.getO2Extra(t)
        gExtra = self.getGlucoseExtra(t)
        Co = Vo * ( oExtra/(oExtra+Ko) )
        Cg = ( (pg*A0/2.0) - (27.0*Co/10.0) ) * ( gExtra/(gExtra+Kg) )
        Pa = ( (2.0 * Cg) + ( (27.0/5.0)*Co ) )
        Ph = Kh*( (29.0*(pg*Vo-Co))/5.0 )
        dxdt[0] = A[0] - D[0] * self.H(oExtra, S[0][1], N, gamma[0][1]) * hif
        dxdt[1] = A[1] * self.H(hif, S[1][2], N, gamma[1][2]) - D[1] * ldh
        dxdt[2] = A[2] * self.H(hif, S[1][3], N, gamma[1][3]) - D[2] * pdk
        dxdt[3] = A[3] * self.H(pdk, S[3][4], N, gamma[3][4]) - D[3] * pdh
        dxdt[4] = Co
        dxdt[5] = Cg
        dxdt[6] = Pa
        dxdt[7] = Ph
        return dxdt

class NewModel(Simu):
    def __init__(self):
        super().__init__()
        self.pOde.update({
            "pg_max" : 50,
            "pg_min" : 1,
            "k" : 100,
            "ldh0" : 2.176675,
            "po_max" : 1,
            "po_min" : 0.0,
            "l" : 100,
            "pdh0" : 0.21245
        })
        self.modelName = "New Model"
        self.description = ""
    
    def model(self, t, x, p):
        dxdt = [None] * len(x)
        hif, ldh, pdk, pdh, oConsumption, gConsumption, atp, h = x
        A, D, N, S, gamma, Vo, Ko, pg, A0,  Kg, Kh, pg_max, pg_min, k, ldh0, po_max, po_min, l, pdh0 = p.values()
        oExtra = self.getO2Extra(t)
        gExtra = self.getGlucoseExtra(t)
        po = (po_max-po_min)/(1+math.exp(-l*(pdh-pdh0))) + po_min
        pg = (pg_max-pg_min)/(1+math.exp(-k*(ldh-ldh0))) + pg_min
        Co = po * Vo * ( oExtra/(oExtra+Ko) )
        Cg = ( (pg*A0/2.0) - (27.0*Co/10.0) ) * ( gExtra/(gExtra+Kg) )
        Pa = ( (2.0 * Cg) + ( (27.0/5.0)*Co ) )
        Ph = Kh * ( (29.0*(pg*Vo-Co))/5.0 )
        dxdt[0] = A[0] - D[0] * self.H(oExtra, S[0][1], N, gamma[0][1]) * hif
        dxdt[1] = A[1] * self.H(hif, S[1][2], N, gamma[1][2]) - D[1] * ldh
        dxdt[2] = A[2] * self.H(hif, S[1][3], N, gamma[1][3]) - D[2] * pdk
        dxdt[3] = A[3] * self.H(pdk, S[3][4], N, gamma[3][4]) - D[3] * pdh
        dxdt[4] = Co
        dxdt[5] = Cg
        dxdt[6] = Pa
        dxdt[7] = Ph
        return dxdt
    
class NewModelWithO2Decreasing(NewModel):

    def __init__(self):
        NewModel.__init__(self)
        self.modelName = "New Model with Varying O2 over time"
        self.description = ""
    
    def getO2Extra(self, t):
        return 0.056/(1+math.exp(0.3*(t-480))) + 0.056/(1+math.exp(-0.3*(t-960)))

class BaseModelWithO2Decreasing(BaseModel):
    def __init__(self):
        BaseModel.__init__(self)
        self.modelName = "Base Model with Varying O2 over time"

    def getO2Extra(self, t):
        return 0.056/(1+math.exp(0.3*(t-480))) + 0.056/(1+math.exp(-0.3*(t-960)))

def run_base_model():
    simulation = BaseModel()
    simulation.setInitialCondition({
            "proliferation" : [ 3.565, 1.726, 3.31, 0.28, 0.056, 5.0, 0.0, 10**(-7.4)/1000]
    })
    simulation.run()
    vars = ["Oxygen Consumption Rate", "Glucose Consumption Rate", "H+ Production Rate", "ATP Production Rate"]
    simulation.plot_vars("proliferation", vars, "Base Model - Proliferation conditions (no nutrient absence)")

def run_new_model():
    simulation = NewModel()
    simulation.setInitialCondition({
            "proliferation" : [ 3.565, 1.726, 3.31, 0.28, 0.056, 5.0, 0.0, 10**(-7.4)/1000]
    })
    simulation.run()
    vars = ["Oxygen Consumption Rate", "Glucose Consumption Rate", "H+ Production Rate", "ATP Production Rate"]
    simulation.plot_vars("proliferation", vars, "New Model - Proliferation conditions (no nutrient absence)")

def run_decreasing_o2_base_model():
    simulation = BaseModelWithO2Decreasing()
    simulation.setInitialCondition({
            "normoxia+normal_glucose" : [ 3.565, 1.726, 3.31, 0.28, 0.056, 5.0, 0.0, 10**(-7.4)/1000]
    })
    simulation.run()
    vars_rates = ["Oxygen Extracellular", "Oxygen Consumption Rate", "Glucose Consumption Rate", "H+ Production Rate", "ATP Production Rate"]
    vars_gene = ["Oxygen Extracellular", "HIF", "LDH", "PDK", "PDH"]
    vars_molecule = ["Oxygen Extracellular", "Glucose Extracellular", "H+ Extracellular", "pH"]
    simulation.plot_vars("normoxia+normal_glucose", vars_rates, "New Model With O2 Decreasing - Reaction Rates")
    simulation.plot_vars("normoxia+normal_glucose", vars_gene, "New Model With O2 Decreasing - Genes level")
    simulation.plot_vars("normoxia+normal_glucose", vars_molecule, "New Model With O2 Decreasing - Nutrient and H+")

def run_decreasing_o2_new_model():
    simulation = NewModelWithO2Decreasing()
    simulation.pOde.update({"k":100, "l":100})
    simulation.setInitialCondition({
            "normoxia+normal_glucose" : [ 3.565, 1.726, 3.31, 0.28, 0.056, 5.0, 0.0, 10**(-7.4)/1000]
    })
    simulation.run()
    simulation.run()
    vars_rates = ["Oxygen Extracellular", "Oxygen Consumption Rate", "Glucose Consumption Rate", "H+ Production Rate", "ATP Production Rate"]
    vars_gene = ["Oxygen Extracellular", "HIF", "LDH", "PDK", "PDH"]
    vars_molecule = ["Oxygen Extracellular", "Glucose Extracellular", "H+ Extracellular", "pH"]
    simulation.plot_vars("normoxia+normal_glucose", vars_rates, "New Model With O2 Decreasing - Reaction Rates")
    simulation.plot_vars("normoxia+normal_glucose", vars_gene, "New Model With O2 Decreasing - Genes level")
    simulation.plot_vars("normoxia+normal_glucose", vars_molecule, "New Model With O2 Decreasing - Nutrient and H+")

def run_physicell_tumor_conditions():
    simulation = BaseModel()
    Vo = 0.01875
    simulation.pOde.update({
            "Vo" : Vo,
            "Ko" : 0.005,
            "pg" : 1.0,
            "A0" : 29.0/5.0 * Vo,
            "Kg" : 0.04,
            "Kh" : 1.0e-5,
            "pg" : 1
    })
    simulation.pSimu.update({"tspan" : [0.0, 1440.0], "dt" : 0.1})
    simulation.setInitialCondition({"PhysiCell" : [3.565, 1.726, 3.31, 0.28, 0.056, 5.0, 0.0, 3.981072e-5, 0.0, 0.0] })
    simulation.run()
    vars = ["Oxygen", "Glucose", "H+", "ATP", "pH"]
    title = '{}, pg = {}'.format(simulation.modelName, simulation.pOde["pg"])
    simulation.plot_vars("PhysiCell", vars, "{} - Substrate concentration".format(title))
    vars = ["Oxygen Consumption Rate", "Glucose Consumption Rate", "H+ Production Rate", "ATP Production Rate"]
    simulation.plot_vars("PhysiCell", vars, "{} Rates".format(title))

if __name__ == "__main__":
    run_base_model()
    run_new_model()
    run_decreasing_o2_base_model()
    run_decreasing_o2_new_model()
    run_physicell_tumor_conditions()