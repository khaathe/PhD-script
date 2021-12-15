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

hypoxic_period = lambda t : 0.056/(1+math.exp(0.3*(t-480))) + 0.056/(1+math.exp(-0.3*(t-960)))
hypoxia_at_8h = lambda t : 0.056/(1+math.exp(0.3*(t-480)))
hypoglycemia_period = lambda t : 5.0/(1+math.exp(0.3*(t-480))) + 5.0/(1+math.exp(-0.3*(t-960)))

class Simu(ABC):
    NORMOXIA = 0.056 # ~5%
    HYPOXIA = 0.010564 # 1%
    ANOXIA = 1.0564e-05 # 0.001%
    GLUCOSE_NORMAL = 5.0
    GLUCOSE_LOW = 1.0
    GLUCOSE_VERY_LOW = 0.1

    CONFIG = {'displaylogo': False}
    VARS_GENES = ["HIF", "LDH", "PDK", "PDH"]
    VARS_METABOLISM = ["Oxygen Consumed", "Glucose Consumed", "ATP", "H+"]
    VARS_RATES = ["Oxygen Consumption Rate", "Glucose Consumption Rate", "H+ Production Rate", "ATP Production Rate"]
    VARS_MICROENVIRONMENT = ["Oxygen Extracellular", "Glucose Extracellular", "H+ Extracellular", "pH"]
    

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
            "Vo" : 0.01875,
            "Ko" : 0.005,
            "pg" : 1.0,
            "A0" : 0.10875,
            "Kg" : 0.04,
            "Kh" : 1e-5
        }
        self.pSimu = {
            "tspan" : [0.0, 1440.0],
            "dt" : 0.1,
            "subRows" : 4,
            "subCols" : 1
        }
        self.microenvironment = {
            "oExtraInitial" : 0.056,
            "gExtraInitial" : 5.0,
            "protonExtraInitial" : 10**(-7.4 + 3),
            "oExtra" : 0.056,
            "gExtra" : 5.0
        }
        self.solutions = None
        self.model = self.model
        self.modelName = ""
        self.description = ""

    @abstractclassmethod
    def model(self, t, x, p):
        pass

    def getO2Extra(self, t):
        return self.microenvironment["oExtra"]

    def getGlucoseExtra(self, t):
        return self.microenvironment["gExtra"]

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
        sol = solve_ivp(self.model, tspan, self.initialCondition, t_eval= t, method='DOP853', args=[self.pOde], rtol = 1e-6, atol=1e-10)
        data = { "Time" : sol.t}
        df_column_names = self.VARS_GENES + self.VARS_METABOLISM
        
        for i in range(len(df_column_names)):
            data[ df_column_names[i] ] = sol.y[i]
        
        for molecule in ["Oxygen", "Glucose"]:
            var_rate = "{} Consumption Rate".format(molecule)
            var_metabolism = "{} Consumed".format(molecule)
            data[var_rate] = self.calculate_rates(data[var_metabolism], self.pSimu["dt"])
        
        for molecule in ["ATP", "H+"]:
            var_rate = "{} Production Rate".format(molecule)
            data[var_rate] = self.calculate_rates(data[molecule], self.pSimu["dt"])

        data["Oxygen Extracellular"] = [self.getO2Extra(t) for t in sol.t ]
        data["Glucose Extracellular"] = [self.getGlucoseExtra(t) for t in sol.t ]
        data["H+ Extracellular"] = [ self.microenvironment["protonExtraInitial"] + x for x in sol.y[7] ]
        data["pH"] = [ -math.log10(1e-3*x) for x in data["H+ Extracellular"]]

        df = pd.DataFrame(data)
        self.solutions = df
            
    def setVo(self, Vo):
        self.pOde["Vo"] = Vo
        self.pOde["A0"] = 29.0/5.0 * Vo

    def plot(self, title, x_title="Time (min)", y_title="Value over Time"):
        data = self.solutions
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

    def plot_vars(self, vars, title, x_title="Time (min)", y_title="Value over Time"):
        i = 1
        sub = make_subplots(rows = len(vars), cols = 1, subplot_titles = vars )
        for v in vars:
            sub.append_trace( go.Scatter(x=self.solutions["Time"], y=self.solutions[v], mode='lines', showlegend=False), row=i, col=1)
            i = i + 1
        sub.update_layout(title = title)
        sub.update_xaxes(title_text = x_title)
        sub.update_yaxes(title_text = y_title, showexponent = 'all', exponentformat = 'e')
        sub.show()
        return sub

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

    def changeGamma(self, i, j, newValue):
        self.pOde["gamma"][i][j] = newValue

    def plot_rates(self, title=None ):
        title = "{} - Reaction Rates".format(self.modelName) if title == None else title
        self.plot_vars(self.VARS_RATES, title)

    def plot_genes(self, title=None ):
        title = "{} - Genes level".format(self.modelName) if title == None else title
        self.plot_vars(self.VARS_GENES, title)

    def plot_microenviroment(self, title=None ):
        title = "{} - Extracellular Concentration".format(self.modelName) if title == None else title
        self.plot_vars(self.VARS_MICROENVIRONMENT, title)

    def plot_metabolism(self, title=None ):
        title = "{} - Metabolism (total Consumption and Production)".format(self.modelName) if title == None else title
        self.plot_vars(self.VARS_METABOLISM, title)

# Genes and Metabolism Separated
class GMSModel(Simu):
    def __init__(self):
        super().__init__()
        self.modelName = "GMS Model"
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

# Genes and Metabolism Connected
class GMCModel(Simu):
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
        self.modelName = "GMC Model"
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

# Proton Linked to Glucose Consumption
class PLGCModel(GMCModel):
    def __init__(self):
        super().__init__()
        self.modelName = "PLGC Model"
    
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
        Ph = Kh * 2 * Cg
        dxdt[0] = A[0] - D[0] * self.H(oExtra, S[0][1], N, gamma[0][1]) * hif
        dxdt[1] = A[1] * self.H(hif, S[1][2], N, gamma[1][2]) - D[1] * ldh
        dxdt[2] = A[2] * self.H(hif, S[1][3], N, gamma[1][3]) - D[2] * pdk
        dxdt[3] = A[3] * self.H(pdk, S[3][4], N, gamma[3][4]) - D[3] * pdh
        dxdt[4] = Co
        dxdt[5] = Cg
        dxdt[6] = Pa
        dxdt[7] = Ph
        return dxdt

def heatmap_h_prod_rate(simulation, title, nbInterval=100):
    oxygen = np.linspace(0.0, 0.056, nbInterval)
    glucose = np.linspace(0.0, 5.0, nbInterval)
    data = np.zeros((nbInterval, nbInterval))
    for i in range(0, nbInterval):
        for j in range(0, nbInterval):
            simulation.microenvironment["gExtra"] = glucose[i]
            simulation.microenvironment["oExtra"] = oxygen[j]
            sol = solve_ivp(simulation.model, simulation.pSimu["tspan"], simulation.initialCondition, method='DOP853', args=[simulation.pOde], rtol = 1e-6, atol=1e-10)
            data[i, j] = sol.y[7, -1] / simulation.pSimu["tspan"][1]
    fig = go.Figure(
        data = go.Heatmap(
            x = oxygen,
            y = glucose,
            z = data,
            zmin = 0.0,
            zmax = data.max(),
            colorbar = {
                "exponentformat" : "e",
                "showexponent" : "all",
                "title" : "H+ Production Rate"
            }
        )
    )
    fig.update_layout(
        title = title,
        xaxis = {
            "title" : "Oxygen"
        },
        yaxis = {
            "title" : "Glucose"
        }
    )
    fig.show()

def run_GMS_model_fixed_condition():
    simulation = GMSModel()
    simulation.initialCondition = [ 3.565, 1.726, 3.31, 0.28, 0.056, 5.0, 0.0, 10**(-7.4)/1000]
    simulation.run()
    simulation.plot_genes()
    simulation.plot_rates()
    simulation.plot_metabolism()
    simulation.plot_microenviroment()

def run_GMC_model_fixed_condition():
    simulation = GMCModel()
    simulation.initialCondition = [ 3.565, 1.726, 3.31, 0.28, 0.056, 5.0, 0.0, 10**(-7.4)/1000]
    simulation.run()
    simulation.plot_genes()
    simulation.plot_rates()
    simulation.plot_metabolism()
    simulation.plot_microenviroment()

def run_GMC_model_fixed_hypoxia():
    simulation = GMCModel()
    simulation.getO2Extra = lambda t : GMCModel.HYPOXIA
    simulation.initialCondition = [ 3.565, 1.726, 3.31, 0.28, 0.056, 5.0, 0.0, 10**(-7.4)/1000]
    simulation.run()
    simulation.plot_genes()
    simulation.plot_rates()
    simulation.plot_metabolism()
    simulation.plot_microenviroment()

def run_GMS_model_hypoxic_period():
    simulation = GMSModel()
    simulation.getO2Extra = hypoxic_period
    simulation.initialCondition = [ 3.565, 1.726, 3.31, 0.28, 0.056, 5.0, 0.0, 10**(-7.4)/1000]
    simulation.modelName = simulation.modelName + " with Hypoxic Period"
    simulation.run()
    simulation.plot_genes()
    simulation.plot_rates()
    simulation.plot_metabolism()
    simulation.plot_microenviroment()

def run_GMC_model_hypoxic_period():
    simulation = GMCModel()
    simulation.getO2Extra = hypoxic_period
    simulation.modelName = simulation.modelName + " with Hypoxic Period"
    simulation.initialCondition = [ 3.565, 1.726, 3.31, 0.28, 0.056, 5.0, 0.0, 10**(-7.4)/1000]
    simulation.run()
    simulation.plot_genes()
    simulation.plot_rates()
    simulation.plot_metabolism()
    simulation.plot_microenviroment()

def run_GMC_model_hypoglycemic_period():
    simulation = GMCModel()
    simulation.getGlucoseExtra = hypoglycemia_period
    simulation.modelName = simulation.modelName + " with Hypoglycemic Period"
    simulation.initialCondition = [ 3.565, 1.726, 3.31, 0.28, 0.056, 5.0, 0.0, 10**(-7.4)/1000]
    simulation.run()
    simulation.plot_genes()
    simulation.plot_rates()
    simulation.plot_metabolism()
    simulation.plot_microenviroment()

def run_GMC_model_physicell_tumor_conditions():
    simulation = GMSModel()
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
    simulation.initialCondition = [ 3.565, 1.726, 3.31, 0.28, 0.056, 5.0, 0.0, 10**(-7.4)/1000]
    simulation.run()
    vars = ["Oxygen Extracellular", "Glucose Extracellular", "H+ Extracellular", "ATP", "pH"]
    title = '{}, pg = {}'.format(simulation.modelName, simulation.pOde["pg"])
    simulation.plot_vars(vars, "{} - Substrate concentration".format(title))
    vars = ["Oxygen Consumption Rate", "Glucose Consumption Rate", "H+ Production Rate", "ATP Production Rate"]
    simulation.plot_vars(vars, "{} Rates".format(title))

def run_GMC_Model_hif_overexpressed():
    simulation = GMCModel()
    simulation.getO2Extra = hypoxic_period
    simulation.modelName = simulation.modelName + " with HIF overexpressed"
    simulation.initialCondition = [ 3.565, 1.726, 3.31, 0.28, 0.056, 5.0, 0.0, 10**(-7.4)/1000]
    simulation.changeGamma(0, 1, 1.0)
    simulation.run()
    simulation.plot_genes()
    simulation.plot_rates()
    simulation.plot_metabolism()
    simulation.plot_microenviroment()

def run_GMC_Model_pdk_overexpressed():
    simulation = GMCModel()
    simulation.getO2Extra = hypoxic_period
    simulation.modelName = simulation.modelName + " with PDK overexpressed"
    simulation.initialCondition = [ 3.565, 1.726, 3.31, 0.28, 0.056, 5.0, 0.0, 10**(-7.4)/1000]
    simulation.changeGamma(1, 3, 8.0)
    simulation.run()
    simulation.plot_genes()
    simulation.plot_rates()
    simulation.plot_metabolism()
    simulation.plot_microenviroment()

def run_GMC_Model_pdk_deactivated():
    simulation = GMCModel()
    simulation.getO2Extra = hypoxic_period
    simulation.modelName = simulation.modelName + " with PDK deactivated"
    simulation.initialCondition = [ 3.565, 1.726, 3.31, 0.28, 0.056, 5.0, 0.0, 10**(-7.4)/1000]
    simulation.changeGamma(3, 4, 1.0)
    simulation.run()
    simulation.plot_genes()
    simulation.plot_rates()
    simulation.plot_metabolism()
    simulation.plot_microenviroment()
    
def run_GMC_Model_hif_overexpressed_and_pdk_deactivated():
    simulation = GMCModel()
    simulation.getO2Extra = hypoxic_period
    simulation.modelName = simulation.modelName + " with HIF overexpressed + PDK deactivated"
    simulation.initialCondition = [ 3.565, 1.726, 3.31, 0.28, 0.0, 0.0, 0.0, 10**(-7.4)/1000]
    simulation.changeGamma(0, 1, 1.0)
    simulation.changeGamma(3, 4, 1.0)
    simulation.run()
    simulation.plot_genes()
    simulation.plot_rates()
    simulation.plot_metabolism()
    simulation.plot_microenviroment()

def run_PLGC_model_fixed_condition():
    simulation = PLGCModel()
    simulation.initialCondition = [ 3.565, 1.726, 3.31, 0.28, 0.0, 0.0, 0.0, 10**(-7.4)/1000]
    simulation.run()
    simulation.plot_genes()
    simulation.plot_rates()
    simulation.plot_metabolism()
    simulation.plot_microenviroment()

def run_PLGC_model_hypoxic_period():
    simulation = PLGCModel()
    simulation.getO2Extra = hypoxic_period
    simulation.modelName = simulation.modelName + " with Hypoxic Period"
    simulation.initialCondition = [ 3.565, 1.726, 3.31, 0.28, 0.0, 0.0, 0.0, 10**(-7.4)/1000]
    simulation.run()
    simulation.plot_genes()
    simulation.plot_rates()
    simulation.plot_metabolism()
    simulation.plot_microenviroment()

def run_PLGC_model_hypoglycemic_period():
    simulation = PLGCModel()
    simulation.getGlucoseExtra = hypoglycemia_period
    simulation.modelName = simulation.modelName + " with Hypoglycemic Period"
    simulation.initialCondition = [ 3.565, 1.726, 3.31, 0.28, 0.0, 0.0, 0.0, 10**(-7.4)/1000]
    simulation.run()
    simulation.plot_genes()
    simulation.plot_rates()
    simulation.plot_metabolism()
    simulation.plot_microenviroment()

def run_PLGC_model_fixed_hypoxia_and_hypoglycemia_period():
    simulation = PLGCModel()
    simulation.getO2Extra = hypoxic_period
    simulation.getGlucoseExtra = hypoglycemia_period
    simulation.modelName = simulation.modelName + " with a Period whithout Nutrient"
    simulation.initialCondition = [ 3.565, 1.726, 3.31, 0.28, 0.0, 0.0, 0.0, 10**(-7.4)/1000]
    simulation.run()
    simulation.plot_genes()
    simulation.plot_rates()
    simulation.plot_metabolism()
    simulation.plot_microenviroment()

def compare_GMC_and_PLGC_model_hypoxic_period():
    run_GMC_model_hypoxic_period()
    run_PLGC_model_hypoxic_period()

def compare_GMC_and_PLGC_model_hypoglycemic_period():
    run_GMC_model_hypoglycemic_period()
    run_PLGC_model_hypoglycemic_period()
   

if __name__ == "__main__":
    pass