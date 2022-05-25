#!/usr/bin/python3

from abc import ABC, abstractclassmethod
from operator import mod
from matplotlib.pyplot import title
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

y0 = [ 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 10**(-7.4)/1000]

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
            "A" : [0.05, 0.005, 0.005, 0.005, 0.005,], # We multiply HIF prod rate by 140 or we should change gamma_O2-> as well
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
                [0.0,10.0,0.0,0.0,0.0],
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
            "Kh" : 2.5e-4
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

    def get_plot_vars(self, vars, title, x_title="Time (min)", y_title="Value over Time"):
        i = 1
        sub = make_subplots(rows = len(vars), cols = 1, subplot_titles = vars )
        for v in vars:
            sub.append_trace( go.Scatter(x=self.solutions["Time"], y=self.solutions[v], mode='lines', showlegend=False), row=i, col=1)
            i = i + 1
        sub.update_layout(title = title)
        sub.update_xaxes(title_text = x_title)
        sub.update_yaxes(title_text = y_title, showexponent = 'all', exponentformat = 'e')
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

    def get_plot_rates(self, title=None ):
        title = "{} - Reaction Rates".format(self.modelName) if title == None else title
        return self.get_plot_vars(self.VARS_RATES, title)

    def get_plot_genes(self, title=None ):
        title = "{} - Genes level".format(self.modelName) if title == None else title
        return self.get_plot_vars(self.VARS_GENES, title)

    def get_plot_microenviroment(self, title=None ):
        title = "{} - Extracellular Concentration".format(self.modelName) if title == None else title
        return self.get_plot_vars(self.VARS_MICROENVIRONMENT, title)

    def get_plot_metabolism(self, title=None ):
        title = "{} - Metabolism (total Consumption and Production)".format(self.modelName) if title == None else title
        return self.get_plot_vars(self.VARS_METABOLISM, title)

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
            "k" : 4,
            "ldh0" : 2.35,
            "po_max" : 1,
            "po_min" : 0.0,
            "l" : 15,
            "pdh0" : 0.575
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

class CustomModel(GMCModel):
    def __init__(self):
        super().__init__()
        self.modelName = "Custom Model"
    
    def model(self, t, x, p):
        dxdt = [None] * len(x)
        hif, ldh, pdk, pdh, oConsumption, gConsumption, atp, h = x
        A, D, N, S, gamma, Vo, Ko, pg, A0,  Kg, Kh, pg_max, pg_min, k, ldh0, po_max, po_min, l, pdh0 = p.values()
        oExtra = self.getO2Extra(t)
        gExtra = self.getGlucoseExtra(t)
        po = (po_max-po_min)/(1+math.exp(-l*(pdh-pdh0))) + po_min
        pg = (pg_max-pg_min)/(1+math.exp(-k*(ldh-ldh0))) + pg_min
        Co = po * Vo * ( oExtra/(oExtra+Ko) )
        Cg = ( (pg*A0/2.0) - (29.0*Co/10.0) ) * ( gExtra/(gExtra+Kg) )
        Pa = ( (2.0 * Cg) + ( (29.0/5.0)*Co ) )
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

def simulate_at_different_extra_concentration(simulation, nbInterval):
    oxygen = np.linspace(0.0, 0.056, nbInterval)
    glucose = np.linspace(0.0, 5.0, nbInterval)
    data = np.zeros((nbInterval, nbInterval))
    for i in range(0, nbInterval):
        for j in range(0, nbInterval):
            simulation.microenvironment["gExtra"] = glucose[i]
            simulation.microenvironment["oExtra"] = oxygen[j]
            sol = solve_ivp(simulation.model, simulation.pSimu["tspan"], simulation.initialCondition, method='DOP853', args=[simulation.pOde], rtol = 1e-6, atol=1e-10)
            data[i, j] = sol.y[7, -1] / simulation.pSimu["tspan"][1]
    return [oxygen, glucose, data]

def get_heatmap_h_prod_rate(simulation, title, nbInterval=10):
    oxygen, glucose, data = simulate_at_different_extra_concentration(simulation, nbInterval)
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
    return fig

def get_contour_h_prod_rate(simulation, title, nbInterval=10):
    oxygen, glucose, data = simulate_at_different_extra_concentration(simulation, nbInterval)
    fig = go.Figure(
        data = go.Contour(
            x = oxygen,
            y = glucose,
            z = data,
            zmin = 0.0,
            zmax = data.max(),
            colorbar = {
                "exponentformat" : "e",
                "showexponent" : "all",
                "title" : "Protons Rate (mmol/L/min)"   
            },
             contours=dict(
                start=0,
                end=0.001,
                size=1e-4
            )
        )
    )
    fig.update_layout(
        title = title,
        xaxis = {
            "title" : "Oxygen concentration (mmol/L)"
        },
        yaxis = {
            "title" : "Glucose concentration (mmol/L)"
        }
    )
    return fig

def run_PLGC_Model_hif_overexpressed():
    simulation = PLGCModel()
    simulation.getO2Extra = hypoxic_period
    simulation.modelName = simulation.modelName + " - gamma_Oxygen->HIF = 1"
    simulation.initialCondition = y0
    simulation.changeGamma(0, 1, 1.0)
    simulation.run()
    simulation.get_plot_genes().show()
    simulation.get_plot_rates().show()
    simulation.get_plot_metabolism().show()
    simulation.get_plot_microenviroment().show()

def run_PLGC_Model_ldh_overexpressed():
    simulation = PLGCModel()
    simulation.getO2Extra = hypoxic_period
    simulation.modelName = simulation.modelName + " - gamma HIF->LDH = 10"
    simulation.initialCondition = y0
    simulation.changeGamma(1, 2, 40.0)
    simulation.run()
    simulation.get_plot_genes().show()
    simulation.get_plot_rates().show()
    simulation.get_plot_metabolism().show()
    simulation.get_plot_microenviroment().show()

def run_PLGC_Model_pdh_overexpressed():
    simulation = PLGCModel()
    simulation.getO2Extra = hypoxic_period
    simulation.modelName = simulation.modelName + " - gamma PDK->PDH = 1"
    simulation.initialCondition = y0
    simulation.changeGamma(3, 4, 1.0)
    simulation.run()
    simulation.get_plot_genes().show()
    simulation.get_plot_rates().show()
    simulation.get_plot_metabolism().show()
    simulation.get_plot_microenviroment().show()
    
def run_PLGC_Model_hif_and_pdh_overexpressed():
    simulation = GMCModel()
    simulation.getO2Extra = hypoxic_period
    simulation.modelName = simulation.modelName + " - gamma Oxygen->HIF = 1, gamma PDK->PDH = 1"
    simulation.initialCondition = y0
    simulation.changeGamma(0, 1, 1.0)
    simulation.changeGamma(3, 4, 1.0)
    simulation.run()
    simulation.get_plot_genes().show()
    simulation.get_plot_rates().show()
    simulation.get_plot_metabolism().show()
    simulation.get_plot_microenviroment().show()

def alter_all_genetic_regulations():
    dir = "/home/spinicck/PhD/Code/output/hif-model/images/contourmap/"
    title = "Protons production rate in various extracellular conditions"
    simulation = PLGCModel()
    simulation.initialCondition = y0
    get_contour_h_prod_rate(simulation, title).write_image(dir + "normal.png", scale=3)

    gammas = [ (0,1), (1,2), (1,3), (3, 4) ]
    activation = [1, 2, 3, 5, 10, 40]
    names = ["O2", "HIF", "LDH", "PDK", "PDH"]
    for (x,y) in gammas:
            for g in activation:
                simulation = PLGCModel()
                simulation.initialCondition = y0
                simulation.changeGamma(x, y, g)
                fig = get_contour_h_prod_rate(
                    simulation, 
                    "{} - γ{}->{} = {}".format(title, names[x], names[y], g)
                )
                file = dir + "{}_{}_{}.png".format(names[x], names[y], g)
                fig.write_image(file, scale=3)

def figure_for_publication():
    nbInterval = 10

    fig = make_subplots(
        rows = 2,
        cols = 2, 
        subplot_titles= (
            "Normal Conditions",
            "Reduced HIF degradation",
            "Increased HIF degration",
            "Reduced LDH upregulation by HIF"
        )
    )
    
    simulation = PLGCModel()
    simulation.initialCondition = y0
    oxygen, glucose, data = simulate_at_different_extra_concentration(simulation, nbInterval)
    fig.add_trace(
        go.Contour(
            x = oxygen,
            y = glucose,
            z = data,
            zmin = 0.0,
            zmax = data.max(),
            contours=dict(
                start=0,
                end=0.001,
                size=1e-4
            ),
            coloraxis="coloraxis"
        ),
        row = 1,
        col = 1
    )

    simulation = PLGCModel()
    simulation.initialCondition = y0
    simulation.changeGamma(0, 1, 1)
    oxygen, glucose, data = simulate_at_different_extra_concentration(simulation, nbInterval)
    fig.add_trace(
        go.Contour(
            x = oxygen,
            y = glucose,
            z = data,
            zmin = 0.0,
            zmax = data.max(),
            contours=dict(
                start=0,
                end=0.001,
                size=1e-4
            ),
            coloraxis="coloraxis"
        ),
        row = 1,
        col = 2
    )
    
    simulation = PLGCModel()
    simulation.initialCondition = y0
    simulation.changeGamma(0, 1, 40)
    oxygen, glucose, data = simulate_at_different_extra_concentration(simulation, nbInterval)
    fig.add_trace(
        go.Contour(
            x = oxygen,
            y = glucose,
            z = data,
            zmin = 0.0,
            zmax = data.max(),
            contours=dict(
                start=0,
                end=0.001,
                size=1e-4
            ),
            coloraxis="coloraxis"
        ),
        row = 2,
        col = 1
    )

    simulation = PLGCModel()
    simulation.initialCondition = y0
    simulation.changeGamma(1, 2, 3)
    oxygen, glucose, data = simulate_at_different_extra_concentration(simulation, nbInterval)
    fig.add_trace(
        go.Contour(
            x = oxygen,
            y = glucose,
            z = data,
            zmin = 0.0,
            zmax = data.max(),
            contours=dict(
                start=0,
                end=0.001,
                size=1e-4
            ),
            coloraxis="coloraxis"
        ),
        row = 2,
        col = 2
    )

    fig.update_xaxes(title_text = "Oxygen concentration (mmol/L)")
    fig.update_yaxes(title_text = "Glucose concentration (mmol/L)")
    fig.update_layout(
        title_text = "Protons production rate in various extracellular conditions",
        title_font_family = "Arial",
        title_font_size = 36,
        coloraxis = dict(
            colorbar = {
                "exponentformat" : "e",
                "showexponent" : "all",
                "title" : "Protons Rate (mmol/L/min)"   
            }
        ),
        font = dict(
            family = "Arial",
            size = 16
        )
    )

    fig.update_annotations(
        font_size = 24
    )
    
    file = "/home/spinicck/PhD/Code/output/hif-model/images/contourmap/publication_ready_plot.png"
    fig.write_image(file, width = 2000, height  = 1000)

def test():
    simulation = PLGCModel()
    simulation.initialCondition = y0
    get_contour_h_prod_rate(
        simulation, 
        "Protons production rate in various\noxygen and glucose extracellular conditions"
    ).show()

if __name__ == "__main__":
    pass