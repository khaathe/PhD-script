#!/usr/bin

import numpy as np
from plotly.subplots import make_subplots
from scipy.integrate import solve_ivp
import plotly.graph_objects as go
import math

NORMOXIA = 0.056 # ~5%
HYPOXIA = 0.010564 # 1%
ANOXIA = 1.0564e-05 # 0.001%
GLUCOSE_NORMAL = 5.0
GLUCOSE_LOW = 1.0
GLUCOSE_VERY_LOW = 0.1

LABELS_PLOT = ["HIF", "LDH", "PDK", "PDH", "Oxygen", "Glucose", "ATP", "H+"]
COLOR_PLOT = ["blue", "green", "goldenrod", "grey", "red", "purple", "cyan", "magenta"]

def get_base_parameters():
    return {
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

def get_base_initial_conditions():
    return {
        "normoxia+normal_glucose" : [ 3.565, 1.726, 3.31, 0.28, NORMOXIA, GLUCOSE_NORMAL, 0.0, 10**(-7.4)/1000],
        "hypoxia+normal_glucose" : [ 3.565, 1.726, 3.31, 0.28, HYPOXIA, GLUCOSE_NORMAL, 0.0, 10**(-7.4)/1000],
        "normoxia+low_glucose" : [ 3.565, 1.726, 3.31, 0.28, NORMOXIA, GLUCOSE_LOW, 0.0, 10**(-7.4)/1000],
        "hypoxia+low_glucose" : [ 3.565, 1.726, 3.31, 0.28, HYPOXIA, GLUCOSE_LOW, 0.0, 10**(-7.4)/1000]
    }

def H(y,s,n,gamma):
    return s**n / ( s**n + y**n ) + gamma * y**n / (s**n + y**n)

def base_model(t, x, p):
    dxdt = [None] * 8
    hif, ldh, pdk, pdh, o, g, atp, h = x
    A, D, N, hif_coeff_prod, S, gamma, Vo, Ko, pg, A0,  Kg, Kh = p.values()
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
    A, D, N, hif_coeff_prod, S, gamma, Vo, Ko, pg, A0,  Kg, Kh = p.values()
    
    L = 50
    beta = 1
    k = 10
    ldh0 = 2.176675
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

def plot_solution(sol, labels=["HIF", "LDH", "PDK", "PDH", "Oxygen", "Glucose", "ATP", "H+"]):
    fig = go.Figure(
        layout=go.Layout(
            title='HIF Model',
            xaxis_title = 'time (min)'
        )
    )
    for i in range(8):
        fig.add_trace( go.Scatter(x = sol.t, y = sol.y[i], mode = "lines", name = labels[i]))
    fig.add_trace( go.Scatter(x = sol.t, y = -np.log10(sol.y[7]*1000), mode = "lines", name = 'pH'))
    return fig

def subplot_solution(sub, row_index, col_index, sol, plot_legend=False):
    for i in range(8):
        sub.append_trace( 
            go.Scatter(
                x = sol.t, 
                y = sol.y[i], 
                mode = "lines", 
                name = LABELS_PLOT[i], 
                legendgroup='group_{}'.format(LABELS_PLOT[i]), 
                showlegend=plot_legend,
                line_color = COLOR_PLOT[i]
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
            line_color = COLOR_PLOT[i]
        ), 
        row=row_index,
        col=col_index
    )
    return sub

def simulate_different_environment(m, initial_conditions, parameter, tspan, dt):
    solutions= dict()
    t = np.arange(tspan[0], tspan[1], dt)
    for condition in initial_conditions:
        solutions[condition] = solve_ivp(m, tspan, initial_conditions[condition], t_eval= t, method='DOP853', args=[parameter], rtol = 1e-6, atol=1e-10)
    return solutions

def run_simulation(model):
    initial_conditions = get_base_initial_conditions()
    p = get_base_parameters()
    solutions = simulate_different_environment(model, initial_conditions, p, [0.0, 1440.0], 0.1)
    sub = make_subplots(rows = 4, cols=1, subplot_titles=list(solutions.keys()))
    subplot_solution(sub, 1, 1, solutions["normoxia+normal_glucose"], plot_legend=True)
    subplot_solution(sub, 2, 1, solutions["hypoxia+normal_glucose"])
    subplot_solution(sub, 3, 1, solutions["normoxia+low_glucose"])
    subplot_solution(sub, 4, 1, solutions["hypoxia+low_glucose"])
    sub.update_layout(title_text="Base Model Result, Vo={}".format(p["Vo"]))
    sub.show()

def run_base_model():
    initial_conditions = get_base_initial_conditions()
    p = get_base_parameters()
    solutions = simulate_different_environment(base_model, initial_conditions, p, [0.0, 1440.0], 0.1)
    sub = make_subplots(rows = 4, cols=1, subplot_titles=list(solutions.keys()))
    subplot_solution(sub, 1, 1, solutions["normoxia+normal_glucose"], plot_legend=True)
    subplot_solution(sub, 2, 1, solutions["hypoxia+normal_glucose"])
    subplot_solution(sub, 3, 1, solutions["normoxia+low_glucose"])
    subplot_solution(sub, 4, 1, solutions["hypoxia+low_glucose"])
    sub.update_layout(title_text="Base Model Result, Vo={}".format(p["Vo"]))
    sub.show()

def run_base_model_increased_o2_consumption():
    initial_conditions = get_base_initial_conditions()
    p = get_base_parameters()
    Vo = 1.16e-4
    p["Vo"] = Vo
    p["A0"] = 29.0/5.0 * Vo
    solutions = simulate_different_environment(base_model, initial_conditions, p, [0.0, 1440.0], 0.1)
    sub = make_subplots(rows = 4, cols=1, subplot_titles=list(solutions.keys()))
    subplot_solution(sub, 1, 1, solutions["normoxia+normal_glucose"], plot_legend=True)
    subplot_solution(sub, 2, 1, solutions["hypoxia+normal_glucose"])
    subplot_solution(sub, 3, 1, solutions["normoxia+low_glucose"])
    subplot_solution(sub, 4, 1, solutions["hypoxia+low_glucose"])
    sub.update_layout(title_text="Base Model Result, Vo={}".format(p["Vo"]))
    sub.show()

def run_base_model_increased_o2_and_glucose_consumption():
    initial_conditions = get_base_initial_conditions()
    p = get_base_parameters()
    Vo = 1.16e-4
    p["Vo"] = Vo
    p["A0"] = 29.0/5.0 * Vo
    p["pg"] = 50
    solutions = simulate_different_environment(base_model, initial_conditions, p, [0.0, 1440.0], 0.1)
    sub = make_subplots(rows = 4, cols=1, subplot_titles=list(solutions.keys()))
    subplot_solution(sub, 1, 1, solutions["normoxia+normal_glucose"], plot_legend=True)
    subplot_solution(sub, 2, 1, solutions["hypoxia+normal_glucose"])
    subplot_solution(sub, 3, 1, solutions["normoxia+low_glucose"])
    subplot_solution(sub, 4, 1, solutions["hypoxia+low_glucose"])
    sub.update_layout(title_text="Base Model Result, Vo = {}, pg = {}".format(p["Vo"], p["pg"]))
    sub.show()

def run_new_model():
    initial_conditions = get_base_initial_conditions()
    p = get_base_parameters()
    solutions = simulate_different_environment(new_model, initial_conditions, p, [0.0, 1440.0], 0.1)
    sub = make_subplots(rows = 4, cols=1, subplot_titles=list(solutions.keys()))
    subplot_solution(sub, 1, 1, solutions["normoxia+normal_glucose"], plot_legend=True)
    subplot_solution(sub, 2, 1, solutions["hypoxia+normal_glucose"])
    subplot_solution(sub, 3, 1, solutions["normoxia+low_glucose"])
    subplot_solution(sub, 4, 1, solutions["hypoxia+low_glucose"])
    sub.update_layout(title_text="New Model Result, Vo={}".format(p["Vo"]))
    sub.show()

def run_new_model_increased_o2_consumption():
    initial_conditions = get_base_initial_conditions()
    p = get_base_parameters()
    Vo = 1.16e-4
    p["Vo"] = Vo
    p["A0"] = 29.0/5.0 * Vo
    solutions = simulate_different_environment(new_model, initial_conditions, p, [0.0, 1440.0], 0.1)
    sub = make_subplots(rows = 4, cols=1, subplot_titles=list(solutions.keys()))
    subplot_solution(sub, 1, 1, solutions["normoxia+normal_glucose"], plot_legend=True)
    subplot_solution(sub, 2, 1, solutions["hypoxia+normal_glucose"])
    subplot_solution(sub, 3, 1, solutions["normoxia+low_glucose"])
    subplot_solution(sub, 4, 1, solutions["hypoxia+low_glucose"])
    sub.update_layout(title_text="New Model Result, Vo={}".format(p["Vo"]))
    sub.show()

if __name__ == "__main__":
    # run_base_model()
    run_base_model_increased_o2_consumption()
    # run_base_model_increased_o2_and_glucose_consumption()
    # run_new_model()
    run_new_model_increased_o2_consumption()