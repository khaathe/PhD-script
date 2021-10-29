using DifferentialEquations
using Plots

plotlyjs()

H(y,s,n,gamma) =  s^n / ( s^n + y^n ) + gamma * y^n / (s^n + y^n)

function model(dxdt, x, p, t)
    # Compute Flux values
    hif, ldh, pdk, pdh, o, g, atp, h = x 
    A, D, N, hif_coeff_prod, S, gamma, Vo, Ko, pg, A0,  Kg, Kh = p

	Co = Vo * ( o/(o+Ko) )
	Cg = ( (pg*A0/2.0) - (27.0*Co/10.0) ) * ( g/(g+Kg) )
	Pa = ( (2.0 * Cg) + ( (27.0/5.0)*Co ) )
	Ph = Kh*( (29.0*(pg*Vo-Co))/5.0 )

    dxdt[1] = A * hif_coeff_prod - D * H(o, S[1][2], N, gamma[1][2]) * hif
    dxdt[2] = A * H(hif, S[2][3], N, gamma[2][3]) - D * ldh
    dxdt[3] = A * H(hif, S[2][4], N, gamma[2][4]) - D * pdk
    dxdt[4] = A * H(pdk, S[4][5], N, gamma[4][5]) - D * pdh
    dxdt[5] = -Co
    dxdt[6] = -Cg
    dxdt[7] = Pa
    dxdt[8] = Ph
end

A = 0.005
D = 0.005
N = 4.0
hif_coeff_prod = 140.0
Vo = 2.1e-11
Ko = 0.005
pg = 1.0
A0 = 29.0/5.0 * Vo
Kg = 0.04
Kh = 2.5e-4
S = [ 
    [0.0,20.85e-3,0.0,0.0,0.0],
    [0.0,0.0,4.64,4.0,0.0],
    [0.0,0.0,0.0,0.0,0.0],
    [0.0,0.0,0.0,0.0,2.2],
    [0.0,0.0,0.0,0.0,0.0]
]
gamma = [
    [0.0,40.0,0.0,0.0,0.0],
    [0.0,0.0,3.81,6.97,0.0],
    [0.0,0.0,0.0,0.0,0.0],
    [0.0,0.0,0.0,0.0,0.14],
    [0.0,0.0,0.0,0.0,0.0]
]
hypoxia = 0.010564
normoxia = 0.056
gluNormal = 5.0
gluLow = 1.0

x0 = [ hif_coeff_prod / gamma[1][2], gamma[2][3], gamma[2][4], gamma[4][5], normoxia, gluNormal, 0.0, 10^(-7.4)/1000]
tspan = (0.0, 480.0)
dt = 0.1
p = [A, D, N, hif_coeff_prod, S, gamma, Vo, Ko, pg, A0,  Kg, Kh]
prob = ODEProblem(model, x0, tspan, p)
sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8, saveat=0.1)

plot(sol, label=["HIF" "LDH" "PDK" "PDH" "O2" "Glucose" "ATP" "H+"])
plot!(sol.t, [ -log10(u[8]*1000) for (u) in sol.u], label="pH")
