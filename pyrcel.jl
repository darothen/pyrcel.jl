## Imports
include("constants.jl")
c = constants

using DifferentialEquations
using Optim
using Plots
using Printf
using Roots
# using StaticArrays  ## DON'T USE THIS.
using Sundials

# Stand-alone run.

## Base model constants 
V = 0.44
T0 = 283.15
S0 = -0.05
# S0 = -0.001
P0 = 85000.0

## Create an initial aerosol distribution
μ_aer = 0.15
N_aer = 1000.0
σ_aer = 1.2
κ_aer = 0.54
n_bins = 250

## Initialize aerosol size distribution
lr = log(μ_aer / 10. / σ_aer)
rr = log(μ_aer * 10. * σ_aer)

base = ℯ
rs = base.^(range(lr, stop=rr, length=n_bins))
mids = sqrt.(rs[1:end-1] .* rs[2:end])
r_drys = mids * 1e-6

function pdf(x::Float64)::Float64
    scaling = N_aer / sqrt(2*π) / log(σ_aer)
    exponent = log(x / μ_aer)^2 / 2 / log(σ_aer)^2
    return (scaling / x) * exp(-exponent)
end
rsl = rs[1:end-1]
rsr = rs[2:end]
Nis = 0.5 * (rsr - rsl) .* (pdf.(rsl) + pdf.(rsr)) * 1e6

## Initialize parcel init conditions
# Water vapor
es(T_c) = 611.2 * exp(17.67 * T_c / (T_c + 243.5))
wv0 = (S0 + 1.0) * (
    c.epsilon * es(T0 - 273.15) / (P0 - es(T0 - 273.15))
)
# Find equilibrium wet particle radius
σ_w(T) = 0.0761 - 1.55e-4 * (T - 273.15)
function Seq(r, r_dry, T, κ)
    A = 2. * c.Mw * σ_w(T) / c.R / T / c.rho_w / r
    B = (r^3 - r_dry^3) / (r^3 - (r_dry^3 * (1.0 - κ)))
    return exp(A) * B - 1.0
end
function kohler_crit(T::Float64, r_dry::Float64, κ::Float64)
    neg_Seq(r) = -1.0 * Seq(r, r_dry, T, κ)
    res = optimize(
        neg_Seq, r_dry, r_dry*1e4, Optim.Brent()
    )
    r_crit = Optim.minimizer(res)
    s_crit = Optim.minimum(res) * -1.0
    return r_crit, s_crit
end
r0s = []
# NOTE: κ is constant here so we don't need the zipped iteration
for r_dry in reverse(r_drys)
    f(r) = Seq(r, r_dry, T0, κ_aer) - S0
    r_b, _ = kohler_crit(T0, r_dry, κ_aer)
    r_a = r_dry
    # @printf "r_b = %g\n" r_b
    r0 = find_zero(f, (r_a, r_b), )
    # @printf "%g | %g | %g | %g \n" r_a r_dry r0 r_b
    # @printf "(%g %g) | (%g %g) | (%g %g)\n" r_a f(r_a) r0 f(r0) r_b f(r_b)
    append!(r0s, r0)
end
r0s = collect(reverse(r0s))

## Total water volume
function water_vol(r0::Float64, r_dry::Float64, Ni::Float64)::Float64
    return (4*π/3.0) * c.rho_w * Ni * (r0^3 - r_dry^3)
end
es(T_c) = 611.2 * exp(17.67 * T_c / (T_c + 243.5))
function rho_air(T, P, RH=1.0)
    qsat = RH * 0.622 * es(T - 273.15) / P
    Tv = T * (1.0 + 0.61 * qsat)
    return P / c.Rd / Tv
end
wc0 = sum(water_vol.(r0s, r_drys, Nis)) / rho_air(T0, P0, 0.0)
wi0 = 0.0

@printf "AEROSOL DISTRIBUTION\n"
for i = 1:length(r_drys)
    @printf "%3.2e %8.1f\n" r_drys[i] Nis[i]
end

@printf "\nInitial Conditions\n"
@printf "------------------\n"
@printf "%9.1f %9.2f %9.1e %9.1e %9.1e %9.3f\n" P0/100. T0 wv0*1e3 wc0*1e3 wi0*1e3 S0

## Set up ODE solver
y0 = [P0, T0, wv0, S0]
append!(y0, r0s)
accom = 1.0
params = [r_drys, Nis, V, κ_aer, accom]

function ka(T, r, rho)
    ka_cont = 1e-3 * (4.39 + 0.071*T)
    denom = 1.0 + (ka_cont / c.at / r / rho / c.Cp) * sqrt(2*π*c.Ma/c.R/T)
    return ka_cont / denom
end
function dv(T, r, P, accom)
    P_atm = P * 1.01325e-5
    dv_cont = 1e-4 * (0.211 / P_atm) * ((T / 273.0)^1.94)
    denom = 1.0 + (dv_cont / accom / r) * sqrt(2*π*c.Mw/c.R/T)
    return dv_cont / denom
end

struct atm_state
    T
    P
    ρ_air
    wv
    S
end

function calc_dr_dt(r, r_dry, κ, atm, accom=accom)
    dv_r = dv(atm.T, r, atm.P, accom)
    ka_r = ka(atm.T, r, atm.ρ_air)

    T_c = atm.T - 273.15
    pv_sat = es(T_c)

    G_a = (c.rho_w * c.R * atm.T) / (pv_sat * dv_r * c.Mw)
    G_b = (c.L * c.rho_w * ((c.L * c.Mw / (c.R * atm.T)) - 1.0)) / (ka_r * atm.T)
    G = 1.0 / (G_a + G_b)

    Seq_r = Seq(r, r_dry, atm.T, κ)
    dS = atm.S - Seq_r

    dr_dt = (G / r) * dS
    return dr_dt
end

function pm_parcel_odes!(du,u,p,t)
  P, T, wv, S = u[1:4]
  rs = u[5:end]
  r_drys, Nis, V, κ, accom = p

  # Compute air densities from current state
  T_c = T - 273.15
  pv_sat = es(T_c)
  wv_sat = wv / (S + 1.0)
  Tv = (1.0 + 0.61 * wv) * T
  e_wv = (1.0 + S) * pv_sat # water vapor pressure
  ρ_air = P / c.Rd / Tv
  ρ_air_dry = (P - e_wv) / c.Rd / T

  # Save state
  atm = atm_state(T, P, ρ_air, wv, S)

  # Tendencies
  dP_dt = -1.0 * ρ_air * c.g * V
  dwc_dt = 0.
  for i = 1:length(rs)
    r = rs[i]
    r_dry = r_drys[i]
    Ni = Nis[i] 

    ## Calc dr/dt here
    # ka_r = ka(T, r, ρ_air)
    # dv_r = dv(T, r, P, accom)

    # G_a = (c.rho_w * c.R * T) / (pv_sat * dv_r * c.Mw)
    # G_b = (c.L * c.rho_w * ((c.L * c.Mw / (c.R * T)) - 1.0)) / (ka_r * T)
    # G = 1.0 / (G_a + G_b)s

    # Seq_r = Seq(r, r_dry, T, κ)
    # dS = S - Seq_r

    # dr_dt = (G / r) * dS

    ## Calc dr/dt in function
    dr_dt = calc_dr_dt(r, r_dry, κ, atm, accom)

    dwc_dt += Ni * r * r * dr_dt
    # drs_dt[i] = dr_dt
    du[4+i] = dr_dt

  end
  dwc_dt *= (4*π * c.rho_w / ρ_air_dry)
  # dwc_dt = sum(Nis.*(rs.^2).*drs_dt) * 4*π * c.rho_w / ρ_air_dry

  dwv_dt = -1. * dwc_dt

  dT_dt = -c.g * V / c.Cp - c.L * dwv_dt / c.Cp

  # α = (c.g * c.Mw * c.L) / (c.Cp * c.R * atm.T^2)
  # α -= (c.g * c.Ma) / (c.R * atm.T)
  # γ = (atm.P * c.Ma) / (c.Mw * pv_sat)
  # γ += (c.Mw * c.L * c.L) / (c.Cp * c.R * atm.T^2)
  α = (c.g * c.Mw * c.L) / (c.Cp * c.R * T * T)
  α -= (c.g * c.Ma) / (c.R * T)
  γ = (P * c.Ma) / (c.Mw * pv_sat)
  γ += (c.Mw * c.L * c.L) / (c.Cp * c.R * T * T)
  dS_dt = (α * V) - (γ * dwc_dt)

  du[1] = dP_dt
  du[2] = dT_dt
  du[3] = dwv_dt
  du[4] = dS_dt
  # du[5:end] = drs_dt

end

## Solver setup
# state_atol = [1e-4, 1e-4, 1e-10, 1e-8]
state_atol = [1e-4, 1e-4, 1e-10, 1e-8]
append!(state_atol, 1e-12*ones(length(rs)))
state_rtol = 1e-7

tspan = (0.0, 280)
output_dt = 1.0
solver_dt = 0.5
# n_out = convert(Integer, solver_dt / output_dt)

prob = ODEProblem(pm_parcel_odes!,y0,tspan,params)
u = y0
p = params

## One-shot solve
@time sol = solve(
    prob, 

    ## SOLVER
    # Rodas4P(), # ~29 seconds
    # Rodas5(), # ~40 seconds on reference problem
    # CVODE_BDF( # ~5 seconds on reference problem
    #     method=:Newton,
    #     max_order=5,
    # ),  
    # RadauIIA3(), # ~14 seconds but incorrect solution
    # Rosenbrock23(), # >60 seconds
    # Kvaerno5(), # ~13 seconds
    # KenCarp4(), # ~9 seconds
    # QNDF(), # ~12 seconds
    # QNDF1(), # ~26 seconds
    QBDF(), # ~4-8 seconds - much faster on second run

    ## USE THESE SETTINGS
    reltol=state_rtol,
    abstol=state_atol,
    # dt=solver_dt,
    tstops=tspan[1]:solver_dt:tspan[2],
    saveat=tspan[1]:output_dt:tspan[2],
)
# plot(sol.t, sol[4,:])
println(maximum(sol[4,:]))

## Integrator solve
integrator = init(
    prob, 
    # QNDF(),
    CVODE_BDF(
        method=:Newton,
        max_order=5,
    ), 
    reltol=state_rtol,
    abstol=state_atol,
    tstops=output_dt,
    saveat=output_dt,
)

ts = tspan[1]:output_dt:tspan[2]
# for (u, t) in TimeChoiceIterator(integrator)
for (u, t) in tuples(integrator)
    @show u[1:7], t
end