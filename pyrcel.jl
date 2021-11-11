## Imports
include("constants.jl")
c = constants
include("thermo.jl")

using CPUTime
using DifferentialEquations
using Plots
using Printf
using Roots
# using StaticArrays  ## DON'T USE THIS.
using Sundials

# Stand-alone run.

## Base model constants 
V = 0.44
T₀ = 283.15
S₀ = -0.05
P₀ = 85000.0

## Create an initial aerosol distribution
μₐ = 0.15
Nₐ = 1000.0
σₐ = 1.2
κₐ = 0.54
n_bins = 251

## Initialize aerosol size distribution
lr = log(μₐ / 10. / σₐ)
rr = log(μₐ * 10. * σₐ)

base = ℯ
rs = base.^(range(lr, stop=rr, length=n_bins))
mids = sqrt.(rs[1:end-1] .* rs[2:end])
r_drys = mids * 1e-6

function pdf(x)
    scaling = Nₐ / sqrt(2*π) / log(σₐ)
    exponent = log(x / μₐ)^2 / 2 / log(σₐ)^2
    return (scaling / x) * exp(-exponent)
end
rsl = rs[1:end-1]
rsr = rs[2:end]
Nis = 0.5 * (rsr - rsl) .* (pdf.(rsl) + pdf.(rsr)) * 1e6

## Initialize parcel init conditions
# Water vapor
es(T_c) = 611.2 * exp(17.67 * T_c / (T_c + 243.5))
wv₀ = (S₀ + 1.0) * (
    c.epsilon * es(T₀ - 273.15) / (P₀ - es(T₀ - 273.15))
)
# Find equilibrium wet particle radius
r0s = []
# NOTE: κ is constant here so we don't need the zipped iteration
for r_dry in reverse(r_drys)
    f(r) = Seq(r, r_dry, T₀, κₐ) - S0
    r_b, _ = kohler_crit(T₀, r_dry, κₐ)
    r_a = r_dry
    r0 = Roots.find_zero(f, (r_a, r_b), )
    # @printf "%g | %g | %g | %g \n" r_a r_dry r0 r_b
    # @printf "(%g %g) | (%g %g) | (%g %g)\n" r_a f(r_a) r0 f(r0) r_b f(r_b)
    append!(r0s, r0)
end
r0s = collect(reverse(r0s))

## Total water volume
𝑉(r, r_dry, Ni) = 
    (4*π/3.0) * c.rho_w * Ni * (r^3 - r_dry^3)
wc₀ = sum(𝑉.(r0s, r_drys, Nis)) / ρ_air(T₀, P₀, 0.0)
wi₀ = 0.0

@printf "AEROSOL DISTRIBUTION\n"
for i = 1:length(r_drys)
    @printf "%3.2e %8.1f\n" r_drys[i] Nis[i]
end

@printf "\nInitial Conditions\n"
@printf "------------------\n"
@printf "%9.1f %9.2f %9.1e %9.1e %9.1e %9.3f\n" P₀/100. T₀ wv₀*1e3 wc₀*1e3 wi₀*1e3 S₀

## Set up ODE solver
y₀ = [P₀, T₀, wv₀, S₀]
append!(y₀, r0s)
accom = 1.0
params = [r_drys, Nis, V, κₐ, accom]

struct atm_state{T}
    T::T
    P::T
    ρ_air::T
    wv::T
    S::T
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

  dwv_dt = -1. * dwc_dt

  dT_dt = -c.g * V / c.Cp - c.L * dwv_dt / c.Cp

  α = (c.g * c.Mw * c.L) / (c.Cp * c.R * T * T)
  α -= (c.g * c.Ma) / (c.R * T)
  γ = (P * c.Ma) / (c.Mw * pv_sat)
  γ += (c.Mw * c.L * c.L) / (c.Cp * c.R * T * T)
  dS_dt = (α * V) - (γ * dwc_dt)

  du[1] = dP_dt
  du[2] = dT_dt
  du[3] = dwv_dt
  du[4] = dS_dt

end

## Solver setup
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

# Warm up ODE RHS
pm_parcel_odes!(zeros(length(y0)), y0, params, 0.)

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
    QNDF();
    # CVODE_BDF( # ~5 seconds on reference problem
    #     method=:Newton,
    #     max_order=5,
    # );  
    dt=output_dt,
    reltol=state_rtol,
    abstol=state_atol,
    tstops=tspan[1]:solver_dt:tspan[2],
    saveat=tspan[1]:output_dt:tspan[2],
)

# step!(integrator, output_dt, true)

@printf "Integration Control\n"
@printf "-------------------\n"
@printf "output_dt: %1.2f\n" output_dt
@printf "solver_dt: %1.2f\n" solver_dt
@printf "\nBEGIN INTEGRATION ->\n\n"
@printf "Integration Loop\n\n"
@printf "  step     time  walltime  Δwalltime |      z       T       S\n"
@printf " ------------------------------------|-----------------------\n"

ts = tspan[1]:output_dt:tspan[2]
# ts = 0:1:10
step = 0
elapsed_time = 0.
CPUtic()
for (u, t) in TimeChoiceIterator(integrator, ts)
    Δt = CPUtoq()
    elapsed_time += Δt
    step += 1
    T = u[2]
    S = u[4]
    z = V * t
    @printf "%6d %7.2fs %8.2fs %9.2fs | %5.1fm %7.2f %6.2f%% \n" step t elapsed_time Δt z T S*100.
    CPUtic()
end
_ = CPUtoq()