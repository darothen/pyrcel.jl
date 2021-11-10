using Optim: Optim

es(T_c) = 611.2 * exp(17.67 * T_c / (T_c + 243.5))
σ_w(T) = 0.0761 - 1.55e-4 * (T - 273.15)

function dv(T, r, P, accom)
    P_atm = P * 1.01325e-5
    dv_cont = 1e-4 * (0.211 / P_atm) * ((T / 273.0)^1.94)
    denom = 1.0 + (dv_cont / accom / r) * sqrt(2*π*c.Mw/c.R/T)
    return dv_cont / denom
end

function ka(T, r, rho)
    ka_cont = 1e-3 * (4.39 + 0.071*T)
    denom = 1.0 + (ka_cont / c.at / r / rho / c.Cp) * sqrt(2*π*c.Ma/c.R/T)
    return ka_cont / denom
end

function kohler_crit(T, r_dry, κ)::Tuple{Float64, Float64}
    neg_Seq(r) = -1.0 * Seq(r, r_dry, T, κ)
    res = Optim.optimize(
        neg_Seq, r_dry, r_dry*1e4, Optim.Brent()
    )
    r_crit = Optim.minimizer(res)
    s_crit = Optim.minimum(res) * -1.0
    return r_crit, s_crit
end

function ρ_air(T, P, RH=1.0)
    qsat = RH * 0.622 * es(T - 273.15) / P
    Tv = T * (1.0 + 0.61 * qsat)
    return P / c.Rd / Tv
end

function Seq(r, r_dry, T, κ)
    A = 2. * c.Mw * σ_w(T) / c.R / T / c.rho_w / r
    B = (r^3 - r_dry^3) / (r^3 - (r_dry^3 * (1.0 - κ)))
    return exp(A) * B - 1.0
end
