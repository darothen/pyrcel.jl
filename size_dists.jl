using SpecialFunctions: erf, erfinv

abstract type SizeDistribution end

struct LogNormal <: SizeDistribution
    μ :: Real
    σ :: Real
    N :: Real
    base :: Real
    function LogNormal(μ::Real, σ::Real, N::Real; base::Real=ℯ)
        @assert(μ > 0, "μ should be positive")
        @assert(σ > 0, "σ should be positive")
        @assert(N ≥ 0, "N should be 0 or positive")
        new(μ, σ, N, base)
    end
end

function Base.show(io::IO, dist::LogNormal)
    print(io, "LogNormal(μ=", dist.μ, ", σ=", dist.σ, ", N=", dist.N, ", base=",
          dist.base, ")")
end

function cdf(x::Real, dist::LogNormal)
    erf_arg = log(dist.base, x / dist.μ) / sqrt(2.0) / log(dist.base, dist.σ)
    return (dist.N / 2.0) * (1.0 + erf(erf_arg))
end

function invcdf(y::Real, dist::LogNormal)
    @assert(0 ≤ y ≤ 1, "y must be between 0 and 1")
    erfinv_arg = 2.0 * y / dist.N - 1.0
    exponent = log(dist.base, dist.σ) * sqrt(2.0) * erfinv(erfinv_arg)
    return dist.μ * exp(exponent)
end

function moment(k::Integer, dist::LogNormal)
    scaling = dist.μ^k * dist.N 
    exponent = (k^2)/2.0 * log(dist.base, dist.σ)^2
    return scaling * exp(exponent)
end

function pdf(x::Real, dist::LogNormal)
    @assert(x > 0, "x must be positive")
    scaling = dist.N / sqrt(2*π) / log(dist.base, dist.σ)
    exponent = log(dist.base, x / dist.μ)^2 / 2 / log(dist.base, dist.σ)^2
    return (scaling / x) * exp(-exponent)
end

# Bind a functor to our distributon so calling it passes through to the PDF func
function (dist::LogNormal)(x::Real)
    pdf(x, dist)
end