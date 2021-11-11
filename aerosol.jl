# include("size_dists.jl")
import .size_dists

abstract type ParticleSpecies end

struct AerosolBin{T}
    r_dry::T
    r₁::T
    r₂::T
    N::T
    κ::T
end

function Base.show(io::IO, bin::AerosolBin)
    print(io, "AerosolBin(r_dry=", bin.r_dry, ", N=", bin.N, ")")
end

struct AerosolSpecies <: ParticleSpecies
    name :: String
    bins :: AbstractVector{AerosolBin}

    function AerosolSpecies(
        name::String, size_dist::SizeDistribution, nbins::Integer, κ::Real
    )
        lr = log(dist.μ / 10. / size_dist.σ)
        rr = log(dist.μ * 10. * size_dist.σ)
        rs = (dist.base).^(range(lr, stop=rr, length=nbins+1))
        mids = sqrt.(rs[1:end-1] .* rs[2:end])
        r_drys = mids * 1e-6

        rsl = rs[1:end-1]
        rsr = rs[2:end]
        Ns = 0.5 * (rsr - rsl) .* (pdf.(rsl, Ref(size_dist)) + pdf.(rsr, Ref(size_dist))) * 1e6

        bins = []
        for i in 1:n_bins
            push!(bins, AerosolBin(r_drys[i], rs[i], rs[i+1], Ns[i], κ))
        end
        
        new(name, bins)
    end
end

function Base.show(io::IO, aer::AerosolSpecies)
    print(io, "AerosolSpecies(", aer.name, ", nbins=", length(aer.bins), ")")
end

# TODO - kill this and just use StructArray?
# id(s::State) = s.ID
# id.(states)
r_dry(bin::AerosolBin) = bin.r_dry
r_drys(aer::AerosolSpecies) = r_dry.(aer.bins)
N(bin::AerosolBin) = bin.N
Ns(aer::AerosolSpecies) = N.(aer.bins)
κ(bin::AerosolBin) = bin.κ
κs(aer::AerosolSpecies) = κ.(aer.bins)
