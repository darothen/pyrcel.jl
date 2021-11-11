
abstract type ParticleSpecies end

struct AerosolSpecies <: ParticleSpecies
    name :: String
    r_drys :: AbstractVector{Real}
    Ns :: AbstractVector{Real}
    κs :: AbstractVector{Real}

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

        κs = κ*ones(eltype(κ), nbins)
        
        new(name, r_drys, Ns, κs)
    end
end

function Base.show(io::IO, aer::AerosolSpecies)
    print(io, "AerosolSpecies(", aer.name, ", nbins=", length(aer.Ns), ")")
end