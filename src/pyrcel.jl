module pyrcel

if VERSION < v"1.5"
    @error "pyrcel requires Julia v1.5 or newer."
end

## Package includes / exports
# include("Package/Package.jl")

## System-level configuration
if get(ENV, "JULIA_PYRCEL_LOAD_PYPLOT", "1")
    # default
    # include("Plotting/Plotting.jl")
else   
    @info "Not eagerly loading plotting functions."
end

end