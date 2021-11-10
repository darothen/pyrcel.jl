"""
This is intended to be a reference of an idealized API for running
the parcel model simulations.
"""

using Pyrcel

aerosols = [
    AerosolSpecies(
        name = 'sulfate',
        size_dist = Lognorm(
            μ = 0.015,
            σ = 1.6,
            N = 85000.
        )
        κ = 0.54,
        bins = 200
    ),
    AerosolSpecies(
        name = 'sea salt',
        size_dist = Lognorm(
            μ = 0.85,
            σ = 1.2,
            N = 10.
        )
        κ = 1.2,
        bins = 40
    )
]

model = AdiabaticParcelModel(
    aerosols = aerosols,
    state = ParcelState(
        T = 283.15,
        P = 1050.,
        S = -0.05
    )
    updraft = ConstantUpdraft(V = 0.5)
)

simulation = Simulation(
    model, 
    Δt_output = 1.,
    Δt_solver = 0.2,
    t_end = 250.,
    terminate_depth = 25.  # terminate 25 seconds after Smax
)

run!(simulation)

# Plot 2-panel plot:
#   Left) Smax / Temperature vs height
#   Right) Aerosol growth traces
plot_summary(
    aerosol_sample_frac = 0.05  # Plot 5% of aerosol bins per mode
)
