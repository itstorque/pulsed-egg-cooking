module PulsedEggSimulations

    using LinearAlgebra
    using DiffEqBase
    using OrdinaryDiffEq: SplitODEProblem, solve, IMEXEuler, KenCarp4, SBDF2, Rosenbrock23
    import SciMLBase
    using Plots
    using Statistics
    using Unitful
    using Unitful: 𝐓, 𝚯, 𝐋

    export thermodynamic_simulation, thermodynamic_simulation_steps, gelation_sim
    export plot_temperature, plot_gelation, plot_temperature_and_gelation

    include("gelation_dynamics.jl")
    include("physical.jl")
    include("plotting.jl")
    include("thermodynamics.jl")

end # module PulsedEggSimulations
