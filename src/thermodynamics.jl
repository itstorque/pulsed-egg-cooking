"""
Perform an FEM thermodynamic simulation of an egg's cross-section to get the temperature profile over time.
Units are in Celsius, mm, and seconds for temperature, length and time respectively.
"""
function thermodynamic_simulation(
        # Drive parameters
        t_total::Float64,
        source_temp::Function,
        condition::Function,
        affect!::Function;
        # Simulation parameters
        Δt::Float64,
        Tegg_init::Real=20.,
        n_cells::Integer=50,
        scaler::Real=0.37,
        scaler_asym::Real=1,
        # Egg parameters
        egg_diameter::Real=6e-2,
        yolk_diameter::Real=3e-2,
        T_init::Union{Nothing, Real}=nothing,
        T_clamp_alpha::Real=200.,
        solver_opts::Dict=Dict(),
        alg=SBDF2(autodiff=false) #IMEXEuler
    )

    @show t_total
    @show source_temp
    @show condition
    @show affect!
    @show Δt
    @show Tegg_init
    @show n_cells
    @show scaler
    @show scaler_asym
    @show egg_diameter
    @show yolk_diameter
    @show T_init
    @show T_clamp_alpha
    @show solver_opts
    @show alg

    if T_init == nothing
        T_init = 24.0
    end

    # renormalize for egg lengthscales
    egg_diameter *= scaler
    yolk_diameter *= scaler

    # Define thermal diffusivity functions using physical parameters
    α_yolk = T -> k_yolk(T) / (ρ_yolk(T) * c_p_yolk * (scaler_asym*egg_diameter)^2) # m^2/s
    α_albumen = T -> k_albumen(T) / (ρ_albumen(T) * c_p_albumen * egg_diameter^2) # m^2/s

    """
    Calculate the dynamic alpha vector for the egg simulation
    """
    function dynamic_alpha_vec(Tvec; n, n_yolk)
        # clamp temperature in Celsius to avoid divergences
        T = clamp.(Tvec, 0, T_clamp_alpha)

        α_vec = α_albumen.(T) #[α_albumen(T) for t in T]
        ishift_test = 1
        i1, i2 = ishift_test + 1 + div(n - 1, 2) - div(n_yolk, 2), ishift_test + 1 + div(n - 1, 2) + div(n_yolk, 2)
        α_vec[i1:i2] .= α_yolk.(T[i1:i2])
        α_vec
    end

    FT = Real; # float type

    a, b, n = 0, 1, n_cells   # zmin, zmax
    n_yolk = Integer(round(yolk_diameter*n_cells/egg_diameter)) #11
    # n̂_min, n̂_max = -1, 1  # Outward facing unit vectors
    β, γ = 0, π;  # source term coefficients
    N_t = Integer(ceil(t_total / Δt));  # number of timesteps to take
    
    Δz = FT(b - a) / FT(n)
    Δz² = Δz^2;
    ∇²_op = [1 / Δz², -2 / Δz², 1 / Δz²];  # interior Laplacian operator

    S(z) = β * sin(γ * z)               # source term, (sin for easy integration)
    zf = range(a, b, length = n + 1);   # coordinates on cell faces

    # pulsing_period = 2*60;

    # α_vec = [α_albumen(80) for i in 1:n+1]
    # α_vec[1 + div(n - 1, 2) - div(n_yolk, 2):1 + div(n - 1, 2) + div(n_yolk, 2)] .= α_yolk(80)

    α_vec = dynamic_alpha_vec(ones(n+1)*80; n=n, n_yolk=n_yolk)

    # plot(α_vec, marker="o")

    # Initialize interior and boundary stencils:
    ∇² = Tridiagonal(ones(FT, n) .* ∇²_op[1],
    ones(FT, n + 1) .* ∇²_op[2],
    ones(FT, n) .* ∇²_op[3]);

    # Modify boundary stencil to account for BCs

    ∇².du[1] = 0  # modified stencil
    ∇².d[1] = 0 # to ensure `∂_t T = 0` at `z=0`
    ∇².dl[1] = 0  # to ensure `∂_t T = 0` at `z=0`

    ∇².du[n] = 0  # modified stencil
    ∇².d[n + 1] = 0 # to ensure `∂_t T = 0` at `z=zmax`
    ∇².dl[n] = 0  # to ensure `∂_t T = 0` at `z=zmax`

    # Calculate 
    D = α_vec .* ∇²

    T0 = zeros(FT, n + 1);
    T0 .= Tegg_init; 
    T0[1] = T_init; # set top BC
    T0[n + 1] = T_init; # set top BC
    
    cb = DiscreteCallback(condition, affect!)

    function rhs!(dT, T, params, t)
        n = params.n
        i = 1:n # interior domain
        
        Tsource = source_temp(params.thermal_state);
        AT_b = Tsource / Δz²;
        
        dT .= 0
        dT[2] = α_albumen(T[2]) * AT_b 
        dT[n] = α_albumen(T[n]) * AT_b

        dT
    end

    params = (; n, n_yolk, thermal_state=0.)

    tspan = (FT(0), N_t * FT(Δt))

    function update_func!(A, u, p, t)
        # @show u[20]
        A .= dynamic_alpha_vec(u; n_yolk=p.n_yolk, n=p.n) .* ∇²
        A
    end

    prob = SplitODEProblem(SciMLBase.DiffEqArrayOperator(D; update_func=update_func!),
        rhs!,
        T0,
        tspan,
        params)

    sol = solve(prob,
        alg;
        dt = Δt,
        saveat = range(FT(0), N_t * FT(Δt), length = N_t),
        # callback = cb,
        progress = true,
        progress_message = (dt, u, p, t) -> t,
        solver_opts...
    )

end

function thermodynamic_simulation_steps(
    # Drive parameters
    t_total::Float64, # seconds
    pulsing_period::Real, # seconds
    pulse_sequence::Array{Real}; # tempratures in Celsius
    solver_opts::Dict=Dict(),
    kwargs...
)

    # initial egg temperature in Celsius
    if T_init == nothing
        T_init = pulse_sequence[1]
    end

    @assert length(pulse_sequence) * pulsing_period >= t_total, "Pulse sequence too short for total simulation time"

    function temp_from_thermal_state(thermal_state)
        pulse_sequence[clamp(Integer(thermal_state)+1, 1, length(pulse_sequence))]
    end

    function condition(u, t, integrator)
        # mod(t, pulsing_period) == 0
        cycle_position = mod(t, pulsing_period)
        cycle_on_time = pulsing_period
        
        (cycle_position == cycle_on_time) || (cycle_position == 0)
    end

    function affect!(integrator)
        n = integrator.p[1]
        n_yolk = integrator.p[2]
        thermal_state = integrator.p[3] + 1
        integrator.p = (; n=n, n_yolk=n_yolk, thermal_state=thermal_state)

        T = temp_from_thermal_state(thermal_state)

        integrator.u[1] = T
        integrator.u[2] = T
        integrator.u[n] = T
        integrator.u[n+1] = T
    end

    solver_opts[:callback] = cb
    solver_opts[:tstops] = (0:pulsing_period:N_t * FT(Δt))

    return thermodynamic_simulation(
        t_total,
        temp_from_thermal_state,
        condition,
        affect!;
        Δt=Δt,
        Tegg_init=Tegg_init,
        n_cells=n_cells,
        scaler=scaler,
        scaler_asym=scaler_asym,
        egg_diameter=egg_diameter,
        yolk_diameter=yolk_diameter,
        T_init=T_init,
        T_clamp_alpha=T_clamp_alpha,
        solver_opts=solver_opts
    )

end

function thermodynamic_simulation_steps(
    # Drive parameters
    t_total::Float64,
    pulsing_period::Real,
    duty_cycle::Real,
    T_hot::Real,
    T_cold::Real;
    Δt::Float64,
    T_init::Union{Nothing, Real}=nothing,
    solver_opts::Dict=Dict(),
    kwargs...
)

    # initial egg temperature in Celsius
    if T_init == nothing
        T_init = T_cold
    end

    function temp_from_thermal_state(thermal_state)
        if thermal_state == 0
            return T_hot
        else
            return T_cold
        end
    end

    function condition(u, t, integrator)
        # mod(t, pulsing_period) == 0
        cycle_position = mod(t, pulsing_period)
        cycle_on_time = pulsing_period * duty_cycle
        
        (cycle_position == cycle_on_time) || (cycle_position == 0)
    end

    function affect!(integrator)
        n = integrator.p[1]
        n_yolk = integrator.p[2]
        thermal_state = mod(integrator.p[3] + 1, 2)
        integrator.p = (; n=n, n_yolk=n_yolk, thermal_state=thermal_state)

        T = temp_from_thermal_state(thermal_state)

        integrator.u[1] = T
        integrator.u[2] = T
        integrator.u[n] = T
        integrator.u[n+1] = T
    end

    if duty_cycle > 0
        cb = DiscreteCallback(condition, affect!)
        solver_opts[:callback] = cb
        solver_opts[:tstops] = (0:pulsing_period * duty_cycle:N_t * FT(Δt))
    else
        solver_opts[:callback] = nothing
    end

    return thermodynamic_simulation(
        t_total,
        temp_from_thermal_state,
        condition,
        affect!;
        Δt=Δt,
        T_init=T_init,
        kwargs...
    )

end

# commented out for now
# """
# Unitful variants of the thermodynamic simulation functions
# """
# function thermodynamic_simulation_steps(
#     t_total::Quantity{<:Real, 𝐓},
#     pulsing_period::Quantity{<:Real, 𝐓},
#     pulse_sequence::Array{<:Quantity{<:Real, 𝚯}};
#     T_init::Union{Nothing, Quantity{<:Real, 𝚯}}=nothing,
#     T_clamp_alpha::Quantity{<:Real, 𝚯}=200.°C,
#     Δt::Quantity{<:Real, 𝐓}=1.0s,
#     egg_diameter::Quantity{<:Real, 𝐋}=6e-2m,
#     yolk_diameter::Quantity{<:Real, 𝐋}=3e-2m,
#     Tegg_init::Quantity{<:Real, 𝚯}=20.°C,
#     kwargs...
# )
#     # Convert all unitful quantities to their underlying numerical values in base units (Celsius)
#     t_total = ustrip(t_total |> s)
#     pulsing_period = ustrip(pulsing_period |> s)
#     pulse_sequence = [ustrip(temp |> °C) for temp in pulse_sequence]
#     Δt = ustrip(Δt |> s)
#     Tegg_init = ustrip(Tegg_init |> °C)
#     egg_diameter = ustrip(egg_diameter |> m)
#     yolk_diameter = ustrip(yolk_diameter |> m)
#     T_clamp_alpha = ustrip(T_clamp_alpha |> °C)

#     if T_init !== nothing
#         T_init = ustrip(T_init |> °C)
#     end

#     # Call the original function with stripped values
#     thermodynamic_simulation(
#         t_total,
#         pulsing_period,
#         pulse_sequence;
#         T_init=T_init,
#         kwargs...
#     )
# end 

# function thermodynamic_simulation_steps(
#     t_total::Quantity{<:Real, 𝐓},
#     pulsing_period::Quantity{<:Real, 𝐓},
#     duty_cycle::Quantity{<:Real, 𝐓},
#     T_hot::Quantity{<:Real, 𝚯},
#     T_cold::Quantity{<:Real, 𝚯};
#     T_clamp_alpha::Quantity{<:Real, 𝚯}=200.°C,
#     Δt::Quantity{<:Real, 𝐓}=1.0s,
#     egg_diameter::Quantity{<:Real, 𝐋}=6e-2m,
#     yolk_diameter::Quantity{<:Real, 𝐋}=3e-2m,
#     Tegg_init::Quantity{<:Real, 𝚯}=20.°C,
#     kwargs...
# )
#     # Convert all unitful quantities to their underlying numerical values in base units (Celsius)
#     t_total = ustrip(t_total |> s)
#     pulsing_period = ustrip(pulsing_period |> s)
#     pulse_sequence = [ustrip(temp |> °C) for temp in pulse_sequence]
#     Δt = ustrip(Δt |> s)
#     Tegg_init = ustrip(Tegg_init |> °C)
#     egg_diameter = ustrip(egg_diameter |> m)
#     yolk_diameter = ustrip(yolk_diameter |> m)
#     T_clamp_alpha = ustrip(T_clamp_alpha |> °C)

#     if T_init !== nothing
#         T_init = ustrip(T_init |> °C)
#     end

#     # Call the original function with stripped values
#     thermodynamic_simulation(
#         t_total,
#         pulsing_period,
#         pulse_sequence;
#         Δt,
#         Tegg_init,
#         n_cells,
#         scaler,
#         scaler_asym,
#         egg_diameter,
#         yolk_diameter,
#         T_init,
#         T_clamp_alpha,
#         kwargs...
#     )
# end