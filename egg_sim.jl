using LinearAlgebra
using DiffEqBase
using OrdinaryDiffEq: SplitODEProblem, solve, IMEXEuler, KenCarp4
import SciMLBase
using Plots
using Statistics

ρ_yolk = T -> 1037.3 - 0.1386 * (T) - 0.0023*(T)^2   # kg/m^3
ρ_albumen = T -> 1043.3 - 0.0115 * (T) - 0.0041*(T)^2  # kg/m^3
c_p_yolk = 3120  # kJ/(kg·°C)
c_p_albumen = 3800  # kJ/(kg·°C)
k_yolk = T -> 0.0008*(T) + 0.395 #0.5  # W/(m·°C)
k_albumen = T -> 0.0013*(T) + 0.5125 #0.6  # W/(m·°C)


# renormalize for egg lengthscales
# 62mm	43mm
egg_diameter = 0.55*43e-3 # 43mm

α_yolk = T -> k_yolk(T) / (ρ_yolk(T) * c_p_yolk * egg_diameter^2) # m^2/s
α_albumen = T -> k_albumen(T) / (ρ_albumen(T) * c_p_albumen * egg_diameter^2) # m^2/s

function dynamic_alpha_vec(Tvec; n, n_yolk)

    T = clamp.(Tvec, 0, 200)

    α_vec = α_albumen.(T) #[α_albumen(T) for t in T]
    i1, i2 = 1 + div(n - 1, 2) - div(n_yolk, 2), 1 + div(n - 1, 2) + div(n_yolk, 2)
    α_vec[i1:i2] .= α_yolk.(T[i1:i2])
    α_vec
end

function run_simulation(t_total, pulsing_period, duty_cycle; T_hot, T_cold, Δt, α_yolk=α_yolk, α_albumen=α_albumen, Tegg_init=20)

    FT = Float64; # float type

    a, b, n = 0, 1, 50   # zmin, zmax, number of cells
    n_yolk = 11
    n̂_min, n̂_max = -1, 1  # Outward facing unit vectors
    β, γ = 0, π;  # source term coefficients
    N_t = Integer(ceil(t_total / Δt));  # number of timesteps to take
    
    Δz = FT(b - a) / FT(n)
    Δz² = Δz^2;
    ∇²_op = [1 / Δz², -2 / Δz², 1 / Δz²];  # interior Laplacian operator

    S(z) = β * sin(γ * z)               # source term, (sin for easy integration)
    zf = range(a, b, length = n + 1);   # coordinates on cell faces

    # pulsing_period = 2*60;

    α_vec = [α_albumen(80) for i in 1:n+1]
    α_vec[1 + div(n - 1, 2) - div(n_yolk, 2):1 + div(n - 1, 2) + div(n_yolk, 2)] .= α_yolk(80)

    # plot(α_vec, marker="o")

    # Initialize interior and boundary stencils:
    ∇² = Tridiagonal(ones(FT, n) .* ∇²_op[1],
    ones(FT, n + 1) .* ∇²_op[2],
    ones(FT, n) .* ∇²_op[3]);

    # Modify boundary stencil to account for BCs

    ∇².du[1] = 0  # modified stencil
    ∇².d[1] = 0 # to ensure `∂_t T = 0` at `z=0`
    ∇².dl[1] = 0  # to ensure `∂_t T = 0` at `z=0`

    # Modify boundary stencil to account for BCs
    ∇².du[n] = 0  # modified stencil
    ∇².d[n + 1] = 0 # to ensure `∂_t T = 0` at `z=zmax`
    ∇².dl[n] = 0  # to ensure `∂_t T = 0` at `z=zmax`

    # Calculate 
    D = α_vec .* ∇²

    T0 = zeros(FT, n + 1);
    T0 .= Tegg_init; 
    T0[1] = T_hot; # set top BC
    T0[n + 1] = T_hot; # set top BC

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

    @assert duty_cycle >= 0 && duty_cycle < 0.99
    
    if duty_cycle >= 0.01
        cb = DiscreteCallback(condition, affect!)
        opts = Dict(:callback => cb, :tstops => (0:pulsing_period * duty_cycle:N_t * FT(Δt)))
    else
        opts = Dict(:callback => nothing)
    end

    function rhs!(dT, T, params, t)
        n = params.n
        i = 1:n # interior domain
        
        source_temp = temp_from_thermal_state(params.thermal_state);
        AT_b = source_temp / Δz²;
        
        dT .= 0
        dT[2] = α_albumen(T[2]) * AT_b 
        dT[n] = α_albumen(T[n]) * AT_b

        dT
    end

    params = (; n, n_yolk, thermal_state=0.)

    tspan = (FT(0), N_t * FT(Δt))


    function update_func!(A, u, p, t)
        A .= dynamic_alpha_vec(u; n_yolk=p.n_yolk, n=p.n) .* ∇²
        A
    end

    prob = SplitODEProblem(SciMLBase.DiffEqArrayOperator(D; update_func=update_func!),
        rhs!,
        T0,
        tspan,
        params)
    alg = IMEXEuler(autodiff=false);
    println("Solving...")
    sol = solve(prob,
        alg;
        dt = Δt,
        saveat = range(FT(0), N_t * FT(Δt), length = N_t),
        # callback = cb,
        progress = true,
        progress_message = (dt, u, p, t) -> t,
        opts...
        # tstops = (duty_cycle != 0) ? (0:pulsing_period * duty_cycle:N_t * FT(Δt)) : nothing,
    )

end


sol = run_simulation(12*60, 2*60, 0.5; T_hot=100., T_cold=30., Δt=1);


plot_in_kelvin = true
temperature_offset = plot_in_kelvin ? 273.15 : 0

# Plot drive
Plots.plot(sol.t/60, temperature_offset .+ hcat(sol.u...)[1, :])
Plots.plot!(sol.t/60,temperature_offset .+ hcat(sol.u...)[n+1, :])


## Plot Time sols 
Plots.plot(; xlabel="Time [min]", 
ylabel = plot_in_kelvin ? "Temperature [K]" : "Temperature [°C]", 
title="Cross-Section Temperature vs. Time")

# Plots.plot!(sol.t/60, hcat(sol.u...)[2, :], label="0mm", color=:grey)
# Plots.plot!(sol.t, hcat(sol.u...)[3, :], label="albumen", color=:grey)
Plots.plot!(sol.t/60, temperature_offset .+ hcat(sol.u...)[7, :], label="albumen", color=:grey)
Plots.plot!(sol.t/60, temperature_offset .+ hcat(sol.u...)[Integer(floor(n/4)), :], label="albumen", color=:grey)
# Plots.plot!(hcat(sol.u...)[5, :], label="albumen", color=:grey)
# Plots.plot!(hcat(sol.u...)[Integer(floor(n/2)) - 2, :], label="yolk", color=:orange)
Plots.plot!(sol.t/60, temperature_offset .+ hcat(sol.u...)[Integer(floor(n/2)) - 1, :], label="yolk", color=:orange)
Plots.plot!(sol.t/60, temperature_offset .+ hcat(sol.u...)[Integer(floor(n/2)) + 0, :], label="yolk", color=:orange)
Plots.plot!(sol.t/60, temperature_offset .+ hcat(sol.u...)[Integer(floor(n/2)) + 1, :], label="yolk", color=:orange)
# Plots.plot!(hcat(sol.u...)[Integer(floor(n/2)) + 2, :], label="yolk", color=:orange)
# Plots.plot!(hcat(sol.u...)[n-3, :], label="albumen", color=:black)
# Plots.plot!(sol.t, hcat(sol.u...)[n-2, :], label="albumen", color=:black)
# Plots.plot!(sol.t, hcat(sol.u...)[n-1, :], label="albumen", color=:black)
# Plots.plot!(hcat(sol.u...)[n, :], label="albumen", color=:black)

ylims!(280,  390)
yticks!(280:20:390)

vline!([4])
hline!([318])

## Plot Agg
plot(zf[2:end-1], median(hcat(sol.u...)[2:end-1, :]; dims=2), ylabel="Aggregate Simulated Temperature", label="Median Temperature", xlabel="Position in Egg")
plot!(zf[2:end-1], mean(hcat(sol.u...)[2:end-1, :]; dims=2), label="Mean Temperature")

vline!(zf[[1 + div(n - 1, 2) - div(n_yolk, 2), 1 + div(n - 1, 2) + div(n_yolk, 2)]], label="Egg Yolk")


# Plot Slices vs X
p1 = Plots.plot(ylabel="Temperature", xlabel="Cross-Section of Egg", legend_title="Time")
for (idx, tvec) in enumerate(sol.u[1:Integer(floor(end//10)):end])
    Plots.plot!(zf[2:end-1], tvec[2:end-1], label = idx, markershape = :diamond)
end
p1