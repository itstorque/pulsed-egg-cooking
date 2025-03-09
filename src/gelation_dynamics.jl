function gelation_sim(sol)
    n, n_yolk, _ = sol.prob.p

    i1, i2 = 1 + div(n - 1, 2) - div(n_yolk, 2), 1 + div(n - 1, 2) + div(n_yolk, 2)

    Ai_vec = Ai_albumen * ones(n+1)
    Ai_vec[i1:i2] .= Ai_yolk

    Eai_vec = Eai_albumen * ones(n+1)
    Eai_vec[i1:i2] .= Eai_yolk

    function cooking!(du, u, p, t)
        temp = sol(t) .+ 273.15
        du .= Ai_vec .* exp.(-Eai_vec ./ (R .* (temp))) .* (1 .- u)
    end

    u0 = zeros(length(sol.u[1]))
    cooking_prob = ODEProblem(cooking!, u0, (sol.t[1], sol.t[end]))
    sol2 = solve(cooking_prob, Rosenbrock23(autodiff = false), dt=0.1)
    
    sol2
end