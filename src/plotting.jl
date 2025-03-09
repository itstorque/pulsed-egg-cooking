function crosssection_variables(sol)
    
    positions = [
        2,
        1 + Integer(round(sol.prob.p.n / 2 * 0.2)),
        1 + Integer(round(sol.prob.p.n / 2 * 0.36)),
        1 + Integer(round(sol.prob.p.n / 2 * 0.44)),

        1 + Integer(round(sol.prob.p.n / 2 * 0.6)),
        1 + Integer(round(sol.prob.p.n / 2 * 0.8)),
        1 + Integer(floor(sol.prob.p.n/2)) + 0
    ];

    labels = [
        "25mm",
        "20mm",
        "16mm",
        "14mm",
        "10mm",
        "5mm",
        "0mm"
    ];

    colors = [
        :black,
        :grey,
        :grey,
        :goldenrod1,
        :orange,
        :darkorange1,
        :orangered
    ];

    styles = [
        :solid,
        :solid,
        :dash,
        :dash,
        :solid,
        :solid,
        :solid
    ];

    positions, labels, colors, styles

end

function plot_temperature(sol; plot_in_kelvin=false)

    temperature_offset = plot_in_kelvin ? 273.15 : 0.0

    positions, labels, colors, styles = crosssection_variables(sol);

    p1 = Plots.plot(; xlabel="Time [min]", 
        ylabel = plot_in_kelvin ? "Temperature [K]" : "Temperature [°C]",
        title="Cross-Section Temperature vs. Time"
    );

    for (i, pos) in enumerate(positions)
        Plots.plot!(sol.t/60, temperature_offset .+ hcat(sol.u...)[pos, :], label=labels[i], color=colors[i], style=styles[i])
    end

    p1

end

function plot_gelation(sol, sol2)

    positions, labels, colors, styles = crosssection_variables(sol);

    p2 = plot(title="Gelation", xlabel="Time [min]", ylabel="Degree of Cooking", legend = false, ylim=(0, 1));

    for (i, pos) in enumerate(positions)
        plot!(sol2.t/60, hcat(sol2.u...)[pos, :], label=labels[i], color=colors[i], style=styles[i])
    end

    p2

end

function plot_temperature_and_gelation(sol, sol2; plot_in_kelvin=false)

    p1 = plot_temperature(sol; plot_in_kelvin=false);
    p2 = plot_gelation(sol, sol2);
    
    plot(p1, p2, layout=(1, 2), size=(800, 300), dpi=300, left_margin = 5Plots.mm, bottom_margin=5Plots.mm)
end