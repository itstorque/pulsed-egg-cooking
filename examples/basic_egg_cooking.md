---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.16.7
  kernelspec:
    display_name: Julia 1.11.1
    language: julia
    name: julia-1.11
---

```julia
] activate ..
```

```julia
using PulsedEggSimulations
using Plots
```

# Hard-Boiling an Egg

```julia
sol = thermodynamic_simulation_steps(12*60., 2*60., 0., 100., 20.; Δt=1., n_cells=50);

sol2 = gelation_sim(sol);

fig = plot_temperature_and_gelation(sol, sol2)
xticks!(0:2:12)

fig
```

# Soft-Boiling an Egg

```julia
sol = thermodynamic_simulation_steps(8*60., 2*60., 0., 100., 20.; Δt=1., n_cells=50);

sol2 = gelation_sim(sol);

fig = plot_temperature_and_gelation(sol, sol2)
xticks!(0:2:12)

fig
```
