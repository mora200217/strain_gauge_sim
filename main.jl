"""
Load Cell Lab – Improved Version

FEM simulation of a load cell with strain gauge.
Average longitudinal strain → R(t)
Exported to ParaView, LTSpice and MATLAB.

Author: Andrés Morales
Date: 18-Dec-2025
"""

using Gridap
using GridapGmsh
using CairoMakie
using LinearAlgebra
using DelimitedFiles

# ============================================================
# 1. Model & FE spaces
# ============================================================

model = GmshDiscreteModel("step/load_cell_beam.msh")

order = 1
reffe = ReferenceFE(lagrangian, VectorValue{3,Float64}, order)

V = TestFESpace(
    model, reffe;
    conformity = :H1,
    dirichlet_tags = ["fixed"]
)

U = TrialFESpace(V)

# ============================================================
# 2. Domains
# ============================================================

Ω   = Triangulation(model)
dΩ  = Measure(Ω, 3)

ΓL  = BoundaryTriangulation(model; tags=["load"])
dΓ  = Measure(ΓL, 2)

Ωg  = Triangulation(model; tags=["strain_gauge"])
dΩg = Measure(Ωg, 3)

# ============================================================
# 3. Material model (Linear Elasticity)
# ============================================================

E, ν = 2.0e11, 0.3

λ = (E*ν)/((1+ν)*(1-2ν))
μ = E/(2*(1+ν))

I = one(TensorValue{3,3,Float64})

ε(u) = 0.5 * (∇(u) + transpose(∇(u)))
σ(ε) = λ*tr(ε)*I + 2μ*ε

# ============================================================
# 4. Strain Gauge Model
# ============================================================

d̂ = normalize(VectorValue(1.0, 0.0, 0.0))
ε_long(u) = d̂ ⋅ ε(u) ⋅ d̂

R0 = 1200.0
GF = 2.0

# ============================================================
# 5. Weak form (time-independent part)
# ============================================================

a(u,v) = ∫( σ(ε(u)) ⊙ ε(v) ) * dΩ

# Load direction
n̂ = VectorValue(0.0, 0.0, -1.0)

# ============================================================
# 6. Time loop
# ============================================================

dt, T = 0.01, 1.0
times = collect(0:dt:T)
Rhist = zeros(length(times))

uh_last = nothing

for (i,t) in enumerate(times)

    F = 1.0e6 * sin(2π*t)
    f(x) = F * n̂

    l(v) = ∫( f ⋅ v ) * dΓ

    op = AffineFEOperator(a, l, U, V)
    uh = solve(op)
    uh_last = uh

    # Average strain in gauge volume
    εg_avg = ∫( ε_long(uh) ) * dΩg / ∫( 1.0 ) * dΩg

    Rhist[i] = R0 * (1 + GF * εg_avg)
end

# ============================================================
# 7. Export to ParaView
# ============================================================

writevtk(
    Ω,
    "load_cell_results";
    cellfields = Dict(
        "eps_long" => ε_long(uh_last),
        "strain"   => ε(uh_last)
    ),
    pointfields = Dict(
        "displacement" => uh_last
    )
)

println("✔ VTK exported")

# ============================================================
# 8. Plot
# ============================================================

fig = Figure(resolution=(800,400))
ax = Axis(
    fig[1,1],
    xlabel="Time [s]",
    ylabel="Resistance [Ω]",
    title="Strain Gauge Response"
)

lines!(ax, times, Rhist, linewidth=2)
fig

# ============================================================
# 9. Export data
# ============================================================

writedlm("resistance.txt", hcat(times, Rhist))
println("✔ Exported: resistance.txt")
