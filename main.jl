"""
Load Cell Lab – Corrected Version

Quasi-static FEM simulation of a load cell with strain gauge.
Average longitudinal strain → R(F)

Author: Andrés Morales
Date: 18-Dec-2025
"""

using Gridap
using GridapGmsh
using CairoMakie
using LinearAlgebra
using DelimitedFiles
using Gridap.Visualization

# ============================================================
# 1. Model & FE spaces
# ============================================================

model = GmshDiscreteModel("mesh/load_cell_beam.msh")

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
# 3. Linear elastic material (isotropic)
# ============================================================

E  = 60e9      # Pa (Almunioum)
ν  = 0.3

λ = (E*ν)/((1+ν)*(1-2ν))
μ = E/(2*(1+ν))

I = one(TensorValue{3,3,Float64})

ε(u) = 0.5*(∇(u) + transpose(∇(u)))
σ(ε) = λ*tr(ε)*I + 2μ*ε

# ============================================================
# 4. Strain gauge model (uniaxial)
# ============================================================

d̂ = normalize(VectorValue(1.0,1.0,0.0))
ε_long(u) = d̂ ⋅ ε(u) ⋅ d̂

R0 = 1200.0
GF = 2.0

# ============================================================
# 5. Weak form (static equilibrium)
# ============================================================

a(u,v) = ∫( σ(ε(u)) ⊙ ε(v) ) * dΩ

n̂ = VectorValue(0.0,0.0,-1.0)

# ============================================================
# 6. Load sweep (quasi-static)
# ============================================================

forces = range(0, 50, length=40)   # 0–2000 N
Rhist  = zeros(length(forces))

uh_last = nothing

for (i,F) in enumerate(forces)

  f(x) = F*n̂
  l(v) = ∫( f ⋅ v ) * dΓ

  op = AffineFEOperator(a, l, U, V)
  uh = solve(op)
  uh_last = uh
 
  εg = sum(∫( ε_long(uh) ) * dΩg) 
  Rhist[i] = R0*(1 + GF*εg)
  @show F εg 
end

# ============================================================
# 7. Export to ParaView
# ============================================================

writevtk(
  Ω,
  "load_cell_results",
  cellfields = Dict(
    "eps_long" => CellField(ε_long(uh_last), Ω),
    "strain"   => CellField(ε(uh_last), Ω)
  )
)

println("✔ VTK exported")

# ============================================================
# 8. Plot sensor response
# ============================================================

fig = Figure(resolution=(800,400))
ax = Axis(
  fig[1,1],
  xlabel = "Force [N]",
  ylabel = "Resistance [Ω]",
  title  = "Load Cell – Strain Gauge Response"
)

scatter!(ax, forces, Rhist)
fig

# ============================================================
# 9. Export data
# ============================================================

writedlm("resistance.txt", hcat(forces, Rhist))
println("✔ Exported resistance.txt")
