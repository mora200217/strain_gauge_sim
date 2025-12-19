using Gridap
using GridapGmsh
using CairoMakie
using LinearAlgebra
using DelimitedFiles
using Gridap.Visualization
using Random

model = GmshDiscreteModel("mesh/load_cell_beam.msh")

order = 1
reffe = ReferenceFE(lagrangian, VectorValue{3,Float64}, order)

V = TestFESpace(model, reffe; conformity=:H1, dirichlet_tags=["fixed"])
U = TrialFESpace(V)

Ω  = Triangulation(model)
dΩ = Measure(Ω, 3)

ΓL = BoundaryTriangulation(model; tags=["load"])
dΓ = Measure(ΓL, 2)

Ωg  = Triangulation(model; tags=["strain_gauge"])
dΩg = Measure(Ωg, 3)

E, ν = 60e9, 0.3
λ = (E*ν)/((1+ν)*(1-2ν))
μ = E/(2*(1+ν))

I = one(TensorValue{3,3,Float64})

ε(u) = 0.5*(∇(u) + transpose(∇(u)))
σ(ε) = λ*tr(ε)*I + 2μ*ε

d̂ = normalize(VectorValue(1.0,1.0,0.0))
ε_long(u) = d̂ ⋅ ε(u) ⋅ d̂

R0, GF = 120.0, 2.0

a(u,v) = ∫( σ(ε(u)) ⊙ ε(v) ) * dΩ
n̂ = VectorValue(0.0,0.0,-1.0)

# Ruido térmico (Johnson)
kB = 1.380649e-23
T  = 300.0
B  = 100.0
σR = sqrt(4*kB*T*R0*B)

Random.seed!(3)

forces = range(0, 50, length=40)
Rhist  = zeros(length(forces))

uh_last = nothing

for (i,F) in enumerate(forces)

  f(x) = F*n̂
  l(v) = ∫( f ⋅ v ) * dΓ

  uh = solve(AffineFEOperator(a, l, U, V))
  uh_last = uh

  εg = sum(∫( ε_long(uh) ) * dΩg)

  noise = σR * randn()
  Rhist[i] = R0*(1 + GF*εg) + noise
end


writevtk(
  Ω, "load_cell_results",
  cellfields = Dict(
    "eps_long" => CellField(ε_long(uh_last), Ω),
    "strain"   => CellField(ε(uh_last), Ω)
  )
)

fig = Figure(resolution=(800,400))
ax = Axis(fig[1,1], xlabel="Fuerza aplicada [N]", ylabel="Resistencia en galga [Ω]")
ax.title = "Simulación de variación resistiva con carga"
scatter!(ax, forces, Rhist)

fig

writedlm("resistance.txt", hcat(forces, Rhist))
# --- Guardar para LTspice ---

# Tomar solo los primeros 10 valores
first10 = Rhist[1:10]

# Convertirlos en string separados por comas
param_list = join(first10, ", ")

# Crear el texto para LTspice
ltspice_txt = """
.STEP PARAM Sensor LIST $param_list
"""

# Guardar en archivo
open("spice/resistance_param.inc", "w") do io
    write(io, ltspice_txt)
end

println("✔ Archivo guardado con los primeros 10 valores")
