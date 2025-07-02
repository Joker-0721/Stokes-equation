using ApproxOperator
using Pardiso
using WriteVTK ,XLSX
using SparseArrays, LinearAlgebra
import ApproxOperator.Stokes:∫∫μ∇u∇vdxdy
import ApproxOperator.Elasticity:∫∫p∇udxdy, ∫vᵢtᵢds, ∫∫vᵢbᵢdxdy, ∫vᵢgᵢds, ∫qpdΩ

include("import_cavity.jl")

ps = MKLPardisoSolver()

ndiv = 4

elements, nodes, nodes_p= import_cavity_RI("Stokes-equation/msh/cav_quad_"*string(ndiv)*".msh", "Stokes-equation/msh/cav_quad_"*string(ndiv)*".msh");

nᵘ = length(nodes)
nᵖ = length(nodes_p)
# nₑ = length(elements["Ω"])

E = 1.0
ν = 1.0
μ = 0.5
b₁ = 0.5
b₂ = 0.5
t₁ = 0.5
t₂ = 0.5
g₁ = 0.0
g₂ = 0.0
P = 80.0

n₁₁(n₁,n₂) = 1.0
n₂₂(n₁,n₂) = 1.0
prescribe!(elements["Ωᵘ"],:μ=>(x,y,z)->μ)
prescribe!(elements["Ωᵘ"],:b₁=>(x,y,z)->b₁)
prescribe!(elements["Ωᵘ"],:b₂=>(x,y,z)->b₂)
prescribe!(elements["Ωᵖ"],:P=>(x,y,z)->P)
prescribe!(elements["Ωᵖ"],:b₁=>(x,y,z)->b₁)
prescribe!(elements["Ωᵖ"],:b₂=>(x,y,z)->b₂)
prescribe!(elements["Ωᵘ"],:E=>(x,y,z)->E)
prescribe!(elements["Ωᵘ"],:ν=>(x,y,z)->ν)
prescribe!(elements["Ωᵖ"],:E=>(x,y,z)->E)
prescribe!(elements["Ωᵖ"],:ν=>(x,y,z)->ν)
prescribe!(elements["Γ₁"],:t₁=>(x,y,z)->t₁)
prescribe!(elements["Γ₁"],:t₂=>(x,y,z)->t₂)
prescribe!(elements["Γ₂"],:t₁=>(x,y,z)->t₁)
prescribe!(elements["Γ₂"],:t₂=>(x,y,z)->t₂)
prescribe!(elements["Γ₃"],:t₁=>(x,y,z)->t₁)
prescribe!(elements["Γ₃"],:t₂=>(x,y,z)->t₂)
prescribe!(elements["Γ₄"],:t₁=>(x,y,z)->t₁)
prescribe!(elements["Γ₄"],:t₂=>(x,y,z)->t₂)
prescribe!(elements["Γ₁"],:g₁=>(x,y,z)->g₁)
prescribe!(elements["Γ₂"],:g₁=>(x,y,z)->g₁)
prescribe!(elements["Γ₃"],:g₁=>(x,y,z)->g₁)
prescribe!(elements["Γ₄"],:g₁=>(x,y,z)->1.0)
prescribe!(elements["Γ₁"],:g₂=>(x,y,z)->g₂)
prescribe!(elements["Γ₂"],:g₂=>(x,y,z)->g₂)
prescribe!(elements["Γ₃"],:g₂=>(x,y,z)->g₂)
prescribe!(elements["Γ₄"],:g₂=>(x,y,z)->1.0)
prescribe!(elements["Γ₁"],:α=>(x,y,z)->1e12*E)
prescribe!(elements["Γ₂"],:α=>(x,y,z)->1e12*E)
prescribe!(elements["Γ₃"],:α=>(x,y,z)->1e12*E)
prescribe!(elements["Γ₄"],:α=>(x,y,z)->1e12*E)
prescribe!(elements["Γ₁"],:n₁₁=>(x,y,z,n₁,n₂)->n₁₁(n₁,n₂))
prescribe!(elements["Γ₁"],:n₂₂=>(x,y,z,n₁,n₂)->n₂₂(n₁,n₂))
prescribe!(elements["Γ₁"],:n₁₂=>(x,y,z)->0.0)
prescribe!(elements["Γ₂"],:n₁₁=>(x,y,z,n₁,n₂)->n₁₁(n₁,n₂))
prescribe!(elements["Γ₂"],:n₂₂=>(x,y,z,n₁,n₂)->n₂₂(n₁,n₂))
prescribe!(elements["Γ₂"],:n₁₂=>(x,y,z)->0.0)
prescribe!(elements["Γ₃"],:n₁₁=>(x,y,z,n₁,n₂)->n₁₁(n₁,n₂))
prescribe!(elements["Γ₃"],:n₂₂=>(x,y,z,n₁,n₂)->0.0)
prescribe!(elements["Γ₃"],:n₁₂=>(x,y,z)->0.0)
prescribe!(elements["Γ₄"],:n₁₁=>(x,y,z,n₁,n₂)->n₁₁(n₁,n₂))
prescribe!(elements["Γ₄"],:n₂₂=>(x,y,z,n₁,n₂)->n₂₂(n₁,n₂))
prescribe!(elements["Γ₄"],:n₁₂=>(x,y,z)->0.0)

aᵘ = ∫∫μ∇u∇vdxdy => elements["Ωᵘ"]
aᵖ = ∫qpdΩ => elements["Ωᵖ"]
bᵖ = ∫∫p∇udxdy => (elements["Ωᵖ"],elements["Ωᵘ"])
# ∫∫vᵢbᵢdxdy => elements["Ω"]
f = [
    ∫∫vᵢbᵢdxdy => elements["Ωᵘ"],
    ∫vᵢtᵢds => elements["Γ₁"],
    ∫vᵢtᵢds => elements["Γ₂"],
    ∫vᵢtᵢds => elements["Γ₃"],
    ∫vᵢtᵢds => elements["Γ₄"]
]
aᵅ = [
    ∫vᵢgᵢds => elements["Γ₁"],
    ∫vᵢgᵢds => elements["Γ₂"],
    ∫vᵢgᵢds => elements["Γ₃"],
    ∫vᵢgᵢds => elements["Γ₄"]
]

kᵘᵘ = zeros(2*nᵘ,2*nᵘ)
kᵖᵘ = zeros(nᵖ,2*nᵘ)
kᵖᵖ = zeros(nᵖ,nᵖ)
fᵖ = zeros(nᵖ)
fᵘ = zeros(2*nᵘ)
d = zeros(2*nᵘ+nᵖ)

aᵘ(kᵘᵘ)
aᵖ(kᵖᵖ)
bᵖ(kᵖᵘ)
f(fᵘ)
aᵅ(kᵘᵘ,fᵘ)

k =[kᵘᵘ kᵖᵘ';kᵖᵘ kᵖᵖ]
f = [fᵘ;fᵖ]

d = k\f

# points = zeros(3,nᵖ)
# for (i,node) in enumerate(nodes)
#     points[1,i] = node.x
#     points[2,i] = node.y
#     points[3,i] = node.d
#     # points[3,i] = us[i]*4
# end
# cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE_STRIP,[x.𝐼 for x in elm.𝓒]) for elm in elements["Ωᵘ"]]
# # vtk_grid("./vtk/hmd_2d/error/non_uniform_Tri3_"*string(ndiv)*".vtu",points,cells) do vtk
# vtk_grid("./png/cav_mix_"*string(ndiv)*".vtu",points,cells) do vtk
#     # vtk["d"] = [node.d for node in nodes]
#     vtk["精确解"] = us
# end


𝑢₁ = d[1:2:2*nᵘ]
𝑢₂ = d[2:2:2*nᵘ]
# 𝑢₃ = d[3:3:3*nᵘ]
𝑝 = d[2*nᵘ+1:2*nᵘ+nᵖ]
push!(nodes,:u₁=>𝑢₁,:u₂=>𝑢₂)
push!(nodes_p,:p=>𝑝)

colors = zeros(nᵘ)
𝗠 = zeros(10)
for (i,node) in enumerate(nodes)
    x = node.x
    y = node.y
    indices = sp(x,y)
    ni = length(indices)
    𝓒 = [nodes_p[i] for i in indices]
    data = Dict([:x=>(2,[x]),:y=>(2,[y]),:z=>(2,[z]),:𝝭=>(4,zeros(ni)),:𝗠=>(0,𝗠)])
    ξ = 𝑿ₛ((𝑔=1,𝐺=1,𝐶=1,𝑠=0), data)
    𝓖 = [ξ]
    a = type(𝓒,𝓖)
    set𝝭!(a)
    p = 0.0
    N = ξ[:𝝭]
    for (k,xₖ) in enumerate(𝓒)
        p += N[k]*xₖ.p
    end
    colors[i] = p
end
α = 1.0
points = [[node.x+α*node.u₁ for node in nodes]';[node.y+α*node.u₂ for node in nodes]']
cells = [MeshCell(VTKCellTypes.VTK_TETRA,[xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements["Ωᵘ"]]
vtk_grid("./png/cav_mix_"*string(ndiv)*"_"*string(nᵖ),points,cells) do vtk
    vtk["u"] = (𝑢₁,𝑢₂,𝑢₃)
    vtk["𝑝"] = colors
end

XLSX.openxlsx("./png/cav_mix.xlsx", mode = "rw") do xf
    sheet = xf[1]
    for (n,u) in zip(indices,uₐ)
        sheet["A"*string(n)] = n
        sheet["B"*string(n)] = u
    end
end