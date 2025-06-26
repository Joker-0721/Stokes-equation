using ApproxOperator, XLSX
using WriteVTK
using CairoMakie
using SparseArrays
import BenchmarkExample: BenchmarkExample
include("import_cavity.jl")
ndiv   = 4
ndivs  = 4
# ndivs2 = 16

elements, nodes, nodes_s= import_cavity_RI("Stokes-equation/msh/cav_quad_"*string(ndiv)*".msh", "Stokes-equation/msh/cav_quad_"*string(ndiv)*".msh");

nᵘ = length(nodes)
nᵖ = length(nodes_s)
nₑ = length(elements["Ω"])
# nₑₛ = length(elements["Ω"])

E = 1.0
ν = 1.0
μ = 0.5
b₁ = 0.5
b₂ = 0.5
t₁ = 0.5
t₂ = 0.5

set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
# set𝝭!(elements["Ωˢ"])
# set∇𝝭!(elements["Ωˢ"])
set𝝭!(elements["Γ₁"])
set𝝭!(elements["Γ₂"])
set𝝭!(elements["Γ₃"])
set𝝭!(elements["Γ₄"])
prescribe!(elements["Ω"],:b₁=>(x,y,z)->b₁)
prescribe!(elements["Ω"],:b₂=>(x,y,z)->b₂)
prescribe!(elements["Γ₁"],:t₁=>(x,y,z)->t₁)
prescribe!(elements["Γ₁"],:t₂=>(x,y,z)->t₂)
prescribe!(elements["Γ₂"],:t₁=>(x,y,z)->t₁)
prescribe!(elements["Γ₂"],:t₂=>(x,y,z)->t₂)
prescribe!(elements["Γ₃"],:t₁=>(x,y,z)->t₁)
prescribe!(elements["Γ₃"],:t₂=>(x,y,z)->t₂)
prescribe!(elements["Γ₄"],:t₁=>(x,y,z)->t₁)
prescribe!(elements["Γ₄"],:t₂=>(x,y,z)->t₂)
prescribe!(elements["Γ₁"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₂"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₃"],:g=>(x,y,z)->0.0)
prescribe!(elements["Γ₄"],:g=>(x,y,z)->0.0)

ops = [
    Operator{:∫∫μ∇u∇vdxdy}(:μ=>μ),
    Operator{:∫∫p∇vdxdy}(),
    Operator{:∫∫vᵢbᵢdxdy}(),
    Operator{:∫vᵢtᵢds}(),
]

kᵘ = zeros(2*nᵘ,2*nᵘ)
kᵘᵖ = zeros(nᵖ,2*nᵘ)
kᵖ = zeros(nᵖ,nᵖ)
f = zeros(2*nᵘ)
# d = zeros(3*nᵇ+2*nˢ)

ops[1](elements["Ω"],kᵘ)
ops[2](elements["Ω"],elements["Ωˢ"],kᵘᵖ)
# ops[3](elements["Ωˢ"],kᵖ)
ops[3](elements["Ω"],f)
ops[4](elements["Γ₁"],f)
ops[4](elements["Γ₂"],f)
ops[4](elements["Γ₃"],f)
ops[4](elements["Γ₄"],f)
# ops[4](elements["Γ₂"],f)
# ops[4](elements["Γ₃"],f)
# ops[4](elements["Γ₄"],f)
# ops[4](elements["Γ₁"],f)
# ops[4](elements["Γ₂"],f)
# ops[4](elements["Γ₃"],f)
# ops[4](elements["Γ₄"],f)


k = [kᵘ kᵘᵖ';kᵘᵖ kᵖ]
f = [f;zeros(nᵖ)]

# k = kʷˢ*inv(kˢˢ)*kʷˢ'
# k = -kʷˢ*(kˢˢ\kʷˢ')
# a = eigvals(k)
# println(log10(a[3*nᵇ-2nˢ+1]))
# println(a[3*nᵇ-2nˢ+1])

# d = k\f
# pardiso(ps,d,k,f)
# d₁ = d[1:3:3*nᵘ]
# d₂ = d[2:3:3*nᵘ] 
# # d₃ = d[3:3:3*nᵇ]
# s₁ = d[3*nᵇ+1:2:3*nᵇ+2*nˢ]
# # s₂ = d[3*nᵇ+2:2:3*nᵇ+2*nˢ]

# push!(nodes,:d₁=>d₁,:d₂=>d₂,:d₃=>d₃)
# push!(nodes_s,:q₁=>s₁,:q₂=>s₂)
# # eval(VTK_mix_pressure)

# set𝝭!(elements["Ωᵍ"])
# set∇𝝭!(elements["Ωᵍ"])
# set𝝭!(elements["Ωᵍˢ"])
# set∇𝝭!(elements["Ωᵍˢ"])

# prescribe!(elements["Ωᵍ"],:u=>(x,y,z)->w(x,y))
# prescribe!(elements["Ωᵍ"],:θ₁=>(x,y,z)->θ₁(x,y))
# prescribe!(elements["Ωᵍ"],:θ₂=>(x,y,z)->θ₂(x,y))
# prescribe!(elements["Ωᵍˢ"],:Q₁=>(x,y,z)->Q₁(x,y))
# prescribe!(elements["Ωᵍˢ"],:Q₂=>(x,y,z)->Q₂(x,y))
# L₂_u = ops[8](elements["Ωᵍ"])
# L₂_q = ops[9](elements["Ωᵍˢ"])
# a = log10(L₂_u)
# b = log10(L₂_q)
# println(a)
# println(b)
# index = 1:20
# XLSX.openxlsx("./xlsx/SquarePlate.xlsx", mode="rw") do xf
#     Sheet = xf[6]
#     ind = findfirst(n->n==ndivs,index)+2
#     Sheet["J"*string(ind)] = nˢ
#     Sheet["K"*string(ind)] = a
#     Sheet["L"*string(ind)] = b
# end

# println(wᶜ)
# e = abs(wᶜ[1]-𝑣)

# index = [8,16,32,64]
# XLSX.openxlsx("cav_patch.xlsx", mode="rw") do xf
#     Sheet = xf[3]
#     ind = findfirst(n->n==ndiv,index)+1
#     Sheet["B"*string(ind)] = log10(1/ndiv)
#     Sheet["C"*string(ind)] = a
# end

# fig = Figure()
# ind = 100
# ax = Axis(fig[1,1], 
#     aspect = DataAspect(), 
#     xticksvisible = false,
#     xticklabelsvisible=false, 
#     yticksvisible = false, 
#     yticklabelsvisible=false,
# )

# hidespines!(ax)
# hidedecorations!(ax)
# xs = LinRange(0, 1, ind)
# ys = LinRange(0, 1, ind)
# zs = zeros(ind,ind)
# 𝗠 = zeros(21)
# # s = 2.1/(ndivs)*ones(length(nodes_s))
# # push!(nodes_s,:s₁=>s,:s₂=>s,:s₃=>s)
# for (i,x) in enumerate(xs)
#     for (j,y) in enumerate(ys)
#         indices = sp(x,y,0.0)
#         ni = length(indices)
#         𝓒 = [nodes_s[i] for i in indices]
#         data = Dict([:x=>(2,[x]),:y=>(2,[y]),:z=>(2,[0.0]),:𝝭=>(4,zeros(ni)),:𝗠=>(0,𝗠)])
#         ξ = 𝑿ₛ((𝑔=1,𝐺=1,𝐶=1,𝑠=0), data)
#         𝓖 = [ξ]
#         a = type(𝓒,𝓖)
#         set𝝭!(a)
#         q = 0.0
#         N = ξ[:𝝭]
#         for (k,xₖ) in enumerate(𝓒)
#             q += N[k]*xₖ.q₁
#             # q += N[k]*xₖ.q₂
#         end
#         zs[i,j] = q
#     end
#  end
# surface!(xs,ys,zeros(ind,ind),color=zs,colorrange=(-0.000025,0.000025),colormap=:lightrainbow)
# contour!(xs,ys,zs,levels=-0.000025:0.00000715:0.000025,color=:azure)
# Colorbar(fig[1,2], limits=(-0.000025,0.000025), colormap=:lightrainbow)
# save("./png/SquarePlate_mix_tri3_q1_"*string(ndiv)*"_"*string(ndivs)*".png",fig, px_per_unit = 3.0)
# save("./png/SquarePlate_mix_tri3_q2_"*string(ndiv)*"_"*string(ndivs)*".png",fig, px_per_unit = 10.0)
# save("./png/SquarePlate_mix_tri6_q1_"*string(ndiv)*"_"*string(ndivs)*".png",fig, px_per_unit = 3.0)
# save("./png/SquarePlate_mix_colorbar.png",fig, px_per_unit = 10.0)
# save("./png/SquarePlate_mix_tri6_q2_"*string(ndiv)*"_"*string(ndivs)*".png",fig, px_per_unit = 10.0)
save("./png/cav_mix_quad4_q1_"*string(ndiv)*"_"*string(ndivs)*".png",fig, px_per_unit = 3.0)
# save("./png/SquarePlate_mix_quad4_q2_"*string(ndiv)*"_"*string(ndivs)*".png",fig, px_per_unit = 10.0)
# save("./png/SquarePlate_mix_quad8_q1_"*string(ndiv)*"_"*string(ndivs)*".png",fig, px_per_unit = 3.0)
# save("./png/SquarePlate_mix_quad8_q2_"*string(ndiv)*"_"*string(ndivs)*".png",fig, px_per_unit = 10.0)

# fig