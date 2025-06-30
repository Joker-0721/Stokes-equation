using ApproxOperator
using Pardiso
using WriteVTK ,XLSX
using SparseArrays, LinearAlgebra
import ApproxOperator.Stokes:∫∫μ∇u∇vdxdy
import ApproxOperator.Elasticity:∫∫p∇udxdy, ∫vᵢtᵢds, ∫∫vᵢbᵢdxdy, ∫vᵢgᵢds

include("import_cavity.jl")

ps = MKLPardisoSolver()

ndiv = 4

elements, nodes, nodes_s= import_cavity_RI("Stokes-equation/msh/cav_quad_"*string(ndiv)*".msh", "Stokes-equation/msh/cav_quad_"*string(ndiv)*".msh");

nᵘ = length(nodes)
nᵖ = length(nodes_s)
nₑ = length(elements["Ω"])

E = 1.0
ν = 1.0
μ = 0.5
b₁ = 0.5
b₂ = 0.5
t₁ = 0.5
t₂ = 0.5

n₁₁(n₁,n₂) = 1.0
n₂₂(n₁,n₂) = 1.0
prescribe!(elements["Ω"],:μ=>(x,y,z)->μ)
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
prescribe!(elements["Γ₁"],:g₁=>(x,y,z)->0.0)
prescribe!(elements["Γ₂"],:g₁=>(x,y,z)->0.0)
prescribe!(elements["Γ₃"],:g₁=>(x,y,z)->0.0)
prescribe!(elements["Γ₄"],:g₁=>(x,y,z)->0.0)
prescribe!(elements["Γ₁"],:g₂=>(x,y,z)->0.0)
prescribe!(elements["Γ₂"],:g₂=>(x,y,z)->0.0)
prescribe!(elements["Γ₃"],:g₂=>(x,y,z)->0.0)
prescribe!(elements["Γ₄"],:g₂=>(x,y,z)->0.0)
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

aᵘ = ∫∫μ∇u∇vdxdy => elements["Ω"]
bᵖ = ∫∫p∇udxdy => (elements["Ω"],elements["Ωˢ"])
# kᵖ = ∫∫vᵢbᵢdxdy => elements["Ω"]
f = [
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
kᵘᵖ = zeros(nᵖ,2*nᵘ)
kᵖᵖ = zeros(nᵖ,nᵖ)
fᵖ = zeros(nᵖ)
fᵘ = zeros(2*nᵘ)
d = zeros(2*nᵘ+nᵖ)

aᵘ(kᵘᵘ)
bᵖ(kᵘᵖ)
# kᵖ(fᵖ)
f(fᵘ)
aᵅ(kᵘᵘ,fᵘ)

k =[kᵘᵘ kᵘᵖ';kᵘᵖ kᵖᵖ]
f = [fᵘ;fᵖ]

# d = k\f
set_matrixtype!(ps, -2)
k = get_matrix(ps,k,:N)
pardiso(ps,d,k,f)
# d₁ = d[1:3:3*nᵘ]
# d₂ = d[2:3:3*nᵘ]

# u = d[1:2nᵘ]
# p = d[2nᵘ+1:end]

# # 创建 VTK 网格
# points = zeros(3, nᵘ)
# for (i, node) in enumerate(nodes)
#     points[1, i] = node.x
#     points[2, i] = node.y
#     points[3, i] = 0.0  # 2D问题，z坐标为0
# end

# # 创建单元连接关系
# cells = []
# for elem in elements["Ω"]
#     # 假设是四边形单元，每个单元有4个节点
#     push!(cells, MeshCell(VTKCellTypes.VTK_QUAD, elem.𝓒))
# end

# # 创建 VTK 文件
# vtk_grid("cavity_flow", points, cells) do vtk
#     # 添加速度场
#     u₁ = u[1:2:end]  # u1分量
#     u₂ = u[2:2:end]  # u2分量
#     vtk["Velocity"] = (u₁, u₂, zeros(nᵘ))  # 2D问题，第三个分量为0
    
#     # 添加压力场（需要插值到节点上）
#     # 这里假设压力节点与速度节点相同，实际情况可能需要调整
#     vtk["Pressure"] = p
# end

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
# s₂ = d[3*nᵇ+2:2:3*nᵇ+2*nˢ]

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
# save("./png/cav_mix_quad4_q1_"*string(ndiv)*"_"*string(ndivs)*".png",fig, px_per_unit = 3.0)
# save("./png/SquarePlate_mix_quad4_q2_"*string(ndiv)*"_"*string(ndivs)*".png",fig, px_per_unit = 10.0)
# save("./png/SquarePlate_mix_quad8_q1_"*string(ndiv)*"_"*string(ndivs)*".png",fig, px_per_unit = 3.0)
# save("./png/SquarePlate_mix_quad8_q2_"*string(ndiv)*"_"*string(ndivs)*".png",fig, px_per_unit = 10.0)

# fig