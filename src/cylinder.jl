using ApproxOperator
using TimerOutputs
using ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
using WriteVTK, XLSX
using SparseArrays, LinearAlgebra

import ApproxOperator.Stokes:∫∫μ∇u∇vdxdy
import ApproxOperator.Elasticity:∫∫p∇udxdy, ∫vᵢtᵢds, ∫∫vᵢbᵢdxdy, ∫vᵢgᵢds, ∫qpdΩ, L₂
import Gmsh: gmsh

# 速度场精确解
function velocity(x, y)
    r = sqrt(x^2 + y^2)
    θ = atan(y, x)
    u_r = U * (1 - a^2 / r^2) * cos(θ)
    u_θ = -U * (1 + a^2 / r^2) * sin(θ)
    u_x = u_r * cos(θ) - u_θ * sin(θ)
    u_y = u_r * sin(θ) + u_θ * cos(θ)
    return u_x, u_y
end

# 压力系数精确解
function pressure(x, y)
    r = sqrt(x^2 + y^2)
    θ = atan(y, x)
    return 1 - 4 * sin(θ)^2
end
const to = TimerOutput()

gmsh.initialize()
# type = "quad"
type = "tri"
ndiv_u = 10
ndiv_p = 2
type_p = :(ReproducingKernel{:Linear2D,:□,:CubicSpline})
integrationOrder = 2
@timeit to "open msh file" gmsh.open("msh/cylinder_"*type*"_"*string(ndiv_p)*".msh")
@timeit to "get nodes_p" nodes_p = get𝑿ᵢ()  
xᵖ = nodes_p.x
yᵖ = nodes_p.y
zᵖ = nodes_p.z
nᵖ = length(nodes_p)
sp = RegularGrid(xᵖ,yᵖ,zᵖ,n = 3,γ = 5)
s = 1/ndiv_p
s₁ = 1.5*s*ones(nᵖ)
s₂ = 1.5*s*ones(nᵖ)
s₃ = 1.5*s*ones(nᵖ)
push!(nodes_p,:s₁=>s₁,:s₂=>s₂,:s₃=>s₃)


@timeit to "open msh file" gmsh.open("msh/cylinder_"*type*"_"*string(ndiv_u)*".msh")
@timeit to "get entities" entities = getPhysicalGroups()
@timeit to "get nodes" nodes = get𝑿ᵢ()
nᵘ = length(nodes)

kᵘᵘ = zeros(2*nᵘ,2*nᵘ)
kᵖᵘ = zeros(nᵖ,2*nᵘ)
kᵖᵖ = zeros(nᵖ,nᵖ)
fᵖ = zeros(nᵖ)
fᵘ = zeros(2*nᵘ)

E = 1.0
ν = 0.3
μ = 0.01
p = 0.0

n₁₁(x,y,z,n₁,n₂) = n₁*n₁
n₁₂(x,y,z,n₁,n₂) = n₁*n₂
n₂₂(x,y,z,n₁,n₂) = n₂*n₂

@timeit to "assembly" begin
    @timeit to "get elements" elements_u = getElements(nodes, entities["Ω"], integrationOrder)
    @timeit to "get elements" elements_p = getElements(nodes_p, entities["Ω"], eval(type_p), integrationOrder, sp)
    prescribe!(elements_u, :μ=>μ)
    prescribe!(elements_p, :E=>E, :ν=>ν)
    @timeit to "calculate shape functions" set∇𝝭!(elements_u)
    @timeit to "calculate shape functions" set𝝭!(elements_p)
    𝑎 = ∫∫μ∇u∇vdxdy=>elements_u
    𝑏 = ∫∫p∇udxdy=>(elements_p, elements_u)
    # 𝑐 = ∫qpdΩ=>elements_p
    𝑓 = ∫∫vᵢbᵢdxdy => elements_u
    @timeit to "assemble" 𝑎(kᵘᵘ)
    # @timeit to "assemble" 𝑐(kᵖᵖ)
    @timeit to "assemble" 𝑏(kᵖᵘ)
end

@timeit to "calculate ∫vᵢgᵢds" begin
    @timeit to "get elements" elements_1 = getElements(nodes, entities["Γ₁"], integrationOrder, normal=true)
    @timeit to "get elements" elements_2 = getElements(nodes, entities["Γ₂"], integrationOrder)
    # @timeit to "get elements" elements_3 = getElements(nodes, entities["Γ₃"], integrationOrder, normal=true)
    prescribe!(elements_1, :g₁=>0.0, :g₂=>0.0, :α=>1e14, :n₁₁=>n₁₁, :n₂₂=>n₂₂, :n₁₂=>n₁₂)
    prescribe!(elements_2, :g₁=>1.0, :g₂=>0.0, :α=>1e14, :n₁₁=>1.0, :n₂₂=>1.0, :n₁₂=>0.0)
    # prescribe!(elements_3, :g₁=>0.0, :g₂=>0.0, :α=>1e14, :n₁₁=>n₁₁, :n₂₂=>n₂₂, :n₁₂=>n₁₂, p=>p)
    @timeit to "calculate shape functions" set𝝭!(elements_1)
    @timeit to "calculate shape functions" set𝝭!(elements_2)
    # @timeit to "calculate shape functions" set𝝭!(elements_3)
    𝑎 = ∫vᵢgᵢds => elements_1∪elements_2
    @timeit to "assemble" 𝑎(kᵘᵘ, fᵘ)
end

k =[kᵘᵘ kᵖᵘ';kᵖᵘ kᵖᵖ]
f = [fᵘ;fᵖ]

@timeit to "solve" d = k\f

push!(nodes, :d₁=>d[1:2:2*nᵘ], :d₂=>d[2:2:2*nᵘ], :d₃ => zeros(nᵘ))
push!(nodes_p, :p=>d[2*nᵘ+1:end])

elements = getElements(nodes, entities["Ω"], 10)
prescribe!(elements, :u₁ => (x,y,z) -> velocity(x,y)[1], :u₂ => (x,y,z) -> velocity(x,y)[2], :u₃ => 0.0)
set∇𝝭!(elements)
L₂error = L₂(elements)
gmsh.finalize()

println(to)
println("L₂ error: ", L₂error)

# pressure = zeros(nᵘ)
# u₁ = zeros(nᵘ)
# u₂ = zeros(nᵘ)
# u₃ = zeros(nᵘ)
# 𝗠 = zeros(10)
# for (i,node) in enumerate(nodes)
#     x = node.x
#     y = node.y
#     z = node.z
#     indices = sp(x,y,z)
#     ni = length(indices)
#     𝓒 = [nodes_p[i] for i in indices]
#     data = Dict([:x=>(2,[x]),:y=>(2,[y]),:z=>(2,[z]),:𝝭=>(4,zeros(ni)),:𝗠=>(0,𝗠)])
#     ξ = 𝑿ₛ((𝑔=1,𝐺=1,𝐶=1,𝑠=0), data)
#     𝓖 = [ξ]
#     a = eval(type_p)(𝓒,𝓖)
#     set𝝭!(a)
#     p = 0.0
#     N = ξ[:𝝭]
#     for (k,xₖ) in enumerate(𝓒)
#         p += N[k]*xₖ.p
#     end
#     pressure[i] = p
#     u₁[i] = node.d₁
#     u₂[i] = node.d₂
# end
# α = 1.0
# points = zeros(3, nᵘ)
# for node in nodes
#     I = node.𝐼
#     points[1, I] = node.x
#     points[2, I] = node.y
#     points[3, I] = node.z
# end

# # 在构建 cells 数组之前，先检查单元节点类型
# println("正在检查单元节点类型...")
# node_counts = Dict{Int, Int}() # 创建一个字典来统计不同节点数的单元数量

# # 首先，遍历所有单元，统计每个单元的节点数
# for elm in elements
#     n_nodes = length(elm.𝓒) # 获取当前单元的节点数量
#     # 统计不同节点数出现的次数
#     node_counts[n_nodes] = get(node_counts, n_nodes, 0) + 1
# end

# # 打印统计信息
# for (n_nodes, count) in node_counts
#     println("具有 $n_nodes 个节点的单元数量: $count")
# end

# # 定义一个函数，根据节点数量确定 VTK 单元类型
# function get_vtk_cell_type(n_nodes)
#     if n_nodes == 3
#         return VTKCellTypes.VTK_TRIANGLE
#     elseif n_nodes == 4
#         # 可能是四边形(QUAD)或四面体(TETRA)，需要根据你的网格类型判断
#         # 这里是二维网格的例子，假设是四边形
#         return VTKCellTypes.VTK_QUAD
#     elseif n_nodes == 8
#         return VTKCellTypes.VTK_HEXAHEDRON
#     # 可以根据需要添加更多的映射关系
#     else
#         error("不支持的单元节点数量: $n_nodes. 无法映射到已知的 VTK 单元类型。")
#     end
# end

# # 现在构建 cells 数组，并对每个单元进行检查
# cells = MeshCell{VTKCellType}[]
# for elm in elements
#     n_nodes = length(elm.𝓒)
#     vtk_cell_type = get_vtk_cell_type(n_nodes)
#     # 也可以在这里打印或记录个别单元的信息（如果怀疑某个单元有问题）
#     # println("单元节点索引: ", [xᵢ.𝐼 for xᵢ in elm.𝓒], " 类型: ", vtk_cell_type)
#     push!(cells, MeshCell(vtk_cell_type, [xᵢ.𝐼 for xᵢ in elm.𝓒]))
# end

# println("单元类型检查完毕，开始写入 VTK...")

# cells = [MeshCell(VTKCellTypes.VTK_QUAD,[xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements]
# cells = [MeshCell(VTKCellTypes.VTK_TETRA,[xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements]
# cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE,[xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements]
# cells = [MeshCell(VTKCellTypes.VTK_HEXAHEDRON,[xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements["Ωᵘ"]]
# vtk_grid("./vtk/cylinder_"*type*"_"*string(ndiv_u)*"_"*string(nᵖ),points,cells) do vtk
#     vtk["u"] = (u₁,u₂,u₃)
#     vtk["p"] = pressure
# end

# println(nodes[5])