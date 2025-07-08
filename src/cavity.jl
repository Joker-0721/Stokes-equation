using ApproxOperator
using Pardiso, TimerOutputs
using ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
using WriteVTK, XLSX
using SparseArrays, LinearAlgebra

import ApproxOperator.Stokes:∫∫μ∇u∇vdxdy
import ApproxOperator.Elasticity:∫∫p∇udxdy, ∫vᵢtᵢds, ∫∫vᵢbᵢdxdy, ∫vᵢgᵢds, ∫qpdΩ
import ApproxOperator.Heat: L₂
import Gmsh: gmsh

𝑢(x,y,z) = x + y
const to = TimerOutput()

gmsh.initialize()
ndiv = 4
@timeit to "open msh file" gmsh.open("msh/cav_quad_"*string(ndiv)*".msh")

@timeit to "get entities" entities = getPhysicalGroups()
@timeit to "get nodes" nodes = get𝑿ᵢ()
@timeit to "get nodes_p" nodes_p = get𝑿ᵢ()  

nᵘ = length(nodes)
nᵖ = length(nodes_p)
kᵘᵘ = zeros(2*nᵘ,2*nᵘ)
kᵖᵘ = zeros(nᵖ,2*nᵘ)
kᵖᵖ = zeros(nᵖ,nᵖ)
fᵖ = zeros(nᵖ)
fᵘ = zeros(2*nᵘ)

E = 1.0
ν = 0.01
μ = 0.5
b₁ = 0.5
b₂ = 0.5
t₁ = 0.5
t₂ = 0.5
g₁ = 0.0
g₂ = 0.0
P = 80.0
n₁₁ = 1.0
n₂₂ = 1.0
n₁₂ = 0.0

# function original_test_system()
#     k = spzeros(11, 11)
#     f = zeros(11)
#     for i in 1:10
#         k[i, i] = 1.0
#         k[i, i+1] = -1.0
#         k[i+1, i] = -1.0
#         k[i+1, i+1] = 1.0
#     end
#     k[1,1] = 1e10
#     k[11,11] = 1e10
#     f[11] = 1.0
#     d = k\f
#     return d
# end

# function fluid_system()

@timeit to "calculate ∫∫μ∇u∇vdxdy" begin
    @timeit to "get elements" elements_u = getElements(nodes, entities["Ω"])
    prescribe!(elements_u, :μ=>0.5, :b₁=>0.5, :b₂=>0.5, :E=>1.0, :ν=>1.0)
    @timeit to "calculate shape functions" set∇𝝭!(elements_u)
    𝑎ᵘ = ∫∫μ∇u∇vdxdy => elements_u
    @timeit to "assemble" 𝑎ᵘ(kᵘᵘ)
end

@timeit to "calculate ∫qpdΩ" begin
    @timeit to "get elements" elements_p = getElements(nodes_p, entities["Ω"])
    prescribe!(elements_p, :P=>80.0, :b₁=>0.5, :b₂=>0.5, :E=>1.0, :ν=>1.0)
    @timeit to "calculate shape functions" set∇𝝭!(elements_p)
    𝑎ᵖ = ∫qpdΩ => elements_p
    @timeit to "assemble" 𝑎ᵖ(kᵖᵖ)
end

@timeit to "calculate ∫∫p∇udxdy" begin
    @timeit to "get elements" elements_p = getElements(nodes_p, entities["Ω"])
    @timeit to "get elements" elements_u = getElements(nodes, entities["Ω"])
    prescribe!(elements_p, :P=>80.0, :b₁=>0.5, :b₂=>0.5, :E=>1.0, :ν=>1.0)
    prescribe!(elements_u, :b₁=>0.5, :b₂=>0.5, :E=>1.0, :ν=>1.0)
    @timeit to "calculate shape functions" set∇𝝭!(elements_p)
    @timeit to "calculate shape functions" set∇𝝭!(elements_u)
    𝑏ᵖ = ∫∫p∇udxdy => (elements_p, elements_u)
    @timeit to "assemble" 𝑏ᵖ(kᵖᵘ)
end

@timeit to "calculate ∫∫vᵢbᵢdxdy" begin
    @timeit to "get elements" elements_u = getElements(nodes, entities["Ω"])
    prescribe!(elements_u, :b₁=>0.5, :b₂=>0.5, :E=>1.0, :ν=>1.0)
    @timeit to "calculate shape functions" set𝝭!(elements_u)
    𝑓 = ∫∫vᵢbᵢdxdy => elements_u
    @timeit to "assemble" 𝑓(fᵘ)
end

@timeit to "calculate ∫vᵢtᵢds" begin
    @timeit to "get elements" elements_1 = getElements(nodes, entities["Γ₁"])
    @timeit to "get elements" elements_2 = getElements(nodes, entities["Γ₂"])
    @timeit to "get elements" elements_3 = getElements(nodes, entities["Γ₃"])
    @timeit to "get elements" elements_4 = getElements(nodes, entities["Γ₄"])
    prescribe!(elements_1, :t₁=>t₁, :t₂=>t₂)
    prescribe!(elements_2, :t₁=>t₁, :t₂=>t₂)
    prescribe!(elements_3, :t₁=>t₁, :t₂=>t₂)
    prescribe!(elements_4, :t₁=>t₁, :t₂=>t₂)
    @timeit to "calculate shape functions" set𝝭!(elements_1)
    @timeit to "calculate shape functions" set𝝭!(elements_2)
    @timeit to "calculate shape functions" set𝝭!(elements_3)
    @timeit to "calculate shape functions" set𝝭!(elements_4)
    𝑓 = ∫vᵢtᵢds => elements_1
    𝑓 = ∫vᵢtᵢds => elements_2
    𝑓 = ∫vᵢtᵢds => elements_3
    𝑓 = ∫vᵢtᵢds => elements_4
    @timeit to "assemble" 𝑓(fᵘ)
end

@timeit to "calculate ∫vᵢgᵢds" begin
    @timeit to "get elements" elements_1 = getElements(nodes, entities["Γ₁"])
    @timeit to "get elements" elements_2 = getElements(nodes, entities["Γ₂"])
    @timeit to "get elements" elements_3 = getElements(nodes, entities["Γ₃"])
    @timeit to "get elements" elements_4 = getElements(nodes, entities["Γ₄"])
    prescribe!(elements_1, :g₁=>g₁, :g₂=>g₂, :α=>1e7*E, :n₁₁=>n₁₁, :n₂₂=>n₂₂, :n₁₂=>n₁₂)
    prescribe!(elements_2, :g₁=>g₁, :g₂=>g₂, :α=>1e7*E, :n₁₁=>n₁₁, :n₂₂=>n₂₂, :n₁₂=>n₁₂)
    prescribe!(elements_3, :g₁=>1.0, :g₂=>1.0, :α=>1e7*E, :n₁₁=>n₁₁, :n₂₂=>0.0, :n₁₂=>0.0)
    prescribe!(elements_4, :g₁=>g₁, :g₂=>g₂, :α=>1e7*E, :n₁₁=>n₁₁, :n₂₂=>n₂₂, :n₁₂=>n₁₂)
    @timeit to "calculate shape functions" set𝝭!(elements_1)
    @timeit to "calculate shape functions" set𝝭!(elements_2)
    @timeit to "calculate shape functions" set𝝭!(elements_3)
    @timeit to "calculate shape functions" set𝝭!(elements_4)
    𝑎ᵅ = ∫vᵢgᵢds => elements_1
    𝑎ᵅ = ∫vᵢgᵢds => elements_2
    𝑎ᵅ = ∫vᵢgᵢds => elements_3
    𝑎ᵅ = ∫vᵢgᵢds => elements_4
    @timeit to "assemble" 𝑎ᵅ(kᵘᵘ, fᵘ)
end

k =[kᵘᵘ kᵖᵘ';kᵖᵘ kᵖᵖ]
f = [fᵘ;fᵖ]

@timeit to "solve" d = k\f

# return d

# end

push!(nodes, :d=>d)

elements = getElements(nodes, entities["Ω"], 10)
prescribe!(elements, :u=>𝑢)
set∇𝝭!(elements)
L₂error = L₂(elements)
gmsh.finalize()

println(to)
println("L₂ error: ", L₂error)

# # 主測試函數
# function test_fluid_code()
#     # 獲取兩個系統的結果
#     d_original = original_test_system()
#     d_fluid = fluid_system()
    
#     # 計算相對誤差 (只比較前11個自由度)
#     n = min(11, length(d_fluid))
#     error = norm(d_fluid[1:n] - d_original[1:n]) / norm(d_original[1:n])
    
#     # 判斷誤差是否在合理範圍內 (這裡設定為1e-5)
#     tolerance = 1e-5
#     if error < tolerance
#         println("測試通過: 誤差值 = ", error, " < ", tolerance)
#         return true
#     else
#         println("測試失敗: 誤差值 = ", error, " >= ", tolerance)
#         return false
#     end
# end

# # 執行測試
# test_fluid_code()