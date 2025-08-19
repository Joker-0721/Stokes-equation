using ApproxOperator
using TimerOutputs
using ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
using WriteVTK, XLSX
using SparseArrays, LinearAlgebra

import ApproxOperator.Stokes:∫∫μ∇u∇vdxdy
import ApproxOperator.Elasticity:∫∫p∇udxdy, ∫vᵢtᵢds, ∫∫vᵢbᵢdxdy, ∫vᵢgᵢds, ∫qpdΩ, L₂
import Gmsh: gmsh

const to = TimerOutput()

gmsh.initialize()
ndiv = 4
ndiv_p = 4
type_p = :(ReproducingKernel{:Linear2D,:□,:CubicSpline})
integrationOrder = 2
@timeit to "open msh file" gmsh.open("msh/cav_quad_"*string(ndiv)*".msh")
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

@timeit to "open msh file" gmsh.open("msh/cav_quad_"*string(ndiv)*".msh")
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
μ = 0.5*E/(1+ν)

@timeit to "assembly" begin
    @timeit to "get elements" elements_u = getElements(nodes, entities["Ω"], integrationOrder)
    @timeit to "get elements" elements_p = getElements(nodes_p, entities["Ω"], eval(type_p), integrationOrder, sp)
    prescribe!(elements_u, :μ=>μ)
    prescribe!(elements_p, :E=>E, :ν=>ν)
    @timeit to "calculate shape functions" set∇𝝭!(elements_u)
    @timeit to "calculate shape functions" set𝝭!(elements_p)
    𝑎 = ∫∫μ∇u∇vdxdy => elements_u
    𝑏 = ∫∫p∇udxdy=>(elements_p, elements_u)
    𝑐 = ∫qpdΩ=>elements_p
    𝑓 = ∫∫vᵢbᵢdxdy => elements_u
    @timeit to "assemble" 𝑎(kᵘᵘ)
    @timeit to "assemble" 𝑐(kᵖᵖ)
    @timeit to "assemble" 𝑏(kᵖᵘ)
end

@timeit to "calculate ∫vᵢgᵢds" begin
    @timeit to "get elements" elements_1 = getElements(nodes, entities["Γ₁"], integrationOrder)
    @timeit to "get elements" elements_2 = getElements(nodes, entities["Γ₂"], integrationOrder)
    @timeit to "get elements" elements_3 = getElements(nodes, entities["Γ₃"], integrationOrder)
    @timeit to "get elements" elements_4 = getElements(nodes, entities["Γ₄"], integrationOrder)
    prescribe!(elements_1, :g₁=>0.0, :g₂=>0.0, :α=>1e14*E, :n₁₁=>1.0, :n₂₂=>1.0, :n₁₂=>0.0)
    prescribe!(elements_2, :g₁=>0.0, :g₂=>0.0, :α=>1e14*E, :n₁₁=>1.0, :n₂₂=>1.0, :n₁₂=>0.0)
    prescribe!(elements_3, :g₁=>1.0, :g₂=>0.0, :α=>1e14*E, :n₁₁=>1.0, :n₂₂=>1.0, :n₁₂=>0.0)
    prescribe!(elements_4, :g₁=>0.0, :g₂=>0.0, :α=>1e14*E, :n₁₁=>1.0, :n₂₂=>1.0, :n₁₂=>0.0)
    @timeit to "calculate shape functions" set𝝭!(elements_1)
    @timeit to "calculate shape functions" set𝝭!(elements_2)
    @timeit to "calculate shape functions" set𝝭!(elements_3)
    @timeit to "calculate shape functions" set𝝭!(elements_4)
    𝑎 = ∫vᵢgᵢds => elements_1∪elements_2∪elements_3∪elements_4
    @timeit to "assemble" 𝑎(kᵘᵘ, fᵘ)
end

k =[kᵘᵘ kᵖᵘ';kᵖᵘ kᵖᵖ]
f = [fᵘ;fᵖ]

@timeit to "solve" d = k\f

push!(nodes, :d₁=>d[1:2:2*nᵘ], :d₂=>d[2:2:2*nᵘ], :d₃=>zeros(nᵘ))

elements = getElements(nodes, entities["Ω"], 10)
prescribe!(elements, :u₁=>𝑢, :u₂=>𝑢, :u₃=>0.0)
set∇𝝭!(elements)
L₂error = L₂(elements)
gmsh.finalize()

println(to)
println("L₂ error: ", L₂error)