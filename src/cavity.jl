using ApproxOperator
using Pardiso, TimerOutputs
using ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ
using WriteVTK, XLSX
using SparseArrays, LinearAlgebra

import ApproxOperator.Stokes:∫∫μ∇u∇vdxdy
import ApproxOperator.Elasticity:∫∫p∇udxdy, ∫vᵢtᵢds, ∫∫vᵢbᵢdxdy, ∫vᵢgᵢds, ∫qpdΩ
import Gmsh: gmsh

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
# d = zeros(2*nᵘ+nᵖ)

@timeit to "calculate ∫∫μ∇u∇vdxdy" begin
    @timeit to "get elements" elements_u = getElements(nodes, nodes, entities["Ω"])
    prescribe!(elements_u, :μ=>0.5, :b₁=>0.5, :b₂=>0.5, :E=>1.0, :ν=>1.0)
    @timeit to "calculate shape functions" set∇𝝭!(elements_u)
    aᵘ = ∫∫μ∇u∇vdxdy=>elements_u
    @timeit to "assemble" aᵘ(kᵘᵘ)
end

# @timeit to "floop k" @floop begin
#     kᵘᵘ = zeros(2*nᵘ,2*nᵘ)
#     for elm in elements
#         ∫∫μ∇u∇vdxdy(elm, kᵘᵘ)
#     end
# end

@timeit to "calculate ∫∫qpdΩ" begin
    @timeit to "get elements" elements = getElements(nodes_p, entities["Ω"])
    prescribe!(elements, :P=>80.0, :b₁=>0.5, :b₂=>0.5, :E=>1.0, :ν=>1.0)
    @timeit to "calculate shape functions" set∇𝝭!(elements)
    aᵖ = ∫∫qpdΩ=>elements
    @timeit to "assemble" aᵖ(kᵖᵖ)
end

# @timeit to "floop k" @floop begin
#     kᵖᵖ = zeros(nᵖ,nᵖ)
#     for elm in elements
#         ∫∫qpdΩ(elm, kᵖᵖ)
#     end
# end

@timeit to "calculate ∫∫p∇udxdy" begin
    @timeit to "get elements" elements = getElements(nodes, entities["Ω"])
    @timeit to "get elements" elements = getElements(nodes_p, entities["Ω"])
    prescribe!(elements, :)
    @timeit to "calculate shape functions" set∇𝝭!(elements)
    bᵖ = ∫∫p∇udxdy=>elements
    @timeit to "assemble" bᵖ(kᵖᵘ)
end

# @timeit to "floop k" @floop begin
#     kᵖᵘ = zeros(nᵖ,2*nᵘ)
#     for elm in elements
#         ∫∫∇v∇udxdy(elm, kᵖᵘ)
#     end
# end

@timeit to "calculate ∫vᵢtᵢds" begin
    @timeit to "get elements" elements = getElements(nodes, entities["Γ₁"], normal=true)
    prescribe!(elements, :t₁=>1.0, :t₂=>1.0, :g₁=>1.0, :g₂=>1.0, :α=>1e12*1.0)
    f = ∫vᵢtᵢds=>elements
    @timeit to "assemble" f(k,fᵘ)
end

# @timeit to "floop f" @floop begin
#     fᵘ = zeros(2*nᵘ)
#     for elm in elements
#         ∫vᵢtᵢds(elm, fᵘ)
#     end
# end

@timeit to "calculate ∫vᵢgᵢds" begin
    @timeit to "get elements" elements = getElements(nodes, entities["Γ₁"], normal=true)
    prescribe!(elements, :k=>1.0)
    aᵅ = ∫vᵢgᵢds=>elements
    @timeit to "assemble" aᵅ(kᵘᵘ, fᵘ)
end

k =[kᵘᵘ kᵖᵘ';kᵖᵘ kᵖᵖ]
f = [fᵘ;fᵖ]

@timeit to "solve" d = k\f

println(to)