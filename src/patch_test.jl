
using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.Heat: ∫vtdΓ, ∫vgdΓ, ∫vbdΩ, L₂, ∫∫∇v∇udxdy

using TimerOutputs 
using SparseArrays
import Gmsh: gmsh

𝑢(x,y,z) = x + y
const to = TimerOutput()

gmsh.initialize()
ndiv = 4
@timeit to "open msh file" gmsh.open("msh/cav_quad_"*string(ndiv)*".msh")
# @timeit to "open msh file" gmsh.open("patchtest_165.msh")

@timeit to "get entities" entities = getPhysicalGroups()
println("entities: ", entities)
@timeit to "get nodes" nodes = get𝑿ᵢ()
println("nodes: ", nodes)

nₚ = length(nodes)
k = zeros(nₚ,nₚ)
f = zeros(nₚ)

@timeit to "calculate ∫∫∇v∇udxdy" begin
    @timeit to "get elements" elements = getElements(nodes, entities["Ω"])
    prescribe!(elements, :k=>1.0)
    @timeit to "calculate shape functions" set∇𝝭!(elements)
    𝑎 = ∫∫∇v∇udxdy=>elements
    @timeit to "assemble" 𝑎(k)
end

# @timeit to "floop k" @floop begin
#     k = zeros(nₚ,nₚ)
#     for elm in elements
#         ∫∫∇v∇udxdy(elm, k)
#     end
# end

@timeit to "calculate ∫vgdΓ" begin
    @timeit to "get elements" elements = getElements(nodes, entities["Γ₁"])
    prescribe!(elements, :g=>𝑢, :α=>1e7)
    @timeit to "calculate shape functions" set𝝭!(elements)
    𝑓 = ∫vgdΓ=>elements
    @timeit to "assemble" 𝑓(k,f)
end

@timeit to "solve" d = k\f

push!(nodes, :d=>d)

elements = getElements(nodes, entities["Ω"], 10)
prescribe!(elements, :u=>𝑢)
set∇𝝭!(elements)
L₂error = L₂(elements)
gmsh.finalize()

println(to)
# println("L₂ error: ", L₂error)
