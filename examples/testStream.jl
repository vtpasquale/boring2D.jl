
using Revise
using boring2D

cd("examples")

Cl = solveStream("n0012.toml")
# println("Cl = $(Cl)")
# println("Cl*180/pi = $(Cl*180/pi)")
# println("∂Cl∂α     = $(∂Cl∂α)")

inputFileName = "n0012.toml"
using SparseArrays
using TOML

# Read input file
inputData = TOML.parsefile(inputFileName)

# Echo input file
TOML.print(inputData)

# Process mesh
mesh = boring2D.readMesh(inputData["mesh"]["fileName"])
nDof = size(mesh.nodes,1)
nEle = size(mesh.triangles,1)
triangleElements = boring2D.TriangleElements(mesh)
airfoilBoundary, farfieldBoundary, f, s = boring2D.processMeshBoundaries(mesh,inputData)
gxγ = boring2D.computeGxγ(mesh,triangleElements,airfoilBoundary,f)

# Assemble coefficent matrix
K = boring2D.assembleConvectionStiffness(mesh,triangleElements,1.0)
g = boring2D.assembleKuttaConditionConstraint(mesh,inputData,triangleElements)
esa = boring2D.assembleEsa(nDof,farfieldBoundary)
A = boring2D.assembleA(K,g,esa,f,s)

# Assemble right-hand side
ψsk = boring2D.assembleψsk(nDof,farfieldBoundary,inputData)
b = boring2D.assembleB(K,ψsk,f,s)

# Solve linear system
x = A\Array(b)

# Recover responses
ψ = boring2D.recoverψ(x,ψsk,esa,f,s)
Cl = (2/inputData["freestream"]["velocity"])* gxγ*x

# Adjoint solution
∂Cl∂x = (2/inputData["freestream"]["velocity"])* gxγ
λ = transpose(A)\transpose(∂Cl∂x)

fx = boring2D.getFxfromF(f)
lamdaOut = zeros(nDof+1,1)
lamdaOut[fx] = λ

# Sensitivity wrt α
∂ψsk_∂α = boring2D.assemble_∂ψsk_∂α(nDof,farfieldBoundary,inputData)
∂b_∂α = boring2D.assemble_∂b_∂α(K,∂ψsk_∂α,f,s)
∂Cl_∂α = transpose(λ)*∂b_∂α # per radian
# ∂Cl_∂α_deg = ∂Cl_∂α*pi/180

pointOutput = Dict("psi"=>ψ,"lamda"=>lamdaOut[1:end-1])
cellOutput = Dict()

# Write output data to file
boring2D.writeSolution("streamSolution.vtu",mesh,pointOutput,cellOutput)