
using Revise
using boring2D

cd("examples")

# solveStream("m03_input.toml")
# solveStream("n0012_input.toml")


inputFileName = "n0012_input.toml"
using SparseArrays
using TOML

# Read input file
inputData = TOML.parsefile(inputFileName)

# Echo input file
TOML.print(inputData)

# Output dictionaries
textOutput =  Dict() # Dict("KuttaJoukowski" => KuttaJoukowskiTextOut,"IntegratedPressure" => IntegratedPressureTextOut, "Adjoint" => AdjointTextOut)
# pointOutput = Dict() # Dict("psi"=>ψ,"lamda"=>lamdaOut[1:end-1])
# cellOutput = Dict() # Dict("vX"=>[velX],"vY"=>[velY],"dVelX"=>[dVelX],"dVelY"=>[dVelY],"velMag"=>[velMag],"Cp"=>[Cp])

# Process mesh
mesh = boring2D.readMesh(inputData["mesh"]["fileName"])
nDof = size(mesh.nodes,1)
nEle = size(mesh.triangles,1)
triangleElements = boring2D.TriangleElements(mesh)
airfoilBoundary, farfieldBoundary, f, s = boring2D.processMeshBoundaries(mesh,inputData)
gxγ = boring2D.computeGxγ(mesh,triangleElements,airfoilBoundary,f)

# Assemble coefficient matrix
K = boring2D.assembleConvectionStiffness(mesh,triangleElements,1.0)
g = boring2D.assembleKuttaConditionConstraint(mesh,inputData,triangleElements)
esa = boring2D.assembleEsa(nDof,farfieldBoundary)
A = boring2D.assembleA(Array(K),g,esa,f,s)

# Assemble right-hand side
ψsk = boring2D.assembleψsk(nDof,farfieldBoundary,inputData)
b = boring2D.assembleB(Array(K),ψsk,f,s)

# Stream function solution
x = convert(Matrix{Float64},A)\convert(Vector{Float64},b)

# Stream function values
ψ = boring2D.recoverψ(x,ψsk,esa,f,s)

# Lift coefficient from Kutta-Joukowski
Cl = (2/inputData["freestream"]["velocity"])* gxγ*x
textOutput["KuttaJoukowski"] = Dict("Cl" => Cl)

# Velocity and Cp
velX, velY = boring2D.recoverVelocities(ψ,mesh,triangleElements)
velMag = sqrt.(velX.^2 + velY.^2)
Cp = 1.0 .- (velMag ./ inputData["freestream"]["velocity"]).^2

# Integrate pressure coefficients
textOutput["IntegratedPressure"] = boring2D.computeCfFromCp(Cp,mesh,airfoilBoundary,inputData)

# Volume output dictionaries
pointOutput = Dict("psi"=>ψ)
cellOutput = Dict("vX"=>[velX],"vY"=>[velY],"velMag"=>[velMag],"Cp"=>[Cp])

# Adjoint solution
∂Cl∂x = (2/inputData["freestream"]["velocity"])* gxγ
λ = transpose(convert(Matrix{Float64},A))\transpose(∂Cl∂x)

# Adjoint contour output
fx = boring2D.getFxfromF(f)
lamdaOut = zeros(nDof+1)
lamdaOut[fx] = λ
cellOutput["λ"] = [lamdaOut]

# Sensitivity wrt α
∂ψsk_∂α = boring2D.assemble_∂ψsk_∂α(nDof,farfieldBoundary,inputData)
∂r_∂α = boring2D.assemble_∂r_∂α(Array(K),∂ψsk_∂α,f,s)
∂Cl_∂α = transpose(λ)*∂r_∂α # per radian | ∂Cl_∂α_deg = ∂Cl_∂α*pi/180
textOutput["Adjoint"] = Dict("∂Cl_∂α_rad" => ∂Cl_∂α,"∂Cl_∂α_deg" => ∂Cl_∂α*pi/180)

# Shape sensitivities
∂g_∂b, ∂r_∂b = boring2D.computeShapePartials(mesh,x,inputData,size(airfoilBoundary.nodeIDs,1))
∂Cl_∂b = ∂g_∂b - transpose( transpose(λ)*∂r_∂b )

# Length-normalized shape sensitivities to surface output file
boring2D.surfaceOutput2Vtk(inputData["output"]["surfaceFile"],airfoilBoundary,∂Cl_∂b[:]./airfoilBoundary.nodeLength)

# Recover element output
dVelX, dVelY = boring2D.recoverVelocities(lamdaOut,mesh,triangleElements)
cellOutput["dVelX"] = [dVelX]
cellOutput["dVelY"] = [dVelY]

# Text output dictionary and output file
TOML.print(textOutput)
open(inputData["output"]["tomlFile"], "w") do io
    TOML.print(io, textOutput)
end

# Volume output dictionaries and output file
boring2D.writeSolution(inputData["output"]["volumeFile"],mesh,pointOutput,cellOutput)
