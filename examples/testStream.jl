
# using Revise
using boring2D

cd("examples")

meshFileName = "n0012.plt";
# meshFileName = "5000NUcav.plt";
# meshFileName = "square.su2";

boundaryConditionType = [2;1] # [farfield;airfoil]
# 0 - Natural BC with zero flux
# 1 - Essential BC with constant value
# 2 - Essential freestream BC

dummyId = 999

freestreamV = 100.
freestreamAlpha = 30. *pi/180.
essentialPhi = [0.;0.]

pointOutput,cellOutput = solveStream(meshFileName,boundaryConditionType,freestreamV,freestreamAlpha,essentialPhi,dummyId)

mesh = readMesh(meshFileName)
boring2D.writeSolution("psi1.vtu",mesh,pointOutput,cellOutput)  


boundaryConditionType = [1;1] # [farfield;airfoil]
# 0 - Natural BC with zero flux
# 1 - Essential BC with constant value
# 2 - Essential freestream BC
essentialPhi = [0.;1.]
pointOutput2,cellOutput2 = solveStream(meshFileName,boundaryConditionType,freestreamV,freestreamAlpha,essentialPhi,dummyId)
boring2D.writeSolution("psi2.vtu",mesh,pointOutput2,cellOutput2)  

boundaryConditionType = [1;1] # [farfield;airfoil]
# 0 - Natural BC with zero flux
# 1 - Essential BC with constant value
# 2 - Essential freestream BC
essentialPhi = [1.;0.]
pointOutput3,cellOutput3 = solvePotentialFlow(meshFileName,boundaryConditionType,freestreamV,freestreamAlpha,essentialPhi,dummyId)
boring2D.writeSolution("psi3.vtu",mesh,pointOutput3,cellOutput3)  


vY1 = cellOutput["vY"][1][219]
vY2 = cellOutput2["vY"][1][219]

phi = pointOutput["phi"] - vY1/vY2*pointOutput2["phi"]
vX = cellOutput["vX"][1] - vY1/vY2*cellOutput2["vX"][1]
vY = cellOutput["vY"][1] - vY1/vY2*cellOutput2["vY"][1]

pointOutput4 = Dict("phi"=>phi)
cellOutput4 = Dict("vX"=>[vX],"vY"=>[vY])

# Write output data to file
boring2D.writeSolution("phi4.vtu",mesh,pointOutput4,cellOutput4) 