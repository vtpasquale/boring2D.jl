
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
essentialPhi = [0.;0.]
zeroVyEle = 219

freestreamV = 100.
freestreamAlpha = 8. *pi/180.
solveStream(meshFileName,boundaryConditionType,freestreamV,freestreamAlpha,essentialPhi,zeroVyEle)