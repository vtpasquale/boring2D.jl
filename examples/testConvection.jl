
using boring2D

cd("examples")

meshFileName = "5000NUcav.plt";
# meshFileName = "square.su2";

k = 1.0 # Conductivity
boundaryTemps = [100.0, 200.0, 300.0, 400.0]

# solve
(Kff,rhs) =solveConvection(meshFileName,k,boundaryTemps)


using SparseArrays
@time Kff\rhs; # 66 seconds - wth?
@time Kff\ sparse(rhs); % 67 seconds
@time Matrix(Kff)\rhs; # .2 seconds



# function profile_test(n,meshFileName,k,boundaryTemps)
#     for i = 1:n
#         solveConvection(meshFileName,k,boundaryTemps)
#     end
# end

# using ProfileView
# ProfileView.@profview profile_test(1,meshFileName,k,boundaryTemps)  # run once to trigger compilation (ignore this one)
# ProfileView.@profview profile_test(1,meshFileName,k,boundaryTemps)

# @time @profview profile_test(1,meshFileName,k,boundaryTemps)  # run once to trigger compilation (ignore this one)
# @time @profview profile_test(4,meshFileName,k,boundaryTemps)



# 
# filename = "5000NUcav.plt";
# mesh = readMesh(filename)

# triangleElements = boring2D.TriangleElements(mesh)

# k = 1.0; # Conductivity

# # Boundary temperatures (at four boundaries)
# boundaryTemps = [100, 200, 300, 400];

# # Number of elements and nodes
# nElements = size(mesh.triangles,1);
# nNodes = size(mesh.nodes,1);
# nEdges = size(mesh.edges,1)

# # One Temperature DOF per node
# nDof = nNodes;

# # # Global matrices (slow method)
# # K = spzeros(nDof,nDof);

# # Global matrices with triplets
# Ks = boring2D.SparseTriplet(nNodes);

# # Assemble stiffness matrix
# for i = 1:nElements
#     dof = mesh.triangles[i,:]
#     dNdX = triangleElements[i].dNdX
#     area = triangleElements[i].area
#     k_e = transpose(dNdX)*k*dNdX*area
#     boring2D.addMatrix!(Ks,k_e,dof)
# end
# K = boring2D.convertToSparseMatrix(Ks,nDof,nDof);

# # Determine fixed and free dof
# uniqueEdges = sort(unique(mesh.edges[:,3]));
# edgeNodes = sort(unique(mesh.edges[:,1:2])); # Fixed dof
# f = trues(nDof);
# f[edgeNodes] .= false;
# s = .!f

# # Assemble boundary temperatures
# u = zeros(nDof,1);
# for i = 1:nEdges
#     # this will overwrite temperatures if nodes are on more than one edge
#     u[mesh.edges[i,1:2]] .= boundaryTemps[mesh.edges[i,3]]; 
# end

# # Partition and solve
# u[f] = K[f,f]\(-K[f,s]*u[s]);

# # Create solution dictionaries
# pointOutput = Dict("Temp"=>u)

# a = Array{Float64}(undef,nElements,1)
# for i in 1:nElements
#     a[i] = triangleElements[i].area
# end
# cellOutput = Dict("area"=>[a])


# # Write output data to file
# boring2D.writeSolution("test.vtu",mesh,pointOutput,cellOutput)
