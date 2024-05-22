"Triangle elements"
struct TriangleElements
    # nodeIDs::Vector{Int32} # [1,3] node numbers
    # x::Vector{Float64} # [1,3] x node locations
    # y::Vector{Float64} # [1,3] y node locations
    # invJ::Matrix{Float64} # [2,2] inverse of Jacobian matrix (constant inside elements)
    area::Float64 # [double] element area
    dNdX::Matrix{Float64} # [2,3] Physical shape function derivatives - constant inside the element
end

"""
    TriangleElements(mesh::Mesh2D)

Construct array of TriangleElements elements from mesh data
"""
function TriangleElements(mesh::Mesh2D)
    
    # Process node IDs
    nTri,mTri = size(mesh.triangles);
    if nTri < 1; error("No elements"); end
    if mTri !=3; error("mTri~=3");     end
    
    # Shape function derivatives are constant inside the element
    dNXi = Matrix{Float64}(undef,2,3);
    dNXi[1,:] = [-1, 1, 0]; # dNdxi  = [-1, 1, 0];
    dNXi[2,:] = [-1, 0, 1]; # dNdeta = [-1, 0, 1];

    # Process elements
    detJinvJ = Array{Float64}(undef,2,2);
    triangleElements = Array{TriangleElements}(undef,nTri); # this doesn't have enough information to allocate memory
    for i = 1:nTri
        nodeIDs = mesh.triangles[i,:];
        x = [ mesh.nodes[nodeIDs[1],1], mesh.nodes[nodeIDs[2],1], mesh.nodes[nodeIDs[3],1]];
        y = [ mesh.nodes[nodeIDs[1],2], mesh.nodes[nodeIDs[2],2], mesh.nodes[nodeIDs[3],2]];
        detJ = (x[2]-x[1])*(y[3]-y[1]) - (x[3]-x[1])*(y[2]-y[1]);
        detJinvJ[1,1] = y[3]-y[1];
        detJinvJ[1,2] = y[1]-y[2];
        detJinvJ[2,1] = x[1]-x[3];
        detJinvJ[2,2] = x[2]-x[1];
        invJ = (1.0/detJ)*detJinvJ;
        area = 0.5*detJ;
        dNdX = invJ*dNXi;

        # triangleElements[i] = TriangleElements(nodeIDs,x,y,invJ,area,dNdX);
        triangleElements[i] = TriangleElements(area,dNdX);
    end
    return triangleElements;
end

# function gaussPointsAndWeights()
#     # 3-Point Gauss integration points & weight factor
#     r = [2/3 1/6 1/6];
#     s = [1/6 1/6 2/3];
#     w3 = 1/3;
#     return r, s, w3
# end
# function [N1,N2,N3]=vectorShapeFunValsAtIntPoints(r,s)
#     # 2D vector shape function values at integration points
#     # shapeFunctions = @(r,s) [1-r-s, r, s];
#     N1 = [1-r(1)-s(1),           0, r(1),    0, s(1),   0   ;
#         0, 1-r(1)-s(1),    0, r(1),    0, s(1) ];
#     N2 = [1-r(2)-s(2),           0, r(2),    0, s(2),   0   ;
#         0, 1-r(2)-s(2),    0, r(2),    0, s(2) ];
#     N3 = [1-r(3)-s(3),           0, r(3),    0, s(3),   0   ;
#         0, 1-r(3)-s(3),    0, r(3),    0, s(3) ];
# end
# function scalarShapeFunValsAtIntPoints(r,s)
#     # Scalar shape function values at integration points
#     # shapeFunctions = @(r,s) [1-r-s, r, s];
#     N1 = [1-r[1]-s[1], r[1], s[1]];
#     N2 = [1-r[2]-s[2], r[2], s[2]];
#     N3 = [1-r[3]-s[3], r[3], s[3]];
#     return N1, N2, N3
# end
