using SparseArrays

mutable struct SparseTriplet
        paddedLength::Int32 # Length of the triplet vectors including padding
        nTriplets::Int32 # Number of triplets
        i::Vector{Int32} # [paddedLength, 1] row indices
        j::Vector{Int32} # [paddedLength, 1] column indices
        s::Vector{Float64} # [paddedLength, 1] matrix values
end

function SparseTriplet(paddedLengthIn::Integer)
    nTriplets = Int32(0);
    i = Vector{Int32}(undef,paddedLengthIn);
    j = Vector{Int32}(undef,paddedLengthIn); 
    s = Vector{Float64}(undef,paddedLengthIn); 
    return SparseTriplet(paddedLengthIn,nTriplets,i,j,s);
end

function padVectors!(sparseTriplet::SparseTriplet)
    # Doubles the size of triplet vectors by padding with zeros
    oldPaddedLength = sparseTriplet.paddedLength;
    newPaddedLength = 2*oldPaddedLength;
    i = Vector{Int32}(undef,newPaddedLength);
    j = Vector{Int32}(undef,newPaddedLength); 
    s = Vector{Float64}(undef,newPaddedLength); 

    i[1:oldPaddedLength] = sparseTriplet.i;
    j[1:oldPaddedLength] = sparseTriplet.j;
    s[1:oldPaddedLength] = sparseTriplet.s;
    
    sparseTriplet.paddedLength = newPaddedLength;
    sparseTriplet.i = i;
    sparseTriplet.j = j;
    sparseTriplet.s = s;
    return 0
end

function addMatrix!(sparseTriplet::SparseTriplet,M::Matrix{Float64},gDof::Vector{Int32})
    # Adds matrix M with global DOF gDof to SparseTriplet
    iM,jM,m=findnz(sparse(M));
    
    # Number management
    lengthM = size(m,1);
    nTripletsOld = sparseTriplet.nTriplets;
    nTripletsNew = nTripletsOld + lengthM;
    sparseTriplet.nTriplets = nTripletsNew;
    while nTripletsNew > sparseTriplet.paddedLength
        padVectors!(sparseTriplet);
    end
    
    # Add values to triplets
    index = nTripletsOld+1:nTripletsNew;
    sparseTriplet.i[index] = gDof[iM];
    sparseTriplet.j[index] = gDof[jM];
    sparseTriplet.s[index] = m;
    return 0
end

function convertToSparseMatrix(sparseTriplet::SparseTriplet,n::Integer,m::Integer)
    # convert triplets to sparse matrix
    return sparse(  sparseTriplet.i[1:sparseTriplet.nTriplets], 
                    sparseTriplet.j[1:sparseTriplet.nTriplets],
                    sparseTriplet.s[1:sparseTriplet.nTriplets],
                    n,m)
end
