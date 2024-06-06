using SparseArrays

mutable struct SparseTriplet
        paddedLength::Int64 # Length of the triplet vectors including padding
        nTriplets::Int64 # Number of triplets
        i::Vector{Int64} # [paddedLength, 1] row indices
        j::Vector{Int64} # [paddedLength, 1] column indices
        s::AbstractVector # [paddedLength, 1] matrix values
end

function SparseTriplet(paddedLengthIn::Int64)
    nTriplets = Int64(0);
    i = Vector{Int64}(undef,paddedLengthIn);
    j = Vector{Int64}(undef,paddedLengthIn); 
    s = Vector{Number}(undef,paddedLengthIn); 
    return SparseTriplet(paddedLengthIn,nTriplets,i,j,s);
end

function padVectors!(sparseTriplet::SparseTriplet)
    # Doubles the size of triplet vectors by padding with zeros
    oldPaddedLength = sparseTriplet.paddedLength;
    newPaddedLength = 2*oldPaddedLength;
    i = Vector{Int64}(undef,newPaddedLength);
    j = Vector{Int64}(undef,newPaddedLength); 
    s = Vector{Number}(undef,newPaddedLength); 

    i[1:oldPaddedLength] = sparseTriplet.i;
    j[1:oldPaddedLength] = sparseTriplet.j;
    s[1:oldPaddedLength] = sparseTriplet.s;
    
    sparseTriplet.paddedLength = newPaddedLength;
    sparseTriplet.i = i;
    sparseTriplet.j = j;
    sparseTriplet.s = s;
    return 0
end

function addMatrix!(sparseTriplet::SparseTriplet,M::AbstractArray,gDof::Vector{Int64})
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

function convertToSparseMatrix(sparseTriplet::SparseTriplet,n::Int64,m::Int64)
    # convert triplets to sparse matrix
    return sparse(  sparseTriplet.i[1:sparseTriplet.nTriplets], 
                    sparseTriplet.j[1:sparseTriplet.nTriplets],
                    sparseTriplet.s[1:sparseTriplet.nTriplets],
                    n,m)
end
