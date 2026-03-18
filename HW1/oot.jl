using LinearAlgebra
"""
The function naive takes in X, Y vectors of length n containing k-tuples, 
these represent the input and output data (points before and after a rigid transformation) 
for calculating the translation vector b and a rotation matrix Q. It returns Q, a kxk matrix, 
and b, a column vector of length k, as a 2-tuple
"""
function naive(X, Y)
    if isempty(X) || isempty(Y)
        throw(ArgumentError("X and Y must not be empty"))
    elseif length(X) != length(Y)
        throw(DimensionMismatch("lengths of the input vectors X and Y arent equal"))
    elseif length(X[1]) != length(Y[1])
        throw(DimensionMismatch("tuples in X and Y must have the same dimension")) 
    end

    n = length(X); #num of points
    k = length(X[1]); #the dimension 
    
    if n < 2*k
        throw(ArgumentError("our assumption is that n >= 2k"))
    end

    X1 = [stack(X)' ones(n)]
    Y1 = stack(Y)'
    
    
    Qb = X1\Y1
    b = Qb[end, :]
    Q = Qb[1:end-1, :]'
    F = qr(Q)

    return (Matrix(F.Q), b)
end


"""
The function kabsch takes two vectors of length n containing k-tuples,
representing points before and after a rigid transformation. It computes
the rotation matrix Q and translation vector b using the Kabsch algorithm
and returns them as a tuple (Q, b). The inputs X and Y are vectors of
n-tuples, Q is a kxk matrix, and b is a column vector of length k. 
"""
function kabsch(X, Y)
    if isempty(X) || isempty(Y)
        throw(ArgumentError("X and Y must not be empty"))
    elseif length(X) != length(Y)
        throw(DimensionMismatch("lengths of the input vectors X and Y arent equal"))
    elseif length(X[1]) != length(Y[1])
        throw(DimensionMismatch("tuples in X and Y must have the same dimension")) 
    end
    
    n = length(X)
    k = length(X[1])
    
    x1 = (1/n).*(reduce(.+, X))
    y1 = (1/n).*(reduce(.+, Y))
    
    X = [X[i] .- x1 for i in 1:n]
    Y = [Y[i] .- y1 for i in 1:n]
    
    Xm = stack(X)
    Ym = stack(Y)

    C = Ym*Xm'
    Cs = svd(C)
    U = Cs.U 
    Vt = Cs.Vt 
    D = diagm(ones(k))
    d = sign(det(U * Vt))
    D[k,k] = d
    Q = U*D*Vt
    
    b = collect(y1) - Q*collect(x1)
    
    return (Q, b)
end
