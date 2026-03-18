function naive(X, Y)
    if isempty(X) || isempty(Y)
        throw(ArgumentError("X and Y must not be empty"))
    elseif length(X) != length(Y)
        throw(DimensionMismatch("lengths of the input vectors X and Y arent equal"))
    elseif length(X[1]) != length(Y[1])
        throw(DimensionMismatch("sizes of tuples arent equal")) 
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