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