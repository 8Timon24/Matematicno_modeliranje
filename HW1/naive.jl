function naive(X, Y)

    n = length(X); #num of points
    k = length(X[1]); #the dimension 


    X1 = [stack(X)' ones(n)]
    Y1 = stack(Y)'
    
    
    Qb = X1\Y1
    b = Qb[end, :]
    Q = Qb[1:end-1, :]'
    F = qr(Q)

    return (Matrix(F.Q), b)
end