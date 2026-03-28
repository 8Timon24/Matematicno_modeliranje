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

    X1 = [stack(X)' ones(n)] #stack column vectors x_i into a matrix as rows and add a column of ones
    Y1 = stack(Y)' #stack column vectors y_i into a matrix
    
    
    Qb = X1\Y1 #Solve the system by least squares method
    b = Qb[end, :] #take b out of the matrix
    Q = Qb[1:end-1, :]' #take Q^T out and transpose it once more so we get Q
    F = qr(Q) #Find the QR decomposition

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
    
    n = length(X) #num of points
    k = length(X[1]) #dimension
    
    x1 = (1/n).*(reduce(.+, X)) #calculate x_
    y1 = (1/n).*(reduce(.+, Y)) #calculate y_
    
    X = [X[i] .- x1 for i in 1:n] #subtract x_ from every x_i to obtain x_i'
    Y = [Y[i] .- y1 for i in 1:n] #subtract y_ from every y_i to obtain y_i'
    
    Xm = stack(X) #order them into a matrix
    Ym = stack(Y) # ^ 

    C = Ym*Xm' #calculate the product Y'X'^T
    Cs = svd(C) #obtain the svd
    U = Cs.U  #take out U from svd
    Vt = Cs.Vt #take out V^T
    D = diagm(ones(k)) #construct an identity of size k
    d = sign(det(U * Vt)) #find d = +-1 
    D[k,k] = d #put d into matrix D
    Q = U*D*Vt #obtain our orthogonal matrix Q
    
    b = collect(y1) - Q*collect(x1) #get b
    
    return (Q, b)
end
