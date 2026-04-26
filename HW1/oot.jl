using LinearAlgebra
"""
The function naive takes in 2 vectors with length n of k-tuples, 
these represent the input and output data (points before and after a rigid transformation) 
for calculating the translation vector b and a rotation matrix Q. It returns Q and b as a 2-tuple
"""
function naive(X, Y)

    n = length(X); #num of points
    k = length(X[1]); #the dimension 
    
    X1 = [stack(X)' ones(n)]
    Y1 = stack(Y)'
    
    
    Qb = X1\Y1
    b = Qb[end, :]
    Q = Qb[1:end-1, :]'
    F = qr(Q)

    #After learning a bit about what qr returns and potentially distorts, it would seem the
    #sign of the columns in Q from qr can differ compared to the ones in the original Q
    #and that happens because we get the original matrix as a product of columns of Q -> q_i with diagonal 
    #elements of R->r_ii, so in order to keep the same orientation we "force" R_ii to be positive, 
    #although we dont explicitly do that to R because it would be a waste, so we just act like it is and
    #fix our orthogonal Q accordingly 
    d = sign.(diag(Matrix(F.R))) #here we collect signs of diagonal elements of R into a vector  
    Q_new = Matrix(F.Q) * Diagonal(d) #we have to multiply Q_new with a diagonal matrix not a vector
    return (Q_new, b)
end

"""
The function kabsch takes two vectors of length n containing k-tuples,
representing points before and after a rigid transformation. It computes
the rotation matrix Q and translation vector b using the Kabsch algorithm
and returns them as a tuple (Q, b). The inputs X and Y are vectors of
n-tuples, Q is a kxk matrix, and b is a column vector of length k. 
"""
function kabsch(X, Y)

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
    d = -sign(det(U * Vt)) #find d = +-1 
    D[k,k] = d #put d into matrix D
    Q = U*D*Vt #obtain our orthogonal matrix Q
    
    b = collect(y1) - Q*collect(x1) #get b
    
    return (Q, b)
end
