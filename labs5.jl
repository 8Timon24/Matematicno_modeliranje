# We'll need the `norm` function from LinearAlgebra
using LinearAlgebra

# the actual `newton` function
"""
    newton(F, JF, X0; kwargs...)

Find one solution of F(X)=0 for a vector-valued function `F` with Jacobian matrix `JF` using Newton's iteration with initial approximation `X0`. Returns a named tuple (X = X, n = n).

Keyword arguments are `tol = 1e-8` and `maxit = 100`.

# Examples
```julia-repl
julia> newton(x -> x^2 - 2, x -> 2x, 1.4).X
1.4142135623730951
```
"""
function newton(F, JF, X0; tol = 1e-8, maxit = 100)
    # set two variables which will be used (also) outside the for loop
    X = X0
    n = 1
    # first evaluation of F
    Y = F(X)
    for outer n = 1:maxit
        # execute one step of Newton's iteration
        X = X0 - JF(X0)\Y
        Y = F(X)
        # check if the result is within prescribed tolerance
        if norm(X-X0) + norm(Y) < tol
            break;
        end
        # otherwise repeat
        X0 = X
    end

    # a warning if maxit was reached
    if n == maxit
        @warn "no convergence after $maxit iterations"
    end
    # let's return a named tuple
    return (X = X, n = n)
end

# the `curve` function
function curve(f, gradf, T0, delta, n)
    # first, a correction to T0  - to obtain a point on the curve
    g = gradf(T0...)
    ng = [-g[2]; g[1]]
    # define the function (and its Jacobi matrix) 
    # for the system which will give us the correction to T0
    F(X) = [f(X...); ng'*(X - collect(T0))]
    JF(X) = [gradf(X...)'; ng']
    # define an empty vector of 2-tuples
    K = Tuple{Float64, Float64}[]
    # add the first point (the correction to T0) to that vector
    push!(K, Tuple(newton(F, JF, collect(T0)).X))

    # obtaining succesive points on the curve
    for k in 2:n
        # define the function (and its Jacobi matrix) 
        local F(X) = [f(X...); (X - collect(K[end]))'*(X - collect(K[end])) - delta^2]
        local JF(X) = [gradf(X...)'; 2*(X - collect(K[end]))']
        # determine the tangent through the current point
        g = gradf(K[end]...)
        ng = [-g[2]; g[1]]
        # evaluate the next point and add it to K
        push!(K, Tuple(newton(F, JF, collect(K[end]) + delta*ng/norm(ng)).X))
        # break the loop once close enough to the starting point (in case of a closed curve)
        if(norm(collect(K[end]) - collect(K[1])) < delta/2)
            break;
        end
    end

    return K
end