using LinearAlgebra
using Test
"""
Checks if the point c is on the line segment ab
"""
function liesOnLineSegment(a, b, c)
    return min(a[1], b[1]) <= c[1] && c[1] <= max(a[1], b[1]) && min(a[2], b[2]) <= c[2] && c[2] <= max(a[2], b[2])
end

"""
Accepts 2 vectors of 4-Tuples of 2 points in the 2D plane and 2 parameters corresponding to them, each representing a polygon 
and returns a vector of 3-tuples of intersections between these 2 polygons 
and the parameters of the curves needed to get the intersection
"""
function poly_intersections(poly_p, poly_q)
    intersections = []
    #p1, p2 are points in R², tp1, tp2 in R are parameters such that p(tp1) = p1 and p(tp2) = p2, same holds for q1, q2, uq1, uq2
    for (p1, p2, tp1, tp2) in poly_p
        for (q1, q2, uq1, uq2) in poly_q
            d1 = det([p2[1]-p1[1] q1[1]-p1[1];
                      p2[2]-p1[2] q1[2]-p1[2]])
            d2 = det([p2[1]-p1[1] q2[1]-p1[1];
                      p2[2]-p1[2] q2[2]-p1[2]])
            d3 = det([q2[1]-q1[1] p1[1]-q1[1];
                      q2[2]-q1[2] p1[2]-q1[2]])
            d4 = det([q2[1]-q1[1] p2[1]-q1[1];
                      q2[2]-q1[2] p2[2]-q1[2]])
            
            if d1*d2 < 0 && d3*d4 < 0
                #we have an intersection, we can calculate it from parametrized line segments
                #solving for p(t) = q(u) where in this case p is a parametrization of the line
                #from points p1 to p2 and q is a parametrization of the line from q2 to q1
                #v1, v2 make it a bit easier to write the system 
                v1 = p2-p1
                v2 = q2-q1
                A = [v1[1] -v2[1]; v1[2] -v2[2]]
                b = [q1[1] - p1[1], q1[2] - p1[2]]

                #backslash operator for solving the linear system
                x = A \ b 
                s, r = x[1], x[2]
                #linear approximation for parameters, will work good if h is small
                t = tp1 + s*(tp2-tp1)
                u = uq1 + r*(uq2-uq1)
                #calculating the actual intersection point
                intersection = p1 + s * v1
                push!(intersections, (intersection, t, u))
            
            elseif iszero(d1) && liesOnLineSegment(p1, p2, q1)
                s = norm(q1 - p1) / norm(p2 - p1)
                t = tp1 + s*(tp2-tp1)
                push!(intersections, (q1, t, uq1))
            
            elseif iszero(d2) && liesOnLineSegment(p1, p2, q2)
                s = norm(q2 - p1) / norm(p2 - p1)
                t = tp1 + s*(tp2-tp1)
                push!(intersections, (q2, t, uq2))
            
            elseif iszero(d3) && liesOnLineSegment(q1, q2, p1)
                r = norm(p1-q1) / norm(q2-q1)
                u = uq1 + r*(uq2-uq1)
                push!(intersections, (p1, tp1, u))
            
            elseif iszero(d4) && liesOnLineSegment(q1, q2, p2)
                r = norm(p2-q1) / norm(q2-q1)
                u = uq1 + r*(uq2-uq1)
                push!(intersections, (p2, tp2, u))
            
            else
                #no intersection
            end
        end
    end
    
    return intersections
end

"""
The function accepts p, q which are collumn vectors and represent the curves of which we are trying to find intersections of, 
pdot, qdot which are also collumn vectors and are their derivatives, intp, intq are 2-tuples and represent the interval at which
we are interested to find intersections at, h is a float and represents the step size for evaluating the curves, so it should be sufficiently small.
polychains is an optional parameter if we want to also return the points used for building the polychains alongside their intersections. 
"""
function intersectionOfCurves(p, pdot, intp, q, qdot, intq, h; polychains = false, maxiter = 1000, tolerance = 1e-12)
    #get the parameters for the polygonal approximation from the given intervals and step size
    n = Int(floor((intp[2]-intp[1])/h))
    m = Int(floor((intq[2]-intq[1])/h))
    p_params = [intp[1]+k*h for k in 0:n]
    q_params = [intq[1]+k*h for k in 0:m]
    
    if p_params[end] < intp[2]
        push!(p_params, intp[2]) #if the length of the interval isnt divisible by h we add the last parameter, which is the right border of the interval, separately
    end 
    
    if q_params[end] < intq[2]    
        push!(q_params, intq[2])
    end
    
    #build the polygonal chain as tuples of points and parameters (p1, p2, tp1, tp2), (p2, p3, tp2, tp3), ..., (pn-1, pn, tpn-1, tpn)
    poly_p = [(p(p_params[i]), p(p_params[i+1]), p_params[i], p_params[i+1]) for i in 1:length(p_params)-1]
    poly_q = [(q(q_params[i]), q(q_params[i+1]), q_params[i], q_params[i+1]) for i in 1:length(q_params)-1]
  
    polychain_intersections = poly_intersections(poly_p, poly_q) #a vector of intersections and parameters at those intersections
    
    improved_intersections = []
    #Run the newton method on all starting approximations of intersections to improve them
    for (I, t0, q0) in polychain_intersections
        
        x_r = [t0, q0]
        i = 1

        for outer i in 1:maxiter
            F = p(x_r[1]) - q(x_r[2])
            JF = [pdot(x_r[1]) -qdot(x_r[2])]
            delta = JF \ F
            x_r = x_r - delta
            if norm(delta) < tolerance
                break
            end
        end

        if i == maxiter
            @warn "no convergence after $maxiter iterations"
        end
        push!(improved_intersections, Tuple(p(x_r[1])))
    end
    
    Q = [Tuple(polychain_intersections[i][1]) for i in 1:length(polychain_intersections)]
    
    #I need this for plotting the polychains
    if polychains
        return (P=improved_intersections, Q=Q, Polys = (poly_p, poly_q))
    end
    
    return (P = improved_intersections, Q = Q)
end

function get_poly_coords(poly_chain)
    # Because we have (p1, p2), (p2, p3), ..., (pn-1, pn) 
    # it is enough to take the first point out from every pair
    x_coords = [pair[1][1] for pair in poly_chain]
    y_coords = [pair[1][2] for pair in poly_chain]

    # And add pn separately in the end 
    push!(x_coords, poly_chain[end][2][1])
    push!(y_coords, poly_chain[end][2][2])
    
    return x_coords, y_coords
end

#   Testing

Test.@testset "Tests for intersection of curves" begin

    # Easy to draw and check 
    Test.@test liesOnLineSegment([0,0], [2,2], [1,1]) == true
    Test.@test liesOnLineSegment([0,0], [1,1], [2,2]) == false

    # Testing the actual algorithm, cos, sin also easy to check by hand
    p_test(t) = [t, cos(t)]
    pdot_test(t) = [1.0, -sin(t)]
    q_test(u) = [u, sin(u)]
    qdot_test(u) = [1.0, cos(u)]

    res = intersectionOfCurves(p_test, pdot_test, (-2π, 2π), q_test, qdot_test, (-2π, 2π), 0.1)
    
    # There should be 4 intersections found
    Test.@test length(res.P) == 4
    
    # (pi/4, cos(pi/4)) is one of the intersections 
    target = (π/4, cos(π/4))
    found = false

    for pt in res.P
        # Check if current point is close enough to the target
        if isapprox(pt[1], target[1], atol=1e-8) && isapprox(pt[2], target[2], atol=1e-8)
            found = true
            break 
        end
    end

    Test.@test found == true

end