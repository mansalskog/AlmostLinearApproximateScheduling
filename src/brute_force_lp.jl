using LinearAlgebra: I, det

"solve the LP max c^T x s.t. A x <= b and x >= 0 where length(b) == 2 by iterating over vertices"
function brute_force_lp(
    A::Matrix{R},
    b::Vector{R},
    c::Vector{R})::Union{Tuple{Vector{R},R},Nothing} where R<:Real
    n = length(c)
    @assert length(b) == 2 "can only handle two constraints"
    @assert size(A) == (2,n) "invalid size for A"

    A_2 = [A  I]
    c_2 = [c ; 0 ; 0]
    obj_star = convert(R, -Inf) # or typemin(R) ?
    kl_star = nothing
    x_star = nothing
    for k in 1:n+2, l in k+1:n+2
        A_kl = A_2[:,[k,l]]
        if det(A_kl) == 0
            # If determinant is zero then x[1],x[2]
            # cannot be a basic solution.
            # Therefore we can skip this case w.l.o.g.
            continue
        end
        x = A_kl \ b
        if x[1] < 0 || x[2] < 0
            # Do not allow negative solutions.
            continue
        end
        obj = c_2[k]*x[1] + c_2[l]*x[2]
        if obj > obj_star
            obj_star = obj
            kl_star = [k,l]
            x_star = x
        end
    end
    if x_star === nothing
        return nothing
    end
    x = zeros(R, n)
    if kl_star[1] <= n
        x[kl_star[1]] = x_star[1]
    end
    if kl_star[2] <= n
        x[kl_star[2]] = x_star[2]
    end
    @assert A * x <= b .+ EXPRESSION_TOL "not feasible!"
    @assert abs(c' * x - obj_star) <= EXPRESSION_TOL "wrong objective!"
    (x,obj_star)
end
