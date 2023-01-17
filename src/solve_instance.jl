using LinearAlgebra: norm

const VALIDATE_SOLUTION = false
const VALIDATE_TRANSFORMED = false

"Compute the total processing times for precedence chains ending at each job."
function longest_prec_chains(inst::Instance, prec_graph::Graph{Vector{Int}}, prec_order::Vector{Int})::Vector{Int}
    q = zeros(Int, inst.n)
    for v in prec_order
        q′ = inst.p[v] + maximum(q[u] for u in Δ_inc(prec_graph, v); init=0)
        @assert q[v] == 0
        q[v] = q′
    end
    q
end

function call_mwu_solver(
    inst::Instance,
    ϵ::R,
    q::Vector{Int},
    τ::Vector{R},
    η::Vector{R},
    d_min::Vector{Int},
    d_max::Vector{Int})::Matrix{R} where R<:Real

    J = 1:inst.n
    D = length(τ)-1
    fix_x = Matrix{Union{Nothing,Int}}(nothing, inst.n, D+1)
    x_to_v = Matrix{Union{Nothing,Int}}(nothing, inst.n, D+1) # x_to_v[j,d] = v
    v_to_x = Tuple{Int,Int}[] # v_to_x[v] = (j,d)
    for j in J, d in 1:D+1
        if d <= d_min[j] || τ[d] < q[j]
            # fixed to 0 by (6)
            fix_x[j,d] = 0
        elseif d_max[j] <= d
            # fixed to 1 by (7)
            fix_x[j,d] = 1
        else
            # add new vertex
            push!(v_to_x, (j,d))
            x_to_v[j,d] = length(v_to_x)
        end
    end
    N̂ = length(v_to_x)

    E = NTuple{2,Int}[]
    for j in J, d in 1:D
        # constraint (3), note upper limit D (not D+1)
        if !isnothing(x_to_v[j,d]) && !isnothing(x_to_v[j,d+1])
            push!(E, (x_to_v[j,d],x_to_v[j,d+1]))
        end
    end
    for (j1,j2) in inst.prec, d in 1:D+1
        # constraint (4), note direction of edge
        if !isnothing(x_to_v[j1,d]) && !isnothing(x_to_v[j2,d])
            push!(E, (x_to_v[j2,d],x_to_v[j1,d]))
        end
    end
    G = Graph{Vector{Int}}(N̂, E)

    P = zeros(R, D+1, N̂)
    for v in 1:N̂
        (j,d) = v_to_x[v]
        P[d,v] = inst.p[j] / inst.m / τ[d]
    end

    a = zeros(R, N̂)
    for v in 1:N̂
        (j,d) = v_to_x[v]
        a[v] = inst.w[j] * η[d]
    end

    ϕ = ϵ * sum(inst.w)
    @assert 0 < ϕ < norm(a, 1) / 2 "bounds on parameter phi"

    ϵ′ = ϵ # TODO later perhaps use ϵ/4 so that (1+ϵ′)²+ϵ′ = 1+3/4*ϵ+ϵ^2/16 <= 1+ϵ
    x_value = mwu_solve_lp(a, P, G, ϵ′, ϕ)

    # convert to solution of original LP
    x_matrix = zeros(R, inst.n, D+1)
    for j in J, d in 1:D+1
        if !isnothing(fix_x[j,d])
            @assert isnothing(x_to_v[j,d])
            x_matrix[j,d] = fix_x[j,d]
        else !isnothing(x_to_v[j,d])
            @assert isnothing(fix_x[j,d])
            x_matrix[j,d] = x_value[x_to_v[j,d]]
        end
    end
    x_matrix
end

"Solve the LP-relaxation, computing the fractional completion times."
function lp_relaxation(inst::Instance, ϵ::R, q::Vector{Int})::Vector{R} where R<:Real
    println("Using real type ", R)
    τ = [zero(R)]
    η = R[]

    # the range of d is [1,D+1] and not [0,D] as in the paper
    d = 1
    p_tot = sum(inst.p)
    while τ[d] < p_tot
        d += 1
        push!(τ, (1+ϵ)^(d-1)) # to τ[d]
        push!(η, τ[d] - τ[d-1]) # to η[d-1]
    end
    D = d-1
    # here we know that τ[D+1] >= sum(inst.p)
    # println("Using ", D+1, " time steps from ", τ[1], " to ", τ[D+1])

    # TODO should be changed later to
    # handle non-polynomial processing times
    # also the weights should be changed
    d_min = ones(Int, inst.n)
    d_max = fill(D+1, inst.n)

    x_value = call_mwu_solver(inst, ϵ, q, τ, η, d_min, d_max)

    C_frac = [τ[D+1] - sum(η[d] * x_value[j,d] for d in 2:D) for j in 1:inst.n]
    C_frac
end

"Round the solution of the LP-relaxation to integer completion times."
function round_completion_times(inst::Instance, C::Vector{<:Real}, prec_order::Vector{Int})::Vector{Int}
    ρ = ordering_to_rank(prec_order)
    J_sorted = sort(collect(1:inst.n), by=j -> (C[j],ρ[j]), alg=MergeSort)

    C = zeros(Int, inst.n)
    for (j,C_j) in zip(J_sorted, cumsum(inst.p[J_sorted]))
        C[j] = C_j
    end
    C
end

"Validate a solution to an instance. This is O(nκ), so it should not be used when run time is measured."
function validate_solution(inst::Instance, C::Vector{Int})
    @assert all(C[j1] <= C[j2]-inst.p[j2] for (j1,j2) in inst.prec) "precedence constraint violated"
    overlap_constr = t -> sum(t in C[j]-inst.p[j]+1:C[j] for j in 1:length(inst.p)) <= inst.m
    @assert all(overlap_constr, 1:sum(inst.p)) "overlap constraint violated"
end

"Find a solution for an instance."
function solve_instance(inst::Instance, ϵ::Real)::Solution
    prec_graph = Graph{Vector{Int}}(inst.n, inst.prec)
    prec_order = topological_ordering(prec_graph)

    q = longest_prec_chains(inst, prec_graph, prec_order)
    C_frac = lp_relaxation(inst, ϵ, q)

    C = round_completion_times(inst, C_frac, prec_order)

    if VALIDATE_SOLUTION
        validate_solution(inst, C)
    end

    wct = sum(inst.w .* C)
    wct_frac = sum(inst.w .* C_frac)
    Solution(ϵ, C_frac, wct_frac, C, wct)
end
