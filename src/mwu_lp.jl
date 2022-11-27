"solve the LP with parameters a, P and G (giving Q) using Algo 1"
function mwu_solve_lp(
    a::Vector{R},
    P::Matrix{R},
    G::Graph,
    ϵ::R,
    ϕ::R)::Vector{R} where R<:Real
    @assert 0 < ϵ < 1
    @assert ϕ > 0

    (m̂,n̂) = size(P)
    @assert length(a) == n̂

    t = zero(R)
    ρ = log(m̂) / ϵ^2
    x = zeros(R, n̂)
    u = ones(R, m̂)' # row vector

    @show ϕ
    @show ϵ
    while t < 1
        b = (u/norm(u, 1) * P)' # row vector but transposed for convenience

        ϵ′ = sqrt(1+ϵ) - 1 # ϵ/4
        y = oracle_solve_lp(a, b, G, ϵ′, ϕ)

        # note that we validate with the original epsilon as we want by ≤ (1+ϵ′)^2 ≤ 1+ϵ
        validate_oracle_solution(b, G, y, ϵ)
        # TODO disable validation

        δ = min(minimum(1/(ρ * P[i,:]'*y) for i in 1:m̂), 1-t)

        for i in 1:m̂
            # P[i,:]'*y == P[i:i,:] * y
            u[i] = u[i] * exp(δ*ϵ*ρ * P[i,:]' * y)
        end
        x = x + δ*y
        t = t + δ
    end
    x
end

function validate_oracle_solution(b::Vector{R}, G::Graph, y::Vector{R}, ϵ::R) where R<:Real
    @assert b'*y <= 1 + ϵ + EXPRESSION_TOL
    # NB. Here we have no numerical problems because of the postprocessing step
    @assert all(0 <= yi <= 1 for yi in y)
    @assert all(y[src_vertex(G, e)] <= y[dst_vertex(G, e)] for e in 1:G.m)
end

function oracle_solve_lp(
    a::Vector{R},
    b::Vector{R},
    G::Graph,
    ϵ::R,
    ϕ::R)::Vector{R} where R<:Real
    @assert 0 < ϵ < 1
    @assert ϕ > 0

    (a,b,G,S,T,fix_y,new_v) = transform_data(a,b,G)
    if VALIDATE_TRANSFORMED
        validate_transformed(a,b,G,S,T)
    end

    γ = ϕ/3
    Γ = [γ]
    while γ < norm(a, 1)
        γ = (1+ϵ) * γ
        push!(Γ, γ)
    end

    ã = zeros(R, length(Γ))
    b̃ = zeros(R, length(Γ))
    S′_list = Vector{Int}[]
    for (i,γ) in enumerate(Γ)
        (S′,f) = find_S′_and_f(G, a, b, S, T, γ, ϵ, ϕ)

        # T(S′) = {t : S′ -> t and t ∈ T}
        TS′ = [t for t in reach_from(G, S′) if b[t] > 0]

        validate_S′_and_f(G, a, b, S, S′, TS′, f, γ, ϵ, ϕ)

        push!(S′_list, S′)

        ã[i] = sum(a[S′])
        # as b[v] == 0 for all v in V \ T we have
        # sum(b[TS′]) == sum(b[reach_from_S])
        b̃[i] = sum(b[TS′])
    end

    y_value = y_from_nfp_gamma(G, ã, b̃, S′_list, ϵ)
    y_value = postprocess_y(y_value, G, a, b)

    # convert back to original y
    y_result = zeros(R, length(fix_y))
    for i in 1:length(fix_y)
        if !isnothing(fix_y[i])
            @assert isnothing(new_v[i])
            y_result[i] = fix_y[i]
        else !isnothing(new_v[i])
            @assert isnothing(fix_y[i])
            y_result[i] = y_value[new_v[i]]
        end
    end
    y_result
end

function validate_S′_and_f(
    G::Graph,
    a::Vector{R},
    b::Vector{R},
    S::Vector{Int},
    S′::Vector{Int},
    TS′::Vector{Int},
    f::Vector{R},
    γ::R,
    ϵ::R,
    ϕ::R) where R<:Real
    # S_diff = S \ S′
    S_diff = [s for s in S if !(s in S′)]
    # S_out = δ_out(S)
    S_out = [e for s in S for e in δ_out(G, s)]
    @assert all(a[src_vertex(G, e)] > 0 for e in S_out)
    @assert all(a[dst_vertex(G, e)] == 0 for e in S_out)
    val_f = sum(f[S_out])
    @assert sum(a[S_diff]) + γ * sum(b[TS′]) / (1+ϵ) <= val_f + ϕ/3 + EXPRESSION_TOL
end

function y_from_nfp_gamma(G::Graph, ã::Vector{R}, b̃::Vector{R}, S′_list::Vector{Vector{Int}}, ϵ::R)::Vector{R} where {R<:Real}
    Γ_length = length(S′_list)

    # solve using our own LP solver
    A = [ones(R, 1, Γ_length) ; reshape(b̃, 1, :)]
    # probably better for accuracy to not divide
    b = [one(R), (1+ϵ)^2]
    c = ã
    (z_value,_) = brute_force_lp(A, b, c)

    # this sets y[v] = sum(z[gamma] where (v is reachable from S′[gamma])
    y = zeros(R, G.n)
    for i in 1:Γ_length
        for v in reach_from(G, S′_list[i])
            y[v] += z_value[i]
        end
    end
    y
end

"transform the output according to the last part of Lemma 3.8"
function postprocess_y(
    y::Vector{R},
    G::Graph,
    a::Vector{R},
    b::Vector{R})::Vector{R} where R<:Real

    # we check a[v] != 0 instead of a in S, because the former is O(1) while the latter is O(|S|)

    G_order = topological_ordering(G)
    # this sets y[v] = max(y[s] for s in S if (v is reachable from s) for v in V - S
    for v in G_order
        if a[v] != 0
            continue
        end
        y[v] = maximum(y[u] for u in Δ_inc(G, v); init=zero(R))
    end
    # this sets y[v] = min(y[t] for t in T if (t is reachable from v) for v in V - T
    for v in Iterators.reverse(G_order)
        if b[v] != 0
            continue
        end
        y[v] = minimum(y[u] for u in Δ_out(G, v); init=one(R))
    end
    y
end

"transform the data according to Lemma 3.8"
function transform_data(
    a::Vector{R},
    b::Vector{R},
    G::Graph)::Tuple{
        Vector{R},
        Vector{R},
        Graph,
        Vector{Int},
        Vector{Int},
        Vector{Union{Nothing,Int}},
        Vector{Union{Nothing,Int}}} where R<:Real
    S = [s for s in 1:G.n if a[s] > 0]
    T = [t for t in 1:G.n if b[t] > 0]

    fix_y = Vector{Union{Nothing,Int}}(nothing, G.n)
    new_v = Vector{Union{Nothing,Int}}(nothing, G.n) # new_v[v] = v′
    old_v = Int[] # old_v[v′] = v
    reach_from_S = falses(G.n)
    for v in reach_from(G, S)
        reach_from_S[v] = true
    end
    reach_to_T = falses(G.n)
    for v in reach_to(G, T)
        reach_to_T[v] = true
    end
    for v in 1:G.n
        if !reach_from_S[v]
            # not reachable from S, remove and fix y[v] to 0
            fix_y[v] = 0
        elseif !reach_to_T[v]
            # not reachable from T, remove and fix y[v] to 1
            fix_y[v] = 1
        else
            # otherwise, keep the vertex
            push!(old_v, v)
            new_v[v] = length(old_v)
        end
    end

    # new_edges = [new_v[v],new_v[u] for (v,u) in edges(G) if v and u is kept]
    new_edges = [(new_v[src_vertex(G, e)],new_v[dst_vertex(G, e)])
                 for e in edges(G)
                     if !isnothing(new_v[src_vertex(G, e)]) && !isnothing(new_v[dst_vertex(G, e)])]

    new_n = length(old_v)
    new_a = zeros(R, new_n)
    new_b = zeros(R, new_n)

    for s in S
        if new_v[s] === nothing
            continue
        end
        s′ = new_n + 1
        new_n += 1
        push!(new_a, a[s])
        push!(new_b, zero(R))
        # push!(new_S, s′) (see below)
        push!(new_edges, (s′,new_v[s]))
    end

    for t in T
        if new_v[t] === nothing
            continue
        end
        t′ = new_n + 1
        new_n += 1
        push!(new_a, zero(R))
        push!(new_b, b[t])
        # push!(new_T, t′) (see below)
        push!(new_edges, (new_v[t],t′))
    end

    new_G = Graph(new_n, new_edges)
    # these could also be created by pushing each s′ and t′, resp. (see above)
    new_S = [s for s in 1:new_n if new_a[s] > 0]
    new_T = [t for t in 1:new_n if new_b[t] > 0]

    (new_a,new_b,new_G,new_S,new_T,fix_y,new_v)
end

"verify that the properties of Lemma 3.8 hold"
function validate_transformed(
    a::Vector{R},
    b::Vector{R},
    G::Graph,
    S::Vector{Int},
    T::Vector{Int}) where R<:Real

    @assert isempty(intersect(S, T)) "intersection of S and T must be empty"

    for s in S
        @assert isempty(δ_inc(G, s)) "sources should have no incoming edges"
    end

    for t in T
        @assert isempty(δ_out(G, t)) "sinks should have no outgoing edges"
    end

    @assert length(reach_from(G, S)) == G.n "all vertices should be reachable from S"
    @assert length(reach_to(G, T)) == G.n "T should be reachable from all vertices"
end
