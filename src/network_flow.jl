using DataStructures: Queue, enqueue!, dequeue!
using DataStructures: BinaryHeap
using LinearAlgebra: norm

# set to true to enable some very expensive asserts
# which are only used for debugging
const VALIDATE_VISITED = false
# if true, print out a debug message on rounding errors
const WARN_ROUNDING_ERRORS = false

"check for rounding errors and check that x > 0.
returns the value x if it is okay, otherwise returns a zero,
so that we can write x = warn_if_negative(x) in order to correct these errors"
function warn_if_negative(x::R) where R<:Real
    if x < 0
        if abs(x) < EXPRESSION_TOL
            WARN_ROUNDING_ERRORS && println(string("ROUNDING ERRORS at ", stacktrace()[2]))
            return zero(R)
        else
            error("Positive variable is negative!")
        end
    end
    return x
end

"find f and S′ satisfying properties of Theorem 3.10, as in algorithm 4 in the paper"
function find_S′_and_f(G::Graph, a::Vector{R}, b::Vector{R}, S::Vector{Int}, T::Vector{Int}, γ::R, ϵ::R, ϕ::R) where R<:Real
    Gl = G
    πl = collect(1:G.n) # Identity mapping of vertices
    π_edgel = collect(1:G.m) # Identity mapping of edges
    fl = zeros(R, Gl.m)
    L = floor(Int, log(1+ϵ, 3*norm(a, 1)/ϕ))

    topo_order = topological_ordering(G)
    ρ = ordering_to_rank(topo_order)

    for ℓ = 0:L
        (Gl,πl,π_edgel,fl) = inc_len(ℓ, G, Gl, πl, π_edgel, a, b, S, T, γ, fl, ρ)
    end
    f′_list = [(π_edgel[e], warn_if_negative(flow)) for (e,flow) in enumerate(fl)]
    S′ = find_S′(Gl, fl, a, b, S, T, γ, ϵ, ϕ)

    # convert flow list to flow vector
    f = zeros(R, G.m)
    for (e,flow) in f′_list
        f[e] += flow
    end
    (S′,f)
end

"find the corresponding edge e′ in a handled sub-graph (G^{i,-},π′),
to e in the graph G"
function get_edge_inv(π′_edge_inv::Vector{Vector{Tuple{Int,Int}}}, e::Int, i::Int)::Int
    e′_list = [e′ for (k,e′) in π′_edge_inv[e] if k == i]
    @assert length(e′_list) == 1
    first(e′_list)
end

"perform the function called inc-len in the paper"
function inc_len(
    l::Int, G::Graph,
    G°::Graph, π°::Vector{Int}, π°_edge::Vector{Int},
    a::Vector{R}, b::Vector{R}, S::Vector{Int}, T::Vector{Int},
    γ::R, f°::Vector{R}, ρ::Vector{Int}) where R<:Real
    (Si,Ti) = find_Si_and_Ti(l-1, G°, a, S, T, f°)
    (adj_Gi_plus,adj_Gi_minus) = find_Gi(G, Si, Ti, l)

    # adjacency lists for supp(f°)
    supp_out = [Int[] for v in 1:G°.n]
    supp_inc = [Int[] for v in 1:G°.n]
    for e in 1:G°.m
        if f°[e] > 0
            push!(supp_out[src_vertex(G°, e)], e)
            push!(supp_inc[dst_vertex(G°, e)], e)
        end
    end

    (G′,π′,π′_edge,π′_edge_inv) = construct_G′(l, G, a, b, adj_Gi_minus)

    f′_plus_list = Vector{Vector{Tuple{Int,R}}}(undef, l+1)
    f′_minus_list = Vector{Vector{Tuple{Int,R}}}(undef, l)
    scratch_space = zeros(R, G°.n)
    for i in 1:l+1
        # f_plus is f∘(i,+) in the paper
        (f_list_plus,visited_V_plus) = find_flow_from_S(G°, π°, f°, supp_out, Si[i], ρ, scratch_space)
        for (e,flow) in f_list_plus
            f°[e] = warn_if_negative(f°[e] - flow)
        end
        update_support!(f°, visited_V_plus, supp_out, supp_inc)
        f′_plus_list[i] = [(π°_edge[e°], warn_if_negative(flow)) for (e°,flow) in f_list_plus]

        if i >= l+1
            break # skip second part on the last iteration
        end

        (f_list_minus,visited_V_minus) = find_flow_to_T(G°, π°, f°, supp_inc, Ti[i], ρ, scratch_space)
        for (e,flow) in f_list_minus
            f°[e] = warn_if_negative(f°[e] - flow)
        end
        update_support!(f°, visited_V_minus, supp_out, supp_inc)
        f′_minus_with_dups = [(π°_edge[e°], warn_if_negative(flow)) for (e°,flow) in f_list_minus]
        f′_minus_list[i] = [(e, flow) for (e, flow) in f′_minus_with_dups
            if !(abs(flow) < EXPRESSION_TOL && count(x -> x[1] == i, π′_edge_inv[e]) != 1)]

    end

    f′ = sum_of_subflows(l, G′.m, f′_plus_list, f′_minus_list, π′_edge_inv)
    (R_graph,C,edge_in_R_graph,s_star,t_star) = construct_R_graph(
        l, G, G′, f′, a, b, γ, adj_Gi_plus, adj_Gi_minus, Si, Ti, π′_edge_inv)
    g = find_blocking_flow(R_graph, C, s_star, t_star)

    for i in 1:l+1
        for (_,out) in adj_Gi_plus[i], e in out
            f′[e] = f′[e] + g[edge_in_R_graph[e]]
        end
    end
    for i in 1:l
        for (_,out) in adj_Gi_minus[i], e in out
            e′ = get_edge_inv(π′_edge_inv, e, i)
            f′[e′] = warn_if_negative(f′[e′] -  g[edge_in_R_graph[e′]])
        end
    end

    # clean up near-zero negative floats
    for e in 1:G′.m
        if f′[e] < VARIABLE_TOL
            f′[e] = zero(R)
        end
    end

    (G′,π′,π′_edge,f′)
end

"update supp_out and supp_inc to represent the support of f°,
by changing them in place. only checks vertices in changed"
function update_support!(
    f°::Vector{R},
    changed_V::Vector{Int},
    supp_out::Vector{Vector{Int}},
    supp_inc::Vector{Vector{Int}}) where R<:Real
    # update the support of f°
    for v in changed_V
        supp_out[v] = [e for e in supp_out[v] if f°[e] > 0]
        supp_inc[v] = [e for e in supp_inc[v] if f°[e] > 0]
    end
end

"construct the graph G′ along with π′_edge and π′_edge_inv, from the graph G and
the subgraph copies G^{i,-} represented by adj_Gi_minus"
function construct_G′(
    l::Int,
    G::Graph,
    a::Vector{R},
    b::Vector{R},
    adj_Gi_minus::Vector{Vector{Tuple{Int,Vector{Int}}}}) where R<:Real
    # construct G′
    # this assumes that none of the edges of G have been removed
    G′_edges = [(src_vertex(G, e),dst_vertex(G, e)) for e in 1:G.m]
    n′ = G.n
    π′ = collect(1:G.n)
    # again assume that no edges have been removed
    π′_edge = collect(1:G.m)
    π′_edge_inv = [Tuple{Int,Int}[] for _ in 1:G.m]
    for i in 1:l
        π′_inv = zeros(Int, G.n)
        for (v,_) in adj_Gi_minus[i]
            if a[v] > 0 || b[v] > 0
                # vertices in S or T are not copied
                π′_inv[v] = v
            else
                n′ += 1 # add new vertex
                push!(π′, v) # record what it's a copy of
                π′_inv[v] = n′ # = the new vertex
            end
        end
        for (v,out) in adj_Gi_minus[i]
            for e in out
                u = dst_vertex(G, e)
                push!(π′_edge, e)

                push!(G′_edges, (π′_inv[v],π′_inv[u]))
                e′ = length(G′_edges)
                push!(π′_edge_inv[e], (i,e′))
            end
        end
    end
    G′ = Graph(n′, G′_edges)
    (G′,π′,π′_edge,π′_edge_inv)
end

function construct_R_graph(
    l::Int,
    G::Graph,
    G′::Graph,
    f′::Vector{R},
    a::Vector{R},
    b::Vector{R},
    γ::R,
    adj_Gi_plus::Vector{Vector{Tuple{Int,Vector{Int}}}},
    adj_Gi_minus::Vector{Vector{Tuple{Int,Vector{Int}}}},
    Si::Vector{Vector{Int}},
    Ti::Vector{Vector{Int}},
    π′_edge_inv::Vector{Vector{Tuple{Int,Int}}}) where R<:Real
    # R = (V′,E_R), C : E_R → [0,∞]
    R_edges = Tuple{Int,Int}[]
    C = R[]
    # number that is larger than any necessary capacity
    inf_cap = 2*sum(a)

    # note that -1 is never a valid edge index
    edge_in_R_graph = fill(-1, G′.m)

    for i in 1:l+1
        for (v,out) in adj_Gi_plus[i]
            for e in out
                u = dst_vertex(G, e)
                push!(R_edges, (v,u))
                edge_in_R_graph[e] = length(R_edges)
                push!(C, inf_cap)
            end
        end
    end
    for i in 1:l
        for (_,out) in adj_Gi_minus[i]
            for e in out
                e′ = get_edge_inv(π′_edge_inv, e, i)

                u = dst_vertex(G′, e′)
                v = src_vertex(G′, e′)
                push!(R_edges, (u,v)) # note reversed direction
                edge_in_R_graph[e′] = length(R_edges)
                push!(C, f′[e′])
            end
        end
    end

    s_star = G′.n+1
    for s in Si[1]
        push!(R_edges, (s_star,s))
        # edge_in_R_graph[e] = 0
        c = a[s] - sum(f′[e] for e in δ_out(G′, s); init=zero(R))
        c = warn_if_negative(c)
        push!(C, c)
    end
    t_star = G′.n+2
    for t in Ti[l+1]
        push!(R_edges, (t,t_star))
        # edge_in_R_graph[e] = 0
        c = γ*b[t] - sum(f′[e] for e in δ_inc(G′, t); init=zero(R))
        c = warn_if_negative(c)
        push!(C, c)
    end

    # clean up near-zero floats (including small positive values)
    for k in 1:length(C)
        if C[k] < 0 && abs(C[k]) < VARIABLE_TOL
            C[k] = zero(R)
        end
    end

    R_graph = Graph(G′.n+2, R_edges)
    (R_graph,C,edge_in_R_graph,s_star,t_star)
end

function sum_of_subflows(
    l::Int, m::Int,
    f′_plus_list::Vector{Vector{Tuple{Int,R}}},
    f′_minus_list::Vector{Vector{Tuple{Int,R}}},
    π′_edge_inv::Vector{Vector{Tuple{Int,Int}}}) where R<:Real
    f′ = zeros(R, m)
    for i in 1:l+1
        for (e,flow) in f′_plus_list[i]
            f′[e] += flow
        end
        if i < l+1
            for (e,flow) in f′_minus_list[i]
                e′ = get_edge_inv(π′_edge_inv, e, i)
                f′[e′] += flow
            end
        end
    end
    f′
end

function general_bfs_fwd(G::Graph, S::Vector{Int}, visited::BitVector, action::Function, pred::Union{Function,Nothing}=nothing)
    Q = Queue{Int}()
    for s in S
        enqueue!(Q, s)
    end
    while !isempty(Q)
        v = dequeue!(Q)
        if !visited[v] && (pred === nothing || pred(v))
            for u in Δ_out(G, v)
                enqueue!(Q, u)
            end
            visited[v] = true
            action(v)
        end
    end
end

function general_bfs_bwd(G::Graph, T::Vector{Int}, visited::BitVector, action::Function, pred::Union{Function,Nothing}=nothing)
    Q = Queue{Int}()
    for t in T
        enqueue!(Q, t)
    end
    while !isempty(Q)
        v = dequeue!(Q)
        if !visited[v] && (pred === nothing || pred(v))
            for u in Δ_inc(G, v)
                enqueue!(Q, u)
            end
            visited[v] = true
            action(v)
        end
    end
end

"find the sets V(G^{i,+}) and V(G^{i,-}) by first finding the smallest i such that Si[i] -> v, for each v,
and then performing BFS backwards from Ti (and Si)"
function find_Gi(G::Graph, Si::Vector{Vector{Int}}, Ti::Vector{Vector{Int}}, l::Int)
    visited = falses(G.n)
    min_i = fill(-1, G.n)
    for i in 1:l+1
        # Here we do not want to reset visited between iterations
        # because we don't want to visit vertices twice
        general_bfs_fwd(G, Si[i], visited, v -> min_i[v] = i)
    end
    # min_i[v] is now the smallest i such that Si[i] -> v

    # clear visited
    fill!(visited, false)

    # BFS to enforce reachability to Ti for V(G^{i,+})
    V_Gi_plus = [Int[] for i in 1:l+1]
    handles_Gi_plus = [Int[] for v in 1:G.n]
    for i in 1:l+1
        @assert !VALIDATE_VISITED || all(!visited[v] for v in 1:G.n)
        general_bfs_bwd(G, Ti[i], visited,
            function(v)
                push!(V_Gi_plus[i], v)
                push!(handles_Gi_plus[v], i)
            end,
            v -> min_i[v] == i)
        for v in V_Gi_plus[i]
            visited[v] = false
        end
    end

    # BFS to enforce reachability to from Si[i+1] for V(G^{i,-})
    handles_from_Si_p1 = [Int[] for v in 1:G.n]
    for i in 1:l
        @assert !VALIDATE_VISITED || all(!visited[v] for v in 1:G.n)
        visited_list = Int[] # used for resetting visited after each search
        general_bfs_fwd(G, Si[i+1], visited,
            function(v)
                push!(handles_from_Si_p1[v], i)
                push!(visited_list, v)
            end,
            v -> i <= min_i[v] <= i+1)
        for v in visited_list
            visited[v] = false
        end
    end

    # BFS to enforce reachability to Ti for V(G^{i,-})
    V_Gi_minus = [Int[] for i in 1:l]
    handles_Gi_minus = [Int[] for v in 1:G.n]
    for i in 1:l
        @assert !VALIDATE_VISITED || all(!visited[v] for v in 1:G.n)
        general_bfs_bwd(G, Ti[i], visited,
            function(v)
                push!(V_Gi_minus[i], v)
                push!(handles_Gi_minus[v], i)
            end,
            v -> i in handles_from_Si_p1[v] && i <= min_i[v] <= i+1)
        for v in V_Gi_minus[i]
            visited[v] = false
        end
    end

    # TODO remove this assert (Observation D.17)
    @assert all(length(handles_Gi_plus[v]) <= 3 for v in 1:G.n) "Observation D.17 should hold"
    @assert all(length(handles_Gi_minus[v]) <= 3 for v in 1:G.n) "Observation D.17 should hold"

    # create adjacency lists from V(G^{i,+}) etc.
    adj_Gi_plus = Vector{Tuple{Int,Vector{Int}}}[]
    for i in 1:l+1
        push!(adj_Gi_plus, Tuple{Int,Vector{Int}}[])
        for v in V_Gi_plus[i]
            push!(adj_Gi_plus[i],
                (v,[e for e in δ_out(G, v) if i in handles_Gi_plus[dst_vertex(G, e)]]))
        end
    end

    adj_Gi_minus = Vector{Tuple{Int,Vector{Int}}}[]
    for i in 1:l
        push!(adj_Gi_minus, Tuple{Int,Vector{Int}}[])
        for v in V_Gi_minus[i]
            push!(adj_Gi_minus[i],
                (v,[e for e in δ_out(G, v) if i in handles_Gi_minus[dst_vertex(G, e)]]))
        end
    end

    # Simpler method of creating adj. list of creating adj. lists:
    #=
     adj_Gi_plus = [
        [(v,[e for e in δ_out(G, v) if i in handles_Gi_plus[dst_vertex(G, e)]) for v in V_Gi_plus[i]]
        [(v,[e for e in δ_out(G, v) if i in handles_Gi_plus[dst_vertex(G, e)]]) for v in V_Gi_plus[i]]
         for i=1:l+1
     ]
     adj_Gi_minus = [
        [(v,[e for e in δ_out(G, v) if i in handles_Gi_minus[dst_vertex(G, e)]) for v in V_Gi_minus[i]]
        [(v,[e for e in δ_out(G, v) if i in handles_Gi_minus[dst_vertex(G, e)]]) for v in V_Gi_minus[i]]
         for i=1:l
     ]
    =#

    (adj_Gi_plus,adj_Gi_minus)
end

"find partitions Si[i] of S s.t. the shortest alt. path has length 2(i-1)
and Ti[i] of T such that the shortest alt. path has length 2(i-1)+1 for i in [1,l+1]
G′ is a handled graph representing the shortcut graph H
Note that Si[i] corresponds to S^{i-1} in the paper"
function find_Si_and_Ti(
    L::Int,
    G′::Graph,
    a::Vector{R},
    S::Vector{Int},
    T::Vector{Int},
    f′::Vector{R})::Tuple{Vector{Vector{Int}},Vector{Vector{Int}}} where R<:Real

    distance = find_distance_from_S_and_T(L, G′, a, S, T, f′)

    # TODO use OffsetArray instead of Vector later
    Si = [Int[] for i=1:L+2]
    Ti = [Int[] for i=1:L+2]

    for s in S
        d = distance[G′.n+s]
        if d <= 2*L
            @assert d >= 0
            @assert rem(d, 2) == 0
            i = div(d, 2)
        else
            i = L + 1
        end
        push!(Si[i+1], s)
    end

    for t in T
        d = distance[t]
        if d <= 2*L+1
            @assert d >= 1
            @assert rem(d, 2) == 1
            i = div(d - 1, 2)
        else
            i = L + 1
        end
        push!(Ti[i+1], t)
    end

    (Si,Ti)
end

function find_distance_from_S_and_T(L::Int, G′::Graph, a::Vector{R}, S::Vector{Int}, T::Vector{Int}, f′::Vector{R}) where R<:Real
    N = G′.n

    # We want the edges which are common to G′ and F to have the same index
    # and according to an implementation detail in Graph, the edge indices in Graph(n, es)
    # will be the same as in enumerate(es)
    # This assumes that no edges are removed from G′
    F_edges = [(src_vertex(G′, e),dst_vertex(G′, e)) for e in 1:G′.m]
    weight = zeros(Int, G′.m)
    for e in 1:G′.m
        if f′[e] > 0
            # the edge (v,u) is in the support of f′, so add the edge (u′,v′) to F
            v = src_vertex(G′, e)
            u = dst_vertex(G′, e)
            push!(F_edges, (N+u,N+v))
            push!(weight, 0)
        end
    end
    for s in S
        # add edge from s′ to s with weight 1
        push!(F_edges, (N+s,s))
        push!(weight, 1)
    end
    for t in T
        # add edge from t to t′ with weight 1
        push!(F_edges, (t,N+t))
        push!(weight, 1)
    end
    F = Graph(2*N, F_edges)

    # construct S[1], ..., S[L+2] and T[1], ..., T[L+2]
    q0 = Queue{Int}()
    q1 = Queue{Int}()
    distance = fill(10*(L+1), F.n) # 10*(L+1) is basically infinity
    visited = falses(F.n)

    for s in S
        out_flow = sum(f′[e] for e in δ_out(G′, s); init=zero(R))
        if VARIABLE_TOL < a[s] - out_flow # s not saturated by flow f′
            enqueue!(q0, N+s)
        end
    end

    d = 0
    while true
        if !isempty(q0)
            v = dequeue!(q0)
            if !visited[v]
                distance[v] = d
                visited[v] = true
                for e in δ_out(F, v)
                    enqueue!(weight[e] == 0 ? q0 : q1, dst_vertex(F, e))
                end
            end
        elseif !isempty(q1)
            (q0,q1) = (q1,q0) # swap q0 and q1
            d += 1
        else
            break
        end
    end
    distance
end

"find S′ according to Lemma D.6"
function find_S′(
    G′::Graph,
    f′::Vector{R},
    a::Vector{R},
    b::Vector{R},
    S::Vector{Int},
    T::Vector{Int},
    γ::R,
    ϵ::R,
    ϕ::R)::Vector{Int} where R<:Real

    L = floor(Int, log(1+ϵ, 3*norm(a, 1) / ϕ))
    (Sl,Tl) = find_Si_and_Ti(L, G′, a, S, T, f′)


    aS = [sum(a[Sl[l]]) for l in 1:L+2]
    bT = [sum(b[Tl[l]]) for l in 1:L+2]

    # compute prefix resp. postfix sums
    aS_geq = reverse(cumsum(reverse(aS)))
    bT_leq = cumsum(bT)
    aS_gt = [aS_geq[i+1] for i=1:L+1]

    # find the l^* that minimises a(S^{> l}) + b(T^{<= l}) / (1+ϵ)
    l_star = argmin([aS_gt[l] + γ * bT_leq[l] / (1+ϵ) for l in 1:L+1])
    S′ = [s for i in 1:l_star for s in Sl[i]]
    S′
end

"find subflow f″ of f′ sent from S according to Lemma D.8.
Returns the flow as a list of Tuples (edge,flow).
Assumes that c is an all-zero verctor of length G.n"
function find_flow_from_S(
    G′::Graph,
    π′::Vector{Int},
    f′::Vector{R},
    supp_out::Vector{Vector{Int}},
    S′::Vector{Int},
    ρ::Vector{Int},
    c::Vector{R})::Tuple{Vector{Tuple{Int,R}},Vector{Int}} where R<:Real

    # NB. we assume that x == 0 for x in c

    # store flow as tuples (edge,flow)
    f″_list = Tuple{Int,R}[]
    # store as tuples (priority,value)
    H = BinaryHeap{Tuple{Int,Int}}(Base.By(first))
    # containts the vertices which we visited
    E″ = Int[]

    for s in S′
        for e in supp_out[s] # δ_out(supp(f′), s)
            c[s] += f′[e]
        end
        if c[s] > 0
            push!(H, (ρ[π′[s]],s))
        end
    end

    while !isempty(H)
        (_,v) = pop!(H)
        push!(E″, v)
        a = c[v]
        for e in supp_out[v] # δ_out(supp(f′), v)
            a′ = min(a, f′[e])
            push!(f″_list, (e, a′))
            u = dst_vertex(G′, e)
            if c[u] <= 0 # TODO more robust than c[u] == 0?
                push!(H, (ρ[π′[u]],u))
            end
            c[u] += a′
            a -= a′
            if a <= 0
                break
            end
        end
    end

    for v in E″
        c[v] = zero(R)
    end

    (f″_list,E″)
end

"find subflow f″ of f′ received by T according to Lemma D.9.
Returns the flow as a list of Tuples (edge,flow).
Assumes that c is an all-zero verctor of length G.n"
function find_flow_to_T(
    G′::Graph,
    π′::Vector{Int},
    f′::Vector{R},
    supp_inc::Vector{Vector{Int}},
    T′::Vector{Int},
    ρ::Vector{Int},
    c::Vector{R})::Tuple{Vector{Tuple{Int,R}},Vector{Int}} where R<:Real

    # NB. we assume that x == 0 for x in c

    # store flow as tuples (edge,flow)
    f″_list = Tuple{Int,R}[]
    # store as tuples (priority,value)
    H = BinaryHeap{Tuple{Int,Int}}(Base.By(first, Base.Reverse))
    E″ = Int[]

    for t in T′
        for e in supp_inc[t] # δ_inc(supp(f′), s)
            c[t] += f′[e]
        end
        if c[t] > 0
            push!(H, (ρ[π′[t]],t))
        end
    end

    while !isempty(H)
        (_,v) = pop!(H)
        push!(E″, v)
        a = c[v]
        for e in supp_inc[v] # δ_inc(supp(f′), v)
            a′ = min(a, f′[e])
            push!(f″_list, (e, a′))
            u = src_vertex(G′, e)
            if c[u] == 0
                push!(H, (ρ[π′[u]],u))
            end
            c[u] += a′
            a -= a′
            if a <= 0
                break
            end
        end
    end

    for v in E″
        c[v] = zero(R)
    end

    (f″_list,E″)
end
