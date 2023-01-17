# TODO: using DataStructures: RBTree

# adjacency list graph representation
mutable struct Graph{AdjT}
    out::Vector{AdjT}
    inc::Vector{AdjT}
    # TODO
    # out::Vector{RBTree{Int}} # out[v] = [e1,e2,...]
    # inc::Vector{RBTree{Int}} # inc[v] = [e1,e2,...]
    src::Vector{Int} # src[e] = v
    dst::Vector{Int} # dst[e] = v
    m::Int # number of edges

    "create a graph from a list of edges on the form (v1,v2)"
    function Graph{AdjT}(n::Int, edges::Vector{Tuple{Int,Int}})::Graph{AdjT} where {AdjT}
        # TODO
        # inc = [RBTree{Int}() for v in 1:n]
        # out = [RBTree{Int}() for v in 1:n]
        inc = [AdjT() for v in 1:n]
        out = [AdjT() for v in 1:n]
        m = length(edges)
        src = zeros(Int, m)
        dst = zeros(Int, m)
        for (e,(v,u)) in enumerate(edges)
            src[e] = v
            dst[e] = u
            push!(out[v], e)
            push!(inc[u], e)
        end
        new(out, inc, src, dst, m)
    end
end

#=
# hack to implement iteration for RBTrees
# TODO change to SplayTree and use tree traversal iterator instead
function Base.iterate(t::RBTree{T})::Union{Nothing,Tuple{T,Int}} where T
    if length(t) == 0
        nothing
    else
        (t[1],1)
    end
end

# hack to implement iteration for RBTrees
function Base.iterate(t::RBTree{T}, i::Int)::Union{Nothing,Tuple{T,Int}} where T
    if length(t) <= i
        nothing
    else
        (t[i+1],i+1)
    end
end
=#

function Base.getproperty(graph::Graph, s::Symbol)
    if s == :n
        return length(getfield(graph, :out)) # number of vertices
    else
        return getfield(graph, s)
    end
end

"get edges as a list of edge indices"
function edges(graph::Graph)::Vector{Int}
    [e for v in 1:graph.n for e in δ_out(graph, v)]
end

"neighbours by outgoing edges"
function Δ_out(graph::Graph, v::Int)::Vector{Int}
    @assert 1 <= v <= graph.n "v must be a valid vertex"
    [graph.dst[e] for e in graph.out[v]]
end

"neighbours by incoming edges"
function Δ_inc(graph::Graph, v::Int)::Vector{Int}
    @assert 1 <= v <= graph.n "v must be a valid vertex"
    [graph.src[e] for e in graph.inc[v]]
end

"outgoing edges from a vertex. this is a O(n) operation"
function δ_out(graph::Graph, v::Int)
    @assert 1 <= v <= graph.n "v must be a valid vertex"
    graph.out[v]
end

"incoming edges to a vertex. this is a O(n) operation"
function δ_inc(graph::Graph, v::Int)
    @assert 1 <= v <= graph.n "v must be a valid vertex"
    graph.inc[v]
end

"find s for an edge e with endpoints (s,d)"
function src_vertex(graph::Graph, e::Int)::Int
    @assert 1 <= e <= length(graph.src) "e must be a valid edge"
    return graph.src[e]
end

"find d for an edge e with endpoints (s,d)"
function dst_vertex(graph::Graph, e::Int)::Int
    @assert 1 <= e <= length(graph.dst) "e must be a valid edge"
    return graph.dst[e]
end

"remove the edge with index e, by only removing e from the adjacency lists of its
endpoints. by storing the adj. lists as balanced binary trees this is done in O(log n)
we assume that the edge e exists and has not yet been removed"
function delete_edge!(graph::Graph, e::Int)
    # TODO remove asserts because of runtime
    v = graph.src[e]
    @assert e in graph.out[v] "e must be a valid out edge"
    delete!(graph.out[v], e)
    u = graph.dst[e]
    @assert e in graph.inc[u] "e must be a valid inc edge"
    delete!(graph.inc[u], e)
    graph.m -= 1

    # set src and dst to invalid values, to prevent using the edge accidentally
    graph.src[e] = -1
    graph.dst[e] = -1
end

"perform DFS by visiting v if mark[v] is false and p(v) is true (or p is nothing). Then we add v to reach and recurse
on the result of the neighborhood function neigh(graph, v)"
function visit!(graph::Graph, p::Union{Function,Nothing}, v::Int, mark::BitVector, reach::Vector{Int}, neigh::Function)
    if mark[v] || (p !== nothing && !p(v))
        return
    end
    mark[v] = true
    push!(reach, v)
    for u in neigh(graph, v)
        visit!(graph, p, u, mark, reach, neigh)
    end
end

"find all vertices reachable from the set S. the optional predicate p
restricts the search to only vertices v such that p(v) is true"
function reach_from(graph::Graph, S::Vector{Int}, p::Union{Function,Nothing}=nothing)::Vector{Int}
    reach = Int[]
    mark = falses(graph.n)
    for s in S
        visit!(graph, p, s, mark, reach, Δ_out)
    end
    reach
end

"find all vertices that the set T is reachable from. the optional predicate p
restricts the search to only vertices v such that p(v) is true"
function reach_to(graph::Graph, T::Vector{Int}, p::Union{Function,Nothing}=nothing)::Vector{Int}
    reach = Int[]
    mark = falses(graph.n)
    for t in T
        visit!(graph, p, t, mark, reach, Δ_inc)
    end
    reach
end

"find a topological ordering of the graph as a sorted list of vertices"
function topological_ordering(graph::Graph)::Vector{Int}
    order = Int[]
    mark = falses(graph.n)
    temp_mark = falses(graph.n)
    function visit_t_o(v)
        if mark[v]
            return
        end
        @assert !temp_mark[v] "graph is not a DAG"
        temp_mark[v] = true
        for u in Δ_out(graph, v)
            visit_t_o(u)
        end
        temp_mark[v] = false
        mark[v] = true
        # we could use pushfirst! here but push! and reverse seems faster
        push!(order, v)
    end
    for v in 1:graph.n
        visit_t_o(v)
    end
    reverse(order)
end

"convert an ordering to a rank function, i.e. if order[i] = v then rank[v] = i"
function ordering_to_rank(order::Vector{Int})::Vector{Int}
    rank = zeros(Int, length(order))
    for i in 1:length(order)
        rank[order[i]] = i
    end
    rank
end
