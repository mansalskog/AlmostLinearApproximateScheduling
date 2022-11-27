using LinkCutTrees
using LinkCutTrees: parent

"find a blocking flow using the algorithm from Sleator and Tarjan, p. 388,
in the graph G from the source s_star and the sink t_star"
function find_blocking_flow(G::Graph, capacity::Vector{R}, s_star::Int, t_star::Int)::Vector{R} where R<:Real
    @assert 1 <= s_star <= G.n "s_star must be a vertex"
    @assert 1 <= t_star <= G.n "t_star must be a vertex"
    trees = [make_tree(Int, Int, R, v) for v in 1:G.n]
    m = length(capacity)
    unused_cap = zeros(R, m) # unused (residual) capacity
    s = trees[s_star]
    t = trees[t_star]
    while true
        while true
            # step 1
            v = find_root(s)
            if v === t break end # go to step 4
            # step 2
            if !isempty(δ_out(G, label(v)))
                e = first(δ_out(G, label(v)))
                link!(v, trees[dst_vertex(G, e)], e, capacity[e])
                continue # go to step 1
            end
            # step 3
            if v === s
                for w_idx in 1:G.n # compute unused capacity
                    w = trees[w_idx]
                    u = parent(w)
                    for e in δ_out(G, w_idx)
                        unused_cap[e] += trees[dst_vertex(G, e)] === u ? cost(w) : capacity[e]
                    end
                end
                return [capacity[e] - unused_cap[e] for e in 1:m] # and stop
            else
                for e in collect(δ_inc(G, label(v)))
                    # record the unused capacity and cut incoming edges
                    u = trees[src_vertex(G, e)]
                    if parent(u) === v
                        unused_cap[e] += cost(u)
                        cut!(u)
                    else
                        unused_cap[e] += capacity[e]
                    end
                    delete_edge!(G, e) # and delete e from the graph
                end
                continue # goto step 1
            end
        end
        # step 4
        (v,c) = find_mincost(s)
        add_cost!(s, -c)
        while true
            # step 5
            (v,c) = find_mincost(s)
            if c > 0
                break # go to step 1
            end
            if parent(v) === nothing
                # v is root, so we cannot remove (v,parent(v)), nor perform cut!(v) 
                break
            end
            delete_edge!(G, edge_label(v))
            cut!(v)
            # record unused capacity of 0 and go to step 5
        end
    end
end