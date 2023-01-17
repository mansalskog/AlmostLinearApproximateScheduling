using Test: @testset, @test

using LinkCutTrees

include("../src/graph.jl")
include("../src/network_flow.jl")
include("../src/blocking_flow.jl")
include("../src/tolerances.jl")

@testset "blocking flow tests" begin
    # test the blocking flow algorithm
    n = 6
    edges = [(1,2), (1,3), (2,3), (2,4), (2,5), (3,5), (5,4), (4,6), (5,6)]
    capacities = [10, 10, 2, 4, 8, 9, 6, 10, 10]
    G = Graph{Set{Int}}(n, edges)
    s = 1
    t = 6
    flow = find_blocking_flow(G, capacities, s, t)
    best_flow = [10, 9, 0, 3, 7, 9, 6, 9, 10]
    for e in 1:length(flow)
        @test flow[e] == best_flow[e]
    end
end

@testset "Find Si and Ti tests" begin
    es = [
        (3,2),
        (4,3), (4,6), (4,7),
        (5,4), (5,6),
        (6,1), (6,2),
        (7,2),
        (8,7),
    ]
    G = Graph{Vector{Int}}(8, es)
    a = [0,0,0,0,5,0,0,5,0,0]
    b = [5,5,0,0,0,0,0,0,0,0]
    S = [v for v in 1:G.n if a[v] > 0]
    @test S == [5,8]
    T = [v for v in 1:G.n if b[v] > 0]
    @test T == [1,2]
    # not really a valid flow, but works for testing
    f = [0,0,0,0,0,0,0,0,1,5]
    # => H has edges (5,1), (5,2), (8,2) and (2,8)
    π = collect(1:G.n) # empty pi
    l = 1
    (Si,Ti) = find_Si_and_Ti(l, G, a, S, T, f)
    @test Si[1] == [5]
    @test Si[2] == [8]
    @test Ti[1] == [1,2]
    @test Ti[2] == []
    (adj_Gi_plus,adj_Gi_minus) = find_Gi(G, Si, Ti, l)
    @test sort(map(first, adj_Gi_plus[1])) == [1,2,3,4,5,6,7]
    @test sort(map(first, adj_Gi_minus[1])) == [2,7,8]
end
