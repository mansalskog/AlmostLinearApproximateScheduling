using Test: @testset, @test

include("../src/graph.jl")

@testset "graph tests" begin
    es = [
        (3,2),
        (4,3), (4,6), (4,7),
        (5,4), (5,6),
        (6,1), (6,2),
        (7,2),
    ]
    dag = Graph(7, es)
    @test dag.n == 7
    
    # check that ordering of edges is preserved (in some sense)
    @test collect(enumerate(es)) == sort([(e,(src_vertex(dag, e),dst_vertex(dag, e))) for e in edges(dag)], by=first)

    @test sort(Δ_out(dag, 4)) == [3, 6, 7]
    @test sort(collect(δ_out(dag, 4))) == [2, 3, 4]
    @test sort(reach_from(dag, [4])) == [1, 2, 3, 4, 6, 7]

    @test sort(Δ_inc(dag, 2)) == [3, 6, 7]
    @test sort(collect(δ_inc(dag, 2))) == [1, 8, 9]
    @test sort(reach_to(dag, [2])) == [2, 3, 4, 5, 6, 7]

    rank = ordering_to_rank(topological_ordering(dag))
    @test all(rank[v1] < rank[v2] for (v1,v2) in es)

    # remove the edge (4,6)
    delete_edge!(dag, 3)

    @test sort(Δ_out(dag, 4)) == [3, 7]
    @test sort(collect(δ_out(dag, 4))) == [2, 4]
    @test sort(reach_from(dag, [4])) == [2, 3, 4, 7]

    @test sort(Δ_inc(dag, 6)) == [5]
    @test sort(collect(δ_inc(dag, 6))) == [6]
    @test sort(reach_to(dag, [6])) == [5, 6]
end
