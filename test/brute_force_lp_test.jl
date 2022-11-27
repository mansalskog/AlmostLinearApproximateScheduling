using Test: @testset, @test

include("../src/brute_force_lp.jl")
include("../src/tolerances.jl")

@testset "simple LP solver tests" begin
    A = Float64[1 1 ; 1.5 1.5]
    b = Float64[2 ; 2]
    c = Float64[3 , 3]
    (x,obj) = brute_force_lp(A, b, c)
    @test x == [4 / 3 , 0]
    @test obj <= 4

    A = Float64[3 2 1 ; 2 5 3]
    b = Float64[10 ; 15]
    c = Float64[2 , 3 , 4]
    (x,obj) = brute_force_lp(A, b, c)
    @test x == [0 , 0 , 5]
    @test obj <= 20

    A = Float64[4 2 ; 2 3]
    b = Float64[32 ; 24]
    c = Float64[5 , 4]
    (x,obj) = brute_force_lp(A, b, c)
    @test x == [6 , 4]
    @test obj == 46

    A = Float64[-2 3 ; -1 3]
    b = Float64[9 ; 12]
    c = Float64[6 , 3]
    (x,obj) = brute_force_lp(A, b, c)
    @test x == [3 , 5]
    @test obj == 33

    # have to use rational here because rounding issues
    A = Rational[1 2 3 ; 5 4 6]
    b = Rational[10 ; 12]
    c = Rational[9 , 8 , 7]
    (x,obj) = brute_force_lp(A, b, c)
    @test x == [0 , 3 , 0]
    @test obj == 24
end
