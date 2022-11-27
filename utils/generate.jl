import Random
import JSON
using Printf: @sprintf

include("../src/data.jl")

# list of all values of n to generate instances for
SIZES = [2^k for k in 3:10]

function generate_instance(n, κ, w_domain, p_domain, m, name)
    job_order = Random.randperm(n)

    max_κ = div(n*(n-1), 2)
    has_edge = Random.shuffle([trues(κ) ; falses(max_κ-κ)])
    @assert length(has_edge) == max_κ
    @assert sum(has_edge) == κ

    edges = NTuple{2,Int}[]
    k = 1
    for i in 1:n, j in i+1:n
        if has_edge[k]
            push!(edges, (job_order[i],job_order[j]))
        end
        k += 1
    end
    @assert length(edges) == κ

    w = Random.rand(w_domain, n)
    p = Random.rand(p_domain, n)

    return Instance(p, w, edges, m, name)
end

function generate_instance_set(seed::Int, out_dir::String)
    Random.seed!(seed)

    for type in ["L","M","S"]
        for n in SIZES
            if type == "L"
                α = 1 # later Random.rand(Float64, (0.10,0.75))
                κ = floor(Int, α * n^1.5)
            elseif type == "M"
                α = 1 # later Random.rand(Float64, (0.8,5.0))
                κ = α * n
            elseif type == "S"
                α = 1 # later Random.rand(Float64, (0.5,5.0))
                κ = floor(Int, α * sqrt(n))
            end
            # TODO later try with larger powers of n, when numerical issues are solved
            w_bounds = 0:n
            p_bounds = 1:n
            m = 1
            println("generating instances with ", n, " jobs and ", κ, " precedence constraints")
            for i in 1:4
                name = @sprintf("%s_%05i_%i_%02i", type, n, κ, i)
                filename = joinpath(out_dir, "$(name).json")
                inst = generate_instance(n, κ, w_bounds, p_bounds, m, name)
                open(filename, "w") do io
                    JSON.print(io, inst)
                end
                println("generated ", filename)
            end
        end
    end
end

seed = parse(Int, ARGS[1])
out_dir = ARGS[2]
generate_instance_set(seed, out_dir)
