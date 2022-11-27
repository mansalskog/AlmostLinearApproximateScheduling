using JuMP

include("../src/graph.jl")
include("../src/data.jl")
include("../src/solve_instance.jl")

# filename extension of the exported files,
# which also determines the export format,
# in this case, gzipped MPS files
FILE_EXT = ".mps.gz"

function construct_ip_model(
    inst::Instance,
    q::Vector{Int},
    ϵ::R,
    int_constr::Bool) where R<:Real

    τ = [zero(R)]
    η = R[]

    # the range of d is [1,D+1] and not [0,D] as in the paper
    d = 1
    while τ[d] < sum(inst.p)
        d += 1
        push!(τ, (one(ϵ)+ϵ)^(d-1)) # to τ[d]
        push!(η, τ[d] - τ[d-1]) # to η[d-1]
    end
    D = d-1
    # here we know that τ[D+1] >= sum(inst.p)
    println("Using ", D+1, " time steps from ", τ[1], " to ", τ[D+1])

    # TODO should be changed later to
    # handle non-polynomial processing times
    # also the weights should be changed
    d_min = ones(Int, inst.n)
    d_max = fill(D+1, inst.n)

    model = JuMP.Model()
    J = 1:inst.n # set of all jobs
    D = length(τ)-1
    @variable(model, x[J,1:D+1])

    if int_constr
        set_integer.(x)
    end

    # (2)
    @objective(model,
               Min,
               sum(inst.w) * τ[D+1] - sum(inst.w[j] * sum(η[2:D] .* x[j,2:D]) for j in J))
    # (3)
    @constraint(model,
                [j in J, d in 1:D],
                x[j,d] <= x[j,d+1])
    # (4)
    @constraint(model,
                [(j1,j2) in inst.prec, d in 1:D+1],
                x[j1,d] >= x[j2,d])
    # (5)
    @constraint(model,
                [d in 2:D+1],
                sum(inst.p .* x[:,d]) <= inst.m * τ[d])
    # (6)
    for j in J, d in 1:D+1
        if d <= d_min[j] || τ[d] < q[j]
            @constraint(model, x[j,d] == 0)
        end
    end

    # (7)
    @constraint(model,
                [j in J, d in d_max[j]:D+1],
                x[j,d] == 1)

    model
end

function export_instance(inst::Instance, ϵ::Real, int_constr::Bool)
    prec_graph = Graph(inst.n, inst.prec)
    prec_order = topological_ordering(prec_graph)

    q = longest_prec_chains(inst, prec_graph, prec_order)
    model = construct_ip_model(inst, q, ϵ, int_constr)

    modelfile = joinpath(out_dir, string(inst.name, FILE_EXT))
    println("writing model to file ", modelfile)
    write_to_file(model, modelfile)
end

inst_dir = ARGS[1]
ϵ = parse(Float64, ARGS[2])
out_dir = ARGS[3]
int_constr = parse(Bool, ARGS[4])

for inst_path in readdir(inst_dir; join=true)
    inst = parse_instance(inst_path)
    export_instance(inst, ϵ, int_constr)
end
