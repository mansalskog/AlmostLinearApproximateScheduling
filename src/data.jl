import JSON

### Instance ###

struct Instance
    p::Vector{Int} # processing times
    w::Vector{Int} # weights
    prec::Vector{NTuple{2,Int}} # contains indices of tasks
    m::Int # number of machines
    name::String
end


function Base.getproperty(inst::Instance, s::Symbol)
    if s == :n
        return length(getfield(inst, :p)) # number of jobs
    else
        return getfield(inst, s)
    end
end

### Solution and Result ###

struct Solution{R <: Real}
    epsilon::R
    C_frac::Vector{R}
    wct_frac::R
    C::Vector{Int}
    wct::Int
end

### JSON parsing ###

function parse_instance(filename::String)::Instance
    data = JSON.parsefile(filename)
    # handle duplicate precedence constraint, i.e.,
    # [(1,2),(2,3),(1,2)] -> [(1,2),(2,3)]
    prec = collect(Set((js[1],js[2]) for js in data["prec"]))
    # handle missing instance name
    name = !("name" in keys(data)) || data["name"] === nothing ? "no_instance_name" : data["name"]
    Instance(data["p"], data["w"], prec, data["m"], name)
end
