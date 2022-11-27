using Printf: @sprintf
using DoubleFloats

include("graph.jl")
include("data.jl")
include("solve_instance.jl")
include("mwu_lp.jl")
include("brute_force_lp.jl")
include("network_flow.jl")
include("blocking_flow.jl")
include("tolerances.jl")

"print results fo a solution to an instance to the IO,
either a file or stdout"
function print_results(io::IO, inst::Instance, stats)
    println(io, "Instance name: ", inst.name)
    println(io, "Total run time: ", stats.time)
    println(io, "Total memory used: ", stats.bytes)
    sol = stats.value
    println(io, "Solution WCT: ", sol.wct)
    println(io, "Fractional WCT: ", @sprintf("%.10e", sol.wct_frac))
end

inst_path = ARGS[1]
ϵ = parse(Double64, ARGS[2])
inst = parse_instance(inst_path)
stats = @timed solve_instance(inst, ϵ)
println("solved instance in ", stats.time, " seconds using ", stats.bytes, " B of memory ")

print_results(stdout, inst, stats)

# pretty ugly way to change the file extension from .json to .sol.txt
out_path = join(vcat(split(inst_path, ".")[1:end-1], "sol", "txt"), ".")
println("Writing results to ", out_path)
if isfile(out_path)
    error("Output file already exists, will not overwrite!")
end
open(out_path, "w") do io
    print_results(io, inst, stats)
end
