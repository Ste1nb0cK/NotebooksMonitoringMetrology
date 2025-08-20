include("envsetup.jl")
include("constants.jl")
include("functions.jl")
addprocs(cfg["nworkers"]; exeflags="--project=" * projectpath) # Call workers with this environment

println(stderr, "Using the configuration:")
println(stderr, cfg)
# # Jump operators and their derivatives

L0 = d -> sqrt(GAMMA) * util.sigma_m
H0 = d -> d * [[0, 0] [0, 1.0 + 0im]] + 0.5 * OMEGA * util.sigma_x
alpha = [0.0 + 0.0im]

@everywhere using BackAction

# @everywhere const tlist = $tlist
println(stderr, "Begin sampling")
@timeerr sim = obtain_ensol(L0, H0, cfg["DELTA"], T, alpha, params, tspan, e_ops, tlist);
println(stderr, "Sampling finished")
using Base.Threads
filename = "currentsample.bin"
ntimes = cfg["ntimes"] - 1

open(filename, "w+") do io
        # Extend file to required size (m*n*element_size)
        seek(io, ntimes * cfg["ntraj"] * sizeof(Float64))
        write(io, UInt8(0))
end

currentsample = open(filename, "r+") do io
        Mmap.mmap(io, Matrix{Float64}, (ntimes, cfg["ntraj"]))
end


# Extract currents
# currentsample = SharedArray{Float64}(ntimes, cfg["ntraj"])
# the idea is to count how many jumps there are in each time bin
# @sync @distributed for trajectoryindex in 1:cfg["ntraj"]
println(stderr, "Writing data")
@timeerr @threads for trajectoryindex in 1:cfg["ntraj"]
        jumptimes = sim[trajectoryindex].prob.kwargs[:callback].continuous_callbacks[1].affect!.jump_times
        njumps = length(jumptimes)
        if njumps == 0
                continue
        end
        jumpindex = 1
        timebinindex = 1

        for jt in jumptimes
                bin = searchsortedlast(tlist, jt)
                if bin <= ntimes
                        currentsample[bin, trajectoryindex] += 1.0
                end

        end
end
println(stderr, "Writing finished, flush and finalize")
Mmap.sync!(currentsample)
finalize(currentsample)

