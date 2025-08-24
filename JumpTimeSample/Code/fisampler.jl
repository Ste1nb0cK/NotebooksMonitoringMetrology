include("envsetup.jl")
include("constants.jl")
include("functions.jl")
addprocs(cfg["nworkers"]; exeflags="--project=" * projectpath) # Call workers with this environment

println(stderr, "Using the configuration:")
println(stderr, cfg)
# # Jump operators and their derivatives
using Base.Threads
@everywhere using BackAction

directoryname = "Sample_$(real(cfield))+i$(imag(cfield))"
try
        mkdir(directoryname)
catch
        @warn "Directory already existed"
end

const deltarange = collect(LinRange(0.0, 1.0, cfg["ndeltas"]))

metadata = Dict(
        "NCHANNELS0" => NCHANNELS0,
        "NLEVELS" => NLEVELS,
        "OMEGA" => OMEGA,
        "GAMMA" => GAMMA,
        "tf" => tf,
        "tspan" => string(tspan),
        "ntimes" => cfg["ntimes"],
        "ntraj" => cfg["ntraj"],
        "nworkers" => cfg["nworkers"],
        "nchannels" => nchannels,
        "nops" => nops,
        "psi0" => string(psi0),
        "u0" => string(u0),
        "cfield" => string(cfield),
        "tlist_first" => first(tlist),
        "tlist_last" => last(tlist),
)
# Write to TOML file
open(directoryname * "/metadata_$(real(cfield))+i$(imag(cfield)).toml", "w") do io
        TOML.print(io, metadata)
end

L0 = d -> sqrt(GAMMA) * util.sigma_m
H0 = d -> d * [[0, 0] [0, 1.0 + 0im]] + 0.5 * OMEGA * util.sigma_x

filenametrsample = directoryname * "/trsample_$(real(cfield))+i$(imag(cfield)).bin"
open(filenametrsample, "w+") do io
        # Extend file to required size (m*n*element_size)
        seek(io, cfg["ndeltas"] * cfg["ntraj"] * sizeof(Float64))
        write(io, UInt8(0))
end

trsample = open(filenametrsample, "r+") do io
        Mmap.mmap(io, Matrix{Float64}, (cfg["ntraj"], cfg["ndeltas"]))
end


Ls_par, H_par, He_par = mo.obtain_parametric_unraveling_operators(L0, H0, T, [cfield], NLEVELS);
cache_aux1 = Vector{ComplexF64}(undef, NLEVELS)
cache_aux2 = Vector{ComplexF64}(undef, NLEVELS)
cache_state = copy(psi0)
psi = Vector{ComplexF64}(undef, NLEVELS)
cache_phi = zeros(ComplexF64, NLEVELS)
# Set the necessary operators
cache_exp = Matrix{ComplexF64}(undef, NLEVELS, NLEVELS)
cache_dexp = Matrix{ComplexF64}(undef, NLEVELS, NLEVELS)

nthreads = Threads.nthreads()

# allocate per thread
nthreads = Threads.nthreads()
psi_pool = [similar(psi0) for _ in 1:nthreads]
cache_state_pool = [copy(psi0) for _ in 1:nthreads]
cache_phi_pool = [zeros(ComplexF64, NLEVELS) for _ in 1:nthreads]
cache_exp_pool = [Matrix{ComplexF64}(undef, NLEVELS, NLEVELS) for _ in 1:nthreads]
cache_dexp_pool = [Matrix{ComplexF64}(undef, NLEVELS, NLEVELS) for _ in 1:nthreads]
cache_aux1_pool = [similar(psi0) for _ in 1:nthreads]
cache_aux2_pool = [similar(psi0) for _ in 1:nthreads]

for i in 1:cfg["ndeltas"]
        delta = deltarange[i]
        L = Ls_par[1](delta)
        dL = ForwardDiff.derivative(Ls_par[1], delta)

        sim = obtain_ensol(L0, H0, delta, T, [cfield], params, tspan, e_ops, tlist)
        @threads for k in 1:cfg["ntraj"]
                tid = Threads.threadid()
                psiphi_finaltime!(sim[k], L, dL, cache_exp_pool[tid], cache_dexp_pool[tid],
                        cache_aux1_pool[tid], cache_aux2_pool[tid], psi_pool[tid],
                        cache_state_pool[tid], cache_phi_pool[tid], delta)
                trsample[k, i] = 2 * real(dot(psi_pool[tid], cache_phi_pool[tid]))
                copyto!(cache_state_pool[tid], psi0)
                cache_phi_pool[tid] .= 0.0 + 0.0im
        end
end

Mmap.sync!(trsample)
finalize(trsample)
