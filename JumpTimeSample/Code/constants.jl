cfg = TOML.parsefile("config.toml")
const NCHANNELS0 = cfg["NCHANNELS0"]
const NLEVELS = cfg["NLEVELS"]
const OMEGA = cfg["OMEGA"]
const GAMMA = cfg["GAMMA"]


const tf = 2 * pi / OMEGA * cfg["tffactor"]
const tspan = (0.0, tf)
const tlist = collect(LinRange(tspan[1], tspan[2], cfg["ntimes"]))

const psi0 = [ComplexF64.(cfg["psi0g_re"], cfg["psi0g_im"]),
        ComplexF64.(cfg["psi0e_re"], cfg["psi0e_im"])]



const params = cs.SimulParameters(
        psi0,
        tf,
        cfg["seed"],
        cfg["ntraj"],
        # Stuff for Gillipsie, unimportant
        1000,
        0.01,
        0.1
)

const T = reshape([1.0 + 0.0im], 1, 1)
const e_ops = [util.sigma_x, util.sigma_y, util.sigma_z]
# const e_ops = []
const nops = length(e_ops)
const u0 = [real(dot(params.psi0, e_ops[1], params.psi0)),
        real(dot(params.psi0, e_ops[2], params.psi0)),
        real(dot(params.psi0, e_ops[3], params.psi0))]
#
# const u0 = [0.0, 0.0, 1.0]
const nchannels = size(T)[1]
# const ndeltas = cfg["ndeltas"]
const nalphas = cfg["nalphas"]

# const delta_range = collect(LinRange(cfg["deltamin"], cfg["deltamax"], ndeltas))
# const alphanorm_range = collect(LinRange(cfg["alphamin"], cfg["alphamax"], nalphas))
# const alphaphase_range = collect(LinRange(0.0, 2 * pi, nalphas))

const nworkers = try
        parse(Int, get(ENV, "SLURM_CPUS_PER_TASK", ""))
catch
        cfg["nworkers"]
end

const pauli_basis = [I(2), util.sigma_x, util.sigma_y, util.sigma_z]
const cfield = 0.0 + 0.0im
const alpha = [cfield]


