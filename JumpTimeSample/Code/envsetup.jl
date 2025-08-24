# projectpath = "./"
projectpath = "../"

if haskey(ENV, "SLURM_JOB_ID")
        # Running inside a Slurm allocation
        project_path = "./"
else
        # Running locally
        project_path = "../"
end

librarypath = "../BackAction.jl"
import Pkg
Pkg.activate(projectpath)
Pkg.develop(path=librarypath)
Pkg.add(["DifferentialEquations", "LinearAlgebra",
        "Distributed", "DelimitedFiles", "TOML", "JLD2", "SharedArrays", "Mmap", "ForwardDiff"])
# Base.Threads
Pkg.instantiate()
Pkg.precompile()
Pkg.resolve()
using DifferentialEquations, LinearAlgebra, Distributed, DelimitedFiles, TOML, JLD2, SharedArrays, Mmap, ForwardDiff
import BackAction.Utilities as util
import BackAction.MonteCarloWaveFunction as mc
import BackAction.CoreStructs as cs
using ForwardDiff, Statistics
import BackAction.MonitoringOperator as mo



