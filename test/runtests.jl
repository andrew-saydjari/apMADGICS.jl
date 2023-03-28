using apMADGICS #Why isn't LowRankOps loaded here?
using Test, Random

using FITSIO, HDF5, LowRankOps, EllipsisNotation, ShiftedArrays, Interpolations, SparseArrays, ParallelDataTransfer
using StatsBase, LinearAlgebra, ProgressMeter, Serialization

src_dir = "../"
include(src_dir*"src/utils.jl")
include(src_dir*"src/gridSearch.jl")
include(src_dir*"src/componentAndPosteriors.jl")
include(src_dir*"src/fileNameHandling.jl")
include(src_dir*"src/ingest.jl")
include(src_dir*"src/lowRankPrescription.jl")
include(src_dir*"src/marginalizeEW.jl")
include(src_dir*"src/spectraInterpolation.jl")
include(src_dir*"src/chi2Wrappers.jl")

include("gridSearch.jl")
include("spectraInterpolation.jl")
include("gaussTests.jl")
