module apMADGICS

    using LowRankOps, EllipsisNotation, ShiftedArrays

    export cache_skyname

    include("utils.jl")
    include("gridSearch.jl")
    include("componentAndPosteriors.jl")
    include("fileNameHandling.jl")
    include("ingest.jl")
    include("lowRankPrescription.jl")
    include("marginalizeEW.jl")
    include("spectraInterpolation.jl")
    include("chi2Wrappers.jl")

end
