module apMADGICS

    using LowRankOps, EllipsisNotation, ShiftedArrays

    export cache_skyname, getSky4visit, cache_starname, stack_out, nanmedian, LowRankMultMatIP
    export sampler_1d_hierarchy_var, update_Ctotinv_Vstarstarlines, deblend_components_all
    export sampler_2d_hierarchy_var, shift_trim_range, sample_chi2_flux_dflux, marginalize_flux_err
    export update_Ctotinv_Vdib, deblend_components_all_asym, nanify

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
