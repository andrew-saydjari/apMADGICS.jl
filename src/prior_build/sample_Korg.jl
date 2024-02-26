## This is script generates a bunch of high resolution model spectra using Korg. These are samples for the starLines prior when on the first pass using pure theory models.
# Author - Andrew Saydjari, CfA

import Pkg; using Dates; t0 = now(); t_then = t0;
using InteractiveUtils; versioninfo()
Pkg.activate("../../"); Pkg.instantiate(); Pkg.precompile()
t_now = now(); dt = Dates.canonicalize(Dates.CompoundPeriod(t_now-t_then)); println("Package activation took $dt"); t_then = t_now; flush(stdout)
using BLISBLAS
using Distributed, SlurmClusterManager, Suppressor, DataFrames
addprocs(SlurmManager(),exeflags=["--project=../../"])
t_now = now(); dt = Dates.canonicalize(Dates.CompoundPeriod(t_now-t_then)); println("Worker allocation took $dt"); t_then = t_now; flush(stdout)

@everywhere begin
    using BLISBLAS
    using LinearAlgebra
    BLAS.set_num_threads(1)
    using FITSIO, Serialization, HDF5, LowRankOps, EllipsisNotation, ShiftedArrays
    using Interpolations, SparseArrays, ParallelDataTransfer, AstroTime, Suppressor
    using ThreadPinning

    prior_dir = "/uufs/chpc.utah.edu/common/home/u6039752/scratch1/working/"
    prior_dir0 = "/uufs/chpc.utah.edu/common/home/u6039752/scratch/working/"
    src_dir = "../../"
    include(src_dir*"src/utils.jl")
    include(src_dir*"src/prior_build/prior_utils.jl")
    
    using StatsBase, ProgressMeter
    using SortFilters, BasisFunctions, Random, DustExtinction, DelimitedFiles
    using Korg, Distributions
end
t_now = now(); dt = Dates.canonicalize(Dates.CompoundPeriod(t_now-t_then)); println("Worker loading took $dt"); t_then = t_now; flush(stdout)

# Task-Affinity CPU Locking in multinode SlurmContext
slurm_cpu_lock()
println(BLAS.get_config()); flush(stdout)
t_now = now(); dt = Dates.canonicalize(Dates.CompoundPeriod(t_now-t_then)); println("CPU locking took $dt"); t_then = t_now; flush(stdout)

using LibGit2; git_branch, git_commit = initalize_git(src_dir); @passobj 1 workers() git_branch; @passobj 1 workers() git_commit

# Set up Korg
@everywhere begin
    runName = "starLine_disk_KnackedKorg";
    wl_lo, wl_hi = 15_000, 17_000 # the wavelength grid is 15_000 : 0.01 : 17_000
    
    apolines = Korg.get_APOGEE_DR17_linelist(include_water=true)
    
    # get the parameter ranges within with the SDSS atmosphere grid exists
    grid_vals = Korg._sdss_marcs_atmospheres[1]
    atm_lbs = first.(grid_vals)
    atm_ubs = last.(grid_vals)

    function query_spectra(intup)
        indx, Teff, logg, A_Xi, vmic = intup
        A_X = collect(A_Xi)
        
        try
            atm = Korg.interpolate_marcs(Teff, logg, A_X; )
            tout = synthesize(atm, apolines, A_X, wl_lo, wl_hi, vmic=vmic, hydrogen_lines=true, hydrogen_line_window_size=300, electron_number_density_warn_threshold=1e10)

            dirindx = lpad(indx ÷ 1000,3,"0")
            specindx = lpad(indx,6,"0")
            dir = "$(runName)/$(dirindx)"
            mkpath(dir)

            writedlm(joinpath(dir, "spectrum_$specindx.dat"), [tout.flux (tout.flux)./tout.cntm])
            return true
        catch
            return false
        end
        GC.gc()
    end
end
Pkg.status("Korg")
t_now = now(); dt = Dates.canonicalize(Dates.CompoundPeriod(t_now-t_then)); println("Korg atm loading took $dt"); t_then = t_now; flush(stdout)

rng = MersenneTwister(368)
# Draw from DR17 Distribution
f = FITS("/uufs/chpc.utah.edu/common/home/sdss/dr17/apogee/spectro/aspcap/dr17/turbo20/allStar-dr17-turbo20.fits")
TEFF = read(f[2],"TEFF")
LOGG = read(f[2],"LOGG")
M_H = read(f[2],"M_H")
ALPHA_M = read(f[2],"ALPHA_M");
close(f)

msk = .!isnan.(TEFF) .& .!isnan.(LOGG) .& .!isnan.(M_H) .& .!isnan.(ALPHA_M)
msk .&= (atm_lbs[1] .< TEFF .< atm_ubs[1])
msk .&= (atm_lbs[2] .< LOGG .< atm_ubs[2])
msk .&= (atm_lbs[3] .< M_H .< atm_ubs[3])
msk .&= (atm_lbs[4] .< ALPHA_M .< atm_ubs[4])
println("(Good Fraction, Num Good) = ",(count(msk)/length(msk), count(msk)))

Nspec_d = 95_000
# A_X[1] is definitionally fixed to 12 and is for H
A_X_d = repeat(Korg.grevesse_2007_solar_abundances', inner=(Nspec_d, 1))
drawinds = rand(rng,1:count(msk),Nspec_d)
m_h_d = M_H[msk][drawinds]
A_X_d[:, 3:Korg.MAX_ATOMIC_NUMBER] .= A_X_d[:, 3:Korg.MAX_ATOMIC_NUMBER] .+ m_h_d
a_m_d = ALPHA_M[msk][drawinds]
for alpha_el in [8, 12, 14, 20, 22] # O Mg Si Ca Ti
    A_X_d[:, alpha_el] = A_X_d[:, alpha_el] .+ a_m_d
end
c_m_d = zeros(Nspec_d)

Teff_d = TEFF[msk][drawinds]
logg_d = LOGG[msk][drawinds]
vmic_d = zeros(Nspec_d)

## Uniform Prior Section
Nspec_u = 5_000
# A_X[1] is definitionally fixed to 12 and is for H
A_X_u = repeat(Korg.grevesse_2007_solar_abundances', inner=(Nspec_u, 1))
A_X_u[:,2:end] .+=  0.01*randn(rng,Nspec_u,Korg.MAX_ATOMIC_NUMBER-1); #0.01 dex sigma

m_h_u = rand(rng,Uniform(atm_lbs[3],atm_ubs[3]), Nspec_u)
A_X_u[:, 3:Korg.MAX_ATOMIC_NUMBER] .= A_X_u[:, 3:Korg.MAX_ATOMIC_NUMBER] .+ m_h_u
# A_X[:,26] .= A_X[:,26] .+ m_h # Fe Only

a_m_u = rand(rng,Uniform(atm_lbs[4],atm_ubs[4]), Nspec_u)
for alpha_el in [8, 12, 14, 20, 22] # O Mg Si Ca Ti
    A_X_u[:, alpha_el] = A_X_u[:, alpha_el] .+ a_m_u
end
c_m_u = zeros(Nspec_u)

Teff_u = rand(rng,Uniform(atm_lbs[1],atm_ubs[1]),Nspec_u)
logg_u = rand(rng,Uniform(atm_lbs[2],atm_ubs[2]),Nspec_u)
vmic_u = zeros(Nspec_u)
# vmic = rand(rng,Uniform(0,5),Nspec);

Nspec = Nspec_d + Nspec_u
Teff = vcat(Teff_d,Teff_u)
logg = vcat(logg_d,logg_u)
vmic = vcat(vmic_d,vmic_u)
m_h = vcat(m_h_d,m_h_u)
a_m = vcat(a_m_d,a_m_u)
c_m = vcat(c_m_d,c_m_u)
A_X = vcat(A_X_d,A_X_u);

mkpath(runName)
h5write("$(runName)/korg_grid_params.h5","Teff",Teff)
h5write("$(runName)/korg_grid_params.h5","logg",logg)
h5write("$(runName)/korg_grid_params.h5","m_h",m_h)
h5write("$(runName)/korg_grid_params.h5","a_m",a_m)
h5write("$(runName)/korg_grid_params.h5","c_m",c_m)
h5write("$(runName)/korg_grid_params.h5","vmic",vmic)
h5write("$(runName)/korg_grid_params.h5","abundances",A_X)

itarg = Iterators.zip(1:Nspec,Teff,logg,eachrow(A_X),vmic);

t_now = now(); dt = Dates.canonicalize(Dates.CompoundPeriod(t_now-t_then)); println("Preparing for and launching pmap $dt"); t_then = t_now; flush(stdout)
pout = @showprogress pmap(query_spectra,itarg);
h5write("$(runName)/korg_grid_params.h5","converged_flag",convert.(Int,pout))
rmprocs(workers())