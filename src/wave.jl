### Wavelength Calibration
# needs to work with pure arclamps and arclamps + FPI
using Optim, FITSIO, Polynomials, StatsBase, Peaks, LsqFit
import FastRunningMedian.running_median
src_dir = abspath("./apMADGICS.jl/") # might need to change if run from pipeline script... but then can use global LibGit2

# Notes
# - build an exclude list into the reader (can still do iterative rejection)
#   but it would help with the cavity fitting on the front end
# - have an exposure pattern recognition function

using LibGit2
git_dir = src_dir
git_commit = LibGit2.head(git_dir)
git_repo = LibGit2.GitRepo(git_dir)
git_head = LibGit2.head(git_repo)
git_branch = LibGit2.shortname(git_head)
println("Running on branch: $git_branch, commit: $git_commit"); flush(stdout)

# arc_group [dither 1, dither 2],[ThAr, UNe], but for now only one dither dimension (removed #1)
arc_grp_tup = [[("sdsswork/mwm", "daily", "apo25m", "59608", 40460007), ("sdsswork/mwm", "daily", "apo25m", "59608", 40460008)]]
fpi_tup = ("sdsswork/mwm", "daily", "apo25m", "59608", 40460010)
function nightly_wavecal(arc_grp_tup, fpi_tup; f2do=1:300, dcav_init=3.73652e7, save_plot_on = true, show_plot_on = false, saveplotspath = "./wavecal_plots/", cache_dir="../local_cache", output_on=false)
    if !ispath(saveplotspath) & save_plot_on
        mkpath(saveplotspath)
    end

    ##### Arc Section #####
    g_flist = []
    for i=1:length(arc_grp_tup)
        push!(g_flist,map(x->build_apLinesPath(x...),arc_grp_tup[i]))
    end
    # Ingest ARC lamp data (consider detector frame corrections, before bad detector replaced)?
    wavelst_arc, pixlst_arc, xpix_arc = ingestArc(g_flist);
    mjd5arc, expidarc = arc_grp_tup[1][1][end-1], string(arc_grp_tup[1][1][end])

    # Fit wavelength solution to ARC Lamp
    fit_pix2wave_arc_partial0(fiber) = fit_pix2wave_arc(wavelst_arc,xpix_arc,fiber)
    outlst_arc0 = fit_pix2wave_arc_partial0.(f2do);
    msklst = map(x->abs.(x[2]).<0.35,outlst_arc0) # residual less than 1 Angstrom cut
    badarcfibers = count.(msklst) .< 75 # this is a fine idea when we have FPI data, but maybe problematic for arcs only
    for badfiber in findall(badarcfibers)
        if length(msklst[badfiber])>0
            msklst[badfiber] .= false ### we should log which fibers this happens to when, rare... but the DRP should fix
        end
    end
    fit_pix2wave_arc_partial(fiber) = fit_pix2wave_arc(wavelst_arc,xpix_arc,fiber,msklst=msklst)
    outlst_arc = fit_pix2wave_arc_partial.(f2do);

    param_lst, offset_param_lst = fit_smoothFibDep(outlst_arc,mjd5arc,expidarc,"Arc";f2do=f2do,save_plot_on=save_plot_on,show_plot_on=show_plot_on,savepath=saveplotspath)
    resid_arc_smooth_partial(fiber) = resid_arc_smooth(wavelst_arc,xpix_arc,param_lst,offset_param_lst,fiber)
    outlst_arc_smooth = resid_arc_smooth_partial.(f2do);
    make_resid_plot2d_arc(outlst_arc,outlst_arc_smooth,mjd5arc,expidarc,save_plot_on=save_plot_on,show_plot_on=show_plot_on,savepath=saveplotspath)
    # save the arc based solution
    wavesavename = cache_wavename(arc_grp_tup[1][1][end-2:end-1]...,cache_dir=cache_dir)
    dirName = splitdir(wavesavename)[1]
    if !ispath(dirName)
        mkpath(dirName)
    end
    save_wavecal(wavesavename,outlst_arc,arc_grp_tup[1][1],"arc",f2do=f2do)
    ## stop if no FPI
    if length(fpi_tup) == 0
        if output_on
            return outlst_arc
        else
            return
        end
    end

    ##### FPI Section #####
    # Ingest FPI data (and transfer wavelength solution to FPI) (consider detector frame corrections, before bad detector replaced)?
    wavelstcombo_FPI0, xpix_FPI0 = ingestFPI(fpi_tup,param_lst,offset_param_lst,f2do=f2do,cache_dir=cache_dir);
    mjd5fpi, expidfpi = fpi_tup[end-1], string(fpi_tup[end])

    # Fit Cavity (derive a "relative" calibration independent of arc lamps based on cavity physics)
    mloclst_FPI, mloclst_msk_FPI, wavelstcombo_FPI = make_mvec(wavelstcombo_FPI0,f2do=f2do,dcav_init=dcav_init);
    xpix_FPI = remove_mask_chip_healper(xpix_FPI0,mloclst_msk_FPI,f2do=f2do)

    cavp = fit_cavity_params(mloclst_FPI,wavelstcombo_FPI,mjd5fpi,expidfpi,f2do=f2do,dcav_init =dcav_init,save_plot_on=false,show_plot_on=false,savepath=saveplotspath)
    ## Fit wavelength solution to FPI; compute/save FPI residuals
    fit_pix2wave_FPI_partial0(fiber) = fit_pix2wave_FPI(mloclst_FPI,xpix_FPI,offset_param_lst,cavp,fiber)
    outlst_FPI = fit_pix2wave_FPI_partial0.(f2do);
    mmsk_lst = make_mmask(mloclst_FPI, outlst_FPI; ccut = 1.0, wid_med = 5, f2do=f2do) # cut outlier FPI modes
    # refit cavity (maybe we need to regenerate wavelstcombo_FPI, this would be an iterative improvement towards a joint cavity and chipgap fit)
    cavp = fit_cavity_params(mloclst_FPI,wavelstcombo_FPI,mjd5fpi,expidfpi,f2do=f2do,dcav_init=dcav_init,msklst=mmsk_lst,save_plot_on=save_plot_on,show_plot_on=show_plot_on,savepath=saveplotspath)
    # refit wavelength solution
    fit_pix2wave_FPI_partial1(fiber) = fit_pix2wave_FPI(mloclst_FPI,xpix_FPI,offset_param_lst,cavp,fiber,msklst=mmsk_lst,scale_on=true)
    outlst_FPI = fit_pix2wave_FPI_partial1.(f2do);
    param_lst_fpi, offset_param_lst_fpi = fit_smoothFibDep(outlst_FPI,mjd5fpi,expidfpi,"FPI",f2do=f2do,save_plot_on=save_plot_on,show_plot_on=show_plot_on,savepath=saveplotspath,scale_on=true)

    # compute/save arclamp residuals
    resid_arc_FPI_partial(fiber) = resid_arc_FPI(wavelst_arc,xpix_arc,outlst_FPI,fiber,scale_on=true)
    outlst_arc_fpi = resid_arc_FPI_partial.(f2do);

    # save plots
    make_resid_plot1d_fpi(outlst_FPI,outlst_arc_fpi,mjd5fpi,expidfpi,fiber=-1,save_plot_on=save_plot_on,show_plot_on=show_plot_on,savepath=saveplotspath)
    make_resid_plot2d_fpi(outlst_FPI,outlst_arc_fpi,mjd5fpi,expidfpi,f2do=f2do,save_plot_on=save_plot_on,show_plot_on=show_plot_on,savepath=saveplotspath)

    # local cache nightly wave cal [complicated since they did arcs every exposure in IV]
    save_wavecal(wavesavename,outlst_FPI,fpi_tup,"fpi",outlst_arc=outlst_arc_fpi,cavp=cavp,write_mode="r+",f2do=f2do)

    if output_on
        return outlst_FPI
    else
        return
    end
end

#### Assorted functions ####

function save_wavecal(wavesavename,outlst,run_tup,wavetype;outlst_arc=nothing,f2do=1:300,cavp=nothing,write_mode="w")
    z = []
    if wavetype == "fpi"
        for fiberiter in f2do
            z = vcat(z,outlst_arc[fiberiter][2])
        end
    else
        for fiberiter in f2do
            z = vcat(z,outlst[fiberiter][2])
        end
    end
    medabs = nanmedian(abs.(z))

    # get the pixel to wavelength polynomial coefficients
    x = extract_nth.(outlst,4)
    binds = findall(length.(x).==0)
    for bind in binds
        x[bind] = NaN*ones(length(x[x.!=0][1]))
    end
    pmat = hcat(x...)

    # get the chip gap polynomial coefficients
    x = extract_nth.(outlst,3)
    binds = findall(length.(x).==0)
    for bind in binds
        x[bind] = NaN*ones(length(x[x.!=0][1]))
    end
    gmat = hcat(x...)

    hdr = FITSHeader(["pipeline","git_branch","git_commit","wave_expid","wavetype"],["apMADGICS.jl",git_branch,git_commit,string(run_tup[end]),wavetype],["","","","",""])

    f = FITS(wavesavename,write_mode)
    write(f,[0],header=hdr,name="$(wavetype)_header_only")
    write(f,pmat,name="$(wavetype)_pix2wave_polycoeff")
    write(f,gmat,name="$(wavetype)_chipgap_polycoeff")
    write(f,[medabs],name="$(wavetype)_medabs_arcresid")
    if wavetype == "fpi"
        # make cavity dictionary
        cavColNames = ["dcav", "m0"]
        cavColVals = [[cavp[1]],[cavp[2]]]
        write(f,cavColNames,cavColVals,name="$(wavetype)_cavity_params")
    end
    close(f)
end

# takes in a list of wave/apLines files and creates full exposure lists of wavelength and pixel
# this wavelength is from the DRP, and we discard it
# can be a list of lists to allow for ThAr/UNe pairs with a dither shift
function ingestArc(glist; chip2do = ["a","b","c"], f2do = 1:300)
    wavelst_arc = []
    pixlst_arc = []
    xpix_arc_g = []
    for sgroupind = 1:length(glist)
        wavelstg = []
        pixlstg = []
        for fname in glist[sgroupind]
            wavelstf = []
            pixlstf = []
            for (chipind, chip) in enumerate(chip2do)
                f = FITS(fname)
                wave = read(f[2],"wave")
                wave_found = read(f[2],"wave_found")
                pixel = read(f[2],"pixel")
                chipv = read(f[2],"chip")
                row = read(f[2],"row");
                close(f)
                fiberlstloc = []
                pixlstloc = []
                for fiber=f2do
                    msk = (row .== fiber-1) .& (chipv.==getChipIndx(chip)) .& (0 .< pixel .< 2049)
                    push!(fiberlstloc,wave[msk])
                    push!(pixlstloc,pixel[msk])
                end
                push!(wavelstf,fiberlstloc)
                push!(pixlstf,pixlstloc)
            end
            push!(wavelstg,wavelstf)
            push!(pixlstg,pixlstf)
        end
        push!(wavelst_arc,wavelstg)
        push!(pixlst_arc,pixlstg)

        xpix_arc = []
        for (findx, fname) in enumerate(glist[sgroupind])
            xpix_arc_loc=[]
            for fiber=f2do
                xloc = []
                for (chipind, chip) in enumerate(chip2do)
                    push!(xloc,pixlstg[findx][chipind][fiber])
                end
                push!(xpix_arc_loc,xloc)
            end
            push!(xpix_arc,xpix_arc_loc)
        end
        push!(xpix_arc_g,xpix_arc)
    end
    return wavelst_arc, pixlst_arc, xpix_arc_g
end

function make_wvec_arc(wavelst_arc,fiber;sgindx=1,chip2do = ["a","b","c"])
    pixlst = []
    for gindx in 1:length(wavelst_arc[sgindx])
        for (chipind, chip) in enumerate(chip2do)
            if length(wavelst_arc[sgindx][gindx][chipind][fiber])>=1
                push!(pixlst,wavelst_arc[sgindx][gindx][chipind][fiber])
            else
                push!(pixlst,[])
            end
        end
    end
    return xvec = vcat(pixlst...)
end

function ingestFPI(fpi_tup, param_lst, offset_param_lst; use_drp=false, peak_wid=5, widx=3, chip2do = ["a","b","c"], f2do = 1:300, cache_dir="./local_cache") 
    wavelst_FPI = []
    pixlst_FPI = []
    for (chipind, chip) in enumerate(chip2do)
        pixcen, row = extractFPIpeaks(fpi_tup,chip,use_drp=use_drp,peak_wid=peak_wid,widx=widx,cache_dir=cache_dir)
        fiberlst = []
        pixlst = []
        for fiber=f2do
            msk = (row .== fiber-1);
            pvec_off = Polynomial.(offset_param_lst)
            offset_params_opt = map(x->x((fiber-150)/150),pvec_off)
            xt_FPI = make_xvec([pixcen[msk]], offset_params_opt, chip2do=[chip])
            Axfpi = positional_poly_mat(xt_FPI,porder=3)
            pvec = Polynomial.(param_lst[2:end])
            if !isnan(param_lst[1][fiber])
                params_opt = vcat(param_lst[1][fiber],map(x->x((fiber-150)/150),pvec))
                wave = Axfpi*params_opt
                push!(pixlst,pixcen[msk])
                push!(fiberlst,wave)
            else
                push!(pixlst,[])
                push!(fiberlst,[])
            end
        end
        push!(pixlst_FPI,pixlst)
        push!(wavelst_FPI,fiberlst)
    end
    
    xpix_FPI=[]
    for fiber=f2do
        xloc = []
        for (chipind, chip) in enumerate(chip2do)
            push!(xloc,pixlst_FPI[chipind][fiber])
        end
        push!(xpix_FPI,xloc)
    end
    
    wavelstcombo_FPI = []
    # pixlstcombo_FPI = []
    for fiber=f2do
        wavelst = []
        # pixlst = []
        for (chipind, chip) in enumerate(chip2do)
            if length(wavelst_FPI[chipind][fiber])>=1
                push!(wavelst,wavelst_FPI[chipind][fiber])
                # push!(pixlst,pixlst_FPI[chipind][fiber])
            else
                push!(wavelst,[])
                # push!(pixlst,[])
            end
        end
        push!(wavelstcombo_FPI,vcat(wavelst...))
        # push!(pixlstcombo_FPI,vcat(pixlst...))
    end

    return wavelstcombo_FPI, xpix_FPI
end

# Takes in wavelength estimates for FPI peaks and converts to integer peak index estimates
function make_mvec(wavelstcombo_FPI; f2do = 1:300, swindow = 41, leadwindow = 10, maxpass = 100, dcav_init = 3.736e7, dvec_cut = 2, dvec_frac_cut = 0.04) 
    mloclst_FPI= []
    mloclst_msk_FPI = []
    for fiber=f2do
        wavelen = length(wavelstcombo_FPI[fiber])
        if wavelen>1
            waveindx = 1:wavelen
            mskg = ones(Bool,wavelen)
            for i=1:maxpass 
                dvec = diff(wavelstcombo_FPI[fiber][mskg])
                mbad = abs.(dvec) .< dvec_cut
                if count(mbad)==0
                    break
                end
                bindx = argmin(abs.(dvec))
                rindx = waveindx[mskg][bindx]:waveindx[mskg][bindx + 1]
                mskg[rindx].=false
                if i==maxpass
                    error("Bad m index removal 1 failed badly, $fiber")
                    # println("Bad m index removal failed badly, $fiber")
                end
            end
            for i=1:maxpass
                dvec = diff(wavelstcombo_FPI[fiber][mskg])
                dm = dvec./running_median(dvec,swindow,:asymmetric_truncated)
                dt = (dm.-roundnan.(dm))./dm
                mbad = (abs.(dt) .> dvec_frac_cut)
                if count(mbad)==0
                    break
                end
                bindx = argmax(abs.(dt))
                rindx = waveindx[mskg][bindx]:waveindx[mskg][bindx + 1]
                mskg[rindx].=false
                if i==maxpass
                    error("Bad m index removal 2 failed badly, $fiber")
                    # println("Bad m index removal failed badly, $fiber")
                end
            end
            mtemp = 2*dcav_init./wavelstcombo_FPI[fiber][mskg] # you don't need to know the cavity all that well to assign deltaM
            mloc = roundnan.(mtemp .- mtemp[1] .+ 1)
            push!(mloclst_FPI,mloc)
            push!(mloclst_msk_FPI,mskg)
        else
            push!(mloclst_FPI,[])
            push!(mloclst_msk_FPI,[])
        end
    end
    
    wavelstcombo_FPI1 = remove_mask_healper(wavelstcombo_FPI,mloclst_msk_FPI,f2do=f2do)
    wvect = map(maximum_empty,wavelstcombo_FPI1)
    mwval, mwind = findmax(filter(!isnan,wvect))
    dstep = median(running_median(diff(wavelstcombo_FPI1[mwind]),swindow,:asymmetric_truncated)[1:leadwindow])
    madd = roundnan.((wvect.-mwval)./dstep);
    
    for fiber=f2do
        mloclst_FPI[fiber].+=madd[fiber]
    end

    return mloclst_FPI, mloclst_msk_FPI, wavelstcombo_FPI1
end

function make_mmask(mloclst, outlst; ccut = 1, wid_med = 5, f2do = 1:300)
    refiqr = fiqr(vcat(map(x->x[2],outlst)...))

    mat = vcat(mloclst...)
    rmat = vcat(map(x->x[2],outlst)...)
    ulst = unique(mat)
    om1_lst = []
    om2_lst = []
    for um in ulst
        msk = (mat.== um)
        mskw = (abs.(mat .- um) .<= wid_med)
        rlqr, rmedian, ruqr = StatsBase.percentile(rmat[mskw],[25,50,75])
        riqr = (ruqr-rlqr)/1.34896
        push!(om1_lst,StatsBase.median(abs.(rmat[msk].-rmedian))/riqr)
        push!(om2_lst,StatsBase.median(abs.(rmat[msk].-rmedian))/refiqr)
    end

    badm = ulst[(om1_lst.>ccut) .| (om2_lst.>ccut)]

    mmsk_lst = []
    for fiber=f2do
        push!(mmsk_lst,.!(mloclst[fiber].∈[badm]))
    end
    return mmsk_lst
end

# x_pixels = pixlst_FPI[chipind][fiber]
# used to make the x-axis for global fit across the chips 
function make_xvec(x_pixels, paramst; chip2do = ["a","b","c"], chip_coord_norm=[1.00032,1.0,0.99991])
    pixlst = []
    for (chipind, chip) in enumerate(chip2do)
        if (length(x_pixels)>=1) & (getChipIndx(chip) == 1)
            push!(pixlst,(x_pixels[chipind])./chip_coord_norm[1].-paramst[1].-1024)
        elseif (length(x_pixels)>=1) & (getChipIndx(chip) == 2)
            push!(pixlst,(x_pixels[chipind]).-1024)
        elseif (length(x_pixels)>=1) & (getChipIndx(chip) == 3)
            push!(pixlst,(x_pixels[chipind])./chip_coord_norm[3].-paramst[2].-1024)
        else
            push!(pixlst,[])
        end
    end
    xvec = vcat(pixlst...)
end

# x_pixels = xpix_arc which is of dimes [dither group][ThAr/UNe][fiber]
function make_xvec_arc_wrap(x_pixels, paramst, fiber; chip2do = ["a","b","c"])
    xvec_arc = []
    for i=1:length(x_pixels[1]) # forcing to only use one dither
        xvec_arc = vcat(xvec_arc,make_xvec(x_pixels[1][i][fiber],paramst,chip2do=chip2do))
    end
    return xvec_arc
end

function positional_poly_mat(x; porder=4)
    xtrans = x./2048; # just transforms x to be in units of "CCD widths" for numerical stability
    if porder == 3
        return [ones(length(xtrans)) xtrans xtrans.^2 xtrans.^3]
    elseif porder == 4
        return [ones(length(xtrans)) xtrans xtrans.^2 xtrans.^3 xtrans.^4]
    end
end

function m2lam(m,p;moff=0)
    # Simple λ = 2D/m model at the moment, but possibility of expanding
    m0 = p[end]
    dparams = p[1:end-1]
    #D(m) = D0 + D1/m0(m-m0-moff) + D2/m0^2(m-m0-moff)^2 + ...
    #madj = (m-moff)
    return 2 ./ (m.+m0).*(dparams[1])# + dparams[2]*madj/m0)# + dparams[3]*(madj/m0)^2)
end

function fit_smoothFibDep(outlst_arc,mjd5,expid,exptype;f2do=1:300,save_plot_on=false,show_plot_on=false,savepath="./",scale_on=false)
    fibaxis = (f2do.-150)/150
    
    # Fit Poly pix2wave Poly Order v Fiber Index
    porder = length(outlst_arc[1][end])
    order_lst = [2*ones(Int,porder-1)...]
    param_lst = []
    fit_xyf = []
    push!(param_lst,map(x->extract_nth(x[end],1),outlst_arc))
    push!(fit_xyf,(fibaxis,param_lst[1],param_lst[1]))
    for i=2:porder
        outP = map(x->extract_nth(x[end],i),outlst_arc);
        msk = (outP .== 0) .| (isnan.(outP))
        xs = fibaxis[.!msk]
        ys = outP[.!msk]
        p = Polynomials.fit(xs,ys,order_lst[i-1])
        push!(param_lst,p.coeffs)
        push!(fit_xyf,(fibaxis,outP,p.(fibaxis)))
    end
    
    # Fit Detector Offsets v Fiber Index
    order_lst = scale_on ? [1,1,1,1,] : [1,1]
    offset_param_lst = []
    fit_offset_xyf = []
    n_params = scale_on ? 4 : 2
    for i=1:n_params
        outP = map(x->extract_nth(x[end-1],i),outlst_arc);
        msk = (outP .== 0) .| (isnan.(outP))
        xs = fibaxis[.!msk]
        ys = outP[.!msk]
        p = Polynomials.fit(xs,ys,order_lst[i])

        push!(offset_param_lst,p.coeffs)
        push!(fit_offset_xyf,(fibaxis,outP,p.(fibaxis)))
    end
    
    if save_plot_on | show_plot_on
        make_fitplots(fit_xyf,fit_offset_xyf,mjd5,expid,exptype,save_plot_on=save_plot_on,show_plot_on=show_plot_on,savepath=savepath,scale_on=scale_on)
    end
    return param_lst, offset_param_lst
end

# refines guesses for the mode spacing and cavity length, then performs curve fitting
# DOES THIS NEED TO BE REFITTING THE CHIP GAPS?
function fit_cavity_params(mloclst_FPI,wavelstcombo_FPI,mjd5,expid; msklst=nothing, f2do=1:300, dcav_init = 3.736e7,save_plot_on=false,show_plot_on=false,savepath="./")
    mvec = if isnothing(msklst)
        vcat(map(fiber->mloclst_FPI[fiber],f2do)...);
    else
        vcat(map(fiber->mloclst_FPI[fiber][msklst[fiber]],f2do)...);
    end
    wvec = if isnothing(msklst)
        convert(Array{Float64},vcat(map(fiber->wavelstcombo_FPI[fiber],f2do)...));
    else
        convert(Array{Float64},vcat(map(fiber->wavelstcombo_FPI[fiber][msklst[fiber]],f2do)...));
    end
    ## refine series of initial guesses
    # first guess for m0
    m0est = vcat(map(fiber->2*dcav_init./wavelstcombo_FPI[fiber].-mloclst_FPI[fiber],f2do)...);
    moff_0 = round(Int,median(filter(!isnan,m0est)))
    # refine dcavity
    A = reshape(2 ./(moff_0.+mvec),(:,1)).* reshape(ones(length(mvec)),:,1);
    d_cav = (A\wvec)[1]
    # final guess for m0
    m0est = vcat(map(fiber->2*d_cav./wavelstcombo_FPI[fiber].-mloclst_FPI[fiber],f2do)...);
    m0final = median(filter(!isnan,m0est))
    # curve fit using LM through LsqFit.jl to obtain cavity parameters
    p0 = [d_cav, m0final]
    # p0 = [dcav_init, moff_0]
    # I guess we need to set these bounds differently for the LCO case
    fit = curve_fit(m2lam, mvec, wvec, p0, autodiff=:forwarddiff, show_trace=false) #, lower=[3.73651e7,1], upper=[3.73653e7,Inf])
    fit_out = fit.param
    if save_plot_on | show_plot_on
        fig = plt.figure(figsize=(8,8),dpi=200)
        ax = fig.add_subplot(1,1,1)
        ax.set_title("MJD5: $(mjd5) \n D Cavity: $(round(fit_out[1],sigdigits=8)) Å \n m0: $(round(fit_out[2],sigdigits=8))")
        ax.hist2d(mvec,wvec.-m2lam(mvec,fit_out),bins=201,cmap="cet_gouldian",range=((0,530),(-0.2,0.2)))
        ax.set_xlabel("FPI Mode (Rel) Index")
        ax.set_ylabel("Residuals for Cavity Fit (Angstroms) \n (FPI Wave given ArcLamp Soln)")
        if save_plot_on
            fig.savefig(savepath*"cavity_fit_fpi_$(mjd5)_$(expid).png", bbox_inches="tight", pad_inches=0.1);
        end
        if !show_plot_on
            plt.close()
        end
    end

    return fit_out
end

function fit_pix2wave_arc(wavelst_arc,xpix_arc,fiber;msklst=nothing)
    wvec_arc = make_wvec_arc(wavelst_arc,fiber,sgindx=1);
    if length(wvec_arc) > 0 
        offset_init = [2191.5, -2202.5]
        msk = if isnothing(msklst)
            ones(Bool,length(wvec_arc))
        else
            msklst[fiber]
        end
        if count(msk) > 0
            partial_loss_arc(offsetv) = loss_arc_NL(xpix_arc,wvec_arc,fiber,offsetv,msk=msk)
            res = optimize(partial_loss_arc, offset_init, LBFGS(), Optim.Options(show_trace=false))
            params_opt = Optim.minimizer(res)
            x_arc = make_xvec_arc_wrap(xpix_arc,params_opt,fiber);
            Axarc = positional_poly_mat(x_arc[msk],porder=3)
            tparam = Axarc\wvec_arc[msk]
            return (x_arc[msk],wvec_arc[msk].-Axarc*tparam,params_opt,tparam)
        else
            return ([],[],[],[])
        end
    else
        return ([],[],[],[])
    end
end


function fit_pix2wave_FPI(mloclst_FPI,xpix_FPI,arc_fit_offset_poly,cavp,fiber;msklst=nothing,scale_on=false)
    if length(mloclst_FPI[fiber]) > 0
        wavec = m2lam(mloclst_FPI[fiber],cavp);
        msk = if isnothing(msklst)
            ones(Bool,length(wavec))
        else
            msklst[fiber]
        end
        pvec_off = Polynomial.(arc_fit_offset_poly)
        offset_init = if scale_on
            [map(x->x((fiber-150)/150),pvec_off)...,1.0001,0.99981] #1.0001,1.0,0.99981 #### 1.00025,0.9999
        else
            map(x->x((fiber-150)/150),pvec_off)
        end
        partial_loss_FPI(offsetv) = loss_FPI_NL(xpix_FPI[fiber],wavec,offsetv,msk=msk,scale_on=scale_on)
        res = optimize(partial_loss_FPI, offset_init, LBFGS(), Optim.Options(show_trace=false))
        # show(res)
        params_opt = Optim.minimizer(res)
        x_FPI = if scale_on 
            make_xvec(xpix_FPI[fiber],params_opt[1:2],chip_coord_norm=[params_opt[3],1.0,params_opt[4]])
        else
            make_xvec(xpix_FPI[fiber],params_opt[1:2])
        end
        Axfpi = positional_poly_mat(x_FPI[msk])
        tparam = Axfpi\wavec[msk]
        return (x_FPI[msk],wavec[msk].-Axfpi*tparam,params_opt,tparam)
    else
        return ([],[],[],[])
    end
end

function resid_arc_smooth(wavelst_arc,xpix_arc,param_lst,offset_param_lst,fiber)
    wvec_arc = make_wvec_arc(wavelst_arc,fiber,sgindx=1);
    if length(wvec_arc) > 1
        xvec_arc = []
        pvec_off = Polynomial.(offset_param_lst)
        offset_params_opt = map(x->x((fiber-150)/150),pvec_off)
        # forcing to only use one dither, here I only have one (last one, but could use better handling)
        for i=1:length(xpix_arc[1]) 
            xvec_arc = vcat(xvec_arc,make_xvec(xpix_arc[1][i][fiber],offset_params_opt))
        end
        Axarc = positional_poly_mat(xvec_arc,porder=3)
        pvec = Polynomial.(param_lst[2:end])
        params_opt = vcat(param_lst[1][fiber],map(x->x((fiber-150)/150),pvec))
        wavepred_arc = Axarc*params_opt
        return  (xvec_arc,wvec_arc.-wavepred_arc)
    else
        return ([],[])
    end
end

function resid_arc_FPI(wavelst_arc,xpix_arc,fpi_out_params,fiber;scale_on=false)
    wvec_arc = make_wvec_arc(wavelst_arc,fiber,sgindx=1);
    if (length(wvec_arc) > 1) & (length(fpi_out_params[fiber][end-1]) > 0)
        xvec_arc = []
        # forcing to only use one dither, here I only have one (last one, but could use better handling)
        for i=1:length(xpix_arc[1])
            params_opt = fpi_out_params[fiber][end-1]
            x_FPI = if scale_on 
                make_xvec(xpix_arc[1][i][fiber],params_opt[1:2],chip_coord_norm=[params_opt[3],1.0,params_opt[4]])
            else
                make_xvec(xpix_arc[1][i][fiber],params_opt[1:2])
            end
            xvec_arc = vcat(xvec_arc,x_FPI)
        end
        Axarc = positional_poly_mat(xvec_arc)
        wavepred_arc = Axarc*fpi_out_params[fiber][end];
        return  (xvec_arc,wvec_arc.-wavepred_arc)
    else
        return ([],[])
    end
end

# I guess I just need to tune Optim's options
function loss_arc_NL(pixlst, wavec, fiber, offsetv; msk=nothing)
    wavecm = wavec[msk]
    xvec = make_xvec_arc_wrap(pixlst, offsetv, fiber)
    Axfpi = positional_poly_mat(xvec[msk],porder=3)
    tparam = Axfpi\wavecm
    return sum((wavecm.-Axfpi*tparam).^2)
end

function loss_FPI_NL(pixlst, wavec, paramvec; msk=nothing, scale_on=false)
    wavecm = wavec[msk]
    offsetv = paramvec[1:2]
    xvec = if scale_on 
        make_xvec(pixlst,offsetv,chip_coord_norm=[paramvec[3],1.0,paramvec[4]])
    else
        make_xvec(pixlst,offsetv)
    end
    Axfpi = positional_poly_mat(xvec[msk])
    tparam = Axfpi\wavecm
    return sum((wavecm.-Axfpi*tparam).^2)
end

function make_fitplots(fit_xyf,fit_offset_xyf,mjd5,expid,exptype;save_plot_on=false,show_plot_on=false,savepath="./",scale_on = false)
    szfit = length(fit_xyf)
    # Poly Fit Plots
    fig = plt.figure(figsize=(24,8),dpi=150)
    plt.suptitle(exptype*"\n MJD5: $(mjd5)",y=0.97,fontsize=14)
    for i=2:szfit
        x,y,f = fit_xyf[i]
        ax = fig.add_subplot(2,szfit-1,i-1)
        ax.scatter(x*150 .+150,y,s=1)
        ax.plot(x*150 .+150,f,color="red",lw=1)
        ax.set_title("Order $(i-1) Poly Term")

        ax = fig.add_subplot(2,szfit-1,szfit-1+i-1)
        ax.scatter(x*150 .+150,y.-f,s=1)
        ax.axhline(0,lw=1,linestyle="--")
        if !scale_on
            ax.set_xlabel("Fiber Index")
        end
        ax.set_ylabel("Residuals")
    end
    if save_plot_on
        fig.savefig(savepath*"poly_fit_$(exptype)_$(mjd5)_$(expid).png", bbox_inches="tight", pad_inches=0.1);
    end
    if !show_plot_on
        plt.close()
    end
    
    # Detector Centering Plots
    n_det_plots = scale_on ? 4 : 2
    y_plot_size = scale_on ? 16 : 8
    fig = plt.figure(figsize=(24,y_plot_size),dpi=150)
    x,y,f = fit_xyf[1]
    plt.suptitle(exptype*"\n MJD5: $(mjd5)",y=0.94,fontsize=14)
    plt.subplots_adjust(hspace=0.2,wspace=0.2)
    ax = fig.add_subplot(n_det_plots,3,1)
    ax.scatter(x*150 .+150,y,s=1)
    ax.set_title("Wavelength at G Detector Center (Å)")

    ax.set_xlabel("Fiber Index")
    title_lst = ["G to B","G to R"]
    for i=1:2
        x,y,f = fit_offset_xyf[i]
        ax = fig.add_subplot(n_det_plots,3,1+i)
        ax.scatter(x*150 .+150,y,s=1)
        ax.plot(x*150 .+150,f,color="red",lw=1)
        ax.set_title("Distance from $(title_lst[i]) Detector Centers (pixels)")

        ax = fig.add_subplot(n_det_plots,3,4+i)
        ax.scatter(x*150 .+150,y.-f,s=1)
        ax.axhline(0,lw=1,linestyle="--")
        if !scale_on
            ax.set_xlabel("Fiber Index")
        end
        ax.set_ylabel("Residuals")
    end
    if scale_on
        title_lst = ["R","B"] # check on the label ordering
        for i=3:4
            x,y,f = fit_offset_xyf[i]
            ax = fig.add_subplot(n_det_plots,3,5+i)
            ax.scatter(x*150 .+150,y,s=1)
            ax.plot(x*150 .+150,f,color="red",lw=1)
            ax.set_title("Scale Factor for $(title_lst[i-2]) Detector")
    
            ax = fig.add_subplot(n_det_plots,3,8+i)
            ax.scatter(x*150 .+150,y.-f,s=1)
            ax.axhline(0,lw=1,linestyle="--")
            ax.set_xlabel("Fiber Index")
            ax.set_ylabel("Residuals")
        end
    end

    if save_plot_on
        fig.savefig(savepath*"detector_fit_arc_$(mjd5)_$(expid).png", bbox_inches="tight", pad_inches=0.1);
    end
    if !show_plot_on
        plt.close()
    end
end

function make_resid_plot2d_arc(outlst_arc,outlst_arc_smooth,mjd5,expid;save_plot_on=false,show_plot_on=false,f2do=1:300,savepath="./")
    fig = plt.figure(figsize=(16,8),dpi=300)
    fig.suptitle("MJD5: $(mjd5) \n EXPID: $(expid)",y=0.97,fontsize=14)

    ax = fig.add_subplot(1,2,1)
    x, y, z = [], [], []
    for fiberiter in f2do
        x = vcat(x,outlst_arc[fiberiter][1])
        y = vcat(y,fiberiter*ones(length(outlst_arc[fiberiter][1])))
        z = vcat(z,outlst_arc[fiberiter][2])
    end
    sc = ax.scatter(x,y,c=z, cmap="cet_bkr", vmin=-0.1, vmax = 0.1, s=1)
    ax.set_xlabel("X Coordinate (pixels)")
    # ax.set_ylabel("Fiber Index")
    ax.set_yticks([])

    divider = mpltk.make_axes_locatable(ax)
    cax = divider.append_axes("top", size="5%", pad=0.05)
    cbar = plt.colorbar(sc, cax=cax, orientation="horizontal")
    cax.xaxis.set_ticks_position("top")
    cax.xaxis.set_label_position("top")
    cax.set_xlabel("Arc Residuals to Arc Fit (Å)",labelpad=10)

    ax = fig.add_subplot(1,2,2)
    x, y, z = [], [], []
    for fiberiter in f2do
        x = vcat(x,outlst_arc_smooth[fiberiter][1])
        y = vcat(y,fiberiter*ones(length(outlst_arc_smooth[fiberiter][1])))
        z = vcat(z,outlst_arc_smooth[fiberiter][2])
    end
    sc = ax.scatter(x,y,c=z, cmap="cet_bkr", vmin=-0.1, vmax = 0.1, s=1)
    ax.set_xlabel("X Coordinate (pixels)")
    # ax.set_ylabel("Fiber Index")
    ax.set_yticks([])

    divider = mpltk.make_axes_locatable(ax)
    cax = divider.append_axes("top", size="5%", pad=0.05)
    cbar = plt.colorbar(sc, cax=cax, orientation="horizontal")
    cax.xaxis.set_ticks_position("top")
    cax.xaxis.set_label_position("top")
    cax.set_xlabel("Poly Arc Residuals to Arc Fit (Å)",labelpad=10)

    if save_plot_on
        fig.savefig(savepath*"resid_plot2d_arc_$(mjd5)_$(expid).png", bbox_inches="tight", pad_inches=0.1);
    end
    if !show_plot_on
        plt.close()
    end
end

function make_resid_plot2d_fpi(outlst_FPI,outlst_arc,mjd5,expid;save_plot_on=false,show_plot_on=false,f2do=1:300,savepath="./")
    fig = plt.figure(figsize=(16,8),dpi=300)
    fig.suptitle("MJD5: $(mjd5) \n EXPID: $(expid)",y=0.97,fontsize=16)
    ax = fig.add_subplot(1,2,1)
    x, y, z = [], [], []
    for fiberiter in f2do
        x = vcat(x,outlst_FPI[fiberiter][1])
        y = vcat(y,fiberiter*ones(length(outlst_FPI[fiberiter][1])))
        z = vcat(z,outlst_FPI[fiberiter][2])
    end
    sc = ax.scatter(x,y,c=z, cmap="cet_bkr", vmin=-0.01, vmax = 0.01, s=1)
    ax.set_xlabel("X Coordinate (pixels)")
    ax.set_ylabel("Fiber Index")

    divider = mpltk.make_axes_locatable(ax)
    cax = divider.append_axes("top", size="5%", pad=0.05)
    cbar = plt.colorbar(sc, cax=cax, orientation="horizontal")
    cax.xaxis.set_ticks_position("top")
    cax.xaxis.set_label_position("top")
    cax.set_xlabel("FPI Residuals to FPI Fit (Å)",labelpad=10)

    ax = fig.add_subplot(1,2,2)
    x, y, z = [], [], []
    for fiberiter in f2do
        x = vcat(x,outlst_arc[fiberiter][1])
        y = vcat(y,fiberiter*ones(length(outlst_arc[fiberiter][1])))
        z = vcat(z,outlst_arc[fiberiter][2])
    end
    sc = ax.scatter(x,y,c=z, cmap="cet_bkr", vmin=-0.1, vmax = 0.1, s=1)
    ax.set_xlabel("X Coordinate (pixels)")
    # ax.set_ylabel("Fiber Index")
    ax.set_yticks([])

    divider = mpltk.make_axes_locatable(ax)
    cax = divider.append_axes("top", size="5%", pad=0.05)
    cbar = plt.colorbar(sc, cax=cax, orientation="horizontal")
    cax.xaxis.set_ticks_position("top")
    cax.xaxis.set_label_position("top")
    cax.set_xlabel("Arc Residuals to FPI Fit (Å)",labelpad=10)

    if save_plot_on
        fig.savefig(savepath*"resid_plot2d_fpi_$(mjd5)_$(expid).png", bbox_inches="tight", pad_inches=0.1);
    end
    if !show_plot_on
        plt.close()
    end
end

function make_resid_plot1d_fpi(outlst_FPI,outlst_arc,mjd5,expid;fiber=-1,save_plot_on=false,show_plot_on=false,savepath="./",f2do=1:300)
    fibertst = if fiber<0
        rng = MersenneTwister(parse(Int,expid))
        rand(rng,f2do)
    else
        fiber
    end

    fig = plt.figure(figsize=(12,5),dpi=200)
    fig.suptitle("Fiber Index: $(fibertst) \n MJD5: $(mjd5)",y=0.97,fontsize=16)
    ax = fig.add_subplot(1,2,1)
    ax.scatter(outlst_FPI[fibertst][1],outlst_FPI[fibertst][2],s=1)
    ax.set_ylim(-0.03,0.03)
    ax.axhline(0,color="red",linestyle="--")
    ax.set_xlabel("X Coordinate (pixels)")
    ax.set_ylabel("Residuals (Å)")
    iqr_fpi = string(round(iqr(outlst_FPI[fibertst][2])/1.34896,digits=4))
    ax.set_title("FPI Residuals to FPI Fit \n "* L"$\sigma_{IQR}$"* " = $iqr_fpi")
    
    ax = fig.add_subplot(1,2,2)
    ax.scatter(outlst_arc[fibertst][1],outlst_arc[fibertst][2],s=1)
    ax.set_ylim(-0.25,0.25)
    ax.axhline(0,color="red",linestyle="--")
    ax.set_xlabel("X Coordinate (pixels)")
    # ax.set_ylabel("Residuals (Å)")
    
    iqr_arc = string(round(iqr(outlst_arc[fibertst][2])/1.34896,digits=4))
    ax.set_title("Arc Residuals to FPI Fit \n "* L"$\sigma_{IQR}$"* " = $iqr_arc")

    if save_plot_on
        fig.savefig(savepath*"resid_plot1d_fpi_$(mjd5)_$(expid).png", bbox_inches="tight", pad_inches=0.1);
    end
    if !show_plot_on
        plt.close()
    end
end

### sort of a local collected utils

#generalize the use of this throughout
function getChipIndx(chip)
    if chip == "a"
        return 1
    elseif chip == "b"
        return 2
    elseif chip=="c"
        return 3
    end
end

function remove_mask_healper(inarray, msk_lst; f2do=1:300)
    outarray = []
    for fiber=f2do
        if length(inarray[fiber])>0
            push!(outarray,inarray[fiber][msk_lst[fiber]])
        else
            push!(outarray,[])
        end
    end
    return outarray
end

function remove_mask_chip_healper(inarray, msk_lst; f2do=1:300)
    outarray = []
    for fiber=f2do
        chiplens = length.(inarray[fiber])
        startind = 1
        loc_outarray = []
        for (chipindx, chiplen) in enumerate(chiplens)
            if chiplen>0
                push!(loc_outarray,inarray[fiber][chipindx][msk_lst[fiber][startind:(startind+chiplen-1)]])
                startind += chiplen
            else
                push!(loc_outarray,[])
            end
        end
        push!(outarray,loc_outarray)
    end
    return outarray
end

function extract_nth(x,n)
    if length(x)>1
        return x[n]
    else
        return NaN
    end
end

function median_empty(x)
    if length(x)>1
        return median(x)
    else
        return NaN
    end
end

function maximum_empty(x)
    if length(x)>1
        return maximum(x)
    else
        return NaN
    end
end

function sum_empty(x)
    if length(x)>1
        return sum(x)
    else
        return NaN
    end
end

function roundnan(x)
    if isnan(x)
        return 0
    else
        return round(Int,x)
    end
end

function fiqr(x)
    return iqr(x)/1.34896
end

####

function extractFPIpeaks(fpi_tup,chip;use_drp=false,peak_wid=5,widx=3,cache_dir="./local_cache")
    if use_drp
        fname = build_apFPILinesPath(fpi_tup[1:end-1]...,chip,fpi_tup[end])
        f = FITS(fname)
        pixcen = read(f[5],"pars")[2,:]
        # chip = read(f[5],"chip")
        row = read(f[5],"row");
        # wavedrp = read(f[5],"wave")
        # linewave = read(f[5],"linewave");
        close(f)
        return pixcen, row
    else
        waveLineFPIsavename = cache_waveLinesFPIname(fpi_tup[end-2:end]...,chip; cache_dir=cache_dir)
        dirName = splitdir(waveLineFPIsavename)[1]
        if !ispath(dirName)
            mkpath(dirName)
        end
        if !isfile(waveLineFPIsavename)
            ap1D_path = build_framepath(fpi_tup[1:end-1]...,fpi_tup[end],chip)
            f = FITS(ap1D_path)
            ref_img = read(f[2]);
            err_img = read(f[3]);
            flg_img = read(f[4]);
            close(f)

            msk_bad = (flg_img .& (2^0 | 2^1 | 2^2 | 2^3 | 2^4 | 2^5 | 2^6 | 2^7 | 2^12 | 2^13 | 2^14) .!=0) #https://www.sdss4.org/dr17/irspec/apogee-bitmasks/
            nan_img = copy(ref_img)
            nan_img[msk_bad].=NaN;

            FPIoutlst = []
            for fibindx in 1:300
                dat_vec = nan_img[:,fibindx]
                err_vec = err_img[:,fibindx]
                pks_i, vals_i = findmaxima(dat_vec,peak_wid,strict=false)
                msk = 0.0 .< vals_i .< 1e5
                pks = pks_i[msk]
                vals = vals_i[msk]
                get_gaussian_indx_partial(gindx) = get_gaussian_indx(gindx,pks,vals,dat_vec,err_vec,fibindx,widx=widx)
                tout = map(get_gaussian_indx_partial,1:length(vals))
                push!(FPIoutlst,tout)
            end
            FPIoutmat = hcat(vcat(FPIoutlst...)...)
            save_FPIpeaks(waveLineFPIsavename,fpi_tup,FPIoutmat)
        end
        f = FITS(waveLineFPIsavename,"r")
        pixcen = read(f["gauss_FPI_params"],"pix_cen").-1 # subtracting 1 to match the DRP python indexing, can remove if I arclamps move to this ecosystem
        row = read(f["gauss_FPI_params"],"fiberindx")
        num_data = read(f["gauss_FPI_params"],"num_data")
        close(f)
        msk = (num_data .== (2*widx + 1)) .& (0 .< pixcen .< 2048) # could loosen or add more conditions, but must be > 4
        return pixcen[msk], row[msk]
    end
end

function get_gaussian_indx(gindx,pks,vals,y,yerr,fibindx;sigma0=1.0,widx=3)
    initial_guess = [vals[gindx], pks[gindx], sigma0, 0]
    slice = (maximum([-widx+pks[gindx],1]):minimum([widx+pks[gindx],length(y)]))
    msknan = .!isnan.(y[slice])
    try
        fit_result = curve_fit(gaussian, slice[msknan], y[slice][msknan], 1 ./yerr[slice][msknan].^2, initial_guess)
        params_opt = coef(fit_result)
        return [params_opt..., sum((fit_result.resid).^2)/count(msknan), count(msknan), fibindx-1] # subtracting 1 to match the DRP python indexing, can remove if I arclamps move to this ecosystem
    catch
        # this catches singular expceptions that I don't fully understand
        return [NaN, NaN, NaN, NaN, NaN, 0, fibindx-1]
    end
end

# Define the Gaussian function with constant offset
function gaussian(x, p)
    return p[1] * exp.(-((x .- p[2]).^2) / (2 * p[3]^2)) .+ p[4]
end

function save_FPIpeaks(fname,run_tup,linesOut;write_mode="w")
    hdr = FITSHeader(["pipeline","git_branch","git_commit","wave_expid"],
        ["apMADGICS.jl",git_branch,git_commit,string(run_tup[end])],
        ["","","",""])
    f = FITS(fname,write_mode)
    try
        write(f,[0],header=hdr,name="header_only")
        colNames = ["amp", "pix_cen", "sigma", "offset", "chi2_dof", "num_data", "fiberindx"]
        colVals = Any[]
        for i=1:7
            if (i==6) | (i==7)
                push!(colVals,Int.(linesOut[i,:]))
            else
                push!(colVals,linesOut[i,:])
            end
        end
        write(f,colNames,colVals,name="gauss_FPI_params")
    finally
        close(f)
    end
end