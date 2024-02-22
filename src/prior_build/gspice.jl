module gspice

# There is a bug in Julia 1.6.4 that requires BLAS.set_num_threads(1)
using LinearAlgebra
# BLAS.set_num_threads(1)

using ShiftedArrays
using Interpolations
using Statistics
using Dates
using FITSIO

export djs_maskinterp
export gspice_standard_scale
export gspice_covar
export gspice_submatrix_inv
export gspice_submatrix_inv_mult
export gspice_gaussian_estimate
export gspice_gp_interp
export gspice_chimask
export gspice_covar_iter_mask


function where(cond)
    ind = findall(cond)
    num = length(ind)
    return ind, num
end


"""Interpolate over a masked, 1-d array.

Parameters
----------
yval : :class:`numpy.ndarray`
    The input values.
mask : :class:`numpy.ndarray`
    The mask (0=good, 1=bad)
xval : :class:`numpy.ndarray`, optional
    If set, use these x values, otherwise use an array.
const : :class:`bool`, optional
    If set to ``True``, bad values around the edges of the array will be
    set to a constant value.  Because of the default behavior of
    :func:`numpy.interp`, this value actually makes no difference in
    the output.

Returns
-------
:class:`numpy.ndarray`
    The `yval` array with masked values replaced by interpolated values.
"""
function djs_maskinterp1(yval, mask, extrap_scheme; xval=nothing)
    igood = findall(mask .==0)
    ngood = length(igood)
    ny    = length(yval)
    ynew  = copy(yval)
    if ngood == ny return ynew end   # if all good
    if ngood == 0 return ynew end    # if none good
    if ngood == 1 return ynew end    # if only one is good (IDL behavior)

# this is Python djs_maskinterp behavior:
#    if ngood == 1
#        igood = findall(mask .==0)
#        return fill(yval[igood[1]], ny)
#    end


    if isnothing(xval)
        ibad  = findall(mask .!= 0)
        itp = LinearInterpolation(igood, ynew[igood]; extrapolation_bc=extrap_scheme)
        ynew[ibad] = itp(ibad)
    else
        ii    = sortperm(xval)
        igood = findall(mask[ii] .==0)
        ibad  = findall(mask[ii] .!= 0)
        itp=LinearInterpolation(xval[ii[igood]], ynew[ii[igood]]; extrapolation_bc=extrap_scheme)
        ynew[ii[ibad]] = itp(xval[ii[ibad]])
    end

    return ynew
end


"""Interpolate over masked pixels in a vector, image or 3-D array.

Parameters
----------
yval : :class:`numpy.ndarray`
    The input values
mask : :class:`numpy.ndarray`
    The mask (0=good, 1=bad)
xval : :class:`numpy.ndarray`, optional
    If set, use these x values, otherwise use an array
axis : :class:`int`, optional
    Must be set if yval has more than one dimension. If set to zero,
    interpolate along the first axis of the array, if set to one,
    interpolate along the second axis of the array, and so on.
const : :class:`bool`, optional
    This value is passed to a helper function, djs_maskinterp1.

Returns
-------
:class:`numpy.ndarray`
    The interpolated array.
"""
function djs_maskinterp(yval, mask; xval=nothing, axis=nothing, constant=false)
    @assert size(mask) == size(yval) "mask must have the same shape as yval."
    if ~isnothing(xval)
        @assert size(mask) == size(yval) "xval must have the same shape as yval."
    end
    sz = size(yval)
    ndim = length(sz)
    ext_scheme = constant ? Flat() : Line()

    if ndim == 1
        ynew = djs_maskinterp1(yval, mask, ext_scheme; xval=xval)
    else
        if isnothing(axis)
            throw("Must set axis if yval has more than one dimension.")
        end
        if (axis < 0) || (axis > ndim-1) || (axis - round(axis)) != 0
            throw("Invalid axis value.")
        end
        ynew = similar(yval)
        if ndim == 2
            if isnothing(xval)
                if axis == 0
                    @views for i in 1:sz[1]
                        ynew[i, :] = djs_maskinterp1(yval[i, :], mask[i, :],
                                                     ext_scheme)
                    end
                else
                    @views for i in 1:sz[2]
                        ynew[:, i] = djs_maskinterp1(yval[:, i], mask[:, i],
                                                     ext_scheme)
                    end
                end
            else
                if axis == 0
                    @views for i in 1:sz[1]
                        ynew[i, :] = djs_maskinterp1(yval[i, :], mask[i, :], ext_scheme,
                                                     xval=xval[i, :])
                    end
                else
                    @views for i in 1:sz[2]
                        ynew[:, i] = djs_maskinterp1(yval[:, i], mask[:, i], ext_scheme,
                                                     xval=xval[:, i])
                    end
                end
            end
        end
    end
    return ynew
end


"""
    gspice_standard_scale(flux, ivar, mask=nothing) -> Dvec, refscale, refmean
Standard scale the data after interpolation over a mask

# Arguments:
- `flux`:      flux array (nspec, npix)
- `ivar`:      ivar array (nspec, npix)

# Keywords:
- `mask`:      mask

# Output:
- `Dvec`:
- `refscale`:  factor each spectrum was multiplied by (nspec)
- `refmean`:   mean of each rescaled wavelength bin (npix)


# Comments:
2021-Oct-21 - Written by Douglas Finkbeiner, CfA
"""
function gspice_standard_scale(flux::Matrix{Float64}, ivar::Matrix{Float64}, mask=nothing)
    # -------- if no mask is passed, use ivar=0 as mask
    pixmask = isnothing(mask) ? ivar .== 0 : mask .!= 0

    # -------- interpolate over masked pixels in the spectral direction
    Dvec = djs_maskinterp(flux, pixmask; axis=0, constant=true)

    # -------- renormalize each spectrum by sqrt(mean ivar)
    wt = .~pixmask
    meanivar = sum(ivar.*wt, dims=2)./sum(wt, dims=2)
    refscale = sqrt.(meanivar)
    Dvec .*= refscale

    Nspec, Npix = size(Dvec)
    refmean = sum(Dvec, dims=1) ./ Nspec
    Dvec = Dvec .- refmean

    return Dvec, refscale, refmean
end


"""
    gspice_covar(spec::Array{Float64,2}; checkmean=false) -> cov, refmean
Compute covariance of spec

# Arguments:
- `spec`:      spectral data array (npix, nspec), must be float64

# Keywords:
- `checkmean`: set to check precision of mean subtraction

# Output:
- `cov`: covariance matrix
- `refmean`: mean of each column (wavelength bin)

# Comments:
All masking of spec must done before calling this. \\
2021-Oct-21 - Written by Douglas Finkbeiner, CfA
"""
function gspice_covar(spec::Array{Float64,2}; checkmean=false)

    Nspec, Npix = size(spec)

    # -------- make columns of spec array mean zero for computation of covariance
    refmean = sum(spec, dims=1) ./ Nspec
    spec0 = spec .- refmean

    # -------- verify mean subtraction worked
    if checkmean
        mnd = sum(spec0, dims=1) ./ Nspec
        println("Min, max, mean ", extrema(mnd), std(mnd))
        flush(stdout)
    end

    # -------- compute covariance
    cov = (spec0'*spec0)./(Nspec-1)
    if !issymmetric(cov)
        println("Warning: Had to symmetrize the matrix!")
        flush(stdout)
        cov .+= cov'
        cov ./=2
    end
    return cov, refmean
end


"""
    gspice_submatrix_inv(M, Minv, imask; bruteforce=false) -> Ainv
Given M and Minv, compute inverse of a submatrix of M

# Arguments:
- `M`:     (Nd,Nd) parent matrix, assumed symmetric and positive semi-definite.
- `Minv`:  inverse of M, usually calculated with Cholesky.
- `imask`: mask of rows/columns to use (1=keep, 0=remove). Contains Nk ones and Nr zeros.

# Keywords:
- `bruteforce`: brute force calculation for testing

# Output:
- `Ainv`: Inverse of the submatrix of M

# Comments:
2021-Oct-21 - Written by Douglas Finkbeiner, CfA
"""
function gspice_submatrix_inv(M, Minv, imask; bruteforce=false)
    # REVISION HISTORY:
    #   2019-Oct-20 - Written by Douglas Finkbeiner, CfA   (At MPIA)
    #   2019-Nov-11 - More compact form for A inverse
    #   2021-Nov   Translated to Julia

    # -------- check dimensions of inputs
    Nd = length(imask)
    @assert Nd == size(M,1)

    # -------- brute force option, for testing
    if bruteforce
        k    = findall(imask)
        Ainv = inv(cholesky(Symmetric(M[k, k])))
        return Ainv
    end

    # -------- index list to remove
    r   = findall(.~imask)
    nr  = length(r)

    # -------- if there is nothing to do, return early
    if nr==0
      println("imask does not remove any rows/columns...")
      flush(stdout)
      return Minv
    end

    # -------- Evaluate  A^{-1} = P - Q U^{-1} Q^T
    k    = findall(imask)
    Uinv = inv(Minv[r, r])
    Q    = Minv[k, r]

    Ainv = Minv[k, k] - (Q * (Uinv * Q'))

    return Ainv
end


"""
    gspice_submatrix_inv_mult(M, Minv, imask, Y, MinvY; irange=nothing, pad=false, bruteforce=false) -> Ainvy
Given M, Minv, and MinvY, compute inverse of submatrix of M, times Y

# Arguments:
- `M`:     (Nd,Nd) parent matrix, assumed symmetric and positive semi-definite.
- `Minv`:  inverse of M, usually calculated with Cholesky.
- `imask`: mask of rows/columns to use (1=keep, 0=remove). Contains Nk ones and Nr zeros.
- `Y`:     Matrix (Nspec,Nd) to multiply Ainv by (Nspec can be 1)
- `MinvY`: Minv times Y

# Keywords:
- `irange`: if imask=0 is contiguous, just give endpoints (faster)
- `pad`: zero-pad removed rows of Ainvy to Nd dimensions
- `bruteforce`: brute force calculation for testing

# Output:
- `Ainvy`: Inverse A (submatrix of M) times Y (Nspec, Nd-Nr). If pad==true, then zero-padded to (Nspec, Nd).

# Comments:
2021-Oct-21 - Written by Douglas Finkbeiner, CfA
"""
function gspice_submatrix_inv_mult(M, Minv, imask, Y, MinvY; irange=nothing, pad=false, bruteforce=false)
    # COMMENTS:
    #   Let M by a block matrix
    #
    #          | A   B |                      | P   Q |
    #     M  = |       |      and    M^{-1} = |       |
    #          | B^T D |                      | Q^T U |
    #
    #   Then the inverse of submatrix A is the Schur complement of U,
    #
    #     A^{-1} = P - Q U^{-1} Q^T
    #
    #   U and M must be invertible and positive semi-definite.
    #   This function returns
    #
    #     A^{-1} Y = P Y - Q U^{-1} Q^T Y
    #

    # -------- catch the case where Y is a row vector
    @assert size(Y,1) !== 1     # 'Y must not be a row vector'
    @assert size(MinvY) == size(Y)

    # -------- check dimensions of inputs
    Nd = length(imask)
    @assert Nd == size(M,1)

    # -------- brute force option (Slow - use only for testing)
    if bruteforce
        k     = findall(imask)
        Ainv  = inv(cholesky(Symmetric(M[k, k])))
        Ainvy = Ainv * Y[k, :]
        return Ainvy
    end

    # -------- index list to remove
    r   = findall(.~imask)
    nr  = length(r)

    # -------- if there is nothing to do, return early
    if nr==0
      println("imask does not remove any rows/columns...")
      return MinvY
    end

    k    = findall(imask)
    Q    = Minv[k, r]
    # -------- Use that Qty = Minv y - U y.
    U    = Symmetric(Minv[r, r])
    Yr   = Y[r, :]
    Qty  = MinvY[r, :] - U * Yr

    # -------- evaluate  A^{-1} Y = P Y - Q U^{-1} Q' Y
    if length(U) == 1
        UinvQtY = Qty/U[1]
    else
        UinvQtY = cholesky(U)\Qty
    end

    # -------- Evaluate   A^{-1} Y = P Y - Q U^{-1} Q^T Y
    #          with some shortcuts, using P Y = Minv Y - Q Y
    if pad
        AinvY0 = copy(MinvY)
        AinvY0[k, :] -= (Q*(UinvQtY+Yr))
        AinvY0[r, :] .= 0
        return AinvY0
    end

    AinvY = MinvY[k, :] - (Q*(UinvQtY+Yr))
    return AinvY
end


"""
    gspice_gaussian_estimate(icond, ipred, covmat, Dvec, covinv; bruteforce=false) -> predoverD, predcovar, predkstar, kstar
Compute Gaussian conditional esimtate

# Arguments:
- `icond`:      array (npix) 1 = pixels to condition on
- `ipred`:      array (npix) 1 = pixels to predict
- `covmat`:     covariance matrix (npix,npix)
- `Dvec`:       data vector (nspec, npix)
- `covinv`:     inverse covariance (to avoid computing again)
- `bruteforce`: brute force calculation for testing

# Output:
- `predoverD`:  matrix to multiply by Dvec to obtain result
- `predcovar`:  predicted covariance
- `predkstar`:  prediction  (bruteforce code only)
- `kstar`:      return kstar index list for debugging

# Comments:
- ipred is 1 for pixels to be predicted, conditional on the reference pixels specified by icond
- In normal operation, this runs with bruteforce=false and \\
   returns predoverD.  The desired result is predoverD*Dvec.  Doing \\
   that as one matrix multiplication is much faster than multiplying \\
   each time in the loop. \\
 \\
2021-Oct-21 - Written by Douglas Finkbeiner, CfA \\
"""
function gspice_gaussian_estimate(icond, ipred, covmat, Dvec, covinv; bruteforce=false)
    # -------- get index lists for reference (k) and interpolation inds
    #          (kstar), following notation of RW Chapter 2 ~ Eq. 2.38
    sz     = size(Dvec)
    single = length(sz) == 1
    predoverD = nothing
    predkstar = nothing

    k, nk         = where(icond .== 1)    # where you have data
    kstar, nkstar = where(ipred .== 1)    # where you want to predict

    if bruteforce  # use old code for testing
        cov_kk         = covmat[k, k]
        # RENAME ???
        cov_kkstar     = covmat[kstar, k] # dim [nkstar,nk]
        cov_kstark     = covmat[k, kstar]
        cov_kstarkstar = covmat[kstar, kstar]

        # -------- Choleksy inversion
        # icov_kk = inv(cholesky(Symmetric(cov_kk)))
        icov_kk = inv(Symmetric(cov_kk))

        # -------- compute the prediction covariance (See RW, Chap. 2)
        predcovar = cov_kstarkstar - (cov_kkstar*(icov_kk*cov_kstark))

        if single       # multiply parts 2 and 3 first
            predkstar = cov_kkstar * (icov_kk * Dvec[k])
        else
            # -------- compute icov_kk * cov_kstark using GSPICE routine
            temp = cov_kkstar * icov_kk
            icov_kk = 0

            if nkstar==1
                temp2 = zeros(sz[2])
                temp2[k] = temp
                predkstar = Dvec * temp2
            else
                predkstar = Dvec[:, k] * temp' # this takes memory
                println("Warning: Using memory intensive code....")
            end
        end
    else             # -------- GSPICE version

        cov_kkstar     = covmat[kstar, k] # dim [nkstar,nk]
        cov_kstark     = covmat[k, kstar]
        cov_kstarkstar = covmat[kstar, kstar]

        # -------- compute icov_kk * cov_kstark using GSPICE routine
        # could set Minvy for a slight speedup
        Y = covmat[:, kstar]
        Minvy = covinv*Y

        Ainvy0 = gspice_submatrix_inv_mult(covmat, covinv, icond, Y, Minvy, pad=true)
                                # Ainvy is icov_kk * cov_kstark
                                # Ainvy0 is that zero padded

        #     predkstar = matrixmult(Dvec, Ainvy0)  ; this takes all the time.
        predoverD = Ainvy0

        # -------- compute the prediction covariance (See RW, Chap. 2)
        predcovar = cov_kstarkstar - Y' *Ainvy0

    end

    return predoverD, predcovar, predkstar, kstar
end


"""
    gspice_gp_interp(Dvec, covmat; irange=nothing, nguard=20, bruteforce=false) -> pred, predvar
Perform GSPICE prediction of each pixel, conditional on others

# Arguments:
- `Dvec`:     data vector, already cleaned and scaled
- `covmat`:   covariance to use for GCE prediction
- `irange`:   range of spectral indices to compute (default to all)
- `nguard`:   number of guard pixels around GCE pixel
- `bruteforce`:  use bruteforce code for debugging

# Output:
- `pred`:     predicted poste rior mean
- `predvar`:  predicted posterior variance

# Comments:
- Computes the Gaussian Conditional Estimate (GCE) of each pixel,
   conditional on the other pixels except those within radius nguard.
 \\
2021-Oct-21 - Written by Douglas Finkbeiner, CfA \\
"""
function gspice_gp_interp(Dvec, covmat; irange=nothing, nguard=20, bruteforce=false, reg_eps=0)

    @assert ~isnothing(covmat)

    # -------- shape of input array
    nspec, npix = size(Dvec)

    # -------- range of spectral bins to operate on
    i1,i2 = isnothing(irange) ? [1,npix] : irange

    # -------- allocate output arrays
    szpred    = (i2-i1)+1
    predvar   = zeros(nspec, szpred)
    if bruteforce
        pred  = zeros(nspec, szpred)
    else
        predoverD = zeros(npix, szpred)
    end

    t0 = now()

    # -------- pre-compute inverse covariance
    covinv = inv(cholesky(Symmetric(covmat+reg_eps*I)))

    # -------- loop over spectral pixels
    for i=i1:i2

        # -------- ipred is 1 for pixels to be predicted
        ipred = falses(npix)
        ipred[i] = true

        # -------- conditional on the reference pixels specified by icond
        icond = trues(npix)
        icond[max(1,i-nguard):min(i+nguard,npix)] .= false

        predoverD0, predcovar, predkstar, kstar =
                gspice_gaussian_estimate(icond, ipred, covmat, Dvec, covinv; bruteforce=bruteforce)
        if bruteforce
            pred[:, kstar.-(i1-1)] = predkstar
        else
            predoverD[:, kstar.-(i1-1)] = predoverD0
        end
        predvar[:, kstar.-(i1-1)] .= diag(predcovar)

        if mod(i,500)==0
            println(i, "  ",(now()-t0).value/1000.0, " sec")
            flush(stdout)
        end
    end

    if ~bruteforce pred = Dvec*predoverD end

    return pred, predvar
end


"""
    gspice_chimask(flux, ivar, mask, nsigma) -> chimask
Compute mask of outliers with respect to GSPICE posterior variance

# Arguments:
- `flux`:      flux array (nspec, npix)
- `ivar`:      ivar array (nspec, npix)
- `mask`:      input mask (0=good).  If not passed, use ivar==0
- `nsigma`:    number of sigma to cut at.

# Output:
- `chimask`:   mask of nonconformant pixels (0=good)

# Comments:
- This routine standard scales the input so each spectrum has the same \\
   mean(ivar), interpolates over bad pixels (from input mask), then \\
   computes a covariance. \\
- The Gaussian pixelwise estimate based on that covariance is compared \\
   to the data, and a Z-score ("chi") is computed using the GSPICE \\
   posterior variance. The output mask is based on that Zscore, and \\
   set to 1 for abs(Z) > nsigma.  \\
- Finally the mask is dilated by marking any pixel on either side of \\
   a bad pixel as also bad.  \\
 \\
2021-Oct-21 - Written by Douglas Finkbeiner, CfA \\
"""
function gspice_chimask(flux, ivar, mask, nsigma; reg_eps=0)

    # -------- interpolate masked pixels in the spectral direction, scale
    Dvec, refscale, refmean = gspice_standard_scale(flux, ivar, mask)

    # -------- obtain the empirical covariance for this Dvec
    covmat, __ = gspice_covar(Dvec)

    # -------- compute GSPICE predicted mean and variance
    pred, predvar = gspice_gp_interp(Dvec, covmat, nguard=20,reg_eps=reg_eps)

    # -------- compute "chi" (really the Z-score)
    chi = (Dvec-pred)./sqrt.(predvar)

    # -------- clip at nsigma, dilate the mask by 1 pixel
    chim = (abs.(chi) .> nsigma)
    chimask = chim .| ShiftedArray(chim,(0,1),default=false) .|
                      ShiftedArray(chim,(0,-1),default=false)

    return chimask
end


"""
    gspice_covar_iter_mask(flux, ivar, mask; nsigma=[20, 8, 6], maxbadpix=64) -> covmat, finalmask
Compute spectral covariance with iterative masking using GSPICE

# Arguments:
- `flux`:      flux array (nspec, npix)
- `ivar`:      ivar array (nspec, npix)
- `mask`:      input mask (0=good).  If not passed, use ivar==0
- `nsigma`:    number of sigma to cut at.
- `maxbadpix`: reject objects with more than maxbadpix pixels masked
- `maxbadpix`: undersampling factor allowed before complain too few samples

# Output:
- `covmat`:    covariance matrix
- `finalmask`: final mask

# Comments:
- Iteratively mask estimates of the covariance matrix
   Pass this routine de-redshifted spectra with mask of definitely bad pixels (1=bad)

2021-Oct-21 - Written by Douglas Finkbeiner, CfA \\
"""
function gspice_covar_iter_mask(flux, ivar, mask; nsigma=[20, 8, 6], maxbadpix=64, usamp_factor=1, reg_eps=0)

    t0 = now()                     # start time
    nspec,npix = size(flux)

    # -------- reject spectra with too many bad pixels
    nbadpix = sum(mask .!= 0, dims=2)[:,1]     # bad pixels for each spectrum
    objmask = nbadpix .<= maxbadpix    # 0=bad (too many pixels masked)
    wmask,nmask = where(objmask)   # index list of spectra to use

    if nmask < npix/usamp_factor throw("Not enough ($(nmask)) good spectra! But you want to model $(npix) pixels!") end

    chimask = false
    thismask = false
    for iter = 1:length(nsigma)
        println("=========================  Pass ", iter, ",   cut at Nsigma = ", nsigma[iter])
        flush(stdout)
        thismask = chimask .| (mask[wmask,:] .!= 0)
        chimask = gspice_chimask(flux[wmask,:], ivar[wmask,:], thismask, nsigma[iter], reg_eps=reg_eps)
        println("mean chimask: ", mean(chimask))
        println("Time: ", (now()-t0).value/1000, " sec")
        flush(stdout)
    end

  finalmask = trues(nspec, npix)    # start from original mask, overwrite mask for good spectra
  finalmask[wmask,:] = thismask .| chimask

  Dvec, refscale, refmean = gspice_standard_scale(flux[wmask,:], ivar[wmask,:], finalmask[wmask,:])
  covmat, __ = gspice_covar(Dvec)

  return covmat, finalmask
end

end  # end of module
