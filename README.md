# apMADGICS

[![Build Status](https://github.com/andrew-saydjari/apMADGICS.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/andrew-saydjari/apMADGICS.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/andrew-saydjari/apMADGICS.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/andrew-saydjari/apMADGICS.jl)


Pipeline for APOGEE spectra using Marginalized Analytic Dataspace Gaussian Inference for Component Separation (MADGICS).

## Installation

This is a pipeline. Not a package. It is not really meant to be installed. The pipeline.jl script is meant to be run. This repo documents the development and versions of the pipeline for transparency and reproducibility.

If you wish to download the code and have the dependencies required to run it installed, you can install directly from GitHub. 

```julia
import Pkg
Pkg.add(url="https://github.com/andrew-saydjari/apMADGICS.jl")
```

## gridSearch Module Flag Bits

There is still a (much smaller dimensional) space that MADGICS needs to sample over (e.g. radial velocity). We have a custom grid-sampler module to implement that sampling. The flag bits from that module are below.

| Value         | Bit         | Meaning     |
| ----------- | ----------- | ----------- |
| 0     | -     | No problems       |
| 1     | 0     | Interpolated minimum not less than minimum (should not occur) |
| 2     | 1     | Minimum index at edge of grid for dimension 1 |
| 4     | 2     | Minimum index at edge of grid for dimension 2 |
| 8     | 3     | Finite difference Hessian beyond grid edge |
| 16    | 4     | Bad curvature of chi2 surface (can't invert full 2d Hessian)|
| 32    | 5     | Very bad curvature of chi2 surface (can't invert diagonal entries)|


