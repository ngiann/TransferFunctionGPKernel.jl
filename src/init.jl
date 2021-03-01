using HCubature, SimpleLinearPiecewiseApproximation, TransferFunctions, Distributions, Printf

@everywhere using SpecialFunctions, Interpolations
@everywhere include("integrals.jl")

include("GKGT.jl")
include("convkernel.jl")
include("convkernel_cuba.jl")
