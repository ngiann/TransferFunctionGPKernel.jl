using HCubature, SimpleLinearPiecewiseApproximation, TransferFunctions, Distributions, Printf

using SpecialFunctions, Interpolations
include("integrals.jl")

include("GKGT.jl")
include("convkernel.jl")
include("convkernel_cuba.jl")
