using HCubature, Cuba, SimpleLinearPiecewiseApproximation, TransferFunctions, Printf

using SpecialFunctions, Interpolations
include("integrals.jl")

include("GKGT.jl")
include("convkernel.jl")
include("convkernel_cuba.jl")
