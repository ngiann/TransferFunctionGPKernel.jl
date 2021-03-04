module TransferFunctionGPKernel

    export symGKGᵀ, testme

    using HCubature, Cuba, SimpleLinearPiecewiseApproximation, TransferFunctions, Printf
    using Distributed

    using SpecialFunctions, Interpolations

    include("integrals.jl")

    include("GKGT.jl")
    include("convkernel.jl")
    include("convkernel_cuba.jl")
    include("testme.jl")

end
