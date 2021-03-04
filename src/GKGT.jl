########################
#   Type Definition    #
########################

abstract type AbstractGKGᵀ end

mutable struct symGKGᵀ <: AbstractGKGᵀ

    support::LinRange{Float64}

    ℓ::Float64

    g::PhysicalTransferFunction

    X::Array{Array{Float64, 1}}

    mass::Float64

    accretion::Float64

    wavelength::Float64

    f::AbstractExtrapolation

    parallelevaluation::Bool

end



##############
#    Info    #
##############

function Base.show(io::IO, k::symGKGᵀ)


    print(io, "Symmetric GKGᵀ product\n",
              @sprintf("parallel evaluation set to %s\n", k.parallelevaluation),
              @sprintf("mass is %e\n", k.mass),
              @sprintf("accretion is %2.2f\n", k.accretion),
              @sprintf("wavelength is %4.1f", k.wavelength))

end



###############
# Constructor #
###############

function symGKGᵀ(mass::Float64=1e8, accretion::Float64=0.5, wavelength::Float64=5000.0, parallelevaluation::Bool=false)

    gₖ = PhysicalTransferFunction(mass=mass, accretion=accretion, wavelength=wavelength)

    numberofsegments = 31 # crucial parameter of approximation

    Xₖ = pw(collect(0.0:0.01:width(gₖ)), gₖ.(collect(0.0:0.01:width(gₖ))), numberofsegments)

    support = LinRange(-30.0, 30.0, 251) # this must also be looked at again

    y = pmap(t -> convkernel(; Xₖ = Xₖ, Xₗ = Xₖ, tᵢ = t, tⱼ = 0.0, ℓ = 1.0), support, distributed = parallelevaluation)

    symGKGᵀ(support, 1.0, gₖ, Xₖ, mass, accretion, wavelength, CubicSplineInterpolation(support, y, extrapolation_bc=0.0), parallelevaluation)

end



##############
# Evaluation #
##############

function (k::symGKGᵀ)(tᵢ, tⱼ, ℓ)

    if ℓ !== k.ℓ # in this case recompute

        k.ℓ = ℓ

        y = pmap(t -> convkernel(; Xₖ = k.X, Xₗ = k.X, tᵢ = t, tⱼ = 0.0, ℓ = k.ℓ), k.support, distributed = k.parallelevaluation)

        k.f = CubicSplineInterpolation(k.support, y, extrapolation_bc=0.0)

    end

    k.f(tᵢ - tⱼ)

end



# ⋆ ⋆ ⋆ ⋆ ⋆ ⋆ ⋆ ⋆ #
#  Verification   #
# ⋆ ⋆ ⋆ ⋆ ⋆ ⋆ ⋆ ⋆ #

function test_GKGT()

    t1, t2    = randn(2)*3
    ℓ         = rand()*5
    w1        = rand()*7000+3000
    mass      = rand()*(9e9-1e6) + 1e6
    accretion = rand()*0.9 + 0.05


    @printf("times: t1 = %.4f, t2 = %.4f\n", t1, t2)
    @printf("ℓ = %.4f\n", ℓ)
    @printf("wavelength = %.4f\n", w1)
    @printf("mass = %.4f\n", mass)
    @printf("accretion = %.4f\n", accretion)


    k1 = symGKGᵀ(mass,accretion, w1, true)
    k2 = symGKGᵀ(mass,accretion, w1, false)

    # all results obtained below must agree with each other

    a1 = k1(t1, t2, ℓ)

    a11 = k2(t2, t1, ℓ)

    a2 = test_approx_physical(t1, t2, ℓ; mass = mass, accretion = accretion, wavelengths=[w1; w1])

    a3 = test_nonsymmetricconvolution(t1, t2, ℓ; mass = mass, accretion = accretion, wavelengths=[w1; w1])

    a4 = test_symmetricconvolution_cuba(t1, t2, ℓ; mass = mass, accretion = accretion, wavelength = w1)

    display([a1; a11; a2; a3; a4])

    tol = 1e-2

    (abs(a1 - a11) < tol) &&
    (abs(a1 - a2)  < tol) &&
    (abs(a1 - a3)  < tol) &&
    (abs(a1 - a4)  < tol) &&
    (abs(a2 - a11) < tol) &&
    (abs(a2 - a3)  < tol) &&
    (abs(a2 - a4)  < tol) &&
    (abs(a3 - a11) < tol) &&
    (abs(a3 -  a4) < tol)

end
