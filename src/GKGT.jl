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

    support = LinRange(-50.0, 50.0, 271) # this must also be looked at again

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
