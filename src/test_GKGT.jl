# ⋆ ⋆ ⋆ ⋆ ⋆ ⋆ ⋆ ⋆ #
#  Verification   #
# ⋆ ⋆ ⋆ ⋆ ⋆ ⋆ ⋆ ⋆ #

function test_GKGT_random( ; tol=1e-1)

    t1, t2    = randn(2)*3
    ℓ         = rand()*(4.95 - 0.05) + 0.05
    w1        = rand()*(10_000 - 3_000) + 3_000.0
    mass      = exp(rand()*(log(1e9) - log(1e6)) + log(1e6))
    accretion = rand()*0.9 + 0.05

    test_GKGT(; t1=t1, t2=t2, ℓ=ℓ, w1=w1, mass=mass,
                accretion=accretion, tol=tol)

end

function test_GKGT(;t1=t1,t2=t2,ℓ=ℓ,w1=w1,mass=mass,accretion=accretion, tol=tol)

    @printf("times: t1 = %.4f, t2 = %.4f\n", t1, t2)
    @printf("ℓ = %.4f\n", ℓ)
    @printf("wavelength = %.4f\n", w1)
    @printf("mass = %.4e\n", mass)
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
