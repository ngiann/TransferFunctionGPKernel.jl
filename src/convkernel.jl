function convkernel(; Xₖ = Xₖ, Xₗ = Xₗ, tᵢ = tᵢ, tⱼ = tⱼ, ℓ = ℓ)

    length(Xₖ) == length(Xₗ) ? nothing : error("Xₖ and Xₗ must be of same length")

    resultapprox = 0.0

    numberofsegments = length(Xₖ)-1

    @inbounds for s in 1:numberofsegments

        # retrieve parameters of s-th segment

        a₁ˢ, a₂ˢ = Xₗ[s][1], Xₗ[s][2]

        b₁ˢ, b₂ˢ = Xₗ[s+1][1], Xₗ[s+1][2]

        γˢₗ = (a₂ˢ - b₂ˢ) / (a₁ˢ - b₁ˢ)

        ψˢₗⱼ = γˢₗ * (tⱼ - a₁ˢ) + a₂ˢ


        for r in 1:numberofsegments

            # retrieve parameters of r-th segment

            a₁ʳ, a₂ʳ = Xₖ[r][1], Xₖ[r][2]

            b₁ʳ, b₂ʳ = Xₖ[r+1][1], Xₖ[r+1][2]

            γʳₖ  = (a₂ʳ - b₂ʳ) / (a₁ʳ - b₁ʳ)

            ψʳₖᵢ = γʳₖ * (tᵢ - a₁ʳ) + a₂ʳ

            # Important: the original function g̃ʳₖ(t) has support [aʳₖ₁, bʳₖ₁].
            # The reflected function h(-t) has support [-bʳₖ₁, -aʳₖ₁].
            # The shifted function in the convolution h(tᵢ-t) has support [-bʳₖ₁+tᵢ, -aʳₖ₁+tᵢ].
            # Hence, the integration bounds are [-bʳₖ₁+tᵢ, -aʳₖ₁+tᵢ]
            # i.e. lowerx is -bʳₖ₁ + tᵢ as opposed to aʳₖ₁
            # and  upperx = -aʳₖ₁ + tᵢ as opposed to bʳₖ₁

            resultapprox += ψˢₗⱼ  * ψʳₖᵢ * integralfunction_1(; lowerx = -b₁ʳ + tᵢ, upperx = -a₁ʳ + tᵢ, lowery = -b₁ˢ + tⱼ, uppery = -a₁ˢ + tⱼ, ℓ = ℓ) -
                             γʳₖ  * ψˢₗⱼ * integralfunction_2(; lowerx = -b₁ʳ + tᵢ, upperx = -a₁ʳ + tᵢ, lowery = -b₁ˢ + tⱼ, uppery = -a₁ˢ + tⱼ, ℓ = ℓ) -
                             γˢₗ  * ψʳₖᵢ * integralfunction_2(; lowerx = -b₁ˢ + tⱼ, upperx = -a₁ˢ + tⱼ, lowery = -b₁ʳ + tᵢ, uppery = -a₁ʳ + tᵢ, ℓ = ℓ) +
                             γˢₗ  * γʳₖ  * integralfunction_3(; lowerx = -b₁ʳ + tᵢ, upperx = -a₁ʳ + tᵢ, lowery = -b₁ˢ + tⱼ, uppery = -a₁ˢ + tⱼ, ℓ = ℓ)

        end
    end

    resultapprox

end




function test_approx_physical(tᵢ, tⱼ, ℓ; mass=1e8, accretion=0.5, wavelengths = [3000.0;5000.0])

    # define transfer function
    gₖ = PhysicalTransferFunction(mass=mass, accretion=accretion, wavelength=wavelengths[1])
    gₗ = PhysicalTransferFunction(mass=mass, accretion=accretion, wavelength=wavelengths[2])

    # discretise transfer functions
    Xₖ = pw(collect(0.0:0.01:50.0), gₖ.(collect(0.0:0.01:50.0)), 117)
    Xₗ = pw(collect(0.0:0.01:50.0), gₗ.(collect(0.0:0.01:50.0)), 117)

    convkernel(; Xₖ = Xₖ, Xₗ = Xₗ, tᵢ = tᵢ, tⱼ = tⱼ, ℓ = ℓ)

end
