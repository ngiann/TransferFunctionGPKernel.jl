using Optim, PyPlot

@everywhere using LinearAlgebra, Cuba, ProgressMeter, TransferFunctions


@everywhere κ(d, r) = exp(- d^2 / (2*r^2))


@everywhere function symmetricconvolution_cuba(; g = g, κ = κ, tᵢ = tᵢ, tⱼ = tⱼ)

    offset = 35

    targetintegral(t) = g(tⱼ-t[2]) * g(tᵢ-t[1]) * κ(t[1]-t[2])

    trans(x) = -offset + (2*offset)*x

    auxiliary(t) = ((2*offset)^2)*targetintegral(trans.(t))

    resultquad = divonne( (t,out) -> out[1] = auxiliary(t), 2, 1, maxevals = 1_000_000, rtol=1e-6)

    resultquad.integral[1]

end



function test_symmetricconvolution_cuba(tᵢ, tⱼ, r; mass=1e8, accretion=0.5, wavelength = 5000.0)

    # define transfer function

    g = PhysicalTransferFunction(mass=mass, accretion=accretion, wavelength=wavelength)

    symmetricconvolution_cuba(; g = g, κ = d -> κ(d, r), tᵢ = tᵢ, tⱼ = tⱼ)

end


#-------------------------------------------------------------------------------


@everywhere function nonsymmetricconvolution(; gₖ = gₖ, gₗ = gₗ, κ = κ, tᵢ = tᵢ, tⱼ = tⱼ)

    offset = 35

    targetintegral(t) = gₖ(tᵢ-t[2]) * gₗ(tⱼ-t[1]) * κ(t[1]-t[2])

    trans(x) = -offset + (2*offset)*x

    auxiliary(t) = ((2*offset)^2)*targetintegral(trans.(t))

    resultquad = divonne( (t,out) -> out[1] = auxiliary(t), 2, 1, maxevals = 1_000_000, rtol=1e-6)

    resultquad.integral[1]

end



function test_nonsymmetricconvolution(tᵢ, tⱼ, r; mass=1e8, accretion=0.5, wavelengths = [3000.0; 5000.0])

    # define transfer functions

    gₖ = PhysicalTransferFunction(mass=mass, accretion=accretion, wavelength=wavelengths[1])

    gₗ = PhysicalTransferFunction(mass=mass, accretion=accretion, wavelength=wavelengths[2])

    nonsymmetricconvolution(; gₖ = gₖ, gₗ = gₗ, κ = d -> κ(d, r) , tᵢ = tᵢ, tⱼ = tⱼ)

end
