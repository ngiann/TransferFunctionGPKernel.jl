
#
# Evaluates definite integral
# ∫ erf(b-y)/(sqrt(2)*ℓ), y∈[c,d]
#
# Wolfram code: Integrate[erf(b-y)/(sqrt(2)*ℓ), {y, c, d}]
#

function sub_1(; b = b, c = c, d = d, ℓ = ℓ)

    @assert(c < d)

    sqrt(2/π) * ℓ * ( exp(-(b - c)^2/(2 * ℓ^2)) - exp(-(b - d)^2/(2 * ℓ^2))) +
        (b - c) * erf((b - c)/(sqrt(2) * ℓ)) +
        (d - b) * erf((b - d)/(sqrt(2) * ℓ))


end


function verify_sub_1(; b = b, c = c, d = d, ℓ = ℓ)

    # integrand
    f(y) = erf((b-y)/(sqrt(2)*ℓ))


    # own numerical verification on 1d regular grid
    dx = 1e-6
    X  = c:dx:d
    num = 0.0
    for x in X
        num += f(x)*dx
    end


    # numerical verification via Cuba.jl library
    h = Cuba.vegas((x, out) -> out[1] = (d-c) * f(c + (d-c)*x[1]), 1, 1, maxevals=1_000_000)

    display(h)

    # numerical verification via HCubature.jl library
    hc = hquadrature(f, c, d)


    # return exact and numerical results
    sub_1(; b = b, c = c, d = d, ℓ = ℓ),
    h.integral[1],
    hc[1],
    num

end




#
# Evaluates definite integral
# ∫∫ exp(- (x-y)² / (2ℓ²) ) dx dy, x∈[a,b], y∈[c,d]
#
# Wolfram code: Integrate[E^(-(x - y)^2/(2 (ℓ^2))), {x, a, b}]
#

function integralfunction_1(; lowerx = lowerx, upperx = upperx, lowery = lowery, uppery = uppery, ℓ = ℓ)

    @assert(lowerx < upperx)
    @assert(lowery < uppery)

    # easier names
    a = lowerx
    b = upperx
    c = lowery
    d = uppery

    # Result of ∫ exp(- (x-y)² / (2ℓ²) ) dx  x∈[a,b] is:
    # sqrt(π/2) r (erf((b - y)/(sqrt(2) ℓ)) - erf((a - y)/(sqrt(2) ℓ)))
    # We need to carry out the second integral over y.
    # The result of this second integral over y is:

    sqrt(π/2) * ℓ * (sub_1(b=b, c=c, d=d, ℓ=ℓ) - sub_1(b=a, c=c, d=d, ℓ=ℓ))

end


#
# Verification: numerical evaluation for definite integral
# ∫∫ exp(- (x-y)² / (2ℓ²) ) dx dy, x∈[a, b], y∈[c, d]
#

function verify_1(; lowerx = lowerx, upperx = upperx, lowery = lowery, uppery = uppery, ℓ = ℓ)

    @assert(lowerx < upperx)
    @assert(lowery < uppery)

    # intergrand
    f(x,y) =  exp(-(x-y)^2 / (2*(ℓ^2)))


    # numerical verification on regular grid
    Δx = 1e-3
    X  = lowerx:Δx:upperx
    Y  = lowery:Δx:uppery

    num = 0.0
    for i in 1:length(X)
        for j in 1:length(Y)
            @inbounds num += f(X[i], Y[j])
        end
    end
    num *= Δx^2


    # numerical verification via Cuba.jl libary
    h = divonne((x, out) -> out[1] = (uppery-lowery)*(upperx-lowerx) * f(lowerx + (upperx-lowerx)*x[1], lowery + (uppery-lowery)*x[2]), 2, 1)

    display(h)


    # return exact and numerical results to verify
    integralfunction_1(lowerx = lowerx, upperx = upperx, lowery = lowery, uppery = uppery, ℓ = ℓ), num, h.integral[1]

end



################################################################################
################################################################################
################################################################################



#
# Evaluates definite integral
# ∫ exp(- (a-y)² / (2ℓ²) ) dy,  y ∈ [c, d]
#
function sub_2a(; a = a, c = c, d = d, ℓ = ℓ)

    # Wolfram code: Integrate[E^(-(a - y)^2/(2 ℓ^2)), {y, c, d}]
    # Solution returned in plain text:
    # integral_c^d e^(-(a - y)^2/(2 ℓ^2)) dy = sqrt(π/2) r (erf((a - c)/(sqrt(2) ℓ)) - erf((a - d)/(sqrt(2) ℓ)))

    sqrt(π/2) * ℓ * (erf((a - c)/(sqrt(2) * ℓ)) - erf((a - d)/(sqrt(2) * ℓ)))

end


#
# Verifies definite integral
# ∫ exp(- (a-y)² / (2ℓ²) ) dy,  y ∈ [c, d]
#
function verify_sub_2a(; a = a, c = c, d = d, ℓ = ℓ)

    f(y) = exp(-(a-y)^2/(2*ℓ^2))

    # own numerical verification on 1d regular grid
    dx = 1e-6
    X  = c:dx:d
    num = 0.0
    for x in X
        num += f(x)*dx
    end


    # numerical verification via Cuba.jl libary
    h = Cuba.vegas((x, out) -> out[1] = (d-c) * f(c + (d-c)*x[1]), 1, 1, maxevals=3_000_000)

    display(h)


    # return exact and numerical results
    sub_2a(; a = a, c = c, d = d, ℓ = ℓ),
    h.integral[1],
    num

end


#
# Evaluates definite integral
# ∫ y erf((a-y)/(sqrt(2)ℓ) dy,  y ∈ [c, d]
#
function sub_2b(; a = a, c = c, d = d, ℓ = ℓ)

    # Wolfram code: Integrate[y Erf[(a - y)/(Sqrt[2] ℓ)], {y, c, d}]
    # Solution returned in plain text:
    # integral_c^d y erf((a - y)/(sqrt(2) ℓ)) dy = 1/2 (a^2 - c^2 + ℓ^2) erf((a - c)/(sqrt(2) ℓ)) - 1/2 (a^2 - d^2 + ℓ^2) erf((a - d)/(sqrt(2) ℓ)) + (ℓ (a + c) e^(-(a - c)^2/(2 ℓ^2)))/sqrt(2 π) - (ℓ (a + d) e^(-(a - d)^2/(2 ℓ^2)))/sqrt(2 π)

    1/2 * (a^2 - c^2 + ℓ^2) * erf((a - c)/(sqrt(2) * ℓ)) -
    1/2 * (a^2 - d^2 + ℓ^2) * erf((a - d)/(sqrt(2) * ℓ)) +
    (ℓ * (a + c) * exp(-(a - c)^2/(2 * ℓ^2)))/sqrt(2*π) - (ℓ*(a + d) * exp(-(a - d)^2/(2*ℓ^2)))/sqrt(2π)


end


#
# Verifies definite integral
# ∫ y erf((a-y)/(sqrt(2)ℓ) dy,  y ∈ [c, d]
#
function verify_sub_2b(; a = a, c = c, d = d, ℓ = ℓ)

    f(y) = y * erf((a-y)/(sqrt(2)*ℓ))

    # own numerical verification on 1d regular grid
    dx = 1e-6
    X  = c:dx:d
    num = 0.0
    for x in X
        num += f(x)*dx
    end


    # numerical verification via Cuba.jl libary
    h = Cuba.vegas((x, out) -> out[1] = (d-c) * f(c + (d-c)*x[1]), 1, 1, maxevals=1_000_000)

    display(h)


    # return exact and numerical results
    sub_2b(; a = a, c = c, d = d, ℓ = ℓ),
    h.integral[1],
    num


end



#
# Evaluates definite integral
# ∫∫ x exp(- (x-y)² / (2ℓ²) ) dx dy, x∈[a, b], y∈[c, d]
#
function integralfunction_2(; lowerx = lowerx, upperx = upperx, lowery = lowery, uppery = uppery, ℓ = ℓ)

    @assert(lowerx < upperx)
    @assert(lowery < uppery)

    # easier names
    a = lowerx
    b = upperx
    c = lowery
    d = uppery




    # Intermediate step:
    # Solution of integral ∫ᵇₐ x exp(-(x-y)^2 / (2ℓ²)) dx is
    # integral_a^b x e^(-(x - y)^2/(2 ℓ^2)) dx = 1/2 ℓ (2 ℓ (e^(-(a - y)^2/(2 ℓ^2)) - e^(-(b - y)^2/(2 ℓ^2))) - sqrt(2 π) y erf((a - y)/(sqrt(2) ℓ)) + sqrt(2 π) y erf((b - y)/(sqrt(2) ℓ)))
    # (Wolfram code: Integrate[x * E^(-((x - y)^2)/(2*ℓ^2)), {x, a, b}]
    #

    #1/2 * ℓ *(2*ℓ * (exp(-(a - y)^2/(2*ℓ^2)) - exp(-(b - y)^2/(2*ℓ^2))) -
    #          sqrt(2π) * y * erf((a - y)/(sqrt(2)*ℓ)) + sqrt(2π) * y * erf((b - y) /(sqrt(2)*ℓ)))

    1/2 * ℓ *(2*ℓ * (sub_2a(a=a, c=c, d=d, ℓ=ℓ) - sub_2a(a=b, c=c, d=d, ℓ=ℓ)) -
              sqrt(2π) * (sub_2b(a = a, c = c, d = d, ℓ = ℓ) -  sub_2b(a = b, c = c, d = d, ℓ = ℓ)))


end



function verify_2(; lowerx = lowerx, upperx = upperx, lowery = lowery, uppery = uppery, ℓ = ℓ)

    f(x,y) =  x * exp(- (x-y)^2 / (2*ℓ^2) )

    # numerical verification on regular grid
    Δx = 1e-3
    X  = lowerx:Δx:upperx
    Y  = lowery:Δx:uppery

    num = 0.0
    for i in 1:length(X)
        for j in 1:length(Y)
            @inbounds num += f(X[i], Y[j])
        end
    end
    num *= Δx^2


    # numerical verification via Cuba.jl libary
    h = divonne((x, out) -> out[1] = (uppery-lowery)*(upperx-lowerx) * f(lowerx + (upperx-lowerx)*x[1], lowery + (uppery-lowery)*x[2]), 2, 1)

    display(h)


    # return exact and numerical results to verify
    integralfunction_2(lowerx = lowerx, upperx = upperx, lowery = lowery, uppery = uppery, ℓ = ℓ), num, h.integral[1]

end



################################################################################
################################################################################
################################################################################

#
# Evaluate integral:
# ∫ y² erf( (a-y) / (sqrt(2) ℓ) ) dy
#

function sub_3a(; a = a, c = c, d = d, ℓ = ℓ)

    # Wolram code: Integrate[y^2 Erf[(a - y)/(Sqrt[2] ℓ)], {y, c, d}]
    #
    # solution:
    # integral_c^d y y erf((a - y)/(sqrt(2) ℓ)) dy = 1/3 (-(a^3 + 3 a ℓ^2) erf((c - a)/(sqrt(2) ℓ)) + (a^3 + 3 a ℓ^2) erf((d - a)/(sqrt(2) ℓ)) + sqrt(2/π) ℓ e^(-(a - c)^2/(2 ℓ^2)) (a^2 + a c + c^2 + 2 ℓ^2) - sqrt(2/π) ℓ e^(-(a - d)^2/(2 ℓ^2)) (a^2 + a d + d^2 + 2 ℓ^2) + c^3 (-erf((a - c)/(sqrt(2) ℓ))) + d^3 erf((a - d)/(sqrt(2) ℓ)))

    1/3 * (-(a^3 + 3*a*ℓ^2) * (erf((c - a)/(sqrt(2)*ℓ)) - erf((d - a)/(sqrt(2)*ℓ))) +
              sqrt(2/π) * ℓ * (exp(-(a - c)^2/(2 * ℓ^2)) * (a^2 + a * c + c^2 + 2 * ℓ^2) -
                               exp(-(a - d)^2/(2 * ℓ^2)) * (a^2 + a * d + d^2 + 2 * ℓ^2)) +
          c^3 * (-erf((a - c)/(sqrt(2) * ℓ))) + d^3 * erf((a - d)/(sqrt(2) * ℓ)))

end

#
# Verify integral:
# ∫ y² erf( (a-y) / (sqrt(2) ℓ) ) dy
#

function verify_sub_3a(; a = a, c = c, d = d, ℓ = ℓ)

    f(y) = y^2 * erf( (a-y) / (sqrt(2) * ℓ) )

    # own numerical verification on 1d regular grid
    dx = 1e-6
    X  = c:dx:d
    num = 0.0
    for x in X
        num += f(x)*dx
    end


    # numerical verification via Cuba.jl libary
    h = Cuba.vegas((x, out) -> out[1] = (d-c) * f(c + (d-c)*x[1]), 1, 1, maxevals=10_000_000)

    display(h)


    # return exact and numerical results
    sub_3a(; a = a, c = c, d = d, ℓ = ℓ),
    h.integral[1],
    num

end




#
# Evaluate integral:
# ∫ y exp( -(a-y)² / (2 ℓ²) ) dy
#

function sub_3b(; a = a, c = c, d = d, ℓ = ℓ)

    # Wolram code: Integrate[y/E^((a - y)^2/(2 ℓ^2)), {y, c, d}]
    #
    # solution:
    # integral_c^d y e^(-(a - y)^2/(2 ℓ^2)) dy = 1/2 ℓ (2 ℓ (e^(-(a - c)^2/(2 ℓ^2)) - e^(-(a - d)^2/(2 ℓ^2))) + sqrt(2 π) a erf((a - c)/(sqrt(2) ℓ)) - sqrt(2 π) a erf((a - d)/(sqrt(2) ℓ)))

    1/2 * ℓ * (2 * ℓ * (exp(-(a - c)^2/(2*ℓ^2)) - exp(-(a - d)^2/(2*ℓ^2))) +
            sqrt(2*π) * a * erf((a - c)/(sqrt(2)*ℓ)) - sqrt(2*π) * a * erf((a - d)/(sqrt(2)*ℓ)))


end

#
# Verify integral:
# ∫ y exp( -(a-y)² / (2 ℓ²) ) dy
#

function verify_sub_3b(; a = a, c = c, d = d, ℓ = ℓ)

    f(y) = y * exp( -(a-y)^2 / (2 * ℓ^2) )

    # own numerical verification on 1d regular grid
    dx = 1e-6
    X  = c:dx:d
    num = 0.0
    for x in X
        num += f(x)*dx
    end


    # numerical verification via Cuba.jl libary
    h = Cuba.vegas((x, out) -> out[1] = (d-c) * f(c + (d-c)*x[1]), 1, 1, maxevals=10_000_000)

    display(h)


    # return exact and numerical results
    sub_3b(; a = a, c = c, d = d, ℓ = ℓ),
    h.integral[1],
    num

end




#
# Evaluate integral:
# ∫∫ y x exp(- (x-y)² / (2ℓ²) ) dx dy, x∈[a, b], y∈[c, d]
#
function integralfunction_3(; lowerx = lowerx, upperx = upperx, lowery = lowery, uppery = uppery, ℓ = ℓ)

    @assert(lowerx < upperx)
    @assert(lowery < uppery)

    # easier names
    a = lowerx
    b = upperx
    c = lowery
    d = uppery


    # Intermediate step: ∫ y x exp(- (x-y)² / (2ℓ²) ) dx dy, x∈[a, b]
    # Wolfram code: Integrate[(x y)/E^((x - y)^2/(2 ℓ^2)), {x, a, b}]
    # solution:
    # integral_a^b y x e^(-(x - y)^2/(2 ℓ^2)) dx = 1/2 ℓ y (2 ℓ (e^(-(a - y)^2/(2 ℓ^2)) - e^(-(b - y)^2/(2 ℓ^2))) - sqrt(2 π) y erf((a - y)/(sqrt(2) ℓ)) + sqrt(2 π) y erf((b - y)/(sqrt(2) ℓ)))
    #
    # We need to integrate the above solution over y y∈[c, d]

    1/2 * ℓ * (2 * ℓ * (sub_3b(; a = a, c = c, d = d, ℓ = ℓ) - sub_3b(; a = b, c = c, d = d, ℓ = ℓ)) -
                    sqrt(2*π) * sub_3a(; a = a, c = c, d = d, ℓ = ℓ) +
                    sqrt(2*π) * sub_3a(; a = b, c = c, d = d, ℓ = ℓ))

end


function verify_3(; lowerx = lowerx, upperx = upperx, lowery = lowery, uppery = uppery, ℓ = ℓ)

    # integrand
    f(x,y) =  y * x * exp(- (x-y)^2 / (2ℓ^2) )

    # numerical verification on regular grid
    Δx = 1e-3
    X  = lowerx:Δx:upperx
    Y  = lowery:Δx:uppery

    num = 0.0
    for i in 1:length(X)
        for j in 1:length(Y)
            @inbounds num += f(X[i], Y[j])
        end
    end
    num *= Δx^2


    # numerical verification via Cuba.jl libary
    h = divonne((x, out) -> out[1] = (uppery-lowery)*(upperx-lowerx) * f(lowerx + (upperx-lowerx)*x[1], lowery + (uppery-lowery)*x[2]), 2, 1)

    display(h)


    # return exact and numerical results to verify
    integralfunction_3(lowerx = lowerx, upperx = upperx, lowery = lowery, uppery = uppery, ℓ = ℓ), num, h.integral[1]

end
