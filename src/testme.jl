using Test

# Verify symmetric kernel calculation by comparing with numerical integration results
for i = 1:10
    @assert(test_GKGT())
end



# Verify integrals
tol = 1e-3

for i = 1:10

    lowery = rand(Uniform(-30.0, 30.0))
    uppery = lowery + rand(Uniform(-30.0, 30.0))

    lowerx = rand(Uniform(-30.0, 30.0))
    upperx = lowerx + rand(Uniform(-30.0, 30.0))

    ℓ = rand(Uniform(0.01, 6.0))

    a = lowerx
    b = upperx
    c = lowery
    d = uppery

    #####
    # 1 #
    #####

    local out = verify_1(lowery = -2.1, uppery = 0.11, lowerx = -2, upperx = 3.1, ℓ=1.6)
    display(out)
    @assert(abs(out[1]-out[3]) < tol)

    out = verify_sub_1(b=-1.1, c=-0.2, d=0.52, ℓ=0.4)
    display(out)
    @assert(abs(out[1]-out[2]) < tol)
    @assert(abs(out[1]-out[3]) < tol)


    #####
    # 2 #
    #####

    out = verify_2(lowery = -2.1, uppery = 0.11, lowerx = -2, upperx = 3.1, ℓ=1.6)
    display(out)
    @assert(abs(out[1]-out[3]) < tol)

    out = verify_sub_2a(a=-1.2,c=-0.25,d=0.9, ℓ=0.9)
    display(out)
    @assert(abs(out[1]-out[2]) < tol)

    out = verify_sub_2b(a=-1.2,c=-0.25,d=0.9, ℓ=0.9)
    display(out)
    @assert(abs(out[1]-out[2]) < tol)


    #####
    # 3 #
    #####

    out = verify_3(lowery = -2.1, uppery = 0.11, lowerx = -2, upperx = 3.1, ℓ=1.6)
    display(out)
    @assert(abs(out[1]-out[3]) < tol)

    out = verify_sub_3a(a=-1.2,c=-0.25,d=0.9, ℓ=0.9)
    display(out)
    @assert(abs(out[1]-out[2]) < tol)

    out = verify_sub_3b(a=-1.2,c=-0.25,d=0.9, ℓ=0.9)
    display(out)
    @assert(abs(out[1]-out[2]) < tol)

end

@info("if you see this message, then tests must be succesful")
