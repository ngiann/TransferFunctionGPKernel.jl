function testme(MAXTRIES=3; tol=0.1)

    # Verify symmetric kernel calculation by comparing with numerical integration results
    for i = 1:MAXTRIES
        @assert(test_GKGT_random(tol=tol))
    end



    # Verify integrals

    for i = 1:MAXTRIES

        lowery = rand() * 60 - 30
        uppery = lowery + rand() * 60

        lowerx = rand() * 60 - 30
        upperx = lowerx + rand() * 60

        ℓ = rand() * (5.0 - 0.01) + 0.01

        a = lowerx
        b = upperx
        c = lowery
        d = uppery

        #####
        # 1 #
        #####
        @printf("verify_1\n")
        local out = verify_1(lowery = lowery, uppery = uppery, lowerx = lowerx, upperx = upperx, ℓ=ℓ)
        display(out)
        @assert(abs(out[1]-out[3]) < tol)

        @printf("\t verify_sub_1\n\n")
        out = verify_sub_1(b=b, c=c, d=d, ℓ=ℓ)
        display(out)
        @assert(abs(out[1]-out[2]) < tol)
        @assert(abs(out[1]-out[3]) < tol)


        #####
        # 2 #
        #####
        @printf("verify_2\n")
        out = verify_2(lowery = lowery, uppery = uppery, lowerx = lowerx, upperx = upperx, ℓ=ℓ)
        display(out)
        @assert(abs(out[1]-out[3]) < tol)

        @printf("\t verify_sub_2a\n")
        out = verify_sub_2a(a=a,c=c,d=d, ℓ=ℓ)
        display(out)
        @assert(abs(out[1]-out[2]) < tol)

        @printf("\t verify_sub_2b\n\n")
        out = verify_sub_2b(a=a,c=c,d=d, ℓ=ℓ)
        display(out)
        @assert(abs(out[1]-out[2]) < tol)


        #####
        # 3 #
        #####
        @printf("verify_3\n")
        out = verify_3(lowery = lowery, uppery = uppery, lowerx = lowerx, upperx = upperx, ℓ=ℓ)
        display(out)
        @assert(abs(out[1]-out[3]) < tol)

        @printf("\t verify_sub_3a\n")
        out = verify_sub_3a(a=a,c=c,d=d, ℓ=ℓ)
        display(out)
        @assert(abs(out[1]-out[2]) < tol)

        @printf("\t verify_sub_3b\n\n")
        out = verify_sub_3b(a=a,c=c,d=d, ℓ=ℓ)
        display(out)
        @assert(abs(out[1]-out[2]) < tol)

    end

    @info("if you see this message, then tests must be succesful")

end
