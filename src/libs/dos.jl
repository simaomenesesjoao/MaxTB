include("aux.jl")

function get_dos(mu, perc, minE_Ha, maxE_Ha, NE, A, B)
    # mu is the array of Chebyshev moments
    # perc is the percentage of moments to keep from the mu array
    # minE_Ha and maxE_Ha (Hartree) are the minimum and maximum values of energy, NE is the size
    # A is the shift and B is the scaling factor for KPM

    M = length(mu)
    Meff = Int(trunc(M*perc/100))

    # Jackson kernel
    jack = jackson_coefs(Meff)

    # Truncate the Chebyshev moment matrix and apply kernel
    mu_ef = mu[1:Meff].*jack

    # Make sure the energy limits are correct
    tol = 0.01
    if (minE_Ha - A)/B <= -1 
        println("Lower limit is too low, updating")
        minE_Ha = -B + A + tol
    end
    if (maxE_Ha - A)/B >= 1
        println("Higher limit is too high, updating")
        minE_Ha = -B + A + tol
        maxE_Ha = B + A - tol
    end

    # Create the energy array
    xlist = LinRange(minE_Ha, maxE_Ha, NE)
    # println(xlist[1], " ", xlist[end])
    xlist = (collect(xlist) .- A)./B
    # xlist = LinRange(-0.99, 0.99, NE)
    # println(xlist[1], " ",xlist[end])

    # Resum the Chebyshev series decomposition
    T = zeros(ComplexF64, 2, NE)
    T[1,:] = xlist.*0.0 .+ 1.0
    T[2,:] = xlist.*1.0

    dos  = zeros(ComplexF64, NE)
    dos += mu_ef[1].*T[1,:]
    dos += mu_ef[2].*T[2,:]

    for n=3:Meff
        i = n%2 + 1
        j = (n+1)%2 + 1
        T[j,:] = xlist.*T[i,:].*2.0 - T[j,:]

        dos += T[j,:].*mu_ef[n]
    end
    
    dos = dos./sqrt.(1 .- xlist.*xlist)./pi.*2

    data = zeros(Float64, NE, 2)
    data[:,1] = (xlist.*B .+ A)*HaeV
    data[:,2] = real.(dos)./B./HaeV

    # println(data[1,1], " " , data[end,1])
    return data
                                
end
