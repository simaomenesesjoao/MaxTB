boltzmann = 8.617333e-5      # [eV/Kelvin] Boltzmann constant
HaeV = 27.211386245988       # [eV]        Hartree in eV


function meshgrid(x,y)
    local Nx = length(x)
    local Ny = length(y)
    xx = zeros(Float64, Nx, Ny)
    yy = zeros(Float64, Nx, Ny)
    for i=1:Nx
        xx[i,:].= x[i]
    end
    for j=1:Ny
        yy[:,j].= y[j]
    end 
    return xx, yy
end     

function jackson_kernel(n,N)
    # Jackson kernel
    a = pi/(N+1)
    b = (N+1-n)*cos(a*n)
    c = sin(a*n)/tan(a)
    return (b+c)/(N+1)
end

function jackson_coefs(N)
    # Build array of N Jackson kernel coefficients
    J = zeros(N,1)
    for i in 1:N
        J[i] = jackson_kernel(i-1, N)
    end
    J[1] = J[1]*0.5
    return J
end


function fermi_dirac(x,beta,mu)
    # Fermi dirac distribution
    arg = beta*(x-mu)
    return 1.0/(1+exp(arg))
end




# function read_positions(filename)
    # # Read the atomic positions 
    
    # fr = open(filename,"r")

    # # Ignore the first 8 lines (header) to output file
    # for i=1:9
        # line = readline(fr)
    # end

    # # Iterate over the positions and save to array
    # flines=readlines(fr)
    # for i in flines
        # j = split(i)
        # x = parse(Float64,j[1])
        # y = parse(Float64,j[2])
        # z = parse(Float64,j[3])
    # end
    # close(fr)

    # return R
# end


function exact(H, Phi, Fermi, beta)
    println("Exact diagonalization")
    # Exact diagonalization
    
    H_dens = Matrix(H)
    Phi_dens = Matrix(Phi)
    a,b = size(H_dens)

    vals, vecs = eigen(H_dens)
    # testing vecs
    # n = 2
    # no = norm(Hdens*vecs[:,n] - vals[n]*vecs[:,n])
    # println("norm", no)
    println(typeof(Phi_dens), typeof(vecs))
    phi_eigen = vecs'*Phi_dens*vecs
    Gamma = zeros(a,a)

    # function fermi(x)
        # return 1.0/(1.0 + exp(beta*(x-Fermi)))
    # end

    for i in 1:a
        Ei = vals[i]
        for j in 1:a
            Ej = vals[j]
            Gamma[i,j] = abs(phi_eigen[j,i])^2#*fermi(Ei)*(1-fermi(Ej))
        end
    end


    # open("matrix.txt", "w") do file
        # write(file, Gamma)
    # end

    writedlm("Phi_nm_jl.dat", Gamma)
    writedlm("eigen_energies_jl.dat", vals)

    

end

