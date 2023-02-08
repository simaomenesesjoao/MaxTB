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

# Important constants
HaeV = 27.211386245988 # [eV]  Hartree in eV

# include("sumcheb.jl")

function fermi_dirac(x,beta,mu)
    # Fermi dirac distribution
    arg = beta*(x-mu)
    return 1.0/(1+exp(arg))
end




function read_positions(filename)
    # Read the atomic positions 
    
    fr = open(filename,"r")

    # Ignore the first 8 lines (header) to output file
    for i=1:9
        line = readline(fr)
    end

    # Iterate over the positions and write to file
    flines=readlines(fr)
    for i in flines
        j = split(i)
        x = parse(Float64,j[1])
        y = parse(Float64,j[2])
        z = parse(Float64,j[3])
    end
    close(fr)

    return R
end


function doscompute_mu!(H,N,Nk)
    # Calculate the Chebyshev moments for the DoS
    a,b = size(H) # a=size of Hilbert space
    list = zeros(ComplexF64,N)
    ring = zeros(ComplexF64,a,2)

    # Run over Nk random vectors
    for k=1:NK

        # Define the random vector
        r = rand(Float64,a).-(0.5)
        r = r./(norm(r)/a^0.5)
        conjr = copy(r) 

        # First two Chebyshev iterations
        ring[1,:] = r
        ring[2,:] = H*r

        # Project back into the original vector
        list[1] = 1
        list[2] = vecdot(conjr, ring[2,:])

        # Remaining Chebyshev iterations
        for n=3:N
            i = n%2 + 1
            j = (n+1)%2 + 1
            ring[i,:] = H*ring[j,:].*2.0 - ring[i,:]

            # Project
            list[n] = vecdot(conjr, ring[i,:])
        end
    end
    return list
end




function compute_mumn!(H,Phi,NNL,NNR,NTx,NTy,Nk)
    # Calculate the Chebyshev moment matrix for the hot carrier generation rate
    # NNL*NTx is the number of Cheb moments

    # Compute a matrix of Chebyshev vectors (starting from n=0)
    # each row n is a vector of the size of the Hilbert space: |n> = Tn(H) |r>
    function init_vec_list!(list,r)
        local a,b = size(list)
        list[1,:] = r 
        list[2,:] = H*r 
        for i=3:a 
            list[i,:] = (H*list[i-1,:]).*2.0-list[i-2,:]
        end 
        return nothing 
    end 

    # Compute a matrix of Chebyshev vectors using the last two 
    # vectors as starting points
    function update_vector_list!(list)
        local a,b = size(list)
        list[1,:] = (H*list[a,:]).*2.0-list[a-1,:]
        list[2,:] = (H*list[1,:]).*2.0-list[a,:]
        for i=3:a 
            list[i,:] = (H*list[i-1,:]).*2.0-list[i-2,:]
        end 
        return nothing 
    end 

    function process_vector_list!(left_list,right_list,phi)
        local a,b = size(right_list) # a=number of polys, b=size of Hilbert space
        local aux_matrix = Array{ComplexF64,2}(undef,b,a)

        # Act with the potential on all the vectors ϕ|n>
        for i=1:a 
            aux_matrix[:,i] = phi*right_list[i,:]
        end 

        # Project on to the left vectors: <m|ϕ|n>
        # Yields a matrix of dimensions a x a
        local mumn = left_list*aux_matrix
        return mumn
    end 

    PhiT = adjoint(Phi) # potential is not hermitian
    a, b = size(H)      # a=size of Hilbert space

    # Lists with Chebyshev moments. 
    # NNL (NNR) is the number of polys per block for the left (right) list
    left_list  = zeros(ComplexF64,NNL,a)
    right_list = zeros(ComplexF64,NNR,a)

    global mumn = zeros(ComplexF64,NNL*NTx,NNR*NTy)

    # Do the calculation for Nk random vectors
    for k=1:Nk
        mumn_k=Array{ComplexF64,2}(undef,NNL*NTx,NNR*NTy)

        # define the random vector.
        # I was initially planning to do complex random vectors, 
        # now I decided to use real random vector, because mumn should be real
        r=rand(Float64,a).-(0.5)
        r=r./(norm(r)/a^0.5)
        conjr=copy(r) 

        # Product of electric potential by the random vector
        Phir=Phi*r

        for i=0:NTx-1
            if i==0
                init_vec_list!(left_list,conjr)
            else
                update_vector_list!(left_list)
            end
        
            for j=0:NTy-1
                if j==0
                    init_vec_list!(right_list,Phir)
                else
                    update_vector_list!(right_list)
                end
                mumn_aux=process_vector_list!(left_list,right_list,PhiT)
                mumn_k[i*NNL+1:(i+1)*NNL,j*NNR+1:(j+1)*NNR]=mumn_aux
            end 
        end

        global mumn = mumn+mumn_k # add to the average
        
    end 
    global mumn=mumn./Nk
    left_list=0
    right_list=0
        
    mumn=real(mumn)
    return mumn
end 


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

