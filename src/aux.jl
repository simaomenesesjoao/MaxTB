

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


function compute_eh_list(mumn,Nx,Ny,hw,A,B,Fermi,beta,kappa,gph)
    # Resum Chebyshev matrix

    function meshgrid(x,y)
        local Nx=length(x)
        local Ny=length(y)
        xx=zeros(Float64,Nx,Ny)
        yy=zeros(Float64,Nx,Ny)
        for i=1:Nx
            xx[i,:].=x[i]
        end
        for j=1:Ny
            yy[:,j].=y[j]
        end 
        return xx,yy
    end     


    xlist=LinRange(-0.99,0.99,Nx)
    ylist=LinRange(-0.99,0.99,Ny)
    NTx,NTy=size(mumn)
    NTxlist=0:1:NTx-1|>collect 
    NTylist=0:1:NTy-1|>collect 

    # going from KPM to eV
    Ei=(xlist.*B).+A 
    Ej=(ylist.*B).+A

    function jackson(n,N)
        a = pi/(N+1)
        b = (N+1-n)*cos(a*n)
        c = sin(a*n)/tan(a)
        return (b+c)/(N+1)
    end

    # apply Jackson kernel to mu matrix (include the 1/2 factor in here)
    Jx = zeros(NTx,1)
    Jy = zeros(NTy,1)
    for i in 1:NTx
        Jx[i] = jackson(i-1, NTx)
    end
    for i in 1:NTy
        Jy[i] = jackson(i-1, NTy)
    end
    Jx[1]=Jx[1]*0.5
    Jy[1]=Jy[1]*0.5
    Jxx,Jyy=meshgrid(Jx,Jy)
    mutmn=mumn.*Jxx.*Jyy 

    # Get the E,E' matrix by summing the Chebyshev series
    NTyy,yy=meshgrid(NTylist,ylist)
    Dny=cos.(NTyy.*acos.(yy))./sqrt.(1.0.-yy.^2)
    xx,NTxx=meshgrid(xlist,NTxlist)
    Dxm=cos.(NTxx.*acos.(xx))./sqrt.(1.0.-xx.^2)

    phi=Dxm*mutmn*Dny 
    factor=4/B^2/pi^2
    phi=phi.*factor


    #writedlm("mutmn.dat", mutmn) 
    # writedlm("phi_cheb.dat", phi) # phi is the E,E' matrix
    # writedlm("energies_jl.dat", Ei)
    
    # Final step of the calculation: final integration of E,E' matrix
    # with a Gaussian broadened Dirac delta

    # The width is just a constant
    # kappa = 0.0
    # gph = 0.06/HaeV
    function getgi(Ei,gph,kappa,Fermi)
        ans=(kappa.*(Ei.-Fermi).^2).+gph
        return ans 
    end 
    gi = getgi(Ei,gph,kappa,Fermi)
    gj = getgi(Ej,gph,kappa,Fermi)

    # Integrating step: kernel
    factor = zeros(Float64,Nx,Ny)
    hole_factor = zeros(Float64,Nx,Ny)
    for i=1:Nx 
        fi = 1.0/(1.0 + exp(beta*(Ei[i]-Fermi)))
        for j=1:Ny 
            fj = 1.0/(1.0 + exp(beta*(Ej[j]-Fermi)))

            gij  = gi[i] + gj[j] 
            amp  = 1.0/sqrt(2.0*pi*gij^2) # amplitude of the gaussian
            d_ij = amp*exp(-(hw+(Ej[j]-Ei[i] ) )^2/2.0/gij^2)
            d_ji = amp*exp(-(hw-(Ej[j]-Ei[i] ) )^2/2.0/gij^2)
            
            factor[i,j]      = d_ji*fi*(1.0-fj)
            hole_factor[i,j] = d_ij*fj*(1.0-fi)
        end 
    end 

    
    #writedlm("mask.dat", factor)
    
    Nelist=zeros(Ny)
    Nhlist=zeros(Ny)
    
    # Integration step
    dx=Ei[2]-Ei[1]
    for i=1:Ny 
        for j=1:Nx
            Nelist[i]=Nelist[i]+phi[j,i]*factor[j,i]*dx
            Nhlist[i]=Nhlist[i]+phi[j,i]*hole_factor[j,i]*dx
        end 
    end 
    Nelist=Nelist.*(2*pi)
    Nhlist=Nhlist.*(2*pi)
    
    return Ej,Nelist,Nhlist
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

