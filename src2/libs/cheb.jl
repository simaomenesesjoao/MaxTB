


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

