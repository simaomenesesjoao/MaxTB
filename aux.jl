

function print_positions(R)
    a,b = size(R)
    f=open("atomic_positions.txt","w")
    for i=1:a
        x,y,z = R[i,:]
        write(f,@sprintf("%f %f %f\n",x,y,z))
    end
    close(f)
    
end



function compute_mumn!(H,Phi,NNL,NNR,NTx,NTy,Nk)
    function init_vec_list!(list,r)
        local a,b=size(list)
        list[1,:]=r 
        list[2,:]=H*r 
        for i=3:a 
            list[i,:]=(H*list[i-1,:]).*2.0-list[i-2,:]
        end 
        return nothing 
    end 

    function update_vector_list!(list)
        local a,b=size(list)
        list[1,:]=(H*list[a,:]).*2.0-list[a-1,:]
        list[2,:]=(H*list[1,:]).*2.0-list[a,:]
        for i=3:a 
            list[i,:]=(H*list[i-1,:]).*2.0-list[i-2,:]
        end 
        return nothing 
    end 

    function process_vector_list!(left_list,right_list,phi)
        local a,b=size(right_list)
        local aux_matrix=Array{ComplexF64,2}(undef,b,a)
        for i=1:a 
            aux_matrix[:,i]=phi*right_list[i,:]
            #project onto the adsorbed states
            #aux_matrix[1:NNA*9,i].=0.0
        end 
        local mumn=left_list*aux_matrix
        return mumn
    end 

    PhiT=adjoint(Phi)
    a,b=size(H)
    left_list=zeros(ComplexF64,NNL,a)
    right_list=zeros(ComplexF64,NNR,a)

    global mumn=zeros(ComplexF64,NNL*NTx,NNR*NTy)

    # Do the calculation for Nk random vectors
    for k=1:Nk
        mumn_k=Array{ComplexF64,2}(undef,NNL*NTx,NNR*NTy)
        #Now set the random seed
        
        r=rand(Float64,a).-(0.5)
        r=r./(norm(r)/a^0.5)
        conjr=copy(r)#I was initially planning to do complex random vectors, 
        #now I decided to use real random vector, because mumn should be real

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
        global mumn=mumn+mumn_k
        
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

    # apply Jackson kernel to mu matrix (include the 1/2 factor in here)
    Jx=1.0/NTx.*((NTx.-NTxlist).*cos.(pi*NTxlist./NTx)+sin.(pi.*NTxlist./NTx).*cos(pi/NTx)./sin(pi/NTx))
    Jy=1.0/NTy.*((NTy.-NTylist).*cos.(pi*NTylist./NTy)+sin.(pi.*NTylist./NTy).*cos(pi/NTy)./sin(pi/NTy))
    Jx[1]=Jx[1]*0.5
    Jy[1]=Jy[1]*0.5
    Jxx,Jyy=meshgrid(Jx,Jy)
    mutmn=mumn.*Jxx.*Jyy 

    # Get the E,E' matrix
    NTyy,yy=meshgrid(NTylist,ylist)
    Dny=cos.(NTyy.*acos.(yy))./sqrt.(1.0.-yy.^2)
    xx,NTxx=meshgrid(xlist,NTxlist)
    Dxm=cos.(NTxx.*acos.(xx))./sqrt.(1.0.-xx.^2)
    phi=Dxm*mutmn*Dny 
    factor=pi^2/B^2/4
    phi=phi.*factor


    writedlm("phi_cheb.dat", phi)
    writedlm("energies_jl.dat", Ei)
    
    # Final step of the calculation: final integration of E,E' matrix
    # with a Gaussian broadened Dirac delta

    # The width is just a constant
    # kappa = 0.0
    # gph = 0.06/HaeV
    function getgi(Ei,gph,kappa,Fermi)
        ans=(kappa.*(Ei.-Fermi).^2).+gph
        return ans 
    end 
    gi=getgi(Ei,gph,kappa,Fermi)
    gj=getgi(Ej,gph,kappa,Fermi)

    # Integrating step
    factor = zeros(Float64,Nx,Ny)
    hole_factor = zeros(Float64,Nx,Ny)
    for i=1:Nx 
        fi = 1.0/(1.0 + exp(beta*(Ei[i]-Fermi)))
        for j=1:Ny 
            fj = 1.0/(1.0 + exp(beta*(Ej[j]-Fermi)))

            gij  = gi[i] + gj[j] 
            amp  = 2*pi*2.0/sqrt(2.0*pi*gij^2) # amplitude of the gaussian
            d_ij = amp*exp(-(hw+(Ej[j]-Ei[i] ) )^2/2.0/gij^2)
            d_ji = amp*exp(-(hw-(Ej[j]-Ei[i] ) )^2/2.0/gij^2)
            
            factor[i,j]      = d_ji*fi*(1.0-fj)
            hole_factor[i,j] = d_ij*fj*(1.0-fi)
        end 
    end 
    

    
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

    writedlm("M_jl.dat", Gamma)
    writedlm("eigen_energies_jl.dat", vals)

    

end

