using DataStructures
using Interpolations
using LinearAlgebra
using SparseArrays
using MKLSparse
using Dates
using Random
using MKL
using Printf
using NearestNeighbors

# KPM shift and scale
goldA=0.5
goldB=0.6

function gold_HV(Elist, Edict)

    # 18 positions of the nearest (and next nearest) neighbors in cartesian coordinates?
    nnlist=[[1,1,0];;[-1,1,0];;[1,-1,0];;[-1,-1,0];;
            [1,0,1];;[-1,0,1];;[1,0,-1];;[-1,0,-1];;
            [0,1,1];;[0,-1,1];;[0,1,-1];;[0,-1,-1];;
            [2,0,0];;[-2,0,0];;[0,2,0];;[0,-2,0];;[0,0,2];;[0,0,-2]]

    # KPM shift and scale
    A=goldA
    B=goldB




        
    # H=getHamiltonian(Elist,Edict,nnlist./2,nnlist,onsite,getVAA,A=A,B=B)
    function getHamiltonian(Elist::Deque{Vector{Int64}}, Edict::Dict{Vector{Int64},Int64},
                            nnlist::Array{Float64,2},    nnlist_for_table::Array{Int64,2},
                            getHopping::Function;      
                            Norbitals::Int64=9,          A::Float64=0.0,      B::Float64=1.0)

        # Elist  - list of all the positions of atoms inside the nanoparticle [(x1,y1,z1), (x2,y2,z2), ...]
        # Edict  - dictionary with all the sites inside the nanoparticle (x,y,z): index
        # nnlist - list of vectors to nearest all nearest neighbors (half of list below)
        # nnlist_for_table - list of vectors to nearest all nearest neighbors 
        # A - KPM shift
        # B - KPM scale

        NNA=length(Elist) # number of atoms

        # Local energies at each orbital
        Es=0.94261/2
        Ep=1.36110/2
        Ed1=0.37285/2
        Ed2=0.36265/2
        onsite=[Es,Ed1,Ed1,Ed1,Ed2,Ed2,Ep,Ep,Ep]
        onsite=(onsite.-A)./B # turning eV into KPM units
        
        function addV!(iidx,jidx,value,V,iptr,jptr)
            # Add the values of the V 9x9 matrix into a sparse matrix format
            # iidx (jidx) is the set of (lines,columns) of the sparse matrix
            # iptr, jptr index the atoms
            local i=0
            local j=0
            local vv=0.0
            local a,b=size(V)
            for i=1:a
                for j=1:b
                    vv=V[i,j]
                    if abs(vv)>1e-10
                        push!(iidx,iptr+i)
                        push!(jidx,jptr+j)
                        push!(value,vv)
                    end 
                end 
            end 
        end 


        # Build a neighbor table for every particle
        # [[n1, n2, ..., n9],     <- indices of neighbors of atom 1
        # [[n1, -1, ..., n9],     <- indices of neighbors of atom 2. -1 means it was not found
        # ...
        global table=Array{Int64,2}(undef,NNA,18)
        for Ri in Elist
            Idx=Edict[Ri]
            for j =1:18
                Rij=nnlist_for_table[:,j]
                Rj=Ri+Rij 
                Idj=get(Edict,Rj,-1)
                table[Idx,j]=Idj 
            end
        end 

        # build the sparse Hamiltonian matrix from the neighbor table and the hoppings
        global iidx=Array{Int64,1}([])
        global jidx=Array{Int64,1}([])
        global value=Array{ComplexF64,1}([])
        for Atomi=1:NNA
            iptr=(Atomi-1)*Norbitals
            for k=1:Norbitals
                push!(iidx,k+iptr)
                push!(jidx,k+iptr)
                push!(value,onsite[k])
            end
            for j=1:18
                Atomj=table[Atomi,j]
                if Atomj>=0
                    #Neighbouring atom is atomA
                    jptr=(Atomj-1)*Norbitals
                    Rij=nnlist[:,j]
                    V=getHopping(Rij)./B
                    addV!(iidx,jidx,value,V,iptr,jptr)
                end
            end
        end
        H=sparse(iidx,jidx,value)
        return H
    end 

                











    function getVAA(R)
        V=Array{Float64,2}(undef,9,9)
        function getHoppingParameters(d)
            if d<0.99999999
                ss_sig=-0.07962/2
                pp_sig=0.17119/2
                pp_pi=-0.00540/2
                dd_sig=-0.05216/2
                dd_pi=0.02878/2
                dd_del=-0.00533/2
                sp_sig=0.11332/2
                sd_sig=-0.04885/2
                pd_sig=-0.06563/2
                pd_pi=0.02124/2
            else
                ss_sig=-0.00105/2
                pp_sig=0.04282/2
                pp_pi=-0.00044/2
                dd_sig=-0.00385/2
                dd_pi=0.00212/2
                dd_del=-0.00026/2
                sp_sig=0.01048/2
                sd_sig=-0.00837/2
                pd_sig=-0.00738/2
                pd_pi=0.00351/2
            end 
            return sp_sig,ss_sig,pp_sig,pp_pi,sd_sig,pd_sig,pd_pi,dd_sig,dd_pi,dd_del
        end
        #compute these elements listed in SK table
        l,m,n=R./norm(R)
        sp_sig,ss_sig,pp_sig,pp_pi,sd_sig,pd_sig,pd_pi,dd_sig,dd_pi,dd_del=getHoppingParameters(norm(R))
        #Same orbital interaction
        V[1,1]=ss_sig
        V[2,2]=( 3*l^2*m^2*dd_sig + ( l^2 + m^2 - 4*l^2*m^2)*dd_pi +(n^2+l^2*m^2)*dd_del)
        V[3,3]=( 3*l^2*n^2*dd_sig + ( l^2 + n^2 - 4*l^2*n^2)*dd_pi +(m^2+l^2*n^2)*dd_del)
        V[4,4]=( 3*m^2*n^2*dd_sig + ( m^2 + n^2 - 4*m^2*n^2)*dd_pi +(l^2+m^2*n^2)*dd_del)
        V[5,5]=(3/4*(l^2-m^2)^2*dd_sig+(l^2+m^2-(l^2-m^2)^2)*dd_pi+(n^2+(l^2-m^2)^2/4)*dd_del)
        V[6,6]=((n^2-.5*(l^2+m^2))^2*dd_sig+3*n^2*(l^2+m^2)*dd_pi +3/4*(l^2+m^2)^2*dd_del)
        V[7,7]=(l^2*pp_sig+(1-l^2)*pp_pi)
        V[8,8]=(m^2*pp_sig+(1-m^2)*pp_pi)
        V[9,9]=(n^2*pp_sig+(1-n^2)*pp_pi)
        #p-p interaction
        V[7,8]=(l*m*pp_sig-l*m*pp_pi)
        V[7,9]=(l*n*pp_sig-l*n*pp_pi)
        V[8,9]=(m*n*pp_sig-m*n*pp_pi)
        #dd interaction
        V[2,3]=(3*l^2*m*n*dd_sig+m*n*(1-4*l^2)*dd_pi+m*n*(l^2-1)*dd_del) #xy-xz
        V[2,4]=(3*l*m^2*n*dd_sig+l*n*(1-4*m^2)*dd_pi+l*n*(m^2-1)*dd_del) #xy-yz
        V[2,5]=(1.5*l*m*(l^2-m^2)*dd_sig+2*l*m*(m^2-l^2)*dd_pi+.5*l*m*(l^2-m^2)*dd_del) #xy - x2-y2
        V[2,6]=(3^.5*l*m*(n^2-.5*(l^2+m^2))*dd_sig-2*3^.5*l*m*n^2*dd_pi+.5*3^.5*l*m*(1+n^2)*dd_del) #xy - 3z^2-r^2
        V[3,4]=(3*n^2*m*l*dd_sig+m*l*(1-4*n^2)*dd_pi+l*m*(n^2-1)*dd_del) #xz-yz
        V[3,5]=(1.5*n*l*(l^2-m^2)*dd_sig+n*l*(1-2*(l^2-m^2))*dd_pi-n*l*(1-.5*(l^2-m^2))*dd_del) #xz ->x^2-y^2
        V[3,6]=(3^.5*l*n*(n^2-.5*(l^2+m^2))*dd_sig+3^.5*l*n*(l^2+m^2-n^2)*dd_pi-.5*3^.5*l*n*(l^2+m^2)*dd_del) #xz ->z^2
        V[4,5]=(1.5*m*n*(l^2-m^2)*dd_sig-m*n*(1+2*(l^2-m^2))*dd_pi+m*n*(1+(l^2-m^2)/2)*dd_del) #yz ->x^2-y^2
        V[4,6]=(3^.5*m*n*(n^2-.5*(l^2+m^2))*dd_sig+3^.5*m*n*(l^2+m^2-n^2)*dd_pi-.5*3^.5*m*n*(l^2+m^2)*dd_del)#yz->z^2
        V[5,6]=(.5*3^.5*(l^2-m^2)*(n^2-.5*(l^2+m^2))*dd_sig+3^.5*n^2*(m^2-l^2)*dd_pi+3^.5*(1+n^2)*(l^2-m^2)/4*dd_del)
        #s-d interaction
        V[1,2]=3^.5 * l*m*sd_sig
        V[1,3]=3^.5*l*n*sd_sig
        V[1,4]=3^.5*n*m*sd_sig
        V[1,5]=3^.5/2*(l^2-m^2)*sd_sig
        V[1,6]=(n^2-.5*(l^2+m^2))*sd_sig
        #sp interaction
        V[1,7]=(l*sp_sig)
        V[1,8]=(m*sp_sig)
        V[1,9]=(n*sp_sig)
        #pd interaction
        V[7,2]=(3^.5*l^2*m*pd_sig+m*(1-2*l^2)*pd_pi)#x->xy 
        V[8,2]=(3^0.5*m^2*l*pd_sig+l*(1-2*m^2)*pd_pi)#y->xy
        V[9,2]=(3^0.5*l*m*n*pd_sig-2*l*m*n*pd_pi)#z->xy
        V[7,3]=(3^.5*l^2*n*pd_sig+n*(1-2*l^2)*pd_pi) #x->xz
        V[8,3]=(3^0.5*l*m*n*pd_sig-2*l*m*n*pd_pi) #xz->y
        V[9,3]=(3^0.5*n^2*l*pd_sig+l*(1-2*n^2)*pd_pi) #xz->z
        V[7,4]=(3^0.5*l*m*n*pd_sig-2*l*m*n*pd_pi) #yz->x
        V[8,4]=(3^0.5*m^2*n*pd_sig+n*(1-2*m^2)*pd_pi) #yz->y
        V[9,4]=(3^0.5*n^2*m*pd_sig+m*(1-2*n^2)*pd_pi) #yz->z
        V[7,5]=(3^0.5/2*l*(l^2-m^2)*pd_sig+l*(1-l^2+m^2)*pd_pi) #x^2-y^2->x
        V[8,5]=(3^0.5/2*m*(l^2-m^2)*pd_sig-m*(1+l^2-m^2)*pd_pi) #x^2-y^2->y
        V[9,5]=(3^0.5/2*n*(l^2-m^2)*pd_sig-n*(l^2-m^2)*pd_pi) #x^2-y^2->z
        V[7,6]=(l*(n^2-(l^2+m^2)/2)*pd_sig-3^0.5*l*n^2*pd_pi) #3z^2-r^2->x
        V[8,6]=(m*(n^2-(l^2+m^2)/2)*pd_sig-3^0.5*m*n^2*pd_pi) #3z^2-r^2->y
        V[9,6]=(n*(n^2-(l^2+m^2)/2)*pd_sig+3^0.5*n*(l^2+m^2)*pd_pi)#3z^2-r^2->z
        #For those hopping not listed in SK table
        #Compute the opposite hopping and use Hermiticity of Hamiltonian to 
        #derive these elements. 
        #Now flip the direction:
        l,m,n=R./(-1*norm(R))
        #dd interaction
        V[3,2]=(3*l^2*m*n*dd_sig+m*n*(1-4*l^2)*dd_pi+m*n*(l^2-1)*dd_del) #xy-xz
        V[4,2]=(3*l*m^2*n*dd_sig+l*n*(1-4*m^2)*dd_pi+l*n*(m^2-1)*dd_del) #xy-yz
        V[5,2]=(1.5*l*m*(l^2-m^2)*dd_sig+2*l*m*(m^2-l^2)*dd_pi+.5*l*m*(l^2-m^2)*dd_del) #xy - x2-y2
        V[6,2]=(3^.5*l*m*(n^2-.5*(l^2+m^2))*dd_sig-2*3^.5*l*m*n^2*dd_pi+.5*3^.5*l*m*(1+n^2)*dd_del) #xy - 3z^2-r^2
        V[4,3]=(3*n^2*m*l*dd_sig+m*l*(1-4*n^2)*dd_pi+l*m*(n^2-1)*dd_del) #xz-yz
        V[5,3]=(1.5*n*l*(l^2-m^2)*dd_sig+n*l*(1-2*(l^2-m^2))*dd_pi-n*l*(1-.5*(l^2-m^2))*dd_del) #xz ->x^2-y^2
        V[6,3]=(3^.5*l*n*(n^2-.5*(l^2+m^2))*dd_sig+3^.5*l*n*(l^2+m^2-n^2)*dd_pi-.5*3^.5*l*n*(l^2+m^2)*dd_del) #xz ->z^2
        V[5,4]=(1.5*m*n*(l^2-m^2)*dd_sig-m*n*(1+2*(l^2-m^2))*dd_pi+m*n*(1+(l^2-m^2)/2)*dd_del) #yz ->x^2-y^2
        V[6,4]=(3^.5*m*n*(n^2-.5*(l^2+m^2))*dd_sig+3^.5*m*n*(l^2+m^2-n^2)*dd_pi-.5*3^.5*m*n*(l^2+m^2)*dd_del)#yz->z^2
        V[6,5]=(.5*3^.5*(l^2-m^2)*(n^2-.5*(l^2+m^2))*dd_sig+3^.5*n^2*(m^2-l^2)*dd_pi+3^.5*(1+n^2)*(l^2-m^2)/4*dd_del)
        #p-p interaction
        V[8,7]=(l*m*pp_sig-l*m*pp_pi)
        V[9,7]=(l*n*pp_sig-l*n*pp_pi)
        V[9,8]=(m*n*pp_sig-m*n*pp_pi)
        #ds interaction
        V[2,1]=3^.5 * l*m*sd_sig
        V[3,1]=3^.5*l*n*sd_sig
        V[4,1]=3^.5*n*m*sd_sig  
        V[5,1]=3^.5/2*(l^2-m^2)*sd_sig
        V[6,1]=(n^2-.5*(l^2+m^2))*sd_sig
        #ps interaction
        V[7,1]=(l*sp_sig)
        V[8,1]=(m*sp_sig)
        V[9,1]=(n*sp_sig)
        #dp interaction
        V[2,7]=(3^.5*l^2*m*pd_sig+m*(1-2*l^2)*pd_pi)#x->xy 
        V[2,8]=(3^0.5*m^2*l*pd_sig+l*(1-2*m^2)*pd_pi)#y->xy
        V[2,9]=(3^0.5*l*m*n*pd_sig-2*l*m*n*pd_pi)#z->xy
        V[3,7]=(3^.5*l^2*n*pd_sig+n*(1-2*l^2)*pd_pi) #x->xz
        V[3,8]=(3^0.5*l*m*n*pd_sig-2*l*m*n*pd_pi) #xz->y
        V[3,9]=(3^0.5*n^2*l*pd_sig+l*(1-2*n^2)*pd_pi) #xz->z
        V[4,7]=(3^0.5*l*m*n*pd_sig-2*l*m*n*pd_pi) #yz->x
        V[4,8]=(3^0.5*m^2*n*pd_sig+n*(1-2*m^2)*pd_pi) #yz->y
        V[4,9]=(3^0.5*n^2*m*pd_sig+m*(1-2*n^2)*pd_pi) #yz->z
        V[5,7]=(3^0.5/2*l*(l^2-m^2)*pd_sig+l*(1-l^2+m^2)*pd_pi) #x^2-y^2->x
        V[5,8]=(3^0.5/2*m*(l^2-m^2)*pd_sig-m*(1+l^2-m^2)*pd_pi) #x^2-y^2->y
        V[5,9]=(3^0.5/2*n*(l^2-m^2)*pd_sig-n*(l^2-m^2)*pd_pi) #x^2-y^2->z
        V[6,7]=(l*(n^2-(l^2+m^2)/2)*pd_sig-3^0.5*l*n^2*pd_pi) #3z^2-r^2->x
        V[6,8]=(m*(n^2-(l^2+m^2)/2)*pd_sig-3^0.5*m*n^2*pd_pi) #3z^2-r^2->y
        V[6,9]=(n*(n^2-(l^2+m^2)/2)*pd_sig+3^0.5*n*(l^2+m^2)*pd_pi)#3z^2-r^2->z
        return V
    end 


    H = getHamiltonian(Elist, Edict, nnlist./2, nnlist, getVAA, A=A, B=B)
    return H
end
