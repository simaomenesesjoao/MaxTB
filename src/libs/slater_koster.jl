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

include("../materials/materials.jl")

function Deque_to_vec(Elist,Edict,vararg...;Transform=nothing)
    NN=length(Elist)
    R=zeros(Float64,NN,3)
    if isnothing(Transform)
        for Ri in Elist
            idx=Edict[Ri]
            R[idx,:]=Ri
        end
        return R
    else
        for Ri in Elist
            idx=Edict[Ri]
            R[idx,:]=Transform(Ri,vararg)
        end
    return R
    end
end

function getVAA(R, hops)
    V=Array{Float64,2}(undef,9,9)

    # Direction cosines
    l,m,n = R./norm(R)

    # basis ordering: s, xy, zx, yz,  x^2-y^2, 3z^2-r, x, y, z

    sp_sig, ss_sig, pp_sig, pp_pi, sd_sig, pd_sig, pd_pi, dd_sig, dd_pi, dd_del = hops

    # s-s interaction
    V[1,1] = ss_sig

    # p-p interaction !
    V[7,8] = l*m*pp_sig -     l*m*pp_pi # x -> y
    V[7,9] = l*n*pp_sig -     l*n*pp_pi # x -> z
    V[8,9] = m*n*pp_sig -     m*n*pp_pi # y -> z
    V[7,7] = l^2*pp_sig + (1-l^2)*pp_pi # x -> x
    V[8,8] = m^2*pp_sig + (1-m^2)*pp_pi # y -> y
    V[9,9] = n^2*pp_sig + (1-n^2)*pp_pi # z -> z
    V[8,7] = V[7,8]                     # y -> x
    V[9,7] = V[7,9]                     # z -> x
    V[9,8] = V[8,9]                     # z -> y

    # d1-d1 interaction !
    V[2,3] = 3*l^2*m*n*dd_sig +            m*n*(1-4*l^2)*dd_pi +   m*n*(l^2-1)*dd_del # xy -> zx
    V[3,4] = 3*n^2*m*l*dd_sig +            m*l*(1-4*n^2)*dd_pi +   l*m*(n^2-1)*dd_del # zx -> yz
    V[2,4] = 3*l*m^2*n*dd_sig +            l*n*(1-4*m^2)*dd_pi +   l*n*(m^2-1)*dd_del # xy -> yz
    V[2,2] = 3*l^2*m^2*dd_sig + ( l^2 + m^2 - 4*l^2*m^2)*dd_pi + (n^2+l^2*m^2)*dd_del # xy -> xy
    V[3,3] = 3*l^2*n^2*dd_sig + ( l^2 + n^2 - 4*l^2*n^2)*dd_pi + (m^2+l^2*n^2)*dd_del # zx -> zx
    V[4,4] = 3*m^2*n^2*dd_sig + ( m^2 + n^2 - 4*m^2*n^2)*dd_pi + (l^2+m^2*n^2)*dd_del # yz -> yz
    V[3,2] = V[2,3]
    V[4,2] = V[2,4]
    V[4,3] = V[3,4]

    # d1-d2 interaction !
    V[2,5] =                   (1.5*l*m*(l^2-m^2)*dd_sig +        2*l*m*(m^2-l^2)*dd_pi +         .5*l*m*(l^2-m^2)*dd_del) #    xy -> x2-y2
    V[2,6] =         (3^.5*l*m*(n^2-.5*(l^2+m^2))*dd_sig -         2*3^.5*l*m*n^2*dd_pi +      .5*3^.5*l*m*(1+n^2)*dd_del) #    xy -> 3z2-r2
    V[3,5] =                   (1.5*n*l*(l^2-m^2)*dd_sig +    n*l*(1-2*(l^2-m^2))*dd_pi -     n*l*(1-.5*(l^2-m^2))*dd_del) #    zx -> x2-y2
    V[3,6] =         (3^.5*l*n*(n^2-.5*(l^2+m^2))*dd_sig + 3^.5*l*n*(l^2+m^2-n^2)*dd_pi -    .5*3^.5*l*n*(l^2+m^2)*dd_del) #    zx -> 3z2-r2
    V[4,5] =                   (1.5*m*n*(l^2-m^2)*dd_sig -    m*n*(1+2*(l^2-m^2))*dd_pi +      m*n*(1+(l^2-m^2)/2)*dd_del) #    yz -> x2-y2
    V[4,6] =         (3^.5*m*n*(n^2-.5*(l^2+m^2))*dd_sig + 3^.5*m*n*(l^2+m^2-n^2)*dd_pi -    .5*3^.5*m*n*(l^2+m^2)*dd_del) #    yz -> 3z2-r2
    V[5:6,2:4] = transpose(V[2:4,5:6])

    # d2-d2 interaction 
    V[5,6] = .5*3^.5*(l^2-m^2)*(n^2-.5*(l^2+m^2))*dd_sig +     3^.5*n^2*(m^2-l^2)*dd_pi + 3^.5*(1+n^2)*(l^2-m^2)/4*dd_del # x2-y2 -> 3z2-r2
    V[5,5] =                      3/4*(l^2-m^2)^2*dd_sig +  (l^2+m^2-(l^2-m^2)^2)*dd_pi +      (n^2+(l^2-m^2)^2/4)*dd_del
    V[6,6] =                 (n^2-.5*(l^2+m^2))^2*dd_sig +        3*n^2*(l^2+m^2)*dd_pi +          3/4*(l^2+m^2)^2*dd_del
    V[6,5] = V[5,6]

    # s-d interaction
    V[1,2] =           3^.5*l*m*sd_sig
    V[1,3] =           3^.5*l*n*sd_sig
    V[1,4] =           3^.5*n*m*sd_sig
    V[1,5] =   3^.5/2*(l^2-m^2)*sd_sig
    V[1,6] = (n^2-.5*(l^2+m^2))*sd_sig
    V[2:6,1] = transpose(V[1,2:6])

    # s-p interaction
    V[1,7] = l*sp_sig
    V[1,8] = m*sp_sig
    V[1,9] = n*sp_sig
    V[7:9,1] = -transpose(V[1,7:9])

    # p-d interaction
    V[7,2] =         3^0.5*l^2*m*pd_sig +         m*(1-2*l^2)*pd_pi # x -> xy 
    V[8,4] =         3^0.5*m^2*n*pd_sig +         n*(1-2*m^2)*pd_pi # y -> yz
    V[9,3] =         3^0.5*n^2*l*pd_sig +         l*(1-2*n^2)*pd_pi # z -> zx

    V[7,4] =         3^0.5*l*m*n*pd_sig -             2*l*m*n*pd_pi # x -> yz
    V[8,3] =         3^0.5*l*m*n*pd_sig -             2*l*m*n*pd_pi # y -> zx
    V[9,2] =         3^0.5*l*m*n*pd_sig -             2*l*m*n*pd_pi # z -> xy

    V[7,3] =         3^0.5*l^2*n*pd_sig +         n*(1-2*l^2)*pd_pi # x -> zx
    V[8,2] =         3^0.5*m^2*l*pd_sig +         l*(1-2*m^2)*pd_pi # y -> xy
    V[9,4] =         3^0.5*n^2*m*pd_sig +         m*(1-2*n^2)*pd_pi # z -> yz

    V[7,5] = 3^0.5/2*l*(l^2-m^2)*pd_sig +       l*(1-l^2+m^2)*pd_pi # x -> x2-y2
    V[8,5] = 3^0.5/2*m*(l^2-m^2)*pd_sig -       m*(1+l^2-m^2)*pd_pi # y -> x2-y2
    V[9,5] = 3^0.5/2*n*(l^2-m^2)*pd_sig -         n*(l^2-m^2)*pd_pi # z -> x2-y2
    V[7,6] = l*(n^2-(l^2+m^2)/2)*pd_sig -         3^0.5*l*n^2*pd_pi # x -> 3z2-r2
    V[8,6] = m*(n^2-(l^2+m^2)/2)*pd_sig -         3^0.5*m*n^2*pd_pi # y -> 3z2-r2
    V[9,6] = n*(n^2-(l^2+m^2)/2)*pd_sig +   3^0.5*n*(l^2+m^2)*pd_pi # z -> 3z2-r2
    V[2:6,7:9] = -transpose(V[7:9,2:6])


    return V
end 

function get_Hk_FCC(kpoint, onsite, nearest, second)
    # 18 positions of the nearest (and next nearest) neighbors in cartesian coordinates?
    nnlist=[[1, 1, 0];;[-1, 1, 0];;[1,-1, 0];;[-1,-1, 0];;
            [1, 0, 1];;[-1, 0, 1];;[1, 0,-1];;[-1, 0,-1];;
            [0, 1, 1];;[ 0,-1, 1];;[0, 1,-1];;[ 0,-1,-1];;
            [2, 0, 0];;[-2, 0, 0];;[0, 2, 0];;[ 0,-2, 0];;[0,0,2];;[0,0,-2]]

    # Hk = Array{Float64,2}(undef,9,9)
    Hk = zeros(ComplexF64, 9, 9)

    for j=1:18
        dR    = nnlist[:,j]
        dist  = norm(dR)
        
        ph = dot(dR,kpoint)
        phase = exp(im*ph*pi)

        if dist < 1.99
            # println("dist<1.99")
            Hk += getVAA(dR, nearest).*phase
        else
            # println("dist>1.99")
            vv = getVAA(dR, second )
            # println("norm", norm(vv - vv'))

            Hk = Hk + vv.*phase
            # println(kpoint," ", j, " ", ph, " ", phase, " ", dR, " ", vv[1,1])
            # println("Diag of matrix:")
            # println(diag(Hk))
        end
    end
    # println("after orbitals:")
    # println(diag(Hk))

    for i=1:9
        Hk[i,i] += onsite[i]
    end

    return Hk

end

function get_bands(onsite, nearest, second)
    # Define the path in k space
    G = [0.0,  0.0,  0.0]
    X = [0.0,  1.0,  0.0]
    W = [0.5,  1.0,  0.0]
    L = [0.5,  0.5,  0.5]
    K = [0.75, 0.75, 0.0]

    path = [G, X, W, L, G, K]
    Npaths = length(path) - 1
    kpp = 30 # number of k points per path
    NK = (kpp-1)*Npaths + 1
    klist = zeros(Float64, NK,  3)
    line  = zeros(Float64, kpp, 3)

    # Create path
    for i in 1:Npaths
        Pi = path[i]
        Pf = path[i+1]
        # println(Pi, Pf)

        for d in 1:3
            line[:,d] = LinRange(Pi[d], Pf[d], kpp)
        end

        # println(line)

        n = (i-1)*(kpp-1) + 1
        klist[n:n+kpp-1,:] = line

    end
    klist[end,:] = path[end]
    # println(klist)

    # calculate eigenenergies for each path
    bands = zeros(Float64, NK, 9+4)
    for k in 1:NK
        kvec = klist[k,:]
        hk = Hermitian(get_Hk_FCC(kvec, onsite, nearest, second))
        # bands[k,4:end-1] = diag(hk)
        bands[k,4:end-1] = eigvals(hk)
        bands[k,end] = norm(hk)
        # println("---", bands[k,:])
    end

    bands[:,1:3] = klist
    return bands
end


function slater_koster_FCC(Elist, Edict, onsite, first_neighbour, second_neighbour, A, B; L1=999999, L2=999999, L3=999999, periodic=false)
    # Use the Slater-Koster parametrization for a set of atomic positions defined in Elist and Edict, 
    # using the parameters onsite, first_neighbour, second_neighbour
    # This is only valid for the FCC lattice

    # Get the tight-binding parameters
    # onsite, first_neighbour, second_neighbour, A, B, fermi, diel_fun = tightbinding(material)
    # println(onsite, first_neighbour, second_neighbour)

    # 18 positions of the nearest (and next nearest) neighbors in cartesian coordinates?
    nnlist=[[1, 1, 0];;[-1, 1, 0];;[1,-1, 0];;[-1,-1, 0];;
            [1, 0, 1];;[-1, 0, 1];;[1, 0,-1];;[-1, 0,-1];;
            [0, 1, 1];;[ 0,-1, 1];;[0, 1,-1];;[ 0,-1,-1];;
            [2, 0, 0];;[-2, 0, 0];;[0, 2, 0];;[ 0,-2, 0];;[0,0,2];;[0,0,-2]]

        
    function getHamiltonian(Elist::Deque{Vector{Int64}}, Edict::Dict{Vector{Int64},Int64},
                            nnlist::Array{Float64,2},    nnlist_for_table::Array{Int64,2},
                            getHopping::Function;      
                            Norbitals::Int64=9,          A::Float64=0.0,      B::Float64=1.0, L1::Int64=999999, L2::Int64=999999, L3::Int64=999999, periodic=false)

        # Elist  - list of all the positions of atoms inside the nanoparticle [(x1,y1,z1), (x2,y2,z2), ...]
        # Edict  - dictionary with all the sites inside the nanoparticle (x,y,z): index
        # nnlist - list of vectors to nearest all nearest neighbors (half of list below)
        # nnlist_for_table - list of vectors to nearest all nearest neighbors 
        # A - KPM shift
        # B - KPM scale

        NNA = length(Elist) # number of atoms
        onsiteKPM = (onsite.-A)./B # turning Ha into KPM units
        
        function addV!(iidx,jidx,value,V,iptr,jptr)
            # Add the values of the V 9x9 matrix into a sparse matrix format
            # iidx (jidx) is the set of (lines,columns) of the sparse matrix
            # iptr, jptr index the atoms
            local i = 0
            local j = 0
            local vv = 0.0
            local a,b = size(V)
            for i = 1:a
                for j = 1:b
                    vv = V[i,j]
                    if abs(vv) > 1e-10
                        push!(iidx, iptr+i)
                        push!(jidx, jptr+j)
                        push!(value, vv)
                    end 
                end 
            end 
        end 


        # Build a neighbor table for every particle
        # [[n1, n2, ..., n9],     <- indices of neighbors of atom 1
        # [[n1, -1, ..., n9],     <- indices of neighbors of atom 2. -1 means it was not found
        # ...
        global table=Array{Int64,2}(undef, NNA, 18)
        for Ri in Elist
            Idx = Edict[Ri]
            for j = 1:18
                Rij = nnlist_for_table[:,j]
                Rj  = Ri + Rij 
                # println(Rj,L1,L2,L3)

                # if it is periodic, then there are no negative positions, by construction
                if periodic
                    Rj[1] = (Rj[1]+L1) % L1
                    Rj[2] = (Rj[2]+L2) % L2
                    Rj[3] = (Rj[3]+L3) % L3
                end

                # println(Rj)
                Idj = get(Edict, Rj, -1)
                table[Idx,j] = Idj 
            end
        end 

        # println("table\n")
        # display(table)

        R = Deque_to_vec(Elist,Edict)

        # build the sparse Hamiltonian matrix from the neighbor table and the hoppings
        global iidx  = Array{Int64,1}([])
        global jidx  = Array{Int64,1}([])
        global value = Array{ComplexF64,1}([]) # hamiltonian

        global iidx_r  = Array{Int64,1}([])
        global value_r = Array{ComplexF64,1}([]) # position operator
        for Atomi=1:NNA
            iptr = (Atomi-1)*Norbitals

            # onsite terms
            for k=1:Norbitals
                push!(iidx, k+iptr)
                push!(jidx, k+iptr)
                push!(value, onsiteKPM[k])

                push!(iidx_r, k+iptr)
                push!(value_r, R[Atomi])
            end


            for j=1:18
                Atomj = table[Atomi,j]
                if Atomj >= 0
                    # Neighbouring atom is atomA
                    jptr = (Atomj-1)*Norbitals
                    Rij  = nnlist[:,j]
                    V    = getHopping(Rij)./B
                    addV!(iidx, jidx, value, V, iptr, jptr)
                end
            end
        end

        # mat = zeros(Float64,length(iidx), 3)
        # mat[:,1] = iidx
        # mat[:,2] = jidx
        # mat[:,3] = value
        # display(mat)
        H = sparse(iidx, jidx, value)
        pos_operator = sparse(iidx_r, iidx_r, value_r)
        v = (H*pos_operator - pos_operator*H).*im
        return H,v
    end 



    function getVAA2(R)
        V=Array{Float64,2}(undef,9,9)
        function getHoppingParameters(d)
            if d<0.9999999
                return first_neighbour
            else
                return second_neighbour
            end
        end

        hops = getHoppingParameters(norm(R))
        return getVAA(R, hops)
    end

    H,v = getHamiltonian(Elist, Edict, nnlist./2, nnlist, getVAA2, A=A, B=B, L1=L1, L2=L2, L3=L3, periodic=periodic)
    return H,v
end
