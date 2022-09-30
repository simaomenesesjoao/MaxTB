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
goldA=5e-1
goldB=6e-1


hw_pd_data=[0.1,0.15,0.2,0.26,0.3,0.36,0.4,0.46,0.5,0.56,0.6,0.72,0.8,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.2,4.4,4.6,4.8,5.0,5.2,5.4,5.6,5.8,6.0,6.2,6.4,6.6,6.8,7.0,7.2,7.4,7.6,7.8,8.0,8.2,8.4,8.6,8.8,9.0,9.5,10.0]./HaeV
e1_pd_data=[-2915.1656,-1273.2755,-697.6032,-396.3504,-285.5793,-191.8077,-157.86,-128.4192,-114.0636,-94.6737,-84.7616,-63.3699,-53.7411,-38.532,-33.8355,-30.1875,-27.1584,-24.7744,-22.5395,-20.1761,-18.2784,-16.9252,-15.6101,-14.4099,-13.452,-12.4944,-11.7216,-10.8712,-10.1223,-9.4127,-8.8201,-8.1328,-7.5933,-7.1307,-6.6825,-6.304,-5.8828,-5.5056,-5.1404,-4.8361,-4.5195,-4.2125,-3.96,-3.7352,-3.2865,-2.9792,-2.7027,-2.5899,-2.538,-2.3941,-2.1675,-1.9683,-1.7632,-1.5249,-1.3293,-1.1312,-0.9457,-0.744,-0.5696,-0.3979,-0.2352,-0.0984,0.0805,0.2415,0.3936,0.4316,0.5775,0.6591,0.7353,0.8319,0.8771]
e2_pd_data=[447.279,224.2332,163.2626,125.333,122.9624,114.7036,113.3258,103.4194,93.808,82.2416,75.696,61.074,54.002,41.2022,36.3052,32.33,28.9,25.74,22.6548,20.592,18.88,17.4336,15.834,14.63,13.4602,12.416,11.475,10.6134,9.8136,9.2064,8.58,8.0754,7.6356,7.2324,6.84,6.4818,6.1104,5.824,5.544,5.292,5.0932,4.8972,4.725,4.5114,4.3472,4.1406,3.9964,3.838,3.5712,3.222,2.89,2.6244,2.4024,2.204,2.0276,1.8834,1.7424,1.6498,1.533,1.47,1.4014,1.343,1.2948,1.2848,1.316,1.344,1.3,1.352,1.3696,1.456,1.482]
hw_au_data=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.1,2.2,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,5.0,5.2,5.4,5.6,5.8,6.0,6.2,6.4,6.6,6.8,7.0,7.2,7.4,7.6,7.8,8.0,8.2,8.4,8.6,8.8,9.0,9.2,9.4,9.6,9.8,10.0]./HaeV
e1_au_data=[-6794.06,-1736.856,-772.9723,-433.5408,-275.74,-189.81,-138.0141,-104.2117,-81.1576,-64.464,-42.7616,-29.5872,-20.7872,-14.5843,-9.9687,-8.0332,-6.394,-3.2096,-1.8557,-0.834,-0.9135,-1.0013,-0.9541,-0.8684,-0.9045,-0.8684,-0.7659,-0.664,-0.4965,-0.3685,-0.3729,-0.5472,-0.748,-0.8736,-0.924,-0.978,-1.0304,-1.1095,-1.2051,-1.1661,-1.0767,-0.9656,-0.8587,-0.7317,-0.4959,-0.3048,-0.1488,0.0241,0.1888,0.3495,0.46,0.5928,0.7491,0.8931,0.9471,0.9828,0.944,0.8424,0.7945,0.8436,0.8979,0.9417,0.9765,1.0011,1.06,1.1605,1.2412,1.2528,1.2369]
e2_au_data=[1353.4422,177.7698,55.0836,24.5794,12.9558,7.7168,5.17,3.6756,2.703,2.0878,1.308,0.8704,0.7296,0.6876,0.8216,1.0224,1.2192,1.86,2.6076,3.8192,4.9192,5.1684,5.37,5.544,5.5748,5.544,5.518,5.4912,5.4668,5.6052,5.74,5.8354,5.7558,5.611,5.4058,5.2688,5.133,4.9632,4.698,4.394,4.1656,3.975,3.7884,3.6356,3.388,3.2186,3.0734,2.904,2.7816,2.7032,2.625,2.5654,2.522,2.546,2.584,2.6496,2.7048,2.673,2.5152,2.392,2.314,2.3056,2.2532,2.158,2.1222,2.0748,2.1216,2.1646,2.192]
eps_m=2.0
e1_pd_interp=LinearInterpolation(hw_pd_data,e1_pd_data)
e2_pd_interp=LinearInterpolation(hw_pd_data,e2_pd_data)
e1_pd=e1_pd_interp(hw)
e2_pd=e2_pd_interp(hw)
eps_pd=e1_pd+1im*e2_pd

e1_au_interp=LinearInterpolation(hw_au_data,e1_au_data)
e2_au_interp=LinearInterpolation(hw_au_data,e2_au_data)
e1_au=e1_au_interp(hw)
e2_au=e2_au_interp(hw)
eps_au=e1_au+1im*e2_au




function compute_potential_sphere(R,Edict)
    # Get the electric potential inside the sphere
    #
    a,b=size(R)
    value_aux=zeros(ComplexF64,a)
    # Potential is proportional to z
    for i=1:a
        RR=R[i,:]
        value_aux[i]=-3*eps_m/(eps_pd+2*eps_m)*RR[3]*a_0/2
    end 

    iidx=1:length(Edict)*9|>collect
    jidx=1:length(Edict)*9|>collect

    value=zeros(ComplexF64,length(Edict)*9)
    for i=1:length(Edict)
        for j=1:9
            value[(i-1)*9+j]=value_aux[i]
        end 
    end 
    Phi=sparse(iidx,jidx,value)
    return Phi
end

function comsol_read(filename;Norbitals=9,Nomean=true,factor=1.0)
    # Read the potential from a COMSOL file
    f=open(filename,"r")

    # Ignore the first 9 lines (header)
    for i=1:9
        😅=readline(f)
    end

    flines=readlines(f)
    φ=Vector{ComplexF64}()
    mean_aux=0+0im
    for i in flines
        if length(i)>3
        j=split(i)
        φi=parse(ComplexF64,j[4])*factor
        push!(φ,φi)
        mean_aux=mean_aux+φi
        end
    end

    # Remove mean, makes simulations better
    if Nomean
        mean=mean_aux/length(φ)
        φ=φ.-mean
    end 

    # Put potential in the form of a sparse matrix
    iidx=zeros(Int64,length(φ)*Norbitals)
    jidx=zeros(Int64,length(φ)*Norbitals)
    value=zeros(ComplexF64,length(φ)*Norbitals)
    for i=1:lastindex(φ)
        for j=1:Norbitals
            iidx[(i-1)*Norbitals+j]=(i-1)*Norbitals+j
            jidx[(i-1)*Norbitals+j]=(i-1)*Norbitals+j
            value[(i-1)*Norbitals+j]=φ[i]
        end 
    end 

    Φ  = sparse(iidx, jidx, value)
    ΦT = sparse(iidx, jidx, conj(value))
    return Φ,ΦT
end

