using DataStructures
using LinearAlgebra
using Printf

function print_positions(R, filename)
    a,b = size(R)
    f=open(filename,"w")
    for i=1:a
        x,y,z = R[i,:]
        write(f,@sprintf("%f %f %f\n",x,y,z))
    end
    close(f)
    
end

function print_positions_comsol(R, filename; rescale::Bool=false)
    # Prints out the positions in a format that is nice to comsol
    # The positions are rescaled so that the shape lies within a cube 
    # of length 1 centered at the origin
    a,b = size(R)
    max_pos = maximum(abs.(R))
    println("Max position is", max_pos)

    header = "% Model:              octa2.mph\n\
    % Version:            COMSOL 5.5.0.359\n\
    % Date:               Sep 29 2022, 16:37\n\
    % Dimension:          3\n\
    % Nodes:              3134\n\
    % Expressions:        1\n\
    % Description:        Electric potential\n\
    % Length unit:        nm\n\
    % x                       y                        z\n"

    f=open(filename,"w")
    write(f, header)
    for i=1:a
        x,y,z = R[i,:]

        # COMSOL expects the nanoparticles to lie within a cube
        # of side 1 centered at the origin, so these points need
        # to be rescaled for that purpose
        if rescale
            x,y,z = R[i,:]/max_pos*0.5
        end
        write(f,@sprintf("%f\t\t%f\t\t%f\n",x,y,z))
    end
    close(f)
    
end



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






function in_shape_planes(R,vararg)
    planelist = vararg[1]
    # Checks if R is within a set of planes
    # planelist is an array of dimensions [N,3] where N is the number of planes
    # the planes are defined through a vector 'p_vector', which is orthogonal to the plane
    # if the projection of R into p_vector is smaller than R, then R is inside

    n,d = size(planelist) # rows (number of vectors), columns

    for i=1:n
        p_vector = planelist[i,:]
        r2 = norm(p_vector)^2

        if dot(R, p_vector) > r2
            return false
        end
    end 
    return true
end 

function in_shape_sphere(R,vararg)
    # Checks if the atom is inside the sphere
    Rmax=vararg[1]
    return (norm(R)-1e-6<Rmax)
end


function explore(nnlist::Array{Int64,2}, in_shape::Function, vararg...; start_pos::Array{Int64,1}=[0,0,0])
    # Finds all the atomic positions belonging to the nanoparticle

    #=
    nnlist: Array{Int64,3,NN} the Neighbouring direction, where NN is the number of neighbouring directions
    in_shape: ans=in_shape(R,vararg) a function that determine whether a given position R is in shape or not, vararg are 
    input variables
    start_point: a vector that defines the starting position
    =#
    N_nearneigh=18 # number of nearest neighbors
    let 
    Udict=Dict{Vector{Int64},Int64}() # dictionary of sites whose neighbors have not been checked yet. 
                                      # Has to be empty by the end of algorithm
    Edict=Dict{Vector{Int64},Int64}() # dictionary with all the sites inside the nanoparticle (x,y,z): index
    
    
    Ulist=Deque{Vector{Int64}}()
    Elist=Deque{Vector{Int64}}() # list of all the sites inside the nanoparticle [(x,y,z), ...]

    idx = 1
    Udict[start_pos] = idx
    push!(Ulist, start_pos) # add the first position to the back of the list
    while (length(Ulist) > 0)

        Ri   = popfirst!(Ulist) # remove element from the front of the list
        idx2 = pop!(Udict,Ri) # remove this element from the dictionary

        Edict[Ri] = idx2 
        push!(Elist, Ri)

        for j=1:N_nearneigh
            Rij = nnlist[:,j] # find all the nearest and next nearest neighbors of this atom
            Rj  = Ri + Rij # calculate their positions
            exists = haskey(Udict,Rj) # check if it already exists in the dictionary
            exists = exists || haskey(Edict,Rj)

            # if it does not exist yet and is inside the nanoparticle, add it
            if ((!exists) && in_shape(Rj,vararg))
                idx = idx + 1
                Udict[Rj] = idx 
                push!(Ulist, Rj)
            end 

        end 
        # println("Edict", Edict, Elist)
        # println("Udict", Udict, Ulist)
        # println("\n")
    end

    return Elist,Edict
    end
end 



function generate_shape_FCC_lat(shape_name, radius, length1=nothing, length2=nothing)
    # Assumes FCC lattice, no units
    # shape_name is one of the following: sphere, cube, octahedron, rhombic, prism, cylinder
    # the sphere, cube, octahedron and rhombic dodecahedron only require one argument (denoted radius)
    # but the prism and the cylinder require two (denoted by radius and length1)

    # 18 positions of the nearest (and next nearest) neighbors in cartesian coordinates?
    nnlist=[[1,1,0];;[-1,1,0];;[1,-1,0];;[-1,-1,0];;
            [1,0,1];;[-1,0,1];;[1,0,-1];;[-1,0,-1];;
            [0,1,1];;[0,-1,1];;[0,1,-1];;[0,-1,-1];;
            [2,0,0];;[-2,0,0];;[0,2,0];;[0,-2,0];;[0,0,2];;[0,0,-2]]

    # Create the atomic positions of the nanoparticle. Sphere by default
    func = in_shape_sphere
    arg  = radius

    if shape_name == "sphere"
        func = in_shape_sphere

    elseif shape_name == "cube"
        func = in_shape_planes
        l = radius/2 # radius is the length of square edges

        # make sure that the floor of l is divisible by 2, so that only the (200) faces appear
        lt = floor(Int, l)
        if lt % 2 != 0
            l += 1.0
        end
        
        # 6 planes defining the cube centered at 0,0,0, with faces parallel to xy, xz and yz
        arg = [l 0 0; -l 0 0 ; 0 l 0 ; 0 -l 0 ; 0 0 l ; 0 0 -l]

    elseif shape_name == "octahedron"
        func = in_shape_planes

        # radius is the length of the octahedron's triangle edges
        s = radius/sqrt(18)

        # 8 planes defining the octahedron centered at (0,0,0)
        arg = [s -s -s; -s s -s; -s -s -s; s s -s; s s s; s -s s; -s s s; -s -s s]

    elseif shape_name == "rhombic"
        func = in_shape_planes

        # radius is the length of the edges of the rhombic dodecahedron
        s = radius/sqrt(3)

        arg = [s s 0; s 0 s; 0 s s; -s -s 0; -s 0 -s; 0 -s -s; 0 -s s; -s 0 s; -s s 0; 0 s -s; s 0 -s; s -s 0]

    elseif shape_name == "periodic"
        func = in_shape_planes
        l = radius
        t = -0.01
        
        # 6 planes defining the cube cornered at 0,0,0, with faces parallel to xy, xz and yz
        arg = [l 0 0; t 0 0 ; 0 l 0 ; 0 t 0 ; 0 0 l ; 0 0 t]

    elseif shape_name == "prism"
	func = in_shape_planes
	l = radius
	r = length1
	h = length2

        # 6 planes defining the prism centered at 0,0,0, with faces parallel to xy, xz and yz
        arg = [l 0 0; -l 0 0 ; 0 r 0 ; 0 -r 0 ; 0 0 h ; 0 0 -h]

    else
        println("Shape "*shape_name*" not supported")
    end

    Elist, Edict = explore(nnlist, func, arg; start_pos = [0,0,0])

    return Elist, Edict
end



function Transform(R,vararg)
    # displacement function
    return R
end 

function generate_shape_FCC(shape_name, a0, radius, length1=nothing, length2=nothing, Transform=nothing)
    # Converts the units to lattice units, obtains positions in those units and
    # then converts back. 
    # a0 is the length of the FCC cube, a is half
    # in units of a, all atomic positions become integers
    # 'Transform' can be used to modify the positions according to some displacement function
    a = a0/2
    if length2 == nothing
	if length1 == nothing
    	    Elist, Edict = generate_shape_FCC_lat(shape_name, radius/a)
	else
    	    Elist, Edict = generate_shape_FCC_lat(shape_name, radius/a, length1/a)
    	end
    else
    	Elist, Edict = generate_shape_FCC_lat(shape_name, radius/a, length1/a, length2/a)
    end

    R = Deque_to_vec(Elist, Edict; Transform=Transform)*a
    return Elist, Edict, R
end


