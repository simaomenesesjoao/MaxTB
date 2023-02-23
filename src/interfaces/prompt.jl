include("nanoparticle.jl")



function run()

    # Properties of the nanoparticle (see documentation for options)
    shape  = "cube"
    mater  = "gold"
    rad_nm = 1.5 # [nm] dimension of nanoparticle

    hb = hamiltonian_builder_predefined(shape, rad_nm, mater)


    # Properties of the external perturbation
    freq_eV = 2.4 # [eV] frequency of optical potential
    eps_m   = 2.0 # relative dielectric constant of the environment (default 1)
    eps     = -200+3im # relative dielectric constant of the nanoparticle (default prebuild)

    pb = potential_builder_sphere(freq_eV, eps=eps, eps_m=eps_m) 
    # pb = potential_builder_import(freq_eV, filename)
    

    # Select the calculation: density of states
    N  = 200       # number of Chebyshev polynomials (default 200)
    NK = 4         # number of random vectors (default 5)
    NE = 1000      # number of energy points (default 1000)
    minE_Ha = -1.0 # minimum energy in Hartree (default -2)
    maxE_Ha = +1.0 # maximum energy in Hartree (default +2)
    dosname = "dos.dat" # output file name (default "dos.dat")

    dosb = calculation_builder_dos(N=N, NE=NE, minE_Ha=minE_Ha, maxE_Ha=maxE_Ha, NK=NK, name=dosname)


    # Select the calculation: band structure
    bandname = "bands.dat" # output file name (default "bands.dat")

    bb = calculation_builder_bands(name=bandname)


    # Select the calculation: hot carrier generation
    T     = 298.0 # [K] temperature (default 298 K)
    gph   = 0.01  # [eV] lifetime (default 10meV)
    conv  = true  # include convergence data in the output (default false)
    N     = 200   # number of Chebyshev polynomials (default 200)
    NR    = 5     # number of random vectors (default 5)
    NB    = 10    # number of Chebyshev polynomial blocks (default 10)
    NE    = 1000  # number of output energies (default 1000)
    hname = "hcg.dat" # hot carrier generation output filename (default hcg.dat)

    gb = calculation_builder_hcg(temp=T, gph=gph, N=N, NR=NR, NB=NB, NE=NE, conv=conv, name=hname)


    # Select the calculation: distribution
    T     = 298.0 # [K] temperature (default 298 K)
    conv  = true  # include convergence data in the output (default false)
    N     = 200   # number of Chebyshev polynomials (default 200)
    NR    = 5     # number of random vectors (default 5)
    NB    = 10    # number of Chebyshev polynomial blocks (default 10)
    NE    = 1000  # number of output energies (default 1000)
    dname = "dist.dat" # hot carrier generation output filename (default dist.dat)

    # Relaxation rates in eV
    function Gamma(E,Ep)
        return 0.1
    end

    db = calculation_builder_dist(temp=T, Gamma=Gamma, N=N, NR=NR, NB=NB, NE=NE, conv=conv, name=dname)

    nanoparticle_run(ham_builder=hb, pot_builder=pb, hcg_builder=gb, dist_builder=db, dos_builder=dosb, band_builder=bb)

end
run()
