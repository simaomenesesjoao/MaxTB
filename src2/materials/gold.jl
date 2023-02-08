using Interpolations

function gold()
    # Gold. Taken from page 298 of Papaconstantopolous's book
    # The factor of two takes the units from Rydberg to Hartree
    
    # Local energies at each orbital
    Es  = 0.56220/2
    Ep  = 1.27897/2
    Ed1 = 0.26097/2
    Ed2 = 0.25309/2
    onsite = [Es,Ed1,Ed1,Ed1,Ed2,Ed2,Ep,Ep,Ep]

    # First neighbour
    ss_sig = -0.06680/2
    pp_sig =  0.17866/2
    pp_pi  = -0.01645/2
    dd_sig = -0.04971/2
    dd_pi  =  0.02624/2
    dd_del = -0.00457/2
    sp_sig =  0.09721/2
    sd_sig = -0.04722/2
    pd_sig = -0.06399/2
    pd_pi  =  0.01896/2
    first_neighbour = [ss_sig, pp_sig, pp_pi, dd_sig, dd_pi, dd_del, sp_sig, sd_sig, pd_sig, pd_pi]

    # Second neighbour
    ss_sig =  0.00277/2
    pp_sig =  0.03707/2
    pp_pi  = -0.01025/2
    dd_sig = -0.00305/2
    dd_pi  =  0.00240/2
    dd_del = -0.00057/2
    sp_sig =  0.00261/2
    sd_sig = -0.00784/2
    pd_sig = -0.00762/2
    pd_pi  =  0.00470/2
    second_neighbour = [ss_sig, pp_sig, pp_pi, dd_sig, dd_pi, dd_del, sp_sig, sd_sig, pd_sig, pd_pi]

    # Frequency gold (eV)
    frequencies = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.1,2.2,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,5.0,5.2,5.4,5.6,5.8,6.0,6.2,6.4,6.6,6.8,7.0,7.2,7.4,7.6,7.8,8.0,8.2,8.4,8.6,8.8,9.0,9.2,9.4,9.6,9.8,10.0]

    # real part of dielectric constant of bulk gold
    eps_real = [-6794.06,-1736.856,-772.9723,-433.5408,-275.74,-189.81,-138.0141,-104.2117,-81.1576,-64.464,-42.7616,-29.5872,-20.7872,-14.5843,-9.9687,-8.0332,-6.394,-3.2096,-1.8557,-0.834,-0.9135,-1.0013,-0.9541,-0.8684,-0.9045,-0.8684,-0.7659,-0.664,-0.4965,-0.3685,-0.3729,-0.5472,-0.748,-0.8736,-0.924,-0.978,-1.0304,-1.1095,-1.2051,-1.1661,-1.0767,-0.9656,-0.8587,-0.7317,-0.4959,-0.3048,-0.1488,0.0241,0.1888,0.3495,0.46,0.5928,0.7491,0.8931,0.9471,0.9828,0.944,0.8424,0.7945,0.8436,0.8979,0.9417,0.9765,1.0011,1.06,1.1605,1.2412,1.2528,1.2369]

    # imaginary part of dielectric constant of bulk gold
    eps_imag = [1353.4422,177.7698,55.0836,24.5794,12.9558,7.7168,5.17,3.6756,2.703,2.0878,1.308,0.8704,0.7296,0.6876,0.8216,1.0224,1.2192,1.86,2.6076,3.8192,4.9192,5.1684,5.37,5.544,5.5748,5.544,5.518,5.4912,5.4668,5.6052,5.74,5.8354,5.7558,5.611,5.4058,5.2688,5.133,4.9632,4.698,4.394,4.1656,3.975,3.7884,3.6356,3.388,3.2186,3.0734,2.904,2.7816,2.7032,2.625,2.5654,2.522,2.546,2.584,2.6496,2.7048,2.673,2.5152,2.392,2.314,2.3056,2.2532,2.158,2.1222,2.0748,2.1216,2.1646,2.192]

    eps_real_interp = LinearInterpolation(frequencies, eps_real)
    eps_imag_interp = LinearInterpolation(frequencies, eps_imag)
    function dielectric(hw)
        return eps_real_interp(hw) + 1im*eps_imag_interp(hw)
    end

    # Lattice constant (length of FCC unit cell in nanometers) - 7.72 Bohr
    a0 = 0.408 

    # KPM shift and scale
    A = 0.5
    B = 1.0
    fermi = 0.5380/2 # Fermi energy in Hartree
    return [onsite, first_neighbour, second_neighbour, A, B, fermi, a0, dielectric] 
end
