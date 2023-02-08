using Interpolations

function palladium()
    # Palladium. Taken from page 237 of Papaconstantopolous's book
    # The factor of two takes the units from Rydberg to Hartree
    
    # Local energies at each orbital
    Es  = 0.94261/2
    Ep  = 1.36110/2
    Ed1 = 0.37285/2
    Ed2 = 0.36265/2
    onsite = [Es,Ed1,Ed1,Ed1,Ed2,Ed2,Ep,Ep,Ep]

    # First neighbour
    ss_sig = -0.07962/2
    pp_sig =  0.17119/2
    pp_pi  = -0.00540/2
    dd_sig = -0.05216/2
    dd_pi  =  0.02878/2
    dd_del = -0.00533/2
    sp_sig =  0.11332/2
    sd_sig = -0.04885/2
    pd_sig = -0.06563/2
    pd_pi  =  0.02124/2
    first_neighbour = [ss_sig, pp_sig, pp_pi, dd_sig, dd_pi, dd_del, sp_sig, sd_sig, pd_sig, pd_pi]

    # Second neighbour
    ss_sig = -0.00105/2
    pp_sig =  0.04282/2
    pp_pi  = -0.00044/2
    dd_sig = -0.00385/2
    dd_pi  =  0.00212/2
    dd_del = -0.00026/2
    sp_sig =  0.01048/2
    sd_sig = -0.00837/2
    pd_sig = -0.00738/2
    pd_pi  =  0.00351/2
    second_neighbour = [ss_sig, pp_sig, pp_pi, dd_sig, dd_pi, dd_del, sp_sig, sd_sig, pd_sig, pd_pi]

    # Frequency paladium (eV)
    frequencies = [0.1,0.15,0.2,0.26,0.3,0.36,0.4,0.46,0.5,0.56,0.6,0.72,0.8,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.2,4.4,4.6,4.8,5.0,5.2,5.4,5.6,5.8,6.0,6.2,6.4,6.6,6.8,7.0,7.2,7.4,7.6,7.8,8.0,8.2,8.4,8.6,8.8,9.0,9.5,10.0]

    # real part of dielectric constant of bulk paladium
    eps_real = [-2915.1656,-1273.2755,-697.6032,-396.3504,-285.5793,-191.8077,-157.86,-128.4192,-114.0636,-94.6737,-84.7616,-63.3699,-53.7411,-38.532,-33.8355,-30.1875,-27.1584,-24.7744,-22.5395,-20.1761,-18.2784,-16.9252,-15.6101,-14.4099,-13.452,-12.4944,-11.7216,-10.8712,-10.1223,-9.4127,-8.8201,-8.1328,-7.5933,-7.1307,-6.6825,-6.304,-5.8828,-5.5056,-5.1404,-4.8361,-4.5195,-4.2125,-3.96,-3.7352,-3.2865,-2.9792,-2.7027,-2.5899,-2.538,-2.3941,-2.1675,-1.9683,-1.7632,-1.5249,-1.3293,-1.1312,-0.9457,-0.744,-0.5696,-0.3979,-0.2352,-0.0984,0.0805,0.2415,0.3936,0.4316,0.5775,0.6591,0.7353,0.8319,0.8771]

    # imaginary part of dielectric constant of bulk paladium
    eps_imag = [447.279,224.2332,163.2626,125.333,122.9624,114.7036,113.3258,103.4194,93.808,82.2416,75.696,61.074,54.002,41.2022,36.3052,32.33,28.9,25.74,22.6548,20.592,18.88,17.4336,15.834,14.63,13.4602,12.416,11.475,10.6134,9.8136,9.2064,8.58,8.0754,7.6356,7.2324,6.84,6.4818,6.1104,5.824,5.544,5.292,5.0932,4.8972,4.725,4.5114,4.3472,4.1406,3.9964,3.838,3.5712,3.222,2.89,2.6244,2.4024,2.204,2.0276,1.8834,1.7424,1.6498,1.533,1.47,1.4014,1.343,1.2948,1.2848,1.316,1.344,1.3,1.352,1.3696,1.456,1.482]

    eps_real_interp = LinearInterpolation(frequencies, eps_real)
    eps_imag_interp = LinearInterpolation(frequencies, eps_imag)
    function dielectric(hw)
        return eps_real_interp(hw) + 1im*eps_imag_interp(hw)
    end

    # Lattice constant (length of FCC unit cell in nanometers) - 7.291 Bohr
    a0 = 0.389

    # KPM shift and scale
    A = 0.5
    B = 1.0
    fermi = 0.5190/2 # Fermi energy in Hartree
    return [onsite, first_neighbour, second_neighbour, A, B, fermi, a0, dielectric]
end
