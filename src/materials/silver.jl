using Interpolations

function silver()
    # Gold 2-center approximation Slater-Koster parameters. Taken from page 242 of Papaconstantopolous's book
    # The factor of two takes the units from Rydberg to Hartree. All units are in Rydberg
    # The dielectric constant is taken from the CRC Handbook of Chemistry and Physics 95th edition, page 12.144
    
    # Local energies at each orbital (Ry)
    Es  = 0.68297
    Ep  = 1.13432
    Ed1 = 0.12249
    Ed2 = 0.12006
    onsite = [Es,Ed1,Ed1,Ed1,Ed2,Ed2,Ep,Ep,Ep]

    # First neighbour (Ry)
    ss_sig = -0.06581
    pp_sig =  0.15752
    pp_pi  =  0.00649
    dd_sig = -0.03151
    dd_pi  =  0.01757
    dd_del = -0.00336
    sp_sig =  0.09781
    sd_sig = -0.03110
    pd_sig = -0.03905
    pd_pi  =  0.01519
    first_neighbour = [sp_sig, ss_sig, pp_sig, pp_pi, sd_sig, pd_sig, pd_pi, dd_sig, dd_pi, dd_del]

    # Second neighbour (Ry)
    ss_sig =  0.00143
    pp_sig =  0.03971
    pp_pi  =  0.00434
    dd_sig = -0.00282
    dd_pi  =  0.00171
    dd_del = -0.00038
    sp_sig =  0.00545
    sd_sig = -0.00462
    pd_sig = -0.00065
    pd_pi  =  0.00172
    second_neighbour = [sp_sig, ss_sig, pp_sig, pp_pi, sd_sig, pd_sig, pd_pi, dd_sig, dd_pi, dd_del]

    # Convert from Rydberg to Hartree
    onsite = onsite./2
    first_neighbour = first_neighbour./2
    second_neighbour = second_neighbour./2

    # Frequency silver (eV)
    frequencies = [0.10, 0.20, 0.30, 0.40, 0.50, 1.00, 1.50, 2.00, 2.50, 3.00, 3.25, 3.50, 3.60, 3.70, 3.77, 3.80, 3.90, 4.00, 4.10, 4.20, 4.30, 4.50, 4.75, 5.00, 5.50, 6.00, 6.50, 7.00, 7.50, 8.00, 9.00, 10.00, 11.00, 12.00, 13.00, 14.00, 14.50, 15.00, 16.00, 17.00, 18.00, 19.00, 20.00, 21.00, 21.50, 22.00, 22.50, 23.00, 23.50, 24.00, 24.50, 25.00, 25.50, 26.00, 26.50, 27.00, 27.50, 28.00, 28.50, 29.00, 30.00, 31.00, 32.00, 33.00, 34.00, 35.00, 36.00, 38.00, 40.00, 42.00, 44.00, 46.00, 48.00, 50.00, 52.00, 54.00, 56.00, 58.00, 60.00, 62.00, 64.00, 66.00, 68.00, 70.00, 72.00, 74.00, 76.00, 78.00, 80.00, 85.00, 90.00, 95.00, 100.00]

    # real part of dielectric constant of bulk silver
    eps_real = [-8050.4648, -2080.4244, -928.8720, -523.1240, -335.1735, -81.4625, -33.4512, -17.3995, -9.4905, -5.1000, -3.4067, -1.9723, -1.2240, -0.5029, 0.1209, 0.4429, 1.5604, 2.2321, 2.2704, 1.9389, 1.7160, 1.2177, 0.7965, 0.5529, 0.3069, 0.1572, 0.1701, 0.2688, 0.4715, 0.7831, 1.4553, 1.8180, 1.9968, 2.2440, 2.3460, 2.3500, 1.9152, 1.5872, 1.1883, 1.0293, 0.9984, 1.0504, 1.1600, 1.2600, 1.2369, 1.0387, 0.7227, 0.4853, 0.3451, 0.2716, 0.2232, 0.2136, 0.2197, 0.2624, 0.3160, 0.3696, 0.4077, 0.4619, 0.5032, 0.5328, 0.5733, 0.5840, 0.5655, 0.5499, 0.5343, 0.5371, 0.5985, 0.6400, 0.6731, 0.6875, 0.7011, 0.7076, 0.6960, 0.6903, 0.7137, 0.7455, 0.6893, 0.6993, 0.7085, 0.7303, 0.7303, 0.7303, 0.7128, 0.6489, 0.6901, 0.6936, 0.6969, 0.7000, 0.7029, 0.7104, 0.7161, 0.7360, 0.7553]

    # imaginary part of dielectric constant of bulk silver
    eps_imag = [1789.1514, 259.5760, 86.0382, 41.6598, 24.5488, 5.0568, 3.1266, 2.2572, 1.4832, 1.0442, 0.8556, 0.5964, 0.5198, 0.4620, 0.4240, 0.4380, 0.9360, 1.9320, 2.9410, 3.7100, 3.9098, 4.3264, 4.3148, 4.2160, 3.8860, 3.4304, 2.9500, 2.5016, 2.0748, 1.7400, 1.4896, 1.6352, 1.7024, 1.8998, 2.1248, 2.6832, 2.8864, 2.8704, 2.5844, 2.2876, 2.0480, 1.9050, 1.8318, 2.0250, 2.1920, 2.3316, 2.3436, 2.1996, 2.0460, 1.8720, 1.7226, 1.5770, 1.4196, 1.3320, 1.2282, 1.1570, 1.1036, 1.0620, 1.0374, 1.0304, 1.0044, 0.9858, 0.9752, 0.9180, 0.8624, 0.7740, 0.7832, 0.6942, 0.6660, 0.6300, 0.5940, 0.5760, 0.5518, 0.5104, 0.4984, 0.2992, 0.4524, 0.4176, 0.3828, 0.3696, 0.3696, 0.3696, 0.3654, 0.3320, 0.3060, 0.2890, 0.2720, 0.2550, 0.2380, 0.1870, 0.1360, 0.1032, 0.0696]

    eps_real_interp = LinearInterpolation(frequencies, eps_real)
    eps_imag_interp = LinearInterpolation(frequencies, eps_imag)
    function dielectric(hw)
        return eps_real_interp(hw) + 1im*eps_imag_interp(hw)
    end

    # Lattice constant (length of FCC unit cell in nanometers) - 7.625 Bohr
    a0 = 0.40349 

    # KPM shift and scale
    A = 0.5
    B = 1.0
    fermi = 0.4635/2 # Fermi energy in Hartree
    return [onsite, first_neighbour, second_neighbour, A, B, fermi, a0, dielectric] 
end
