import matplotlib.pyplot as plt
import numpy as np

file = "dist.dat"
data = np.loadtxt(file)

E = data[:,0]

dist = data[:,1]
dist_3Ne = data[:,2] # triple the number of integration points
dist_fer = data[:,3] # change the fermi energy
dist_hw  = data[:,4] # change the frequency
dist_kbT = data[:,5] # change the temperature
dist_gam = data[:,6] # change the broadening
dist_ph2 = data[:,7] # change the number of polynomials
dist_ph3 = data[:,8] # change the number of polynomials again

hcg     = data[:,9]
hcg_ph2 = data[:,10]
hcg_ph3 = data[:,11]

fig, axs = plt.subplots(2,2,figsize=(10,5))
ax = axs[0,0]
ax.set_title("dist parameters")
ax.plot(E,dist,     label = "orig")
ax.plot(E,dist_3Ne, label = "triple int")
ax.plot(E,dist_fer,'--', label = "fermi")
ax.plot(E,dist_kbT, label = "kbT")
ax.plot(E,dist_hw, '--', label = "hw")
ax.legend()

ax = axs[0,1]
ax.set_title("dist polynomials")
ax.plot(E, dist, label="100%")
ax.plot(E, dist_ph2, label="70%")
ax.plot(E, dist_ph3, label="50%")
ax.legend()

ax = axs[1,0]
ax.set_title("dist broadening")
ax.plot(E, dist, label="original")
ax.plot(E, dist_gam, label="double")
ax.legend()

ax = axs[1,1]
ax.set_title("hcg polynomials")
ax.plot(E, hcg, label="100%")
ax.plot(E, hcg_ph2, label="70%")
ax.plot(E, hcg_ph3, label="50%")
ax.legend()


plt.show()


