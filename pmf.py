import numpy as np
import matplotlib.pyplot as plt

"""
script used in this twitter thread discussion: 
https://twitter.com/agrossfield/status/1410252908248145923
"""

## generate collective variable r (CV1) and shifted CV rr=r+shift (CV2)
shift=1
r=np.linspace(0.5,4,100) # CV1
rr=r+shift # CV2

## use LJ as PMF
epsilon=1
sigma=1
x=sigma/r
PMF = 4*epsilon*(x**12 - x**6)

## provide min and max for integration
r_c = 3
rr_c = r_c+shift
r_min = 1
rr_min = r_min+shift

## povide dr
dr = r[2]-r[1]

## indexes of values we include in integration (size of binding site)
idx = np.where((r >= r_min) & (r <= r_c))
# note: indices are the same for CV1 and CV2

## use kT=1 for simplicity
kT=1

## calculate dissociation const and free energy for CV 1
K_r = 4*np.pi*sum(np.exp(-PMF[idx]/kT)*r[idx]**2)*dr
G_r = kT*np.log(K_r)

## calculate dissociation const and free energy for CV 2
K_rr = 4*np.pi*sum(np.exp(-PMF[idx]/kT)*rr[idx]**2)*dr
G_rr = kT*np.log(K_rr)

## note: difference in calculated G from choice of CV
print(G_r)
print(G_rr)

## plotting
plt.plot(r,PMF,color='blue',label='PMF with CV1 = r, K1 = %1.1f, G1 = %1.1f' % (K_r,G_r))
plt.plot(rr,PMF,color='red',label='PMF with CV2 = r+1, K2 = %1.1f, G2 = %1.1f' % (K_rr,G_rr))
plt.plot([0,4],[0,0],linestyle='--',color='grey')
plt.plot([0,4],[-epsilon,-epsilon],linestyle='--',color='grey')
plt.ylim(-1.25,0.5)
plt.legend()
plt.xlabel('CV1 (blue) or CV2 (red). CV2 = CV1+1')
plt.ylabel('Energy')
plt.show()
