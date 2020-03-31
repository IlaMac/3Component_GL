import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import sys
import os
import math
from astropy.stats import jackknife_resampling

nbeta=4
beta_low=0.244
beta_high=0.247
beta=np.zeros((nbeta))
DS=np.zeros((nbeta))
DS_dev=np.zeros((nbeta))
L=8
h=5
N=L*L*L
V=N*h*h*h

for b in range(nbeta):
    print(b)
    beta[b]=beta_low +b*((beta_high-beta_low)/(nbeta-1))
    print(beta[b])
    file_DS=("beta_%d/Dual_Stiffness.txt" %b)
    DualS=np.loadtxt(file_DS, usecols=1, unpack=True)
    DS[b]=np.mean(DualS[-2000:])
    DS_dev[b]=np.std(DualS[-2000:])/np.sqrt(len(DualS[-2000:])-1)

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size="15")
plt.rc('text.latex', preamble=r'\usepackage{bm}')
#fig, ((ax1, ax2))= plt.subplots(2, 1)
fig, ((ax1))= plt.subplots(1, 1)
ax1.grid(True)
ax1.set_xlabel(r'$\beta$')
ax1.set_ylabel(r'$h\rho$')
ax1.plot(beta,h*DS, 'o-')
#ax2.set_xlabel(r'$\beta$')
#ax2.set_ylabel(r'$dev(\rho)$')
#ax2.plot(beta, DS_dev, 'o-')
plt.tight_layout()
plt.show()
