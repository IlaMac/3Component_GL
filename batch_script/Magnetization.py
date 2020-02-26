import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import sys
import os
import math
from astropy.stats import jackknife_resampling

nbeta=64
beta_low=0.2245
beta_high=0.229
beta=np.zeros((nbeta))
M=np.zeros((nbeta))
U=np.zeros((nbeta))
L=np.array([8, 10, 12, 16])
h=5.4
V=(L*h)**3

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size='18')
plt.rc('text.latex', preamble=r'\usepackage{bm}')
fig, ((ax1, ax2))= plt.subplots(2, 1)

for l in range(len(L)):
    BASEDIR=("/home/ilaria/Desktop/MultiComponents_SC/Output/L%d_a0_b1_eta1_e0.5_h5.4_bmin0.2245_bmax0.229" %L[l])

    for b in range(nbeta):
        beta[b]=beta_low +b*((beta_high-beta_low)/(nbeta-1))
        file_M=("%s/beta_%d/Magnetization.txt" %(BASEDIR, b))
        Magn, Magn2, Magn4=np.loadtxt(file_M, usecols=(1,2,3), unpack=True)
        Half=int(0.5*len(Magn))
        M[b]=np.mean(Magn[:Half])
        U[b]=np.mean(Magn4[:Half])/(3*np.mean(Magn2[:Half])**2)

    ax1.plot(beta, M, '-', label=str(L[l]))
    ax2.plot(beta, U, '-', label=str(L[l]))

ax1.legend(loc='best')
ax1.grid(True)
ax1.set_xlabel(r'$\beta$')
ax1.set_ylabel(r'$m$')
ax2.grid(True)
ax2.set_xlabel(r'$\beta$')
ax2.set_ylabel(r'$U$')
plt.tight_layout()
plt.show()
