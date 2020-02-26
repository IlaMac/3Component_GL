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
Psi1=np.zeros((nbeta))
Psi2=np.zeros((nbeta))
Psi3=np.zeros((nbeta))

L=8
h=5
N=L*L*L
V=N*h*h*h

for b in range(nbeta):
    print(b)
    beta[b]=beta_low +b*((beta_high-beta_low)/(nbeta-1))
    print(beta[b])
    file_Psi=("beta_%d/Psi_density.txt" %b)
    p1, p2, p3=np.loadtxt(file_Psi, usecols=(1,2,3), unpack=True)
    Psi1[b]=np.mean(p1)
    Psi2[b]=np.mean(p2)
    Psi3[b]=np.mean(p3)

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size='15')
plt.rc('text.latex', preamble=r'\usepackage{bm}')
fig, ((ax1, ax2, ax3))= plt.subplots(1, 3)
ax1.set_xlabel(r'$\beta$')
ax1.set_ylabel(r'$\Psi_1$')
ax1.plot(beta, Psi1, 'o-')
ax2.set_xlabel(r'$\beta$')
ax2.set_ylabel(r'$\Psi_2$')
ax2.plot(beta, Psi2, 'o-')
ax3.set_xlabel(r'$\beta$')
ax3.set_ylabel(r'$\Psi_3$')
ax3.plot(beta, Psi3, 'o-')
plt.subplots_adjust(top=0.964,bottom=0.115,left=0.076,right=0.98,hspace=0.14,wspace=0.4)
plt.tight_layout()
plt.show()
