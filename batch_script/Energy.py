import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import sys
import os
import math
from astropy.stats import jackknife_resampling

nbeta=64
beta_low=0.15
beta_high=0.25
beta=np.zeros((nbeta))

L=np.array([8])
E=np.zeros((nbeta))
E_var=np.zeros((nbeta))

print(L)
e=0.5
h=5.4
V=(L*h)**3

fig, ((ax1, ax2))= plt.subplots(2, 1)
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size='18')
plt.rc('text.latex', preamble=r'\usepackage{bm}')

for l in range(len(L)):

    BASEDIR=("/home/ilaria/Desktop/MultiComponents_SC/Output_2C/L%d_e0.5_h5.4_bmin%s_bmax%s" %(L[l], beta_low, beta_high))

    for b in range(nbeta):
        beta[b]=beta_low +b*((beta_high-beta_low)/(nbeta-1))
        file_E=("%s/beta_%d/Energy.txt" %(BASEDIR, b))
        En=np.loadtxt(file_E, usecols=1, unpack=True)
        Half=int(0.5*len(En))
        E[b]=np.mean(En[Half:])*beta[b]/V[l]
        E_var[b]=np.var(En[Half:]*beta[b]/V[l])

    ax1.plot(beta, E, '-', label=str(L[l]))
    ax2.plot(beta, E_var, '-', label=str(L[l]))

ax1.legend(loc='best')
ax1.grid(True)
ax1.set_xlabel(r'$\beta$')
ax1.set_ylabel(r'$\beta E/V$')
ax2.set_xlabel(r'$\beta$')
ax2.grid(True)
ax2.set_ylabel(r'$var(\beta E/V)$')
plt.tight_layout()

plt.show()
