import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import sys
import os
import math
from statsmodels.graphics.tsaplots import plot_acf
import statsmodels.api as sm
from statsmodels.tsa.stattools import acf
import scipy.integrate as integrate
import random
import h5py


beta_low=float(sys.argv[1])
beta_high=float(sys.argv[2])
nbeta=int(sys.argv[3])
e=float(sys.argv[4])
h=float(sys.argv[5])
nu=float(sys.argv[6])


beta=np.zeros((nbeta))
if( (h).is_integer()): h=int(h)
if( (nu).is_integer()): nu=int(nu)

transient_time=float(sys.argv[7])
tau_max=float(sys.argv[8])

transient_time=int(transient_time)
tau_max=int(tau_max)


L=[]
for ind in range(9, len(sys.argv)):
    L.append(int(sys.argv[ind]))

block_size=20*tau_max

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=18)
plt.rc('text.latex', preamble=r'\usepackage{bm}')
fig, ax1 = plt.subplots(2, 1, figsize=(9,12))
ax1[0].set_title(r"$h=%s$; $e=%s$; $\nu=%s$" %(h, e, nu))
ax1[0].set_xlabel(r"$\beta$")
ax1[0].set_ylabel(r"$E/V$")
ax1[1].set_xlabel(r"$\beta$")
ax1[1].set_ylabel(r"$C_{v}$")

for l in range(len(L)):
    BASEDIR=("/home/ilaria/Desktop/MultiComponents_SC/Output_2C/L%d_e%s_h%s_nu%s_bmin%s_bmax%s" %(L[l], e,  h, nu, beta_low, beta_high))

    Cv_mean=np.zeros((nbeta))
    Cv_err=np.zeros((nbeta))
    E_mean=np.zeros((nbeta))
    E_err=np.zeros((nbeta))

    for b in range(nbeta):
        beta[b]=beta_low +b*(beta_high -beta_low)/(nbeta-1)
        
        #fileE=("%s/beta_%d/Energy.npy" %(BASEDIR, b))
        #E=np.load(fileE)

        file=h5py.File('%s/beta_%d/Output.h5' %(BASEDIR, b), 'r')
        E=np.asarray(file['Measurements']['E'])

        #cut of the transient regime:
        E=E[transient_time:]
#        print(len(E))
        
        E_mean[b]=np.mean(E)/(L[l]**3)
        E_err[b]=np.sqrt(np.var(E/(L[l]**3))/(len(E)-1)) 
        
        nblocks=int(len(E)/block_size)
#        print(nblocks)

        varE_resampling=np.zeros((nblocks))

        for block in range(nblocks):
                # <E²> - <E>²
                varE_resampling[block]=np.var(np.random.choice(E, size=block_size))
        
        Cv_mean[b]=beta[b]*np.var(E)/(L[l]**3)
        Cv_err[b]= (beta[b]/(L[l]**3))*np.sqrt(np.var(varE_resampling)/(nblocks-1))

    ax1[0].plot(beta, E_mean, '-')
    ax1[0].errorbar(beta, E_mean, yerr=E_err, capsize=2,label="L=%s" %L[l])
    ax1[1].plot(beta, Cv_mean, '-')
    ax1[1].errorbar(beta, Cv_mean, yerr=Cv_err, capsize=2)

ax1[0].legend(loc="best")
plt.tight_layout()
plt.show()



