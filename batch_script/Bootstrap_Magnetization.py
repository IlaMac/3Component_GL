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
ax1[0].set_ylabel(r"$M/V$")
ax1[1].set_xlabel(r"$\beta$")
ax1[1].set_ylabel(r"$U$")

for l in range(len(L)):
    BASEDIR=("/home/ilaria/Desktop/MultiComponents_SC/Output_2C/L%d_e%s_h%s_nu%s_bmin%s_bmax%s" %(L[l], e,  h, nu,  beta_low, beta_high))

    U_mean=np.zeros((nbeta))
    U_err=np.zeros((nbeta))
    M_mean=np.zeros((nbeta))
    M_err=np.zeros((nbeta))

    for b in range(nbeta):
        beta[b]=beta_low +b*(beta_high -beta_low)/(nbeta-1)

#        fileM=("%s/beta_%d/Magnetization.npy" %(BASEDIR, b))
#        M=np.load(fileM)

        file=h5py.File('%s/beta_%d/Output.h5' %(BASEDIR, b), 'r')
        M=np.asarray(file['Measurements']['m'])

        #cut of the transient regime:
        M=M[transient_time:]
#        print(len(E))
        
        M_mean[b]=np.mean(M)/(L[l]**3)
        M_err[b]=np.sqrt(np.var(M/(L[l]**3))/(len(M)-1)) 
        
        nblocks=int(len(M)/block_size)
#        print(nblocks)

        U_resampling=np.zeros((nblocks))

        for block in range(nblocks):
                # U=<m⁴>/(3<m²>²)
                Mblock=np.random.choice(M, size=block_size)
                U_resampling[block]=np.mean(np.power(Mblock,4))/(3*np.power(np.mean(np.power(Mblock,2)),2))
        
        U_mean[b]=np.mean(U_resampling)
        U_err[b]= np.sqrt(np.var(U_resampling)/(nblocks-1))

    ax1[0].plot(beta, M_mean, '-')
    ax1[0].errorbar(beta, M_mean, yerr=M_err, capsize=2,label="L=%s" %L[l])
    ax1[1].plot(beta, U_mean, '-')
    ax1[1].errorbar(beta, U_mean, yerr=U_err, capsize=2)

ax1[0].legend(loc="best")
plt.tight_layout()
plt.show()



