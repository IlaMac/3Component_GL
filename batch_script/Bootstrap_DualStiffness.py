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
fig, ax1 = plt.subplots(1, 1, figsize=(9,6))
ax1.set_title(r"$h=%s$; $e=%s$; $\nu=%s$" %(h, e, nu))
ax1.set_xlabel(r"$\beta$")
ax1.set_ylabel(r"$L\rho$")

for l in range(len(L)):
    BASEDIR=("/home/ilaria/Desktop/MultiComponents_SC/Output_2C/L%d_e%s_h%s_nu%s_bmin%s_bmax%s" %(L[l], e,  h, nu, beta_low, beta_high))

    Ds_mean=np.zeros((nbeta))
    Ds_err=np.zeros((nbeta))

    for b in range(nbeta):
        beta[b]=beta_low +b*(beta_high -beta_low)/(nbeta-1)

#        fileDs=("%s/beta_%d/Dual_Stiffness.npy" %(BASEDIR, b))
#        Ds=np.load(fileDs)

        file=h5py.File('%s/beta_%d/Output.h5' %(BASEDIR, b), 'r')
        Ds=np.asarray(file['Measurements']['ds'])

        #cut of the transient regime:
        Ds=Ds[transient_time:]
        
        Ds_mean[b]=np.mean(Ds)
        Ds_err[b]=np.sqrt(np.var(Ds)/(len(Ds)-1)) 
        

    ax1.plot(beta, Ds_mean, '-')
    ax1.errorbar(beta, Ds_mean, yerr=Ds_err, capsize=2,label="L=%s" %L[l])

ax1.legend(loc="best")
plt.tight_layout()
plt.show()



