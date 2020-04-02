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
fig, ax1 = plt.subplots(1, 3, constrained_layout=True, figsize=(18,9))
fig.suptitle(r"$h=%s$; $e=%s$; $\nu=%s$" %(h, e, nu))
ax1[0].set_xlabel(r"$\beta$")
ax1[0].set_ylabel(r"$|\Psi_1|$")
ax1[1].set_xlabel(r"$\beta$")
ax1[1].set_ylabel(r"$|\Psi_2|$")
ax1[2].set_xlabel(r"$\beta$")
ax1[2].set_ylabel(r"$|\Psi_3|$")

for l in range(len(L)):
    BASEDIR=("/home/ilaria/Desktop/MultiComponents_SC/Output_3C/L%d_a0_b1_eta1_e%s_h%s_bmin%s_bmax%s" %(L[l], e,  h, beta_low, beta_high))

    Psi1_mean=np.zeros((nbeta))
    Psi1_err=np.zeros((nbeta))
    Psi2_mean=np.zeros((nbeta))
    Psi2_err=np.zeros((nbeta))
    Psi3_mean=np.zeros((nbeta))
    Psi3_err=np.zeros((nbeta))

    for b in range(nbeta):
        beta[b]=beta_low +b*(beta_high -beta_low)/(nbeta-1)
       
        file=h5py.File('%s/beta_%d/Output.h5' %(BASEDIR, b), 'r')
        Psi=np.asarray(file['Measurements']['rho'])
        print(len(Psi))
        
        if( (len(Psi)/3).is_integer() ):     
            Psi=np.reshape(Psi, ((int(len(Psi)/3, 3))))
        else: 
            print("Error len(Psi)/3 not integer!") 
            sys.exit()

        Psi1=Psi[:,0]
        Psi2=Psi[:,1]
        Psi3=Psi[:,2]

        #cut of the transient regime:
        Psi1=Psi1[transient_time:]
        Psi2=Psi2[transient_time:]
        Psi3=Psi3[transient_time:]

        Psi1_mean[b]=np.mean(Psi1)
        Psi1_err[b]=np.sqrt(np.var(Psi1)/(len(Psi1)-1)) 
        Psi2_mean[b]=np.mean(Psi2)
        Psi2_err[b]=np.sqrt(np.var(Psi2)/(len(Psi2)-1))
        Psi3_mean[b]=np.mean(Psi3)
        Psi3_err[b]=np.sqrt(np.var(Psi3)/(len(Psi3)-1))

    ax1[0].plot(beta, Psi1_mean, '-')
    ax1[0].errorbar(beta, Psi1_mean, yerr=Psi1_err, capsize=2,label="L=%s" %L[l])
    ax1[1].plot(beta, Psi2_mean, '-')
    ax1[1].errorbar(beta, Psi2_mean, yerr=Psi2_err, capsize=2,label="L=%s" %L[l])
    ax1[2].plot(beta, Psi3_mean, '-')
    ax1[2].errorbar(beta, Psi3_mean, yerr=Psi3_err, capsize=2,label="L=%s" %L[l])


ax1[0].legend(loc="best")
plt.show()



