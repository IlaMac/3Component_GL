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
import h5py

beta_low=float(sys.argv[1])
beta_high=float(sys.argv[2])
nbeta=int(sys.argv[3])
e=sys.argv[4]
h=float(sys.argv[5])
nu=float(sys.argv[6])

beta=np.zeros((nbeta))
if( (h).is_integer()): h=int(h)
if( (nu).is_integer()): nu=int(nu)


L=[]
for ind in range(7, len(sys.argv)):
    L.append(int(sys.argv[ind]))


Observables=["E", "m", "ds"]

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('text.latex', preamble=r'\usepackage{bm}')

tau_max=0
for l in range(len(L)):
    BASEDIR=("/Users/ilaria/Desktop/MultiComponents_SC/Output_3C/L%d_a0_b1_eta1_e%s_h%s_nu%s_bmin%s_bmax%s" %(L[l], e,  h, nu, beta_low, beta_high))
#    fig, ax1 = plt.subplots(1, 1)
#    ax1.set_title(r"$L=%s$" %L[l] )
#    ax1.set_xlabel(r"$\beta$")
#    ax1.set_ylabel(r"$\tau$")
    tau=np.zeros((nbeta, 3))
    for b in range(nbeta):
        beta[b]=beta_low +b*(beta_high -beta_low)/(nbeta-1)
        for name in range(len(Observables)):
            Obs_mean=np.zeros((nbeta))
            Obs_var=np.zeros((nbeta))
            file=h5py.File('%s/beta_%d/Output.h5' %(BASEDIR, b), 'r')
            print(file, b, l, Observables[name])
            A=np.asarray(file['Measurements']['%s' %(Observables[name])])
            #file=h5py.File('%s/beta_%d/Output.h5' %(BASEDIR, b), 'r')
            #Obs=np.asarray(file['Measurements']['%s' %(Observables[name])])
            #A_Obs=acf(Obs, fft=True)
            print(A)
            fig, ax1 = plt.subplots(1, 1)
            ax1.set_title(r"$L=%s; beta=%s$" %(L[l], beta[b]) )
            ax1.set_xlabel(r"$t$")
            ax1.set_ylabel(r"$Autocorr$")
            ax1.plot(A[:], "-")
            ax1.grid()
            plt.show()
            
#            A_Obs=acf(A, fft=True)
#            temp=np.where(A_Obs[:]<0.1)
#            print(temp)
#            tmax_int=10*temp[0][0]
#            temp_tau=[]
#            time_int=1000
#            tmax_int=max(time_int, tmax_int)
#      
#            if(len(A_Obs)<tmax_int): 
#                tmax_int= len(A_Obs)
#            if(len(A_Obs)<time_int):
#                time_int=len(A_Obs)
#
##            if(len(A_Obs[:time_int]) == len(A_Obs)): print("Before while", len(A_Obs[:time_int]), len(A_Obs), time_int)
#
#            while(time_int<= tmax_int):
##                temp_tau=np.append(temp_tau, np.trapz(A_Obs[:time_int], x=np.arange(time_int)))
#                temp_tau=np.append(temp_tau, np.sum(A_Obs[:time_int]))
#                if( np.sum(A_Obs[:time_int]) >1000): 
#                    print(beta[b], time_int)
#                    tstop=time_int
#                time_int=time_int+1000
#
#            tau[b, name]=np.amax(temp_tau)
#            temp_taumax=np.amax(tau)
#            if(temp_taumax > tau_max): tau_max=temp_taumax
##    ax1.plot(beta, tau[:,0], "-", label="$%s$" %Observables[0])
##    ax1.plot(beta, tau[:,1], "-", label="$%s$" %Observables[1])
##    ax1.plot(beta, tau[:,2], "-", label="$%s$" %Observables[2])
##
##    ax1.annotate(r' $\tau_{MAX}=%s$' %np.amax(tau), xy=(0.05, 0.85), xycoords='axes fraction', bbox=dict(boxstyle="round", edgecolor='orange', fc="w"))
##    ax1.legend(loc="best")
##    plt.tight_layout()
##    plt.show()
##    plt.savefig('%s/tau_L%s.png' %(BASEDIR, L[l]))
#
#print(tau_max)
