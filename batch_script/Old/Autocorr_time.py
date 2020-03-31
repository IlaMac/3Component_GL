
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


beta_low=float(sys.argv[1])
beta_high=float(sys.argv[2])
nbeta=int(sys.argv[3])
h=float(sys.argv[4])
e=sys.argv[5]

beta=np.zeros((nbeta))
if( (h).is_integer()): h=int(h)

#L=np.array([8, 10, 12, 16])

L=np.array([8])

Observables=["Energy", "Magnetization", "Dual_Stiffness"]

Observables=["Energy"]

#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')
#plt.rc('text.latex', preamble=r'\usepackage{bm}')

tau_max=0
for l in range(len(L)):
    BASEDIR=("/home/ilaria/Desktop/MultiComponents_SC/Output_3C/L%d_a0_b1_eta1_e%s_h%s_bmin%s_bmax%s" %(L[l], e,  h, beta_low, beta_high))

    fig, ax1 = plt.subplots(1, 1)
    ax1.set_title(r"$L=%s$" %L[l] )
    ax1.set_xlabel(r"$\beta$")
    ax1.set_ylabel(r"$\tau$")

    tau=np.zeros((nbeta, 3))

    for b in range(1):
        beta[b]=beta_low +b*(beta_high -beta_low)/(nbeta-1)

        for name in range(len(Observables)):
            Obs_mean=np.zeros((nbeta))
            Obs_var=np.zeros((nbeta))
            fileO=("%s/beta_%d/%s.npy" %(BASEDIR, b, Observables[name]))
            Obs=np.load(fileO)
            #A_Obs=acf(Obs, nlags=int(len(Obs)*0.5), fft=True)
            A_Obs=acf(Obs, nlags=(len(Obs)//10), unbiased=True, fft=True)
            A_Obs1=acf(Obs, nlags=(len(Obs)//10), fft=True)
            #temp=np.where(A_Obs[:]<0)
            #tmax_int=10*temp[0][0]
            #temp_tau=[]
            #time_int=1000
            #tmax_int=max(time_int, tmax_int)
            #while(time_int<= tmax_int):
            #    temp_tau=np.append(temp_tau, np.trapz(A_Obs[:time_int]))
            #    time_int=time_int+1000
            print(np.trapz(A_Obs))
            print(np.trapz(A_Obs1))
    
#            temp_tau=np.sum(A_Obs[:len(A_Obs)//2])

            ax1.plot(A_Obs)
            ax1.plot(A_Obs1)
            plt.show()
            #tau[b, name]=np.amax(temp_tau)
            #temp_taumax=np.amax(tau)
            #if(temp_taumax > tau_max): tau_max=temp_taumax
#    ax1.plot(beta, tau[:,0], "-", label="$%s$" %Observables[0])
#    ax1.plot(beta, tau[:,1], "-", label="$%s$" %Observables[1])
#    ax1.plot(beta, tau[:,2], "-", label="$%s$" %Observables[2])
#
#    ax1.annotate(r' $\tau_{MAX}=%s$' %np.amax(tau), xy=(0.05, 0.85), xycoords='axes fraction', bbox=dict(boxstyle="round", edgecolor='orange', fc="w"))
#    ax1.legend(loc="best")
#    plt.tight_layout()
#    plt.savefig('%s/tau_L%s.png' %(BASEDIR, L[l]))

print(tau_max)
