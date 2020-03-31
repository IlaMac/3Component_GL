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
h=float(sys.argv[4])
e=sys.argv[5]

#tau_max=sys.argv[6]

tau_max=100

beta=np.zeros((nbeta))
if( (h).is_integer()): h=int(h)

L=np.array([8, 10, 12, 16])

block_size=30*tau_max
nblocks=30

varE_resampling=np.zeros((nblocks))

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('text.latex', preamble=r'\usepackage{bm}')
fig, ax1 = plt.subplots(1, 1)


for l in range(len(L)):
    BASEDIR=("/home/ilaria/Desktop/MultiComponents_SC/Output_3C/L%d_a0_b1_eta1_e%s_h%s_bmin%s_bmax%s" %(L[l], e,  h, beta_low, beta_high))
    Cv=np.zeros((nbeta))
    Cv_err=np.zeros((nbeta))

    for b in range(nbeta):
        beta[b]=beta_low +b*(beta_high -beta_low)/(nbeta-1)
        
        fileE=("%s/beta_%d/Energy.npy" %(BASEDIR, b))
        E=np.load(fileE)

        for block in range(nblocks):
                # <E²> - <E>²
                varE_resampling[block]=np.var(np.random.choice(E, size=block_size))
        
        Cv[b]=np.var(E)
        Cv_err[b]= np.sqrt(np.var(varE_resampling)/(nblocks-1))

    ax1.plot(beta, Cv, '-')
    ax1.errorbar(beta, Cv, yerr=Cv_err, capsize=2,label="L=%s" %L[l])

plt.legend(loc="best")
plt.show()



