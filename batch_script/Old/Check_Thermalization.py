import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import sys
import os
import math
from statsmodels.graphics.tsaplots import plot_acf
import statsmodels.api as sm
from statsmodels.tsa.stattools import acf


#L, BASEDIR, nb= input("Insert the size L, the working folder and the number of temperatures: ").split()

L=10
BASEDIR="/home/ilaria/Desktop/MultiComponents_SC/Output/L8_a0_b1_eta1_e0.5_h5.0_bmin0.224_bmax0.2265"
nb=1


#Create an array of length nbeta for each observable A

BASEDIR=("/home/ilaria/Desktop/MultiComponents_SC/Output/L%d_a0_b1_eta1_e0.5_h5.4_bmin0.2245_bmax0.229")


#Files of the observables measured for each temperature
for b in range(nbeta):
    beta=H_blow +b*((H_bhigh -H_blow)/nbeta)
    fileM=("%s/beta_%d/Magnetization.txt" %(BASEDIR, beta))
    M, M2, M4=np.loadtxt(fileM, usecols=(1,2,3), unpack=True)
    print(len(M))
    N=M[:10000]
    print(len(N))

    hist, bin_edges = np.histogram(M)
    PM = sm.nonparametric.KDEUnivariate(M)
    PM.fit()

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.rc('text.latex', preamble=r'\usepackage{bm}')
    fig, ((ax1, ax2, ax3)) = plt.subplots(3, 1)
    ax1.set_title(r"$\beta=%lf$" %beta)
    ax1.set_xlabel("t")
    ax1.set_ylabel("M(t)")
    ax1.plot(N)
    ax2.set_xlabel("t")
    ax2.set_ylabel("$A_M(t)$")
    ax2.plot(acf(N,  nlags=100, fft=False), 'o-') 
    ax3.set_xlabel("M")
    ax3.set_ylabel("P(M)")
    ax3.hist(M, bins=20, density=True)
    ax3.plot(PM.support, PM.density)
    fig.tight_layout()
    plt.savefig('%s/beta_%d/Check_M.png' %(BASEDIR, b))
    plt.close()
