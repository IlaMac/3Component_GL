import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import sys
import os
import math

BASEDIR=str(sys.argv[1])
A_name=str(sys.argv[2])
nbeta=int(sys.argv[3])

#the generic observable is called here A
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size='18')
plt.rc('text.latex', preamble=r'\usepackage{bm}')

for b in range(nbeta):

    base=2
    exp=2
    box_length=base**exp

    file_A=("%s/beta_%d/%s.txt" %(BASEDIR, b, A_name))
    A=np.loadtxt(file_A, usecols=1, unpack=True)
    
    tot_length=box_length
    start=0
    A_mean=[]
    A_std=[]
    bins=[]
    while (tot_length< len(A)):
        
        print(start)
        print(tot_length)
        print(len(A))
        A_mean.append(np.mean(A[start:tot_length]))
        A_std.append(np.std(A[start:tot_length]))
        bins.append(tot_length)

        start=tot_length
        exp+=1
        box_length=base**exp
        tot_length+=box_length
    A_mean=np.array(A_mean)
    A_std=np.array(A_std)
    bins=np.array(bins)
    print(A_mean)
    file_Aout=("%s/beta_%d/Thermalization_%s.txt" %(BASEDIR, b, A_name))
    np.savetxt(file_Aout, np.transpose( [bins, A_mean, A_std]), fmt='%19.12e', delimiter=' ',  newline=os.linesep)
    
    fig, ((ax1))= plt.subplots(1, 1)
    ax1.set_title("beta "+str(b))
    ax1.errorbar(bins, A_mean,yerr=A_std, fmt='o')
    ax1.grid(True)
    ax1.set_xlabel(r'$t_{MC}$')
    ax1.set_xscale("log")
    ax1.set_ylabel("$%s$" %A_name)
    plt.tight_layout()
    plt.show()

  


