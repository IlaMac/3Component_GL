import numpy as np
#import matplotlib.pyplot as plt
#import matplotlib.colors as mcolors
import sys
import os
import math

BASEDIR=str(sys.argv[1])
A_name=str(sys.argv[2])
nbeta=int(sys.argv[3])

print(nbeta)

base=2
exp=9

box_length=base**exp
print(box_length)
#the generic observable is called here A


for b in range(nbeta):

    file_A=("%s/beta_%d/%s.txt" %(BASEDIR, b, A_name))
    A=np.loadtxt(file_A, usecols=1, unpack=True)
    
    tot_length=box_length
    start=0
    A_mean=[]
    A_std=[]
    bins=[]
    while (tot_length< len(A)):
    
        A_mean.append=np.mean(A[start:tot_length])
        A_std.append=np.std(A[start:tot_length)], axes=0)/np.sqrt(box_length -1)
        bins.append=tot_length

        print(box_length)
        print(tot_length -start)

        start+=box_length
        exp+=1
        box_length=base**exp
        tot_length+=box_length
    A_mean=np.array(A_mean)
    A_std=np.array(A_std)
    bins=np.array(bins)
    file_Aout=("%s/beta_%d/Thermalization_%s.txt" %(BASEDIR, b, A_name))
    np.savetxt(file_Aout, np.transpose( [bins, A_mean, A_std]), fmt='%19.12e', delimiter=' ',  newline=os.linesep)

       


