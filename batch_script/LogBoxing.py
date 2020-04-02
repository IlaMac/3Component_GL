import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import sys
import os
import math
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


#L=sys.argv[6].split(',')
#print(L)
#L=L.astype(np.int)

#the generic observable is called here A
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size='18')
plt.rc('text.latex', preamble=r'\usepackage{bm}')

#this are the tag used in writing the h5 file
Observables=["E", "m", "ds"]

transient_max=0

for name in range(len(Observables)):
    A_name=Observables[name]
    transient_list=np.zeros((nbeta, len(L)))
    for l in range(len(L)):
        BASEDIR=("/home/ilaria/Desktop/MultiComponents_SC/Output_3C/L%d_a0_b1_eta1_e%s_h%s_bmin%s_bmax%s" %(L[l], e,  h, beta_low, beta_high))

        for b in range(nbeta):
            beta[b]=beta_low +b*(beta_high -beta_low)/(nbeta-1)
            base=2
            exp=2
            box_length=base**exp
#            file=("%s/beta_%d/%s.npy" %(BASEDIR, b, A_name))
#            A=np.load(file)

            file=h5py.File('%s/beta_%d/Output.h5' %(BASEDIR, b), 'r')
            A=np.asarray(file['Measurements']['%s' %(Observables[name])])

            tot_length=box_length
            start=0
            A_mean=[]
            A_std=[]
            bins=[]
            while (tot_length< len(A)):
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
            ##### Find where the plateau starts by taking the minimum values of the derivative of the binned funciton ####
            A_diff=np.diff(A_mean)
            A_diffmin=np.min(np.sqrt(A_diff*A_diff))
            index=np.where(np.sqrt(A_diff*A_diff)==A_diffmin)
            #### The information to be extracted is the time up to the end of the bin where the plateau is observed: bin[index] ####            
            transient_list[b, l]=bins[index]

            file_Aout=("%s/beta_%d/Thermalization_%s.txt" %(BASEDIR, b, A_name))
            np.savetxt(file_Aout, np.transpose( [bins, A_mean, A_std]), fmt='%19.12e', delimiter=' ',  newline=os.linesep)

#            fig, ((ax1))= plt.subplots(1, 1)
#            ax1.set_title("beta=%s" %beta[b])
#            #ax1.plot(A_mean)
#            ax1.axvline(x=bins[index], color='red', linestyle='--')
#            ax1.errorbar(bins, A_mean,yerr=A_std, fmt='o')
#            ax1.grid(True)
#            ax1.set_xlabel(r'$t_{MC}$')
#            ax1.set_xscale("log")
#            ax1.set_ylabel("$%s$" %A_name)
#            plt.tight_layout()
#            plt.show()

    temp_transient_max=np.amax(transient_list)
    if(transient_max < temp_transient_max): transient_max=temp_transient_max

print(transient_max)

 


