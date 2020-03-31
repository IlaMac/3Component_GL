
import numpy as np
import sys
import os
import math

BASEDIR=("/home/cluster_users/x_ilaria/MultiComponents_SC")

SCRIPT_DIR=("%s/Multi_Components_GL/batch_script" %(BASEDIR))

L=8

############# Parameters of the Hamiltonian ---> HP_init.txt in a directory whose name contains the main parameters values##################
H_a=0
H_b=1
H_eta=1
H_e=0.5
H_h=5.0
H_blow=0.224
H_bhigh=0.2265

############Creation of the output folder and of the two files of initialization####################

OUTPUT=("%s/Output/L%d_a%d_b%d_eta%d_e%.1lf_h%.1lf_bmin%.3lf_bmax%.4lf" %(BASEDIR, L, H_a, H_b, H_eta, H_e, H_h, H_blow, H_bhigh))

nrank=32
files_name=['Energy', ' Magnetization', 'Dual_Stiffness', 'Psi_Density']

beta=np.zeros((nrank))

delta_beta=(H_bhigh-H_blow)/(nrank-1);

for b in range(nrank):

    beta[b]= H_blow +b*delta_beta
    BETA_DIR = ("%s/beta_%d" %(OUTPUT, b))
    if not (os.path.isdir(BETA_DIR)):
        os.makedirs(BETA_DIR)

    for n in range(nrank):
        RANK_DIR=("%s/rank_%d" %(OUTPUT, n))
        file_E=("%s/Energy.txt" %(RANK_DIR))
        En=np.loadtxt(file_E, usecols=(0,1,2,3,4,5,6), unpack=True)

        



