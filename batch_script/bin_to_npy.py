import numpy as np
import sys
import os



beta_low=float(sys.argv[1])
beta_high=float(sys.argv[2])
nbeta=int(sys.argv[3])
h=float(sys.argv[4])
e=sys.argv[5]

beta=np.zeros((nbeta))
if( (h).is_integer()): h=int(h)


L=np.array([8, 10, 12])


############### Case of different output binary files ###############

#for l in range(len(L)):
#    BASEDIR=("/home/ilaria/Desktop/MultiComponents_SC/Output_3C/L%d_a0_b1_eta1_e%s_h%s_bmin%s_bmax%s" %(L[l], e,  h, beta_low, beta_high))
#    for b in range(nbeta):
#        
#        beta[b]=beta_low +b*((beta_high-beta_low)/(nbeta-1))
#        
#        file_DS=("%s/beta_%d/Dual_Stiffness" %(BASEDIR, b))
#        Ds=np.fromfile(file_DS + ".bin", dtype=np.float64)
#        np.save(file_DS + ".npy", Ds)
#
#        file_M=("%s/beta_%d/Magnetization" %(BASEDIR, b))
#        M=np.fromfile(file_M + ".bin", dtype=np.float64)
#        np.save(file_M + ".npy", M)
#        
#        file_E=("%s/beta_%d/Energy" %(BASEDIR, b))
#        E=np.fromfile(file_E + ".bin", dtype=np.float64)
#        np.save(file_E + ".npy", E)
#
#        file_Psi=("%s/beta_%d/Psi_density" %(BASEDIR, b))
#        Psi=np.fromfile(file_Psi + ".bin", dtype=np.float64)
#        Psi=np.reshape(Psi, (len(Psi)/3, 3))
#        np.save(file_Psi + ".npy", Psi)


############### Case of just one output binary file ###############

for l in range(len(L)):
    BASEDIR=("/home/ilaria/Desktop/MultiComponents_SC/Output_2C/L%d_e%s_h%s_bmin%s_bmax%s" %(L[l], e,  h, beta_low, beta_high))
    for b in range(nbeta):
        
        beta[b]=beta_low +b*((beta_high-beta_low)/(nbeta-1))

        file_Out=("%s/beta_%d/Output.bin" %(BASEDIR, b))
        Out=np.fromfile(file_Out, dtype=np.float64)
#        open(file_Out, "rb")
#        Out=file_Out.read()
        print(len(Out)/5)

        # E M DS Psi1 Psi2 -> 5 columns
        Out=np.reshape(Out, (int(len(Out)/5),5))
#        E=Out[:,0]
#        M=Out[:,1]
#        Ds=Out[:,2]
#        Psi1=Out[:,3]
#        Psi2=Out[:,4]
#        Psi=np.column_stack((Psi1, Psi2))
#
#        file_M=("%s/beta_%d/Magnetization" %(BASEDIR, b))
#        np.save(file_M + ".npy", M)
#        
#        file_E=("%s/beta_%d/Energy" %(BASEDIR, b))
#        np.save(file_E + ".npy", E)
#
#        file_DS=("%s/beta_%d/Dual_Stiffness" %(BASEDIR, b))
#        np.save(file_DS + ".npy", Ds)
#
#        file_Psi=("%s/beta_%d/Psi_density" %(BASEDIR, b))
#        np.save(file_Psi + ".npy", Psi)


