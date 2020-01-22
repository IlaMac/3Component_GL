#!/bin/bash

BASEDIR=${HOME}/MultiComponents_SC
EXECUTE_DIR=${BASEDIR}/Multi_Components_GL/build/Release
SCRIPT_DIR=${BASEDIR}/Multi_Components_GL/batch_script

L=8

############# Parameters of the Hamiltonian ---> HP_init.txt in a directory whose name contains the main parameters values##################
H_a=0
H_b=1
H_eta=1 
H_e=0.5
H_h=5.0
H_blow=0.224
H_bhigh=0.2265

############ Parameters for the Monte Carlo simulations --> MC_init.txt#####################

Nmisu=1000000
ntau=1
nautosave=100000
l_box=1.0
rho_box=0.5
theta_box=3.141592653
A_box=0.1


############Creation of the output folder and of the two files of initialization####################

cd ${BASEDIR}/Output

if [ ! -d ./SL${L}_a${H_a}_b${H_b}_eta${H_eta}_e${H_e}_h${H_h}_bmin${H_blow}_bmax${H_bhigh} ]; then
   mkdir -p L${L}_a${H_a}_b${H_b}_eta${H_eta}_e${H_e}_h${H_h}_bmin${H_blow}_bmax${H_bhigh}
fi

OUTPUT=${BASEDIR}/Output/L${L}_a${H_a}_b${H_b}_eta${H_eta}_e${H_e}_h${H_h}_bmin${H_blow}_bmax${H_bhigh}

cd ${OUTPUT}

#THE ORDER OF WRITING DOES MATTER
echo $H_a > HP_init.txt
echo $H_b >> HP_init.txt
echo $H_eta >> HP_init.txt
echo $H_e >> HP_init.txt
echo $H_h >> HP_init.txt
echo $H_blow >> HP_init.txt
echo $H_bhigh >> HP_init.txt

#THE ORDER OF WRITING DOES MATTER
echo $Nmisu > MC_init.txt
echo $ntau >> MC_init.txt
echo $nautosave >> MC_init.txt
echo $l_box >> MC_init.txt
echo $rho_box >> MC_init.txt
echo $theta_box >> MC_init.txt
echo $A_box >> MC_init.txt

#################Creation of the submit_runs script#########################

jobname="L${L}_a${H_a}_b${H_b}_eta${H_eta}_e${H_e}_h${H_h}_bmin${H_blow}_bmax${H_bhigh}"
nnodes=1
ntasks=32 #parallel tempering over 32 temperatures

#I create ntasks folder: one for each rank.

cd ${OUTPUT}

for ((rank=0; rank<${ntasks}; rank++)); do

if [ ! -d ./Srank_${rank} ]; then
   mkdir -p rank_${rank}
fi

done

DIR_PAR=${OUTPUT}
#SEED= If I want to repeat exactly a simulation I could initialize the random number generator exactly at the same way

echo "#!/bin/bash
#SBATCH --job-name=${jobname}          # Name of the job
#SBATCH --time=7-00:00:00               # Allocation time
#SBATCH --mem-per-cpu=2000              # Memory per allocated cpu
#SBATCH --nodes=${nnodes}               # Number of nodes
#SBATCH --ntasks=${ntasks}
#SBATCH --output=${SCRIPT_DIR}/logs/log_${jobname}.o
#SBATCH --error=${SCRIPT_DIR}/logs/log_${jobname}.e

mkdir -p ${SCRIPT_DIR}/logs

###srun ${EXECUTE_DIR}/GL_3components ${DIR_PAR} &> ${SCRIPT_DIR}/logs/log_${jobname}.o


" >  ${SCRIPT_DIR}/submit_run


#Submission of the work --> sbatch submit_runs

sbatch ${SCRIPT_DIR}/submit_run
