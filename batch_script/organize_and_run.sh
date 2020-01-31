#!/bin/bash

BASEDIR=${HOME}/MultiComponents_SC
SCRIPT_DIR=${BASEDIR}/Multi_Components_GL/batch_script

LLIST="10 12 16"

############# Parameters of the Hamiltonian ---> HP_init.txt in a directory whose name contains the main parameters values##################
H_a=0
H_b=1
H_eta=1 
H_e=0.5
H_h=5.4
H_blow=0.2245
H_bhigh=0.229

############ Parameters for the Monte Carlo simulations --> MC_init.txt#####################

Nmisu=1000000
ntau=20
nautosave=100000
l_box=1.0
rho_box=0.5
theta_box=3.141592653
A_box=0.1

for L in $LLIST; do

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
nnodes=2
ntasks=64 #parallel tempering over ntasks temperatures

#I create ntasks folder: one for each rank.

cd ${OUTPUT}

for ((rank=0; rank<${ntasks}; rank++)); do

if [ ! -d ./Sbeta_${rank} ]; then
   mkdir -p beta_${rank}
fi

done

cd ${SCRIPT_DIR}
DIR_PAR="${OUTPUT}"

#SEED= If I want to repeat exactly a simulation I could initialize the random number generator exactly at the same way

EXECUTE_DIR="../build/Release"

echo "#!/bin/bash
#SBATCH --job-name=${jobname}          # Name of the job
#SBATCH --time=7-00:00:00               # Allocation time
#SBATCH --mem-per-cpu=2000              # Memory per allocated cpu
#SBATCH --nodes=${nnodes}               # Number of nodes
#SBATCH --ntasks=${ntasks}
#SBATCH --clusters=draken
#SBATCH --output=${DIR_PAR}/logs/log_${jobname}.o
#SBATCH --error=${DIR_PAR}/logs/log_${jobname}.e


srun ${EXECUTE_DIR}/GL_3components ${L} ${DIR_PAR} &> ${DIR_PAR}/logs/log_${jobname}.o


" >  submit_run

#Submission of the work --> sbatch submit_runs

mkdir -p ${DIR_PAR}/logs

sbatch submit_run

done
