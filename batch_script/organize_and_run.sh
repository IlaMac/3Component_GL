#!/bin/bash

BASEDIR=${HOME}/MultiComponents_SC
SCRIPT_DIR=${BASEDIR}/3Component_GL//batch_script

cd /tmp/

if [ ! -d ./SOutput_x_ilaria ]; then
   mkdir -p Output_x_ilaria
fi

#RESTART=0-> Restart from scratch
#RESTART=1-> Restart from interrupted run
#RESTART=2-> Restart from the previois final scenario

RESTART=0

LLIST="8 10 12 16 20"
#LLIST="8 10 12"
############ Parameters of the Hamiltonian ---> HP_init.txt in a directory whose name contains the main parameters values##################
H_a=0
H_b=1
H_eta=1 
H_e=0
H_h=1
H_nu=5
H_blow=10.5
H_bhigh=13.5
H_init=0 #If H_init=0: phases initialized to zero; H_init=1: phase initialized to m=1; H_init=2: phases initialized randomly

############ Parameters for the Monte Carlo simulations --> MC_init.txt#####################

Nmisu=200000
ntau=32
nautosave=100000
l_box=1.0
rho_box=0.5
theta_box=3.141592653
A_box=0.1

for L in $LLIST; do

############Creation of the output folder and of the two files of initialization####################

cd ${BASEDIR}/Output_3C

if [ ! -d ./Se_${H_e} ]; then
   mkdir -p e_${H_e}
fi

cd e_${H_e}

if [ ! -d ./Snu_${H_nu} ]; then
   mkdir -p nu_${H_nu}
fi

cd nu_${H_nu}

if [ ! -d ./Sh_${H_h} ]; then
   mkdir -p h_${H_h}
fi

cd h_${H_h}

if [ ! -d ./SL${L}_a${H_a}_b${H_b}_eta${H_eta}_e${H_e}_h${H_h}_nu${H_nu}_bmin${H_blow}_bmax${H_bhigh}_init${H_init} ]; then
   mkdir -p L${L}_a${H_a}_b${H_b}_eta${H_eta}_e${H_e}_h${H_h}_nu${H_nu}_bmin${H_blow}_bmax${H_bhigh}_init${H_init}
fi

OUTPUT=${BASEDIR}/Output_3C/e_${H_e}/nu_${H_nu}/h_${H_h}/L${L}_a${H_a}_b${H_b}_eta${H_eta}_e${H_e}_h${H_h}_nu${H_nu}_bmin${H_blow}_bmax${H_bhigh}_init${H_init}

cd /tmp/Output_x_ilaria

if [ ! -d ./Se_${H_e} ]; then
   mkdir -p e_${H_e}
fi

cd e_${H_e}

if [ ! -d ./Snu_${H_nu} ]; then
   mkdir -p nu_${H_nu}
fi

cd nu_${H_nu}

if [ ! -d ./Sh_${H_h} ]; then
   mkdir -p h_${H_h}
fi

cd h_${H_h}

if [ ! -d ./SL${L}_a${H_a}_b${H_b}_eta${H_eta}_e${H_e}_h${H_h}_nu${H_nu}_bmin${H_blow}_bmax${H_bhigh}_init${H_init} ]; then
   mkdir -p L${L}_a${H_a}_b${H_b}_eta${H_eta}_e${H_e}_h${H_h}_nu${H_nu}_bmin${H_blow}_bmax${H_bhigh}_init${H_init}
fi

OUTPUT_TEMP=/tmp/Output_x_ilaria/e_${H_e}/nu_${H_nu}/h_${H_h}/L${L}_a${H_a}_b${H_b}_eta${H_eta}_e${H_e}_h${H_h}_nu${H_nu}_bmin${H_blow}_bmax${H_bhigh}_init${H_init}

cd ${OUTPUT}

#THE ORDER OF WRITING DOES MATTER
echo $H_a > HP_init.txt
echo $H_b >> HP_init.txt
echo $H_eta >> HP_init.txt
echo $H_e >> HP_init.txt
echo $H_h >> HP_init.txt
echo $H_nu >> HP_init.txt
echo $H_blow >> HP_init.txt
echo $H_bhigh >> HP_init.txt
echo $H_init >> HP_init.txt

#THE ORDER OF WRITING DOES MATTER
echo $Nmisu > MC_init.txt
echo $ntau >> MC_init.txt
echo $nautosave >> MC_init.txt
echo $l_box >> MC_init.txt
echo $rho_box >> MC_init.txt
echo $theta_box >> MC_init.txt
echo $A_box >> MC_init.txt

#################Creation of the submit_runs script#########################

jobname="L${L}_a${H_a}_b${H_b}_eta${H_eta}_e${H_e}_h${H_h}_nu${H_nu}_bmin${H_blow}_bmax${H_bhigh}_init${H_init}"
nnodes=2
ntasks=64 #parallel tempering over ntasks temperatures
time_limit=7-00:00:00

#I create ntasks folder: one for each rank.

cd ${OUTPUT}

for ((rank=0; rank<${ntasks}; rank++)); do

if [ ! -d ./Sbeta_${rank} ]; then
   mkdir -p beta_${rank}
fi

done

cd ${OUTPUT_TEMP}

for ((rank=0; rank<${ntasks}; rank++)); do

if [ ! -d ./Sbeta_${rank} ]; then
   mkdir -p beta_${rank}
fi

done

cd ${SCRIPT_DIR}
DIR_PAR="${OUTPUT}"
DIR_PAR_TEMP="${OUTPUT}"

#SEED= If I want to repeat exactly a simulation I could initialize the random number generator exactly at the same way

EXECUTE_DIR="../build/Release"

#SBATCH --nodes=${nnodes}               # Number of nodes

echo "#!/bin/bash
#SBATCH --job-name=${jobname}          # Name of the job
#SBATCH --time=${time_limit}            # Allocation time
#SBATCH --mem-per-cpu=2000              # Memory per allocated cpu
#SBATCH --nodes=${nnodes}               # Number of nodes
#SBATCH --ntasks=${ntasks}
#SBATCH --output=${DIR_PAR}/logs/log_${jobname}.o
#SBATCH --error=${DIR_PAR}/logs/log_${jobname}.e

srun ${EXECUTE_DIR}/GL_3component ${L} ${DIR_PAR} ${DIR_PAR} ${RESTART} &> ${DIR_PAR}/logs/log_${jobname}.o

" >  submit_run

#Submission of the work --> sbatch submit_runs

mkdir -p ${DIR_PAR}/logs

sbatch submit_run

echo ${L} ${DIR_PAR} ${DIR_PAR} ${RESTART}
done
