//
// Created by ilaria on 2019-11-13.
//
#include "initialization.h"
#include <filesystem>

void initialize_Hparameters(struct H_parameters &Hp, const fs::path & directory_parameters){

    fs::path hp_init_file = directory_parameters / "HP_init.txt";
    if(fs::exists(hp_init_file)){
        FILE *fin= nullptr;
        if((fin=fopen(hp_init_file.c_str(), "r"))) {
            fscanf(fin, "%d" , &Hp.a);
            fscanf(fin, "%d" , &Hp.b);
            fscanf(fin, "%d" , &Hp.eta);
            fscanf(fin, "%lf" , &Hp.e);
            fscanf(fin, "%lf" , &Hp.h);
            fscanf(fin, "%lf" , &Hp.b_low);
            fscanf(fin, "%lf" , &Hp.b_high);
            fclose(fin);
            //With this modification Hp.beta in not anymore part of the Hamiltonian parameters list
        }
    }else{
        Hp.a=0;
        Hp.b=1;
        Hp.eta=1;
        Hp.e=0.5;
        Hp.h= 5.0;
        Hp.b_low=0.244;
        Hp.b_high=0.247;
    }

}

void initialize_MCparameters(struct MC_parameters &MCp, const fs::path & directory_parameters){

    fs::path mc_init_file = directory_parameters / "MC_init.txt";
    if(fs::exists(mc_init_file)){
        FILE *fin= nullptr;
        if((fin=fopen(mc_init_file.c_str(), "r"))) {
            fscanf(fin, "%d", &MCp.nmisu);
            fscanf(fin, "%d", &MCp.tau);
            fscanf(fin, "%d", &MCp.n_autosave);
            fscanf(fin, "%lf", &MCp.lbox_l);
            fscanf(fin, "%lf", &MCp.lbox_rho);
	        fscanf(fin, "%lf", &MCp.lbox_theta);
            fscanf(fin, "%lf", &MCp.lbox_A);
            fclose(fin);
        }
    }else{
        MCp.nmisu=1000;
        MCp.tau=10;
        MCp.n_autosave=20000;
        MCp.lbox_l=1.0;
        MCp.lbox_rho=0.5;
        MCp.lbox_theta=C_PI;
        MCp.lbox_A=0.1;
    }

}

void initialize_lattice(struct Node* Site, const fs::path & directory_read){

    fs::path psi_init_file = directory_read / "Psi_final.txt";
    fs::path a_init_file = directory_read / "A_final.txt";
    unsigned int i=0;
    unsigned int alpha=0;

    if(fs::exists(psi_init_file)){
        FILE *fPsi= nullptr;
        if((fPsi=fopen(psi_init_file.c_str(), "r"))) {
            for (i = 0; i < N; i++) {
                fread(Site[i].Psi, sizeof(struct O2), 3, fPsi);
            }
            fclose(fPsi);
        }
    }

    if(fs::exists(a_init_file)){
        FILE *fA= nullptr;
        if((fA=fopen(a_init_file.c_str(), "r"))) {
            for (i = 0; i < N; i++) {
                fread(Site[i].A, sizeof(struct O2), 3, fA);
            }
            fclose(fA);
        }
    }

 /*
    for(i=0; i<N; i++){
        for(alpha=0; alpha<3; alpha++){
            Site[i].Psi[alpha].r=sqrt(1./3.);
            Site[i].Psi[alpha].t=0.;
            polar_to_cartesian(Site[i].Psi[alpha]);
        }
    }
*/

}

void initialize_PTarrays(struct PT_parameters &PTp, struct PTroot_parameters &PTroot, struct H_parameters &Hp){
    int p;
    double  beta_low, beta_high, delta_beta;

    if(Hp.b_high>Hp.b_low){ //Paranoic check
        beta_low=Hp.b_low;
        beta_high=Hp.b_high;
    }else{
        beta_low=Hp.b_high;
        beta_high=Hp.b_low;
    }
    PTroot.beta.resize(PTp.np, 0.0);
    PTroot.All_Energies.resize(PTp.np, 0.0);
    PTroot.ind_to_rank.resize(PTp.np, 0);
    PTroot.rank_to_ind.resize(PTp.np, 0);
    delta_beta=(beta_high-beta_low)/(PTp.np-1);
    for(p=0; p<PTp.np; p++){
        PTroot.rank_to_ind[p]=p;
        PTroot.ind_to_rank[p]=p;
        PTroot.beta[p]=beta_low + p*delta_beta;
    }
};