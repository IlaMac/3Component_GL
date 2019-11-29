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
	    fscanf(fin, "%lf" , &Hp.beta);
            fclose(fin);
            printf("Initialization:\n a %d \n b %d \n eta %d \n e %lf \n h %lf \n beta %lf", Hp.a, Hp.b, Hp.eta, Hp.e, Hp.h, Hp.beta);
        }
    }else{
        Hp.a=0;
        Hp.b=1;
        Hp.eta=1;
        Hp.e=0.5;
        Hp.h= 5.0;
        Hp.beta=0.244;
    }

}

void initialize_MCparameters(struct MC_parameters &MCp, const fs::path & directory_parameters){

    fs::path mc_init_file = directory_parameters / "MC_init.txt";
    if(fs::exists(mc_init_file)){
        FILE *fin= nullptr;
        if((fin=fopen(mc_init_file.c_str(), "r"))) {
            fscanf(fin, "%d", &MCp.rnd_seed);
	    fscanf(fin, "%d", &MCp.tau);
            fscanf(fin, "%d", &MCp.nmisu);
            fscanf(fin, "%d", &MCp.n_autosave);
            fscanf(fin, "%lf", &MCp.lbox_l);
            fscanf(fin, "%lf", &MCp.lbox_rho);
	    fscanf(fin, "%lf", &MCp.lbox_theta);
            fscanf(fin, "%lf", &MCp.lbox_A);
            fclose(fin);
            printf("Initialization:\n %d \n %d \n %d \n %d \n %lf \n %lf \n %lf \n %lf", MCp.rnd_seed, MCp.tau, MCp.nmisu, MCp.n_autosave, MCp.lbox_l, MCp.lbox_rho, MCp.lbox_theta, MCp.lbox_A);
        }
    }else{
        MCp.rnd_seed=2;
        MCp.tau=10;
        MCp.nmisu=100000;
        MCp.n_autosave=20000;
        MCp.lbox_l=1.0;
        MCp.lbox_rho=0.5;
        MCp.lbox_theta=C_PI;
        MCp.lbox_A=0.1;
    }

    MCp.rng.seed(MCp.rnd_seed);

}

void initialize_lattice(struct Node* Site, const fs::path & directory_read){

    fs::path psi_init_file = directory_read / "Psi_final.txt";
    fs::path a_init_file = directory_read / "A_final.txt";
    unsigned int i=0;

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

}

