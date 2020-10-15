//
// Created by ilaria on 2019-11-13.
//
#include "initialization.h"
#include "robust_filesystem.h"

void initialize_Hparameters(struct H_parameters &Hp, const fs::path & directory_parameters){

    fs::path hp_init_file = directory_parameters / "HP_init.txt";
    Hp.a=0;
    Hp.b=1;
    Hp.eta=1;
    Hp.e=0.5;
    Hp.h= 5.4;
    Hp.nu=0;
    Hp.b_low=0.2245;
    Hp.b_high=0.2265;
    Hp.init=0;

    if(fs::exists(hp_init_file)){
        FILE *fin= nullptr;
        if((fin=fopen(hp_init_file.c_str(), "r"))) {
            fscanf(fin, "%d" , &Hp.a);
            fscanf(fin, "%d" , &Hp.b);
            fscanf(fin, "%d" , &Hp.eta);
            fscanf(fin, "%lf" , &Hp.e);
            fscanf(fin, "%lf" , &Hp.h);
            fscanf(fin, "%lf" , &Hp.nu);
            fscanf(fin, "%lf" , &Hp.b_low);
            fscanf(fin, "%lf" , &Hp.b_high);
            fscanf(fin, "%d", &Hp.init);
            fclose(fin);
            //With this modification Hp.beta in not anymore part of the Hamiltonian parameters list
        }
    }else{
        Hp.a=0;
        Hp.b=1;
        Hp.eta=1;
        Hp.e=0.5;
        Hp.h= 5.0;
        Hp.nu=0.1;
        Hp.b_low=0.2245;
        Hp.b_high=0.2265;
        Hp.init=1;
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
        MCp.tau=1;
        MCp.n_autosave=20000; //not used now
        MCp.lbox_l=1.0;
        MCp.lbox_rho=0.5;
        MCp.lbox_theta=C_PI*0.25;
        MCp.lbox_A=0.1;
    }

}

void initialize_lattice(struct Node* Site, const fs::path & directory_read, int RESTART, struct H_parameters &Hp){

    unsigned int i, alpha;
    fs::path psi_init_file = directory_read / "Psi_restart.bin";
    fs::path a_init_file = directory_read / "A_restart.bin";

    if(RESTART==1){
        fs::path psi_init_file = directory_read / "Psi_restart.bin";
        fs::path a_init_file = directory_read / "A_restart.bin";
    }
    else if (RESTART==2){
        fs::path psi_init_file = directory_read / "Psi_final.bin";
        fs::path a_init_file = directory_read / "A_final.bin";
    }

    if((fs::exists(psi_init_file)) and (fs::exists(a_init_file)) and (RESTART!=0)){

        FILE *fPsi= nullptr;
        FILE *fA= nullptr;
        if((fPsi=fopen(psi_init_file.c_str(), "r")) and (fA=fopen(a_init_file.c_str(), "r")) ) {
            for (i = 0; i < N; i++) {
                fread(Site[i].Psi, sizeof(struct O2), NC, fPsi);
                fread(Site[i].A, sizeof(double), 3, fA);
            }
            fclose(fA);
            fclose(fPsi);
        }
    }else{
        if(Hp.init==0) {
            for (i = 0; i < N; i++) {
                for (alpha = 0; alpha < 3; alpha++) {
                    Site[i].Psi[0].t = 0.;
                    Site[i].Psi[alpha].r =sqrt(1. / 2.);
                    polar_to_cartesian(Site[i].Psi[alpha]);
                }
            }
        }
        else if(Hp.init==1) {
            for (i = 0; i < N; i++) {
                Site[i].Psi[0].t = 0.;
                Site[i].Psi[1].t = -2.*C_PI/3;
                Site[i].Psi[2].t = -4.*C_PI/3;
                for (alpha = 0; alpha < 3; alpha++) {
                    Site[i].Psi[alpha].r =sqrt(1. / 2.);
                    polar_to_cartesian(Site[i].Psi[alpha]);
                }
            }
        }
        else if((Hp.init!=0) && (Hp.init!=1)) {
            for (i = 0; i < N; i++) {
                for (alpha = 0; alpha < 3; alpha++) {
                    Site[i].Psi[alpha].r = sqrt(rn::uniform_real_box(0, 1));
                    Site[i].Psi[alpha].t = rn::uniform_real_box(0, C_TWO_PI);
                    polar_to_cartesian(Site[i].Psi[alpha]);
                }
            }
        }

    }

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
