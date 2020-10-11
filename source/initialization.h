//
// Created by ilaria on 2019-11-13.
//

#ifndef INITIALIZATION_H
#define INITIALIZATION_H

#include "main.h"
#include "robust_filesystem.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#define C_TWO_PI (6.2831853071795864769252867665590058L)
#define C_PI (3.1415926535897932384626433832795029L)
#define Annealing (0) //To be implemented

//static constexpr double C_TWO_PI = 6.2831853071795864769252867665590058;


struct Node{
    /*three spatial dimensions*/
    //std::array<double,3> A;
    double* A;
    /*three SC components*/
    //std::array<O2, NC> Psi; 
    struct O2* Psi;
};

struct H_parameters{
    /*These values are externally given by an input file*/
    int a;
    int b;
    int eta;
    double e;
    double h;
    double nu;
    double b_low; //lowest beta for the parallel tempering procedure
    double b_high; // highest beta for the parallel tempering procedure
    int init; //Initial conditions
};

struct PT_parameters{
    /*Parallel Tempering parameters*/
    int np;
    int rank;
    int root=0;
};

struct PTroot_parameters{
    /*Arrays root needs to handle the swaps*/
    std::vector <double> beta;
    std::vector <double> All_Energies;
    std::vector <int> ind_to_rank;
    std::vector <int> rank_to_ind;
};

struct MC_parameters{
    /*These values are externally given by an input file*/
    int tau; // estimate of the auto-correlation time
    int nmisu; //total number of independent measures
    int n_autosave; //frequency at which intermediate configuration are saved
    double lbox_l; //length of the box for the uniform distribution of l (cartesian transformation of Psi)
    double lbox_rho; //length of the box for the uniform distribution of rho (polar transformation of Psi --> modulus)
    double lbox_theta; //length of the box for the uniform distribution of theta (polar transformation of Psi --> phase)
    double lbox_A; //length of the box for the uniform distribution of dA (transformation of the vector potential)
};

void initialize_lattice(struct Node* Site, const fs::path & directory_read, int RESTART, struct H_parameters &Hp);
void initialize_Hparameters(struct H_parameters &Hp, const fs::path & directory_parameters);
void initialize_MCparameters(struct MC_parameters &MCp, const fs::path & directory_parameters);
void initialize_PTarrays(struct PT_parameters &PTp, struct PTroot_parameters &PTroot, struct H_parameters &Hp);

#endif //INITIALIZATION_H
