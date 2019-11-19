//
// Created by ilaria on 2019-11-13.
//

#ifndef INITIALIZATION_H
#define INITIALIZATION_H

#include<cstdio>
#include<cmath>
#include<cstdlib>
#include <iostream>
#include <fstream>
#include "main.h"
#include <filesystem>
#define C_TWO_PI (6.2831853071795864769252867665590058L)
#define C_PI (3.1415926535897932384626433832795029L)

namespace fs = std::filesystem;

struct Node{
    double* A; /*three spatial dimensions*/
    struct O2* Psi; /*three SC components*/
};

struct H_parameters{
    /*These values are externally given by an input file*/
    int a;
    int b;
    int eta;
    double e;
    double h;
    double beta;
};

struct MC_parameters{
    /*These values are externally given by an input file*/
    std::mt19937 rng; //Standard mersenne_twister_engine seeded with
    int rnd_seed; //seed to initialize the random number generator
    int tau; // estimate of the auto-correlation time
    int nmisu; //total number of independent measures
    int n_autosave; //frequency at which intermediate configuration are saved
    double lbox_l; //length of the box for the uniform distribution of l (cartesian trasformafion of Psi)
    double lbox_rho; //length of the box for the uniform distribution of rho (polar trasformafion of Psi --> modulus)
    double lbox_theta; //length of the box for the uniform distribution of theta (polar trasformafion of Psi --> phase)
    double lbox_A; //length of the box for the uniform distribution of dA (trasformafion of the vector potential)

};

void initialize_lattice(struct Node* Site, const fs::path & directory_read);
void initialize_Hparameters(struct H_parameters &Hp, const fs::path & directory_parameters);
void initialize_MCparameters(struct MC_parameters &MCp, const fs::path & directory_parameters);

#endif //INITIALIZATION_H
