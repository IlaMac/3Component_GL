#ifndef MONTECARLO_H
#define MONTECARLO_H
#include <iostream>
#include <cstring>
#include "o2.h"

void metropolis(struct Node* Site, struct MC_parameters &MCp, struct H_parameters &Hp, double my_beta);
double local_HPsi(struct O2 Psi, unsigned int ix, unsigned int iy, unsigned int iz, unsigned int alpha,  struct H_parameters &Hp, struct Node* Site);
double local_Htheta(struct O2 Psi, unsigned int ix, unsigned int iy, unsigned int iz, unsigned int alpha,  struct H_parameters &Hp, struct Node* Site);
double local_HA(double A, unsigned int ix, unsigned int iy, unsigned int iz, unsigned int alpha,  struct H_parameters &Hp, struct Node* Site);
double F_2(double newA, unsigned int vec, unsigned int ix, unsigned int iy, unsigned int iz, struct Node* Site);
#endif