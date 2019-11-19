//
// Created by ilaria on 2019-11-16.
//

#ifndef MEASURES_H
#define MEASURES_H

#include "main.h"
#include <fstream>

struct Measures{

    double m; //magnetization (for the phase chirality of the three components
    //Binder cumulant U=<m⁴>/(3*<m²>²)

    double E; //Energy
    double E_pot;
    double E_kin;
    double E_Josephson;
    double E_B;

    double d_rhoz; //Dual stiffness along z

    double *density_psi;
};


void measures_init(struct Measures &mis);
void measures_reset(struct Measures &mis);
void dual_stiffness(struct Measures &mis, struct H_parameters &Hp, struct Node* Site);
void magnetization(struct Measures &mis, struct Node* Site);
void energy(struct Measures &mis, struct H_parameters &Hp, struct Node* Site);
void density_psi(struct Measures &mis, struct Node* Site);


#endif //MEASURES_H
