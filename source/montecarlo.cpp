#include "montecarlo.h"
#include "main.h"
#include "rng.h"

void metropolis( struct Node* Site, struct NN_Node* NN_Site, struct MC_parameters &MCp, struct H_parameters &Hp,  double my_beta){

    double l, phi, d_theta, d_A, rand;
    unsigned int ix=0, iy=0, iz=0, alpha=0, i, vec=0;
    double acc_rate=0.5, acc_l=0., acc_A=0., acc_theta=0.; // acc_rho=0.,
    struct O2 NewPsi{};
    struct O2 OldPsi{};
    double NewA, OldA;
    double newE, oldE, minus_deltaE;
    double h3=(Hp.h*Hp.h*Hp.h);

    for (ix= 0; ix < Lx; ix++) {
        for (iy = 0; iy < Ly; iy++) {
            for (iz = 0; iz < Lz; iz++) {
                i=ix + Lx * (iy + iz * Ly);
                /*************PSI UPDATE: move in the plane ImPsi, RePsi**********/
                for (alpha = 0; alpha < 3; alpha++) {
                    OldPsi = Site[i].Psi[alpha];
                    oldE = local_HPsi(OldPsi, ix, iy, iz, alpha, Hp, Site, NN_Site);
                    l=rn::uniform_real_box(0, MCp.lbox_l);
                    phi=rn::uniform_real_box(0, C_TWO_PI);
                    NewPsi.x = OldPsi.x + (l * cos(phi));
                    NewPsi.y = OldPsi.y + (l * sin(phi));
                    cartesian_to_polar(NewPsi);
                    newE = local_HPsi(NewPsi, ix, iy, iz, alpha, Hp, Site, NN_Site);
                    minus_deltaE = h3*(oldE - newE);
                    if (minus_deltaE > 0) {
                        Site[i].Psi[alpha] = NewPsi;
                        acc_l++;
                    } else {
                        rand= rn::uniform_real_box(0,1);
                        //Boltzmann weight: exp(-\beta \Delta E) E= h³ \sum_i E(i)
                        if (rand < exp(my_beta * minus_deltaE)) {
                            Site[i].Psi[alpha] = NewPsi;
                            acc_l++;
                            for (vec = 0; vec < 3; vec++) {//Update of the auxiliary struct
                                NN_Site[nn(i, vec, 1)].Psi_minusk[alpha +3*vec].r=NewPsi.r;
                                NN_Site[nn(i, vec, 1)].Psi_minusk[alpha +3*vec].t=NewPsi.t - Hp.e*Hp.h*Site[i].A[vec];
                                polar_to_cartesian(NN_Site[nn(i, vec, 1)].Psi_minusk[alpha +3*vec]);
                                NN_Site[nn(i, vec, -1)].Psi_plusk[alpha +3*vec].r=NewPsi.r;
                                NN_Site[nn(i, vec, -1)].Psi_plusk[alpha +3*vec].t=NewPsi.t + Hp.e*Hp.h*Site[nn(i, vec, -1)].A[vec];
                                polar_to_cartesian(NN_Site[nn(i, vec, -1)].Psi_plusk[alpha +3*vec]);
                            }
                        }
                    }
                }
            }
        }
    }

    for (ix= 0; ix < Lx; ix++) {
        for (iy = 0; iy < Ly; iy++) {
            for (iz = 0; iz < Lz; iz++) {
                i=ix + Lx * (iy + iz * Ly);
                /*******PHASE ONLY UPDATE**************/
                for (alpha = 0; alpha < 3; alpha++) {
                    OldPsi = Site[i].Psi[alpha];
                    oldE = local_Htheta(OldPsi, ix, iy, iz, alpha, Hp, Site, NN_Site);
                    d_theta=rn::uniform_real_box(-MCp.lbox_theta, MCp.lbox_theta);
                    NewPsi.t = fmod(OldPsi.t + d_theta, C_TWO_PI);
                    NewPsi.r= OldPsi.r;
                    polar_to_cartesian(NewPsi);
                    newE = local_Htheta(NewPsi, ix, iy, iz, alpha, Hp, Site, NN_Site);
                    minus_deltaE = h3*(oldE - newE);
                    if (minus_deltaE > 0) {
                        Site[i].Psi[alpha] = NewPsi;
                        acc_theta++;
                    } else {
                        rand= rn::uniform_real_box(0,1);
                        //Boltzmann weight: exp(-\beta E) E= h³ \sum_i E(i)
                        if (rand < exp(my_beta * minus_deltaE)) {
                            Site[i].Psi[alpha] = NewPsi;
                            acc_theta++;
                            for (vec = 0; vec < 3; vec++) {//Update of the auxiliary struct
                                NN_Site[nn(i, vec, 1)].Psi_minusk[alpha +3*vec].r=NewPsi.r;
                                NN_Site[nn(i, vec, 1)].Psi_minusk[alpha +3*vec].t=NewPsi.t - Hp.e*Hp.h*Site[i].A[vec];
                                polar_to_cartesian(NN_Site[nn(i, vec, 1)].Psi_minusk[alpha +3*vec]);
                                NN_Site[nn(i, vec, -1)].Psi_plusk[alpha +3*vec].r=NewPsi.r;
                                NN_Site[nn(i, vec, -1)].Psi_plusk[alpha +3*vec].t=NewPsi.t + Hp.e*Hp.h*Site[nn(i, vec, -1)].A[vec];
                                polar_to_cartesian(NN_Site[nn(i, vec, -1)].Psi_plusk[alpha +3*vec]);
                            }
                        }
                    }
                }
            }
        }
    }

    if(Hp.e!=0) {
        for (ix = 0; ix < Lx; ix++) {
            for (iy = 0; iy < Ly; iy++) {
                for (iz = 0; iz < Lz; iz++) {
                    i = ix + Lx * (iy + iz * Ly);
                    /**********VECTOR POTENTIAL UPDATE********/
                    for (vec = 0; vec < 3; vec++) {
                        //Update of A
                        OldA = Site[i].A[vec];
                        oldE = local_HA(OldA, ix, iy, iz, vec, Hp, Site, NN_Site);
                        d_A = rn::uniform_real_box(-MCp.lbox_A, MCp.lbox_A);
                        NewA = OldA + d_A;
                        newE = local_HA(NewA, ix, iy, iz, vec, Hp, Site, NN_Site);
                        minus_deltaE = h3 * (oldE - newE);
                        if (minus_deltaE > 0.) {
                            Site[i].A[vec] = NewA;
                            acc_A++;
                        } else {
                            rand = rn::uniform_real_box(0, 1);
                            //Boltzmann weight: exp(-\beta E) E= h³ \sum_i E(i)
                            if (rand < exp(my_beta * minus_deltaE)) {
                                Site[i].A[vec] = NewA;
                                acc_A++;
                                for (alpha = 0; alpha < 3; alpha++) {//Update of the auxiliary struct
                                    NN_Site[nn(i, vec, 1)].Psi_minusk[alpha + 3 * vec].r = Site[i].Psi[alpha].r;
                                    NN_Site[nn(i, vec, 1)].Psi_minusk[alpha + 3 * vec].t =
                                            Site[i].Psi[alpha].t - Hp.e * Hp.h * Site[i].A[vec];
                                    polar_to_cartesian(NN_Site[nn(i, vec, 1)].Psi_minusk[alpha + 3 * vec]);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    acc_l=(double) acc_l/(3*N);
    acc_theta=(double) acc_theta/(3*N);
    acc_A=(double) acc_A/(3*N);
    MCp.lbox_l= MCp.lbox_l*(acc_l/acc_rate);
    MCp.lbox_theta= MCp.lbox_theta*(acc_theta/acc_rate);
    MCp.lbox_A= MCp.lbox_A*(acc_A/acc_rate);
}


double local_HPsi(struct O2 Psi, unsigned int ix, unsigned int iy, unsigned int iz, unsigned int alpha,  struct H_parameters &Hp, struct Node* Site, struct NN_Node* NN_Site){

    double h_Potential, h_Kinetic=0., h_Josephson=0., h_AB=0., h_tot;
    double h2=(Hp.h*Hp.h);
    double J_alpha1=0, J_beta1=0, J_alpha2=0, J_beta2=0;
    double gauge_phase1, gauge_phase2;
    unsigned int beta=0, vec=0, i;

    i=ix +Lx*(iy+Ly*iz);
    //Compute the local Energy respect to a given component (alpha) of the matter field Psi and a given spatial position (r=(ix, iy, iz))
    //We need to compute just the part of the Hamiltonian involving Psi

    //Potential= (a+ 3/h²)*|Psi_{alpha}(r)|² + b/2*|Psi_{alpha}(r)|⁴
    h_Potential= O2norm2(Psi)*((Hp.a + (3./h2))+  (0.5*Hp.b*O2norm2(Psi)));

    //Kinetic= -(1/h²)*\sum_k=1,2,3 (|Psi_{alpha}(r)||Psi_{alpha}(r+k)|* cos(theta_{alpha}(r+k) - theta_{alpha}(r) +h*e*A_k(r))) + (|Psi_{alpha}(r-k)||Psi_{alpha}(r)|* cos(theta_{alpha}(r) - theta_{alpha}(r-k) +h*e*A_k(r-k)))
    for(vec=0; vec<3; vec++){

        h_Kinetic-=(1./h2)*(O2prod(Psi, NN_Site[i].Psi_plusk[alpha+3*vec]) + O2prod(NN_Site[i].Psi_minusk[alpha+3*vec], Psi));

        //Andreev-Bashkin term = \sum_beta!=alpha \sum_k=1,2,3 nu*(J^k_alpha - J^k_beta)^2;
        // with J^k_alpha= |Psi_{alpha}(r)||Psi_{alpha}(r+k)|* sin(theta_{alpha}(r+k) - theta_{alpha}(r) +h*e*A_k(r))) 
        if(Hp.nu !=0 ) {
            //+k
            J_alpha1=(1./Hp.h)*O2vprod(Psi, NN_Site[i].Psi_plusk[alpha+3*vec]); //this order corresponds to gauge_phase1
            //-k
            J_alpha2=(1./Hp.h)*O2vprod(NN_Site[i].Psi_minusk[alpha+3*vec], Psi); //this order corresponds to gauge_phase2

            for (beta = 0; beta < 3; beta++) {
            	if (beta != alpha) {
                	//+k
                    J_beta1=O2vprod(Site[i].Psi[beta], NN_Site[i].Psi_plusk[beta+3*vec]);
                	h_AB += Hp.nu * (pow((J_alpha1 -J_beta1),2));
                	//-k
                    J_beta2=O2vprod(NN_Site[i].Psi_minusk[beta+3*vec], Site[i].Psi[beta]);
                	h_AB += Hp.nu * (pow((J_alpha2 -J_beta2),2));
			}
		}
	}
    }

    //Josephson= eta* \sum_beta!=alpha |Psi_{alpha}(r)||Psi_{beta}(r)|* cos(theta_{alpha}(r) - theta_{beta}(r))
    for(beta=0; beta<3; beta++){
        if(beta != alpha) {
            h_Josephson += (Hp.eta * O2prod(Psi, Site[i].Psi[beta]));
        }
    }

    h_tot= h_Potential + h_Kinetic + h_Josephson + h_AB;
    return h_tot;
}

double local_Htheta(struct O2 Psi, unsigned int ix, unsigned int iy, unsigned int iz, unsigned int alpha,  struct H_parameters &Hp, struct Node* Site, struct NN_Node* NN_Site){

    double h_Kinetic=0., h_Josephson=0., h_AB=0., h_tot;
    double h2=(Hp.h*Hp.h);
    double J_alpha1=0, J_beta1=0, J_alpha2=0, J_beta2=0;
    double gauge_phase1=0, gauge_phase2=0;
    unsigned int beta=0, vec=0, i;
    struct O2 Psi_aux_plus;
    struct O2 Psi_aux_minus;
    struct O2 Psib_aux_plus;
    struct O2 Psib_aux_minus;

    i=ix +Lx*(iy+Ly*iz);
    //We need to compute just the part of the Hamiltonian involving Psi.t

    //Kinetic= -(1/h²)*\sum_k=1,2,3 (|Psi_{alpha}(r)||Psi_{alpha}(r+k)|* cos(theta_{alpha}(r+k) - theta_{alpha}(r) +h*e*A_k(r))) + (|Psi_{alpha}(r-k)||Psi_{alpha}(r)|* cos(theta_{alpha}(r) - theta_{alpha}(r-k) +h*e*A_k(r-k)))
    for(vec=0; vec<3; vec++){

        h_Kinetic-=(1./h2)*(O2prod(Psi, NN_Site[i].Psi_plusk[alpha+3*vec]) + O2prod(NN_Site[i].Psi_minusk[alpha+3*vec], Psi));

        //Andreev-Bashkin term = \sum_beta!=alpha \sum_k=1,2,3 nu*(J^k_alpha - J^k_beta)^2;
        // with J^k_alpha= |Psi_{alpha}(r)||Psi_{alpha}(r+k)|* sin(theta_{alpha}(r+k) - theta_{alpha}(r) +h*e*A_k(r)))
        if(Hp.nu !=0 ) {
            //+k
            J_alpha1=(1./Hp.h)*O2vprod(Psi, Psi_aux_plus); //this order corresponds to gauge_phase1

            //-k
            J_alpha2=(1./Hp.h)*O2vprod(Psi_aux_minus, Psi); //this order corresponds to gauge_phase2

            for (beta = 0; beta < 3; beta++) {
                if (beta != alpha) {
		            //+k
                    J_beta1=O2vprod(Site[i].Psi[beta], NN_Site[i].Psi_plusk[beta+3*vec]);
                    h_AB += Hp.nu * (pow((J_alpha1 -J_beta1),2));
                    //-k
                    J_beta2=O2vprod(NN_Site[i].Psi_minusk[beta+3*vec], Site[i].Psi[beta]);
                    h_AB += Hp.nu * (pow((J_alpha2 -J_beta2),2));
                }
            }
        }
    }

    //Josephson= eta* \sum_beta!=alpha |Psi_{alpha}(r)||Psi_{beta}(r)|* cos(theta_{alpha}(r) - theta_{beta}(r))
    for(beta=0; beta<3; beta++){
        if(beta != alpha) {
            h_Josephson += (Hp.eta * O2prod(Psi, Site[i].Psi[beta]));
        }
    }

    h_tot= h_Kinetic + h_Josephson + h_AB;
    return h_tot;
}

double local_HA(double A, unsigned int ix, unsigned int iy, unsigned int iz,  unsigned int vec,  struct H_parameters &Hp, struct Node* Site, struct NN_Node* NN_Site){

    double h_Kinetic=0., h_B, h_AB=0., h_tot=0.;
    double h2=(Hp.h*Hp.h);
    double J_alpha, J_beta, J_gamma;
    unsigned int alpha=0, i;

    i=ix +Lx*(iy+Ly*iz);
    //Compute the local Energy respect to a given component (alpha) of the vector potential A and a given spatial position (r=(ix, iy, iz))
    //We need to compute just the part of the Hamiltonian involving A

    //Kinetic= -(1/h²)*\sum_k=1,2,3 |Psi_{alpha}(r)||Psi_{alpha}(r+k)|* cos(theta_{alpha}(r+k) - theta_{alpha}(r) +h*e*A_k(r))
    for(alpha=0; alpha<3; alpha++) {
        h_Kinetic -= (1. / h2) *( O2prod(Site[i].Psi[alpha], NN_Site[i].Psi_plusk[alpha+3*vec]));
    }
    //Andreev-Bashkin term = \sum_beta!=alpha \sum_k=1,2,3 nu*(J^k_alpha - J^k_beta)^2;
    // with J^k_alpha= |Psi_{alpha}(r)||Psi_{alpha}(r+k)|* sin(theta_{alpha}(r+k) - theta_{alpha}(r) +h*e*A_k(r)))
    if(Hp.nu !=0 ) {
        /*In this case is more convenient to use the sine here*/
        J_alpha = (1./Hp.h) * O2vprod(Site[i].Psi[0], NN_Site[i].Psi_plusk[0+3*vec]);
        J_beta = (1./Hp.h) * O2vprod(Site[i].Psi[1], NN_Site[i].Psi_plusk[1+3*vec]);
        J_gamma = (1./Hp.h) * O2vprod(Site[i].Psi[2], NN_Site[i].Psi_plusk[2+3*vec]);

        h_AB += Hp.nu * (pow((J_alpha - J_beta), 2)  + pow((J_alpha - J_gamma), 2) + pow((J_gamma - J_beta), 2));
    }

    h_B=(0.5/h2)*F_2(A, vec, ix, iy, iz, Site);

    h_tot= h_Kinetic + h_AB+ h_B;
    return h_tot;
}

double F_2(double newA, unsigned int k, unsigned int ix, unsigned int iy, unsigned int iz, struct Node* Site){

    unsigned int l;
    unsigned int i;
    double F2_A=0., F_A;
    i=ix +Lx*(iy+Ly*iz);
    //All the plaquettes involving A_vec(i)
    for(l=0; l<3; l++){
        if(l!= k){
            F_A=(newA + Site[nn(i, k, 1)].A[l] - Site[nn(i, l, 1)].A[k] - Site[i].A[l]);
            F2_A+=(F_A*F_A);
        }
    }
    for(l=0; l<3; l++){
        if(l!= k){
            F_A=(Site[nn(i, l, -1)].A[k]+ Site[nn(nn(i, l, -1), k, 1)].A[l] - newA - Site[nn(i, l, -1)].A[l]);
            F2_A+=(F_A*F_A);
        }
    }

    return F2_A;
}

