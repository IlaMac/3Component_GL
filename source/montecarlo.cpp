#include "montecarlo.h"
#include "main.h"

void metropolis( struct Node* Site, struct MC_parameters &MCp, struct H_parameters &Hp){

    std::uniform_real_distribution<double> rand_pos(0, 1.);
    std::uniform_real_distribution<double> rand_sym(-1., 1.);
    double l=0., phi=0., d_theta=0., d_rho=0., d_A=0.,rand =0.;
    unsigned int ix=0, iy=0, iz=0, alpha=0;
    double acc_rate=0.5, acc_l=0., acc_A=0.; //acc_theta=0., acc_rho=0.,
    struct O2 NewPsi;
    struct O2 OldPsi;
    double NewA=0.;
    double newE=0., oldE=0., minus_deltaE=0.;

    for (ix= 0; ix < Lx; ix++) {
        for (iy= 0; iy < Ly; iy++) {
            for (iz = 0; iz < Lz; iz++) {
                for (alpha = 0; alpha < 3; alpha++) {
                    OldPsi=Site[ix+ Lx*(iy+ iz*Ly)].Psi[alpha];
                    oldE = local_HPsi(OldPsi, ix, iy, iz, alpha, Hp, Site);
                    //Update of Psi: move in the plane ImPsi, RePsi
                    l = MCp.lbox_l * rand_pos(MCp.rng);
                    phi = C_TWO_PI *rand_pos(MCp.rng);
                    NewPsi.x = Site[ix+ Lx*(iy+ iz*Ly)].Psi[alpha].x + (l*cos(phi));
                    NewPsi.y = Site[ix+ Lx*(iy+ iz*Ly)].Psi[alpha].x + (l*sin(phi));
                    cartesian_to_polar(NewPsi);
                    newE = local_HPsi(NewPsi, ix, iy, iz, alpha, Hp, Site);
                    minus_deltaE= (oldE- newE);
                    if(minus_deltaE>0){
                        Site[ix+ Lx*(iy+ iz*Ly)].Psi[alpha]=NewPsi;
                        acc_l++;
                    }
                    else{
                        rand=rand_pos(MCp.rng);
                        if( rand < exp(Hp.beta*minus_deltaE) ) {
                            Site[ix + Lx * (iy + iz * Ly)].Psi[alpha] = NewPsi;
                            acc_l++;
                        }else{
                            Site[ix + Lx * (iy + iz * Ly)].Psi[alpha] = OldPsi;
                        }
                    }

                    //Update of A
                    rand=rand_sym(MCp.rng);
                    d_A=((MCp.lbox_A)*rand);
                    NewA= Site[ix + Lx * (iy + iz * Ly)].A[alpha] + d_A;
                    newE = local_HA(NewA, ix, iy, iz, alpha, Hp, Site);
                    oldE = local_HA(Site[ix+ Lx*(iy+ iz*Ly)].A[alpha], ix, iy, iz, alpha, Hp, Site);
                    minus_deltaE= (oldE -newE);
                    if(minus_deltaE>0.){
                        Site[ix+ Lx*(iy+ iz*Ly)].A[alpha]=NewA;
                        acc_A++;
                    }else{
                        rand=rand_pos(MCp.rng);
                        if( rand < exp(Hp.beta*minus_deltaE) ) {
                            Site[ix + Lx * (iy + iz * Ly)].A[alpha] = NewA;
                            acc_A++;
                        }
                    }
                }
            }
        }
    }

    acc_l= (double) acc_l/(3*N);
    acc_A= (double) acc_A/(3*N);
    MCp.lbox_l= MCp.lbox_l*(acc_l/acc_rate);
    MCp.lbox_A= MCp.lbox_A*(acc_A/acc_rate);
}

double local_HPsi(struct O2 Psi, unsigned int ix, unsigned int iy, unsigned int iz, unsigned int alpha,  struct H_parameters &Hp, struct Node* Site){

    double h_Potential=0., h_Kinetic=0., h_Josephson=0., h_tot=0.;
    double h2=(Hp.h*Hp.h);
    unsigned int beta=0, vec=0;
    struct O2 S_gauge;

    //Compute the local Energy respect to a given component (alpha) of the matter field Psi and a given spatial position (r=(ix, iy, iz))
    //We need to compute just the part of the Hamiltonian involving Psi

    //Potential= (a+ 3/h²)*|Psi_{alpha}(r)|² + b/2*|Psi_{alpha}(r)|⁴
    h_Potential+= (O2norm2(Psi)*((Hp.a + (3./h2))+  0.5*Hp.b*O2norm2(Psi)));

    //Kinetic= -(1/h²)*\sum_k=1,2,3 |Psi_{alpha}(r)||Psi_{alpha}(r+k)|* cos(theta_{alpha}(r+k) - theta_{alpha}(r) +h*e*A_k(r))
    for(vec=0; vec<3; vec++){
        //I can rewrite the kinetic energy in terms of the scalar product between the vector in (r) and the vector in (r+k) rotated according to A
        S_gauge.t= (Site[nn(ix, iy, iz, vec, 1)].Psi[alpha].t + Hp.h*Hp.e*Site[ix+ Lx*(iy+ iz*Ly)].A[vec]);
        S_gauge.r= Site[nn(ix, iy, iz, vec, 1)].Psi[alpha].r;
        polar_to_cartesian(S_gauge);
        h_Kinetic-=((1./h2)*O2prod(Psi, S_gauge));
    }

    //Josephson= eta* \sum_beta!=alpha |Psi_{alpha}(r)||Psi_{beys}(r)|* cos(theta_{alpha}(r) - theta_{beta}(r))
    for(beta=0; beta<3; beta++){
        if(beta != alpha) {
            h_Josephson += (Hp.eta * O2prod(Psi, Site[ix + Lx * (iy + iz * Ly)].Psi[beta]));
        }
    }

    h_tot= h_Potential + h_Kinetic + h_Josephson;
    return h_tot;
}

double local_HA(double A, unsigned int ix, unsigned int iy, unsigned int iz, unsigned int vec,  struct H_parameters &Hp, struct Node* Site){

    double h_Kinetic=0., h_B=0., h_tot=0.;
    double h2=(Hp.h*Hp.h);
    unsigned int alpha=0;
    struct O2 S_gauge;

    //Compute the local Energy respect to a given component (alpha) of the vector potential A and a given spatial position (r=(ix, iy, iz))
    //We need to compute just the part of the Hamiltonian involving A

    //Kinetic= -(1/h²)*\sum_k=1,2,3 |Psi_{alpha}(r)||Psi_{alpha}(r+k)|* cos(theta_{alpha}(r+k) - theta_{alpha}(r) +h*e*A_k(r))
    for(alpha=0; alpha<3; alpha++){
        //I can rewrite the kinetic energy in terms of the scalar product between the vector in (r) and the vector in (r+k) rotated according to A
        S_gauge.t= (Site[nn(ix, iy, iz, vec, 1)].Psi[alpha].t + Hp.h*Hp.e*A);
        S_gauge.r= Site[nn(ix, iy, iz, vec, 1)].Psi[alpha].r;
        polar_to_cartesian(S_gauge);
        h_Kinetic-=((1./h2)*O2prod(Site[ix + Lx * (iy + iz * Ly)].Psi[alpha], S_gauge));
    }

    h_B +=(0.5/h2)*F_2(A, vec, ix, iy, iz, Site);
    h_tot= h_Kinetic + h_B ;

    return h_tot;
}

double F_2(double newA, unsigned int vec, unsigned int ix, unsigned int iy, unsigned int iz, struct Node* Site){

    unsigned int m;
    double F2_A=0., F_A=0.;
    for(m=0; m<3; m++){
        if(m != vec){
            F_A=(Site[ix+ Lx*(iy+ iz*Ly)].A[m] + Site[nn(ix, iy, iz, m, 1)].A[vec] - Site[nn(ix, iy, iz, vec, 1)].A[m] - newA);
            F2_A+=(F_A*F_A);
        }
    }
    return F2_A;
}

