#include "montecarlo.h"
#include "main.h"
#include "rng.h"
#include "class_tic_toc.h"

void metropolis( struct Node* Site, struct MC_parameters &MCp, struct H_parameters &Hp,  double my_beta){

    double l, phi, d_theta, d_A, rand;
    unsigned int ix=0, iy=0, iz=0, alpha=0, vec=0, i;
    double acc_rate=0.5, acc_l=0., acc_A=0., acc_theta=0.; // acc_rho=0.,
    struct O2 NewPsi{};
    struct O2 OldPsi{};
    double NewA, OldA;
    double newE, oldE, minus_deltaE;
    double h3=(Hp.h*Hp.h*Hp.h);
    class_tic_toc t_localHpsi(true,5,"local_Hpsi");
    class_tic_toc t_localHtheta(true,5,"local_Htheta");
    class_tic_toc t_localHA(true,5,"local_HA");

    //t_localHpsi.tic();
    for (iz= 0; iz < Lz; iz++) {
        for (iy = 0; iy < Ly; iy++) {
            for (ix = 0; ix < Lx; ix++) {
                i = ix + Lx * (iy + iz * Ly);
                /*************PSI UPDATE: move in the plane ImPsi, RePsi**********/
                for (alpha = 0; alpha < 3; alpha++) {
                    OldPsi = Site[i].Psi[alpha];
                    oldE = local_HPsi(OldPsi, ix, iy, iz, alpha, Hp, Site);
                    l = rn::uniform_real_box(0, MCp.lbox_l);
                    phi = rn::uniform_real_box(0, C_TWO_PI);
                    NewPsi.x = OldPsi.x + (l * cos(phi));
                    NewPsi.y = OldPsi.y + (l * sin(phi));
                    cartesian_to_polar(NewPsi);
                    newE = local_HPsi(NewPsi, ix, iy, iz, alpha, Hp, Site);
                    minus_deltaE = h3 * (oldE - newE);
                    if (minus_deltaE > 0) {
                        Site[i].Psi[alpha] = NewPsi;
                        acc_l++;
                    } else {
                        rand = rn::uniform_real_box(0, 1);
                        //Boltzmann weight: exp(-\beta \Delta E) E= h³ \sum_i E(i)
                        if (rand < exp(my_beta * minus_deltaE)) {
                            Site[i].Psi[alpha] = NewPsi;
                            acc_l++;
                        }
                    }
                }
            }
        }
    }
    //t_localHpsi.toc();

    //t_localHtheta.tic();
    for (iz= 0; iz < Lz; iz++) {
        for (iy = 0; iy < Ly; iy++) {
            for (ix = 0; ix < Lx; ix++) {
                i = ix + Lx * (iy + iz * Ly);

                /*******PHASE ONLY UPDATE**************/
                for (alpha = 0; alpha < 3; alpha++) {
                    OldPsi = Site[i].Psi[alpha];
                    oldE = local_Htheta(OldPsi, ix, iy, iz, alpha, Hp, Site);
                    d_theta = rn::uniform_real_box(-MCp.lbox_theta, MCp.lbox_theta);
                    NewPsi.t = fmod(OldPsi.t + d_theta, C_TWO_PI);
                    NewPsi.r = OldPsi.r;
                    polar_to_cartesian(NewPsi);
                    newE = local_Htheta(NewPsi, ix, iy, iz, alpha, Hp, Site);
                    minus_deltaE = h3 * (oldE - newE);
                    if (minus_deltaE > 0) {
                        Site[i].Psi[alpha] = NewPsi;
                        acc_theta++;
                    } else {
                        rand = rn::uniform_real_box(0, 1);
                        //Boltzmann weight: exp(-\beta E) E= h³ \sum_i E(i)
                        if (rand < exp(my_beta * minus_deltaE)) {
                            Site[i].Psi[alpha] = NewPsi;
                            acc_theta++;
                        }
                    }
                }
            }
        }
    }
    //t_localHtheta.toc();

    //t_localHA.tic();
    if (Hp.e != 0) {
    /**********VECTOR POTENTIAL UPDATE********/
        for (iz= 0; iz < Lz; iz++) {
            for (iy = 0; iy < Ly; iy++) {
                for (ix = 0; ix < Lx; ix++) {
                    i = ix + Lx * (iy + iz * Ly);
                        for (vec = 0; vec < 3; vec++) {
                            //Update of A
                            OldA = Site[i].A[vec];
                             oldE = local_HA(OldA, ix, iy, iz, vec, Hp, Site);
                             d_A = rn::uniform_real_box(-MCp.lbox_A, MCp.lbox_A);
                             NewA = OldA + d_A;
                             newE = local_HA(NewA, ix, iy, iz, vec, Hp, Site);
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
                                }
                            }
                        }
                    }
                }
            }
        }
    //t_localHA.tic();

    acc_l=(double) acc_l/(3*N);
    acc_theta=(double) acc_theta/(3*N);
    acc_A=(double) acc_A/(3*N);
    MCp.lbox_l= MCp.lbox_l*(0.5*(acc_l/acc_rate)+0.5);
    MCp.lbox_theta= MCp.lbox_theta*((0.5*acc_theta/acc_rate)+0.5);
    MCp.lbox_A= MCp.lbox_A*((0.5*acc_A/acc_rate)+0.5);
/*
    std::cout << "local_Hpsi" <<std::endl;
    t_localHpsi.print_measured_time_w_percent();
    std::cout << "local_Htheta"<<std::endl;
    t_localHtheta.print_measured_time_w_percent();
    std::cout << "local_HA"<<std::endl;
    t_localHA.print_measured_time_w_percent();
*/
}



double local_HPsi(struct O2 Psi, unsigned int ix, unsigned int iy, unsigned int iz, unsigned int alpha,  struct H_parameters &Hp, struct Node* Site){

    double h_Potential, h_Kinetic=0., h_Josephson=0., h_AB=0., h_tot;
    double inv_h2=1./(Hp.h*Hp.h);
    double inv_h=1./(Hp.h);
    double J_alpha1=0, J_beta1=0, J_alpha2=0, J_beta2=0;
    double gauge_phase1, gauge_phase2;
    unsigned int beta=0, vec=0, i;
    unsigned int nn_ip, nn_im;
    i=ix +Lx*(iy+Ly*iz);

    unsigned int ip=(ix == Lx-1 ? 0: ix+1);
    unsigned int ipx= ip+Lx*(iy+Ly*iz);
    unsigned int jp=(iy == Ly-1 ? 0: iy+1);
    unsigned int ipy= ix+(Lx*jp)+(Lx*Ly*iz);
    unsigned int kp=(iz == Lz-1 ? 0: iz+1);
    unsigned int ipz= ix+Lx*(iy+Ly*kp);

    unsigned int imx= (ix == 0 ? Lx-1: ix-1)+Lx*(iy+Ly*iz);
    unsigned int imy= ix+Lx*((iy == 0 ? Ly-1: iy-1)+Ly*iz);
    unsigned int imz= ix+Lx*(iy+Ly*(iz == 0 ? Lz-1: iz-1));

    //Compute the local Energy respect to a given component (alpha) of the matter field Psi and a given spatial position (r=(ix, iy, iz))
    //We need to compute just the part of the Hamiltonian involving Psi

    //Potential= (a+ 3/h²)*|Psi_{alpha}(r)|² + b/2*|Psi_{alpha}(r)|⁴
    h_Potential= O2norm2(Psi)*((Hp.a + (3.*inv_h2))+  (0.5*Hp.b*O2norm2(Psi)));

    //Kinetic= -(1/h²)*\sum_k=1,2,3 (|Psi_{alpha}(r)||Psi_{alpha}(r+k)|* cos(theta_{alpha}(r+k) - theta_{alpha}(r) +h*e*A_k(r))) + (|Psi_{alpha}(r-k)||Psi_{alpha}(r)|* cos(theta_{alpha}(r) - theta_{alpha}(r-k) +h*e*A_k(r-k)))

    for(vec=0; vec<3; vec++){
        if(vec==0){
            nn_ip=ipx;
            nn_im=imx;
        }else if(vec==1){
            nn_ip=ipy;
            nn_im=imy;
        }else if(vec==2) {
            nn_ip = ipz;
            nn_im = imz;
        }

        gauge_phase1=Site[nn_ip].Psi[alpha].t - Psi.t + Hp.h*Hp.e*Site[i].A[vec];
        gauge_phase2=Psi.t -Site[nn_im].Psi[alpha].t + Hp.h*Hp.e*Site[nn_im].A[vec];
        h_Kinetic-=inv_h2*(Psi.r*Site[nn_ip].Psi[alpha].r)*cos(gauge_phase1);
        h_Kinetic-=inv_h2*(Psi.r*Site[nn_im].Psi[alpha].r)*cos(gauge_phase2);

        //Andreev-Bashkin term = \sum_beta!=alpha \sum_k=1,2,3 nu*(J^k_alpha - J^k_beta)^2;
        // with J^k_alpha= |Psi_{alpha}(r)||Psi_{alpha}(r+k)|* sin(theta_{alpha}(r+k) - theta_{alpha}(r) +h*e*A_k(r)))
        if(Hp.nu !=0 ) {
            J_alpha1=inv_h*(Psi.r*Site[nn_ip].Psi[alpha].r)*sin(gauge_phase1);
            J_alpha2= inv_h*(Psi.r*Site[nn_im].Psi[alpha].r)*sin(gauge_phase2);
            for (beta = 0; beta < 3; beta++) {
                if (beta != alpha) {
                    //+k
                    J_beta1=inv_h*(Site[i].Psi[beta].r*Site[nn_ip].Psi[beta].r)*sin(Site[nn_ip].Psi[beta].t - Site[i].Psi[beta].t + Hp.h*Hp.e*Site[i].A[vec]);
                    h_AB += Hp.nu * (pow((J_alpha1 -J_beta1),2));
                    //-k
                    J_beta2=inv_h*(Site[i].Psi[beta].r*Site[nn_im].Psi[beta].r)*sin( Site[i].Psi[beta].t -Site[nn_im].Psi[beta].t + Hp.h*Hp.e*Site[nn_im].A[vec]);
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

double local_Htheta(struct O2 Psi, unsigned int ix, unsigned int iy, unsigned int iz, unsigned int alpha,  struct H_parameters &Hp, struct Node* Site){

    double h_Kinetic=0., h_Josephson=0., h_AB=0., h_tot;
    double inv_h2=1./(Hp.h*Hp.h);
    double inv_h=1./(Hp.h);
    double J_alpha1=0, J_beta1=0, J_alpha2=0, J_beta2=0;
    double gauge_phase1=0, gauge_phase2=0;
    unsigned int beta=0, vec=0, i;
    unsigned int nn_ip, nn_im;
    i=ix +Lx*(iy+Ly*iz);

    unsigned int ip=(ix == Lx-1 ? 0: ix+1);
    unsigned int ipx= ip+Lx*(iy+Ly*iz);
    unsigned int jp=(iy == Ly-1 ? 0: iy+1);
    unsigned int ipy= ix+(Lx*jp)+(Lx*Ly*iz);
    unsigned int kp=(iz == Lz-1 ? 0: iz+1);
    unsigned int ipz= ix+Lx*(iy+Ly*kp);

    unsigned int imx= (ix == 0 ? Lx-1: ix-1)+Lx*(iy+Ly*iz);
    unsigned int imy= ix+Lx*((iy == 0 ? Ly-1: iy-1)+Ly*iz);
    unsigned int imz= ix+Lx*(iy+Ly*(iz == 0 ? Lz-1: iz-1));
    //We need to compute just the part of the Hamiltonian involving Psi.t

    //Kinetic= -(1/h²)*\sum_k=1,2,3 (|Psi_{alpha}(r)||Psi_{alpha}(r+k)|* cos(theta_{alpha}(r+k) - theta_{alpha}(r) +h*e*A_k(r))) + (|Psi_{alpha}(r-k)||Psi_{alpha}(r)|* cos(theta_{alpha}(r) - theta_{alpha}(r-k) +h*e*A_k(r-k)))
    for(vec=0; vec<3; vec++){
        if(vec==0){
            nn_ip=ipx;
            nn_im=imx;
        }else if(vec==1){
            nn_ip=ipy;
            nn_im=imy;
        }else if(vec==2) {
            nn_ip = ipz;
            nn_im = imz;
        }
        gauge_phase1=Site[nn_ip].Psi[alpha].t - Psi.t + Hp.h*Hp.e*Site[i].A[vec];
        gauge_phase2=Psi.t -Site[nn_im].Psi[alpha].t + Hp.h*Hp.e*Site[nn_im].A[vec];
        h_Kinetic-=inv_h2*(Psi.r*Site[nn_ip].Psi[alpha].r)*cos(gauge_phase1);
        h_Kinetic-=inv_h2*(Psi.r*Site[nn_im].Psi[alpha].r)*cos(gauge_phase2);

       //Andreev-Bashkin term = \sum_beta!=alpha \sum_k=1,2,3 nu*(J^k_alpha - J^k_beta)^2;
        // with J^k_alpha= |Psi_{alpha}(r)||Psi_{alpha}(r+k)|* sin(theta_{alpha}(r+k) - theta_{alpha}(r) +h*e*A_k(r)))
        if(Hp.nu !=0 ) {

        J_alpha1= inv_h*(Psi.r*Site[nn_ip].Psi[alpha].r)*sin(gauge_phase1);
        J_alpha2= inv_h*(Psi.r*Site[nn_im].Psi[alpha].r)*sin(gauge_phase2);

        for (beta = 0; beta < 3; beta++) {
                if (beta != alpha) {
                J_beta1=inv_h*(Site[i].Psi[beta].r*Site[nn_ip].Psi[beta].r)*sin(Site[nn_ip].Psi[beta].t - Site[i].Psi[beta].t + Hp.h*Hp.e*Site[i].A[vec]);
                h_AB += Hp.nu * (pow((J_alpha1 -J_beta1),2));
                J_beta2=inv_h*(Site[i].Psi[beta].r*Site[nn_im].Psi[beta].r)*sin( Site[i].Psi[beta].t -Site[nn_im].Psi[beta].t + Hp.h*Hp.e*Site[nn_im].A[vec]);
                h_AB += Hp.nu * (pow((J_alpha2 -J_beta2),2));
                	}
        	}
        }
    }

    //Josephson= eta* \sum_beta!=alpha |Psi_{alpha}(r)||Psi_{beta}(r)|* cos(theta_{alpha}(r) - theta_{beta}(r))
    for(beta=0; beta<3; beta++){
	if(beta != alpha) {
            h_Josephson += (Hp.eta * Site[i].Psi[beta].r*Site[i].Psi[alpha].r*cos(Psi.t - Site[i].Psi[beta].t));
        }
    }

    h_tot= h_Kinetic + h_Josephson + h_AB;
    return h_tot;
}

double local_HA(double A, unsigned int ix, unsigned int iy, unsigned int iz,  unsigned int vec,  struct H_parameters &Hp, struct Node* Site){

    double h_Kinetic=0., h_B, h_AB=0., h_tot=0.;
    double inv_h2=1./(Hp.h*Hp.h);
    double inv_h=1./(Hp.h);
    double F2_A=0., F_A;
    double J_alpha, J_beta, J_gamma;
    double gauge_phase1;
    unsigned int alpha=0, i, l;

    unsigned int nn_ip, nn_im, nn_ipl, nn_iml;
    i=ix +Lx*(iy+Ly*iz);

    unsigned int ip=(ix == Lx-1 ? 0: ix+1);
    unsigned int ipx= ip+Lx*(iy+Ly*iz);
    unsigned int jp=(iy == Ly-1 ? 0: iy+1);
    unsigned int ipy= ix+(Lx*jp)+(Lx*Ly*iz);
    unsigned int kp=(iz == Lz-1 ? 0: iz+1);
    unsigned int ipz= ix+Lx*(iy+Ly*kp);

    unsigned int imx= (ix == 0 ? Lx-1: ix-1)+Lx*(iy+Ly*iz);
    unsigned int imy= ix+Lx*((iy == 0 ? Ly-1: iy-1)+Ly*iz);
    unsigned int imz= ix+Lx*(iy+Ly*(iz == 0 ? Lz-1: iz-1));

    if(vec==0){
        nn_ip=ipx;
        nn_im=imx;
    }else if(vec==1){
        nn_ip=ipy;
        nn_im=imy;
    }else if(vec==2) {
        nn_ip = ipz;
        nn_im = imz;
    }


    //Compute the local Energy respect to a given component (alpha) of the vector potential A and a given spatial position (r=(ix, iy, iz))
    //We need to compute just the part of the Hamiltonian involving A

    //Kinetic= -(1/h²)*\sum_k=1,2,3 |Psi_{alpha}(r)||Psi_{alpha}(r+k)|* cos(theta_{alpha}(r+k) - theta_{alpha}(r) +h*e*A_k(r))
    for(alpha=0; alpha<3; alpha++) {
        gauge_phase1=Site[nn_ip].Psi[alpha].t - Site[i].Psi[alpha].t + Hp.h*Hp.e*A;
        h_Kinetic -= inv_h2* Site[nn_ip].Psi[alpha].r*Site[i].Psi[alpha].r*cos(gauge_phase1);
	 }
    //Andreev-Bashkin term = \sum_beta!=alpha \sum_k=1,2,3 nu*(J^k_alpha - J^k_beta)^2;
    // with J^k_alpha= |Psi_{alpha}(r)||Psi_{alpha}(r+k)|* sin(theta_{alpha}(r+k) - theta_{alpha}(r) +h*e*A_k(r)))
    if(Hp.nu !=0 ) {

        J_alpha= inv_h * Site[nn_ip].Psi[0].r*Site[i].Psi[0].r*sin(Site[nn_ip].Psi[0].t - Site[i].Psi[0].t + Hp.h*Hp.e*A);
        J_beta= inv_h * Site[nn_ip].Psi[1].r*Site[i].Psi[1].r*sin(Site[nn_ip].Psi[1].t - Site[i].Psi[1].t + Hp.h*Hp.e*A);
        J_gamma= inv_h * Site[nn_ip].Psi[2].r*Site[i].Psi[2].r*sin(Site[nn_ip].Psi[2].t - Site[i].Psi[2].t + Hp.h*Hp.e*A);
//        J_beta= inv_h* O2vprod(Site[i].Psi[1], aux_field[1]);
//        J_gamma= inv_h * O2vprod(Site[i].Psi[2], aux_field[2]);

            h_AB += Hp.nu * (pow((J_alpha - J_beta), 2)  + pow((J_alpha - J_gamma), 2) + pow((J_gamma - J_beta), 2));
    }

    //All the plaquettes involving A_vec(i)
    for(l=0; l<3; l++){
        if(l==0){
            nn_ipl=ipx;
        }else if(l==1){
            nn_ipl=ipy;
        }else if(l==2) {
            nn_ipl = ipz;
        }
        if(l!= vec){
            F_A=(A + Site[nn_ip].A[l] - Site[nn_ipl].A[vec] - Site[i].A[l]);
            F2_A+=(F_A*F_A);
        }
    }
    for(l=0; l<3; l++){
        if(l==0){
            nn_iml=imx;
        }else if(l==1){
            nn_iml=imy;
        }else if(l==2) {
            nn_iml= imz;
        }
        if(l!= vec){
            F_A=(Site[nn_iml].A[vec]+ Site[nn(nn_iml, vec, 1)].A[l] - A - Site[nn_iml].A[l]);
            F2_A+=(F_A*F_A);
        }
    }

    h_B=(0.5*inv_h2)*F2_A;

    h_tot= h_Kinetic + h_AB+ h_B;
    return h_tot;
}

/*
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
*/
