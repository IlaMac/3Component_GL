#include <cstring>
#include <string>
#include "main.h"
#include "initialization.h"
#include "measures.h"
#include "rng.h"
#include "memory_check.h"


int main(int argc, char *argv[]){

    struct Node* Lattice;
    struct H_parameters Hp;
    struct MC_parameters MCp;
    struct PT_parameters PTp;
    struct PTroot_parameters PTroot;
    unsigned int i;
    int p;
    long int seednumber=-1; /*by default it is a negative number which means that rng will use random_device*/

    double my_beta=0.226;

    std::string directory_read;
    std::string directory_write;
    std::string directory_parameters;

    if(argc > 3 ){
        printf("Too many arguments!");
        myhelp(argc, argv);
    }
    else if(argc < 2){
        printf("Not enough arguments --> Default Initialization. \n");
        myhelp(argc, argv);
    }
    else if(argc ==2) {
        directory_parameters = argv[1];
    }
    else if(argc == 3){
        directory_parameters = argv[1];
        seednumber= reinterpret_cast<long> (argv[2]);
    }

    //initialization of the random number generator
    rn::seed(seednumber);

    //Declaration of structure Lattice
    Lattice=(struct Node*)calloc(N,sizeof(struct Node));
    for(i=0; i<N; i++) {
        Lattice[i].A = (double *) calloc(3, sizeof(double));
        Lattice[i].Psi = (struct O2 *) calloc(3, sizeof(struct O2));
    }

    //Initialize H_parameters: file "H_init.txt"
    initialize_Hparameters(Hp, directory_parameters);
    //Initialize MC_parameters: file "MC_init.txt"
    initialize_MCparameters(MCp, directory_parameters);

    MPI_Status status;
    MPI_Init(NULL, NULL); /* START MPI */
/*DETERMINE RANK OF THIS PROCESSOR*/
    MPI_Comm_rank(MPI_COMM_WORLD, &PTp.rank);
/*DETERMINE TOTAL NUMBER OF PROCESSORS*/
    MPI_Comm_size(MPI_COMM_WORLD, &PTp.np);


    if(PTp.rank == PTp.root) {
        //Initialization ranks arrays
        initialize_PTarrays( PTp, PTroot, Hp);
    }
    MPI_Scatter(PTroot.beta.data(), 1, MPI_DOUBLE, &my_beta, 1, MPI_DOUBLE, PTp.root, MPI_COMM_WORLD);

    printf("I'm rank %d and this is my beta %lf\n", PTp.rank, my_beta);

    directory_write=directory_parameters+"/rank_"+std::to_string(PTp.rank);
    directory_read=directory_parameters+"/rank_"+std::to_string(PTp.rank);

    //Initialize Lattice: files "Psi_init.txt" and "A_init.txt" in the directory DIRECTORY_READ  which depends on the rank!!!!!I HAVE TO INITIALIZE ALSO MY_BETA FROM THE SAME FOLDER!!!!!!!
    initialize_lattice(Lattice, directory_read);
    //Mainloop
    mainloop(Lattice, MCp, Hp, my_beta, PTp, PTroot, directory_write);

    std::cout << "Proccess current resident ram usage: " << process_memory_in_mb("VmRSS") << " MB" << std::endl;
    std::cout << "Proccess maximum resident ram usage: " << process_memory_in_mb("VmHWM") << " MB" << std::endl;
    std::cout << "Proccess maximum virtual  ram usage: " << process_memory_in_mb("VmPeak") << " MB" << std::endl;

    return 0;
}

void mainloop(struct Node* Site, struct MC_parameters &MCp, struct H_parameters &Hp, double my_beta, struct PT_parameters PTp, struct PTroot_parameters PTroot, const fs::path & directory_write) {

    int n = 0, t = 0;
    /*Measurements*/
    struct Measures mis;
    std::ofstream Energy_file;
    std::ofstream Magnetization_file;
    std::ofstream DualStiff_file;
    std::ofstream DensityPsi_file;

    //  fs::path energy_file = directory_write / "Energy.txt"; Use this sintax to send the path of the file to the measures function so that the writing to a file would occur there.
    char Efile_name[256];
    char Mfile_name[256];
    char DSfile_name[256];
    char DPsi_name[256];

    sprintf(Efile_name, "%s/Energy.txt", directory_write.c_str());
    sprintf(Mfile_name, "%s/Magnetization.txt", directory_write.c_str());
    sprintf(DSfile_name, "%s/Dual_Stiffness.txt", directory_write.c_str());
    sprintf(DPsi_name, "%s/Psi_density.txt", directory_write.c_str());

    Energy_file.open(Efile_name, std::ios::out);
    Energy_file.close();
    Energy_file.open(Efile_name, std::ios::app);
    DualStiff_file.open(DSfile_name, std::ios::out);
    DualStiff_file.close();
    DualStiff_file.open(DSfile_name, std::ios::app);
    Magnetization_file.open(Mfile_name, std::ios::out);
    Magnetization_file.close();
    Magnetization_file.open(Mfile_name, std::ios::app);
    DensityPsi_file.open(DPsi_name, std::ios::out);
    DensityPsi_file.close();
    DensityPsi_file.open(DPsi_name, std::ios::app);

    measures_init(mis);

    for (n = 0; n<MCp.nmisu; n++) {
        for (t = 0; t < MCp.tau; t++) {
            metropolis(Site, MCp, Hp,  my_beta);
        }
        //Measures
        measures_reset(mis);
        energy(mis, Hp, my_beta, Site);
        Energy_file<<n<<"\t"<<my_beta<<"\t"<< mis.E<< "\t"<< mis.E_pot<< "\t"<< mis.E_kin<< "\t"<< mis.E_Josephson<< "\t"<< mis.E_B<<  std::endl;
        dual_stiffness(mis, Hp, Site);
        DualStiff_file<<n<<"\t"<<my_beta<<"\t"<<mis.d_rhoz<<std::endl;
        magnetization(mis, Site);
        Magnetization_file<<n<<"\t"<<my_beta<<"\t"<<mis.m<<"\t"<<(mis.m*mis.m)<<"\t"<<(mis.m*mis.m*mis.m*mis.m)<<std::endl;
        density_psi(mis, Site);
        DensityPsi_file<<n<<"\t"<<my_beta<<"\t"<<mis.density_psi[0]<<"\t"<<mis.density_psi[1]<<"\t"<<mis.density_psi[2]<<std::endl;
        //Parallel Tempering swap
        parallel_temp(mis.E, my_beta, PTp, PTroot);

        if ((n % MCp.n_autosave) == 0) {
            //Save a configuration for the restarting
            save_lattice(Site, directory_write, std::string("n") + std::to_string(n) );
            }
    }
    save_lattice(Site, directory_write, std::string("final"));
    Energy_file.close();
    Magnetization_file.close();
    DualStiff_file.close();
}

void parallel_temp(double my_E , double my_beta, struct PT_parameters PTp, struct PTroot_parameters PTroot){

    double coin;
    double n_rand=0., delta_E, delta_beta;
    double beta_lower, beta_higher;
    int p, i=0, nn=0;

    MPI_Gather(&my_E, 1, MPI_DOUBLE, PTroot.All_Energies.data(), 1, MPI_DOUBLE, PTp.root, MPI_COMM_WORLD);
    if (PTp.rank == PTp.root) { //Root forms the pairs and decides (given the energies and the betas) which pairs will swap
        //Pair Formation
        coin = rn::uniform_real_box(0,1);
        if(coin < 0.5) { //each even rank wil be paired with its right neighbour
            nn= +1;
        }else if(coin >= 0.5){ //each even rank wil be paired with its left neighbour
            nn=-1;
        }
        while (i < PTp.np) {
            n_rand=rn::uniform_real_box(0,1);
            delta_E = PTroot.All_Energies[PTroot.ind_to_rank[i]] - PTroot.All_Energies[PTroot.ind_to_rank[(PTp.np + i + nn) % PTp.np]];
//            printf("rank %d b1 %lf rank2 %d b2 %lf \n", PTroot.ind_to_rank[i], PTroot.beta[PTroot.ind_to_rank[i]],
//            PTroot.ind_to_rank[(PTp.np + i + nn) % PTp.np], PTroot.beta[PTroot.ind_to_rank[(PTp.np+ i + nn) % PTp.np]]);
            delta_beta = PTroot.beta[PTroot.ind_to_rank[i]] - PTroot.beta[PTroot.ind_to_rank[(PTp.np + i + nn) % PTp.np]];
//            printf("index %d rand %lf delta_b %lf delta_E %lf exp %lf\n", i, n_rand, delta_beta, delta_E,exp(-delta_beta * delta_E));
            //swapping condition
            if (n_rand < exp(-delta_beta * delta_E)) {
//                printf("SWAP i %d and j %d\n", i, (PTp.np + i + nn) % PTp.np);
                //swap indices
                PTroot.rank_to_ind[PTroot.ind_to_rank[i]] = (PTp.np + i + nn) % PTp.np;
                PTroot.rank_to_ind[PTroot.ind_to_rank[(PTp.np + i + nn) % PTp.np]] = i;
                //swap beta
                beta_lower= PTroot.beta[PTroot.ind_to_rank[i]];
                beta_higher= PTroot.beta[PTroot.ind_to_rank[(PTp.np + i + nn) % PTp.np]];
                PTroot.beta[PTroot.ind_to_rank[i]] = beta_higher;
                PTroot.beta[PTroot.ind_to_rank[(PTp.np + i + nn) % PTp.np]] = beta_lower;
                }
                i+= 2;
            }
            for (p = 0; p < PTp.np; p++) {
                PTroot.ind_to_rank[PTroot.rank_to_ind[p]] = p;
            }
        }
    MPI_Scatter(PTroot.beta.data(), 1, MPI_DOUBLE, &my_beta, 1, MPI_DOUBLE, PTp.root, MPI_COMM_WORLD);
}


unsigned int nn(unsigned int i, unsigned int coord, int dir){

    unsigned int ix, iy, iz;
    int ix_new=0, iy_new=0, iz_new=0;

    ix=i%Lx;
    iy=(i/Lx)%Ly;
    iz=(i/(Lx*Ly));

    if(coord==0){
        ix_new= ix + (int)dir/sqrt(dir*dir);
        if(ix_new==Lx){ ix_new=0;}
        if(ix_new < 0){ ix_new=(Lx-1);}
        iy_new=iy;
        iz_new=iz;
    }
    if(coord==1){
        iy_new= iy + dir/sqrt(dir*dir);
        if(iy_new==Ly){ iy_new=0;}
        if(iy_new<0){ iy_new=(Ly-1);}
        ix_new=ix;
        iz_new=iz;
    }
    if(coord==2){
        iz_new= iz + dir/sqrt(dir*dir);
        if(iz_new==Lz){ iz_new=0;}
        if(iz_new<0){ iz_new=(Lz-1);}
        ix_new=ix;
        iy_new=iy;
    }
    return (ix_new+ Lx*(iy_new+ iz_new*Ly));
}

void myhelp(int argd, char** argu) {
    int i;
    fprintf(stderr,"Errore nei parametri su linea di comando; hai scritto:\n");
    for (i=0;i<argd;i++) fprintf(stderr," %s",argu[i]);
    fprintf(stderr,"\n");
    fprintf(stderr,"%s <DIRECTORY_PARAMETERS> <SEED> \n",argu[0]);
    exit (EXIT_FAILURE);
    return;
}
