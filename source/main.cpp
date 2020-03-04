#include <cstring>
#include <string>
#include "main.h"
#include "initialization.h"
#include "measures.h"
#include "rng.h"
#include "memory_check.h"
#include <h5pp/h5pp.h>


unsigned int Lx, Ly, Lz, N;

int main(int argc, char *argv[]){

    struct Node* Lattice;
    struct H_parameters Hp;
    struct MC_parameters MCp;
    struct PT_parameters PTp;
    struct PTroot_parameters PTroot;
    unsigned int i;
    long int seednumber=-1; /*by default it is a negative number which means that rng will use random_device*/
    double my_beta=0.244;
    int my_ind=0;
    time_t tic, toc;

    std::string directory_read;
    std::string directory_parameters;

    if(argc > 4 ){
        printf("Too many arguments!");
        myhelp(argc, argv);
    }
    else if(argc < 3){
        printf("Not enough arguments --> Default Initialization. \n");
        myhelp(argc, argv);
    }
    else if(argc ==3) {
        /*Rude way*/
        Lx=Ly=Lz=std::atoi(argv[1]);
        N=Lx*Ly*Lz;
        directory_parameters = argv[2];
    }
    else if(argc == 4){
        Lx=Ly=Lz=std::atoi(argv[1]);
        N=Lx*Ly*Lz;
        directory_parameters = argv[2];
        seednumber= reinterpret_cast<long> (argv[2]);
    }

    //initialization of the random number generator
    rn::seed(seednumber);

    //Declaration of structure Lattice
    Lattice=(struct Node*)calloc(N,sizeof(struct Node));
    for(i=0; i<N; i++) {
        Lattice[i].A = (double *) calloc(3, sizeof(double));
        Lattice[i].Psi = (struct O2 *) calloc(NC, sizeof(struct O2));
    }

    //Initialize H_parameters: file "H_init.txt"
    initialize_Hparameters(Hp, directory_parameters);
    //Initialize MC_parameters: file "MC_init.txt"
    initialize_MCparameters(MCp, directory_parameters);

    MPI_Init(NULL, NULL); /* START MPI */
/*DETERMINE RANK OF THIS PROCESSOR*/
    MPI_Comm_rank(MPI_COMM_WORLD, &PTp.rank);
/*DETERMINE TOTAL NUMBER OF PROCESSORS*/
    MPI_Comm_size(MPI_COMM_WORLD, &PTp.np);

    if(PTp.rank == PTp.root) {
        //Initial time
        tic = time(0);
    }

    if(PTp.rank == PTp.root) {
        //Initialization ranks arrays
        initialize_PTarrays( PTp, PTroot, Hp);
    }
    MPI_Scatter(PTroot.beta.data(), 1, MPI_DOUBLE, &my_beta, 1, MPI_DOUBLE, PTp.root, MPI_COMM_WORLD);
    MPI_Scatter(PTroot.rank_to_ind.data(), 1, MPI_INT, &my_ind, 1, MPI_INT, PTp.root, MPI_COMM_WORLD);

    printf("I'm rank %d and this is my beta %lf\n", PTp.rank, my_beta);

    directory_read=directory_parameters+"/beta_"+std::to_string(my_ind);

    //Initialize Lattice: files "Psi_init.txt" and "A_init.txt" in the directory DIRECTORY_READ!!! I HAVE TO INITIALIZE PROPERLY!!!!!!!
    initialize_lattice(Lattice, directory_read);
    //Mainloop
    mainloop(Lattice, MCp, Hp, my_beta, my_ind, PTp, PTroot, directory_parameters);

    if(PTp.rank == PTp.root) {
        //Final time
        toc = time(0);
    }


    std::cout << "Proccess current resident ram usage: " << process_memory_in_mb("VmRSS") << " MB" << std::endl;
    std::cout << "Proccess maximum resident ram usage: " << process_memory_in_mb("VmHWM") << " MB" << std::endl;
    std::cout << "Proccess maximum virtual  ram usage: " << process_memory_in_mb("VmPeak") << " MB" << std::endl;

    if(PTp.rank == PTp.root) {
        std::cout << "Total runtime: "
                  << ((toc - tic) / 3600 > 0 ? std::to_string((toc - tic) / 3600) + " h " : "")
                  << ((toc - tic) / 60 > 0 ? std::to_string(((toc - tic) % 3600) / 60) + " min " : "")
                  << (toc - tic) % 60 << " sec" << "\n" << std::endl;
    }

    return 0;
}

void mainloop(struct Node* Site, struct MC_parameters &MCp, struct H_parameters &Hp, double &my_beta, int &my_ind, struct PT_parameters PTp, struct PTroot_parameters PTroot, std::string directory_parameters) {

    int n = 0, t = 0;

    /*Measurements*/
    Measures mis;

    std::string directory_write;
    directory_write=directory_parameters+"/beta_"+std::to_string(my_ind);

//    // Initialize a file
    h5pp::File file(directory_write+"/Output.h5", h5pp::AccessMode::READWRITE, h5pp::CreateMode::TRUNCATE);
//    // Register the compound type
    h5pp::hid::h5t MY_HDF5_MEASURES_TYPE = H5Tcreate(H5T_COMPOUND, sizeof(Measures));
    H5Tinsert(MY_HDF5_MEASURES_TYPE, "E", HOFFSET(Measures, E), H5T_NATIVE_DOUBLE);
    H5Tinsert(MY_HDF5_MEASURES_TYPE, "E_pot", HOFFSET(Measures, E_pot), H5T_NATIVE_DOUBLE);
    H5Tinsert(MY_HDF5_MEASURES_TYPE, "E_kin", HOFFSET(Measures, E_kin), H5T_NATIVE_DOUBLE);
    H5Tinsert(MY_HDF5_MEASURES_TYPE, "E_Josephson", HOFFSET(Measures, E_Josephson), H5T_NATIVE_DOUBLE);
    H5Tinsert(MY_HDF5_MEASURES_TYPE, "E_B", HOFFSET(Measures, E_B), H5T_NATIVE_DOUBLE);
    H5Tinsert(MY_HDF5_MEASURES_TYPE, "m", HOFFSET(Measures, m), H5T_NATIVE_DOUBLE);
    H5Tinsert(MY_HDF5_MEASURES_TYPE, "ds", HOFFSET(Measures, d_rhoz), H5T_NATIVE_DOUBLE);
    H5Tinsert(MY_HDF5_MEASURES_TYPE, "rho", HOFFSET(Measures, density_psi), H5T_NATIVE_DOUBLE);
    H5Tinsert(MY_HDF5_MEASURES_TYPE, "rank", HOFFSET(Measures, my_rank), H5T_NATIVE_INT);

    file.createTable(MY_HDF5_MEASURES_TYPE, "Measurements", "Measures");

    for (n = 0; n<MCp.nmisu; n++) {
        for (t = 0; t < MCp.tau; t++) {
            metropolis(Site, MCp, Hp,  my_beta);
        }

        //Measures
        mis.reset();
        energy(mis, Hp, my_beta, Site);
        dual_stiffness(mis, Hp, Site);
        magnetization(mis, Site);
        density_psi(mis, Site);
        mis.my_rank=PTp.rank;

        file.appendTableEntries(mis, "Measurements");

        if ((n % MCp.n_autosave) == 0) {
            //Save a configuration for the restarting
            save_lattice(Site, directory_write, std::string("n") + std::to_string(n) );
            }

        //Parallel Tempering swap
        parallel_temp(mis.E, my_beta, my_ind, PTp, PTroot);

        //Files and directory
        directory_write=directory_parameters+"/beta_"+std::to_string(my_ind);
        file = h5pp::File(directory_write+"/Output.h5", h5pp::AccessMode::READWRITE, h5pp::CreateMode::OPEN,0);

    }
    save_lattice(Site, directory_write, std::string("final"));


}

void parallel_temp(double &my_E , double &my_beta, int &my_ind, struct PT_parameters &PTp, struct PTroot_parameters &PTroot){

    double coin;
    double n_rand, delta_E, delta_beta;
    double oldbeta_i, oldbeta_nn;
    int i=0, nn=0, ind_nn=0;
    int oldrank_i=0, oldrank_nn=0;
    int newrank_i=0, newrank_nn=0;


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
            ind_nn=(PTp.np + i + nn) % PTp.np;
            oldrank_i=PTroot.ind_to_rank[i];
            oldrank_nn=PTroot.ind_to_rank[ind_nn];
            delta_E = PTroot.All_Energies[oldrank_i] - PTroot.All_Energies[oldrank_nn];
            delta_beta = PTroot.beta[oldrank_i] - PTroot.beta[oldrank_nn];
            //swapping condition
            //Boltzmann weight: exp(-\beta E) E= h³ \sum_i E(i)
            if (n_rand < exp(delta_beta * delta_E)) {
                //swap indices in the rank_to_ind array

                PTroot.rank_to_ind[oldrank_i] = ind_nn;
                PTroot.rank_to_ind[oldrank_nn] = i;

                //swap indices in the ind_to_rank array
                newrank_i= oldrank_nn;
                PTroot.ind_to_rank[i]= newrank_i;
                newrank_nn=oldrank_i;
                PTroot.ind_to_rank[ind_nn] =newrank_nn;
                //swap beta

                oldbeta_i= PTroot.beta[oldrank_i];
                oldbeta_nn= PTroot.beta[oldrank_nn];
                PTroot.beta[oldrank_i] = oldbeta_nn;
                PTroot.beta[oldrank_nn] = oldbeta_i;
            }
                i+= 2;
        }
    }
    MPI_Scatter(PTroot.beta.data(), 1, MPI_DOUBLE, &my_beta, 1, MPI_DOUBLE, PTp.root, MPI_COMM_WORLD);
    MPI_Scatter(PTroot.rank_to_ind.data(), 1, MPI_INT, &my_ind, 1, MPI_INT, PTp.root, MPI_COMM_WORLD);

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
