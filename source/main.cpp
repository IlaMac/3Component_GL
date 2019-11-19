#include <cstring>
#include <string>
#include "main.h"
#include "initialization.h"
#include "measures.h"

int main(int argc, char *argv[]){

    struct Node* Lattice;
    struct H_parameters Hp;
    struct MC_parameters MCp;
    unsigned int i;
    std::string directory_read;
    std::string directory_write;
    std::string directory_parameters;

    if(argc > 4 ){
        printf("Too many arguments!");
        myhelp(argc, argv);
    }
    else if(argc < 4){
        printf("Not enough arguments --> Default Initialization. \n");
        directory_write = "/home/ilaria/Desktop/Multi_Components_GL/source";
    }
    else if(argc ==4) {
        directory_read  = argv[1];
        directory_write = argv[2];
        directory_parameters = argv[3];

    }

    //Declaration of structure Lattice
    Lattice=(struct Node*)calloc(N,sizeof(struct Node));
    for(i=0; i<N; i++) {
        Lattice[i].A = (double *) calloc(3, sizeof(double));
        Lattice[i].Psi = (struct O2 *) calloc(3, sizeof(struct O2));
    }

    //Initialize H_parameters: file "H_init.txt"
    initialize_Hparameters(Hp, directory_parameters);
    printf("beta %lf T %lf \n", Hp.beta, 1./Hp.beta);
    printf("a %d b %d eta %d h %lf e %lf \n", Hp.a, Hp.b, Hp.eta, Hp.h, Hp.e);

    //Initialize MC_parameters: file "MC_init.txt"
    initialize_MCparameters(MCp, directory_parameters);
    //Initialize Lattice: files "Psi_init.txt" and "A_init.txt" in the directory DIRECTORY_READ
    initialize_lattice(Lattice, directory_read);
    //Mainloop
    mainloop(Lattice, MCp, Hp, directory_write);

    return 0;
}

void mainloop(struct Node* Site, struct MC_parameters &MCp, struct H_parameters &Hp,  const fs::path & directory_write) {

    int n = 0, t = 0;
    struct Measures mis;
    std::ofstream Energy_file;
    std::ofstream Magnetization_file;
    std::ofstream DualStiff_file;
    char Efile_name[256];
    char Mfile_name[256];
    char DSfile_name[256];

    sprintf(Efile_name, "%s/Energy.txt", directory_write.c_str());
    sprintf(Mfile_name, "%s/Magnetization.txt", directory_write.c_str());
    sprintf(DSfile_name, "%s/Dual_Stiffness.txt", directory_write.c_str());

    Energy_file.open(Efile_name, std::ios::out);
    Energy_file.close();
    Energy_file.open(Efile_name, std::ios::app);
    Magnetization_file.open(Mfile_name, std::ios::app);
    DualStiff_file.open(DSfile_name, std::ios::app);

    for (n = 0; n<MCp.nmisu; n++) {

        for (t = 0; t < MCp.tau; t++) {
            metropolis(Site, MCp, Hp);
        }
        //Measures
        measures_init(mis);
        energy(mis, Hp, Site);
        Energy_file<<n<<"\t"<< mis.E<< "\t"<< mis.E_pot<< "\t"<< mis.E_kin<< "\t"<< mis.E_Josephson<< "\t"<< mis.E_B<<  std::endl;


        if ((n % MCp.n_autosave) == 0) {
            //Save a configuration for the restarting
        }

    }
    Energy_file.close();
    Magnetization_file.close();
    DualStiff_file.close();
}

unsigned int nn(unsigned int ix, unsigned int iy, unsigned int iz, unsigned int coord, int dir){

    int ix_new=0, iy_new=0, iz_new=0;
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
    fprintf(stderr,"%s <DIRECTORY_READ> <DIRECTORY_WRTE> <DIRECTORY_PARAMETERS> \n",argu[0]);
    exit(1);
    return;

}