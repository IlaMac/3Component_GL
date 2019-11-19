
#ifndef MAIN_H
#define MAIN_H
#include<cstdio>
#include<cmath>
#include<cstdlib>
#include <iostream>
#include <fstream>
#include <random>
#include "o2.h"
#include "montecarlo.h"
#include "initialization.h"
#include <filesystem>

namespace fs = std::filesystem;


/*Temporary defined as global variables here*/
#define Lx 8
#define Ly 8
#define Lz 8

#define N (Lx*Ly*Lz)

void mainloop(struct Node* Site, struct MC_parameters &MCp, struct H_parameters &Hp,  const fs::path & directory_write);
unsigned int nn(unsigned int ix, unsigned int iy, unsigned int iz, unsigned int coord, int dir);
void myhelp(int argd, char** argu);
#endif
