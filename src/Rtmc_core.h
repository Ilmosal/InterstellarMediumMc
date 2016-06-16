#ifndef RTMC_CORE
#define RTMC_CORE

//Includes//
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <stdlib.h>

//Dust grid variables//
#define SIZE_OF_GRID 32
#define DUST_RHO 0.2
#define DUST_KSCA 1.0
#define DUST_ALBEDO 0.99
#define DUST_GRID_BETA -2

//Simulation variables//
#define PHOTON_PACKS 500000 
#define THETA_RESULT_COEF 0.02
#define PHI_RESULT_COEF 0.062832
#define PHOTON_INTENSITY_RESULT_COEF 0.01

//Photon Variables//
#define RANDOM_PHOTON_DIR 0
#define PHOTON_START_X 0.0
#define PHOTON_START_Y 16.5
#define PHOTON_START_Z 16.5

//Other variables//
#define PI 3.14159265359

//Density Modes//
enum
{
	UNIFORM_DENSITY,
	NON_UNIFORM_DENSITY
};

//Photon struct//
struct Photon {
	float pos[3];
	float dir[2];
	float intensity;
};

//Run function 
void run_RTMC();

//Grid Functions
void init_grid(float dustGrid[SIZE_OF_GRID][SIZE_OF_GRID][SIZE_OF_GRID], int mode);
void uniform_grid(float dustGrid[SIZE_OF_GRID][SIZE_OF_GRID][SIZE_OF_GRID]);
void non_uniform_grid(float dustGrid[SIZE_OF_GRID][SIZE_OF_GRID][SIZE_OF_GRID]);
void fast_fourier_shift(float GRID[SIZE_OF_GRID][SIZE_OF_GRID][SIZE_OF_GRID]);
void inverse_fast_fourier_shift(float GRID[SIZE_OF_GRID][SIZE_OF_GRID][SIZE_OF_GRID]);
void ifftn(complex *v, int n, complex *tmp);;

//Photon functions
void init_phot(struct Photon *phot);
void scatter_photon(struct Photon *phot);
int isPhotonOutOfSystem(struct Photon phot);
float cellDensity(struct Photon phot, float dustGrid[SIZE_OF_GRID][SIZE_OF_GRID][SIZE_OF_GRID]);
void print_photon(struct Photon phot);

//Simulation functions
struct Photon photon_run(float dustGrid[SIZE_OF_GRID][SIZE_OF_GRID][SIZE_OF_GRID]);

//Math functions
float distanceBetween(float vec1[3], float vec2[3]);
float distanceToNextCell(struct Photon phot);
int isRayInTheSameCell(struct Photon phot, float rayPos[3]);
float randFloat();

#endif
