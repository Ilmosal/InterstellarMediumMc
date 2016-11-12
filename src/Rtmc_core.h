#ifndef RTMC_CORE
#define RTMC_CORE

//Includes//
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

//Dust grid variables//
#define SIZE_OF_GRID 32
#define DUST_RHO 0.6
#define DUST_KSCA 1.0
#define DUST_ALBEDO 0.99
#define DUST_GRID_BETA -2

//Simulation variables//
#define PHOTON_PACKS 5000000 
#define THETA_RESULT_COEF 0.02
#define PHI_RESULT_COEF 0.062832
#define PHOTON_INTENSITY_RESULT_COEF 0.1

//Observer Variables//
#define OBS_PHI 0.0
#define OBS_THETA 0.0
#define OBS_DOTP_LIMIT 0.05

//Photon Variables//
#define RANDOM_PHOTON_DIR 0
#define PHOTON_START_X 0.0
#define PHOTON_START_Y 16.5
#define PHOTON_START_Z 16.5

//Other variables//
#define PI 3.14159265359

//Photon struct//
struct Photon {
	float pos[3];
	float dir[2];
	float intensity;
};

//Run function//
void run_RTMC();

//Grid Functions//
void init_grid(float dustGrid[SIZE_OF_GRID][SIZE_OF_GRID][SIZE_OF_GRID], int mode);
void uniform_grid(float dustGrid[SIZE_OF_GRID][SIZE_OF_GRID][SIZE_OF_GRID]);
void non_uniform_grid(float dustGrid[SIZE_OF_GRID][SIZE_OF_GRID][SIZE_OF_GRID]);
void cube2sphere(float dustGrid[SIZE_OF_GRID][SIZE_OF_GRID][SIZE_OF_GRID]);
void fast_fourier_shift(float GRID[SIZE_OF_GRID][SIZE_OF_GRID][SIZE_OF_GRID]);
void inverse_fast_fourier_shift(float GRID[SIZE_OF_GRID][SIZE_OF_GRID][SIZE_OF_GRID]);
void ifftn(complex *v, int n, complex *tmp);;

//Photon functions//
void init_phot(struct Photon *phot);
void scatter_photon(struct Photon *phot);
int isPhotonOutOfSystem(struct Photon phot);
float cellDensity(struct Photon phot, float dustGrid[SIZE_OF_GRID][SIZE_OF_GRID][SIZE_OF_GRID]);
void print_photon(struct Photon phot);
void photDirFloat(struct Photon phot, float dir[3]);

//Simulation functions//
struct Photon photon_run(float dustGrid[SIZE_OF_GRID][SIZE_OF_GRID][SIZE_OF_GRID]);
void photon_run_viewFile(float result[SIZE_OF_GRID][SIZE_OF_GRID],
						 float dustGrid[SIZE_OF_GRID][SIZE_OF_GRID][SIZE_OF_GRID],
						 char *progressionS);

//Math functions//
float distanceBetween(float vec1[3], float vec2[3]);
float distanceToNextCell(struct Photon phot);
int isRayInTheSameCell(struct Photon phot, float rayPos[3]);
float randFloat();
float dotProduct(float vec1[3], float vec2[3]);

#endif
