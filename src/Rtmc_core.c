#include "Rtmc_core.h"

//Main loop for simulation//
void run_RTMC(int mode)
{
	//Declaring necessary variables//
	float results[SIZE_OF_GRID][SIZE_OF_GRID], dustGrid[SIZE_OF_GRID][SIZE_OF_GRID][SIZE_OF_GRID], dir[3], norm[3],
	res1[SIZE_OF_GRID][SIZE_OF_GRID], res2[SIZE_OF_GRID][SIZE_OF_GRID], res3[SIZE_OF_GRID][SIZE_OF_GRID], res4[SIZE_OF_GRID][SIZE_OF_GRID];
	struct Photon phot;
	int i, j, photonInGrid;
	FILE *fp;
	time_t t;

	omp_set_num_threads(4);

	printf("\n-Interstellar photon scattering simulation-\n");
	printf("-------------------------------------------\n");
	printf("MODE: %d - PHOTONS: %d - RHO: %.2f\n", mode, PHOTON_PACKS, DUST_RHO);
	printf("-------------------------------------------\n");

	//Initialising random seed//
	srand((unsigned) time(&t));

	//Initialising dust grid//
	init_grid(dustGrid, mode);

	//Opening output file 
	if (mode == 0)
		fp = fopen("../run/uniform-outputFile", "w+");
	else if (mode == 1)
		fp = fopen("../run/non-uniform-outputFile", "w+");
	else if (mode == 2)
		fp = fopen("../run/uniform-viewFile", "w+");	
	else if (mode == 3)
		fp = fopen("../run/non-uniform-viewFile", "w+");
	
	//Printing the amount of photon packs in the simulation
	fprintf(fp, "%d\n", PHOTON_PACKS);

	//Running the simulation//	
	if (mode == 0 || mode== 1)
	{
		//Printing the amount of photon packs in the simulation
		fprintf(fp, "%d\n", PHOTON_PACKS);

		for (i = 0; i < PHOTON_PACKS; i++)
		{ 
			if (i % (PHOTON_PACKS / 20) == 0)
				printf("Progress: %d %c \n", ((i*100)+100)/PHOTON_PACKS, '%');	
	
			//Running the route for one photon
			phot = photon_run(dustGrid);
		
			//printing the directional information of the photon to an output file
			fprintf(fp, "%f\n%f\n%f\n", phot.dir[0] / PI, phot.dir[1], phot.intensity);	
		}
	} else if (mode == 2 || mode == 3)
	{
		//Initilising the result grid
		for (i = 0; i < SIZE_OF_GRID; i++)
		{
			for (j = 0; j < SIZE_OF_GRID; j++)
			{
				results[i][j] = 0; res1[i][j] = 0; res2[i][j] = 0; res3[i][j] = 0; res4[i][j] = 0;
			}
		}

		#pragma omp parallel sections
		{
			#pragma omp section
			photon_run_viewFile(res1, dustGrid);

			#pragma omp section
			photon_run_viewFile(res2, dustGrid);

			#pragma omp section
			photon_run_viewFile(res3, dustGrid);

			#pragma omp section
			photon_run_viewFile(res1, dustGrid);
		}	

		for (i = 0; i < SIZE_OF_GRID; i++)
			for (j = 0; j < SIZE_OF_GRID; j++)
				results[i][j] += res1[i][j] + res2[i][j] + res3[i][j] + res4[i][j];

/*
		#pragma omp parallel shared(results) private(dustGrid, res1, i, j) for
		{
			photon_run_viewFile(res1, dustGrid);

			#pragma omp critical
			{	
				float test = 0;

				for (i = 0; i < SIZE_OF_GRID; i++)
					for (j = 0; j < SIZE_OF_GRID; j++)
						results[i][j] += res1[i][j];
		
				for (i = 0; i < SIZE_OF_GRID; i++)
					for (j = 0; j < SIZE_OF_GRID; j++)
						test += res1[i][j];
	
				printf("%f\n", test);
			}
		}

		photonInGrid = 0;

		for (i = 0; i < PHOTON_PACKS; i++)
		{ 
			if (i % (PHOTON_PACKS / 20) == 0)
				printf("Progress: %d %c \n", ((i*100)+100)/PHOTON_PACKS, '%');	
	
			//Running the route for one photon
			phot = photon_run(dustGrid);

			//Calculating if the photon hits the grid
			if (phot.pos[0] > 32)
			{
				if (phot.pos[1] > 0 && phot.pos[1] < 32 && phot.pos[2] > 0 && phot.pos[2] < 32)
				{
					photDirFloat(phot, dir); 
					norm[0] = 1; norm[1] = 0; norm[2] = 0;
	
					if (OBS_DOTP_LIMIT > 1 - dir[0] * norm[0])
					{
						results[((int) phot.pos[1])][((int) phot.pos[2])] += phot.intensity;
						photonInGrid++;
					}
				} 
			} 
		}	
*/
		
		//Writing the data to a file
		fprintf(fp, "%d\n", SIZE_OF_GRID);

		for (i = 0; i < SIZE_OF_GRID; i++)
			for (j = 0; j < SIZE_OF_GRID; j++)
				fprintf(fp, "%f\n", results[i][j]);
	}	
	
	fclose(fp);
}

//Grid initialisation function//
void init_grid(float dustGrid[SIZE_OF_GRID][SIZE_OF_GRID][SIZE_OF_GRID], int mode)
{
	if (mode == 0 || mode == 2)
		uniform_grid(dustGrid);
	else if (mode == 1 || mode == 3)
		non_uniform_grid(dustGrid);

}

//function for creating an uniform grid//
void uniform_grid(float dustGrid[SIZE_OF_GRID][SIZE_OF_GRID][SIZE_OF_GRID])
{
	int i, j, k;
	for (i = 0; i < SIZE_OF_GRID; i++)
		for (j = 0; j < SIZE_OF_GRID; j++)
			for (k = 0; k < SIZE_OF_GRID; k++)
				dustGrid[i][j][k] = DUST_RHO;
	
	cube2sphere(dustGrid);
}

//Function for creating a non-uniform grid//
void non_uniform_grid(float dustGrid[SIZE_OF_GRID][SIZE_OF_GRID][SIZE_OF_GRID])
{
	int i, j, k;
	float I_GRID[SIZE_OF_GRID][SIZE_OF_GRID][SIZE_OF_GRID],
		  J_GRID[SIZE_OF_GRID][SIZE_OF_GRID][SIZE_OF_GRID],
		  K_GRID[SIZE_OF_GRID][SIZE_OF_GRID][SIZE_OF_GRID],
		  R_GRID[SIZE_OF_GRID][SIZE_OF_GRID][SIZE_OF_GRID],
		  A_GRID[SIZE_OF_GRID][SIZE_OF_GRID][SIZE_OF_GRID],
		  F_GRID[SIZE_OF_GRID][SIZE_OF_GRID][SIZE_OF_GRID],
		  NORM_GRID[SIZE_OF_GRID][SIZE_OF_GRID][SIZE_OF_GRID];

	float sum, A, norm_sum;

	fftw_complex in[SIZE_OF_GRID][SIZE_OF_GRID][SIZE_OF_GRID],
		    out[SIZE_OF_GRID][SIZE_OF_GRID][SIZE_OF_GRID];
	fftw_plan p;

	//Initialising plan/
	p = fftw_plan_dft_3d(SIZE_OF_GRID, SIZE_OF_GRID, SIZE_OF_GRID, &in[0][0][0], &out[0][0][0], 
			     FFTW_BACKWARD, FFTW_PATIENT);
	
	//Starting generation//
	for (i = 0; i < SIZE_OF_GRID; i++)
	{
		for (j = 0; j < SIZE_OF_GRID; j++)
		{
			for (k = 0; k < SIZE_OF_GRID; k++)
			{	
				//Fill the arrays with values//
				I_GRID[i][j][k] = i - 0.5 * ( SIZE_OF_GRID - 1);
				J_GRID[i][j][k] = j - 0.5 * ( SIZE_OF_GRID - 1);
				K_GRID[i][j][k] = k - 0.5 * ( SIZE_OF_GRID - 1);

				//Power function//
				R_GRID[i][j][k] = sqrt(I_GRID[i][j][k] * I_GRID[i][j][k] + 
						  J_GRID[i][j][k] * J_GRID[i][j][k] + 
						  K_GRID[i][j][k] * K_GRID[i][j][k]);
				
				//Amplitude from power//
				A_GRID[i][j][k] = pow(R_GRID[i][j][k], DUST_GRID_BETA / 2);

				//Making the normalization Grid for later//
				NORM_GRID[i][j][k] = DUST_RHO;
			}
		}
	}

	fast_fourier_shift(A_GRID);

	//Remove DC = offset //
 	for (i = 0; i < SIZE_OF_GRID; i++)
		A_GRID[0][0][i] = 0.0;


	//Setting A into a array of randomized values between [0,1]
	for (i = 0; i < SIZE_OF_GRID; i++)
		for (j = 0; j < SIZE_OF_GRID; j++)
			for (k = 0; k < SIZE_OF_GRID; k++)
				F_GRID[i][j][k] = randFloat();

	//Combining the grids
	for (i = 0; i < SIZE_OF_GRID; i++)
		for (j = 0; j < SIZE_OF_GRID; j++)
			for (k = 0; k < SIZE_OF_GRID; k++)
				in[i][j][k] = (fftw_complex) A_GRID[i][j][k] * cexp(I*2*PI*F_GRID[i][j][k]);					
			
	//Using fftw3 library for solving the inverse fourier transform//	
	p = fftw_plan_dft_3d(SIZE_OF_GRID, SIZE_OF_GRID, SIZE_OF_GRID, &in[0][0][0], &out[0][0][0], 
			     FFTW_BACKWARD, FFTW_PATIENT);

	fftw_execute(p);

	//Building the dust Grid
	for (i = 0; i < SIZE_OF_GRID; i++)
		for (j = 0; j < SIZE_OF_GRID; j++)
			for (k = 0; k < SIZE_OF_GRID; k++)
				dustGrid[i][j][k] = fabs(crealf((complex) in[i][j][k])); 

	//Converting the cubes into spheres
	cube2sphere(dustGrid);
	cube2sphere(NORM_GRID);

	//Getting a factor A for normalization purposes and also calculating the total mass of normalization grid//
	sum = 0; norm_sum = 0;
	for (i = 0; i < SIZE_OF_GRID; i++)
	{
		for (j = 0; j < SIZE_OF_GRID; j++)
		{
			for (k = 0; k < SIZE_OF_GRID; k++)
			{
				sum += fabs(crealf((complex) in[i][j][k]));	
				norm_sum += NORM_GRID[i][j][k];
			}
		}
	}

	//Normalization factor
	A = norm_sum / sum;

	//Normalizing the dust grid
	for (i = 0; i < SIZE_OF_GRID; i++)
		for (j = 0; j < SIZE_OF_GRID; j++)
			for (k = 0; k < SIZE_OF_GRID; k++)
				dustGrid[i][j][k] = A * dustGrid[i][j][k];

	//fftw cleanup
	fftw_destroy_plan(p);
	fftw_cleanup();
}

//Function that corrects the cube into a spheere//
void cube2sphere(float dustGrid[SIZE_OF_GRID][SIZE_OF_GRID][SIZE_OF_GRID])
{
	int i, j, k;
	float center[3], sphereC[] = { (float)SIZE_OF_GRID/2 + 0.5, (float)SIZE_OF_GRID/2 + 0.5, (float) SIZE_OF_GRID/2 +0.5 }, radius = (float) SIZE_OF_GRID;

	for (i = 0; i < SIZE_OF_GRID; i++)
	{
		for (j = 0; j < SIZE_OF_GRID; j++)
		{
			for (k = 0; k < SIZE_OF_GRID; k++)
			{
				center[0] = (float) i + 0.5; center[1] = (float) j + 0.5; center[2] = (float) k + 0.5;

				if (distanceBetween(center, sphereC) > radius)
					dustGrid[i][j][k] = 0.0;
				
			}
		}
	}
}

//Function for fourier shift//
void fast_fourier_shift(float GRID[SIZE_OF_GRID][SIZE_OF_GRID][SIZE_OF_GRID])
{
	int i, j, k;
	float F_GRID[SIZE_OF_GRID][SIZE_OF_GRID][SIZE_OF_GRID];

	//Fourier Shift for K//
	for (i = 0; i < SIZE_OF_GRID; i++)
	{
		for (j = 0; j < SIZE_OF_GRID; j++)
		{
			for (k = 0; k < SIZE_OF_GRID; k++)
			{	
				if (i < SIZE_OF_GRID / 2)
					F_GRID[i][j][k] = GRID[i + SIZE_OF_GRID/2][j][k];
				else
					F_GRID[i][j][k] = GRID[i - SIZE_OF_GRID/2][j][k];
			}
		}
	}

	//Fourier Shift for J//
	for (i = 0; i < SIZE_OF_GRID; i++)
	{
		for (j = 0; j < SIZE_OF_GRID; j++)
		{
			for (k = 0; k < SIZE_OF_GRID; k++)
			{	
				if (j < SIZE_OF_GRID / 2)
					GRID[i][j][k] = F_GRID[i][j + SIZE_OF_GRID/2][k];
				else
					GRID[i][j][k] = F_GRID[i][j - SIZE_OF_GRID/2][k];
			}
		}
	}

	//Fourier Shift for I//
	for (i = 0; i < SIZE_OF_GRID; i++)
	{
		for (j = 0; j < SIZE_OF_GRID; j++)
		{
			for (k = 0; k < SIZE_OF_GRID; k++)
			{	
				if (k < SIZE_OF_GRID / 2)
					F_GRID[i][j][k] = GRID[i][j][k + SIZE_OF_GRID/2];
				else
					F_GRID[i][j][k] = GRID[i][j][k - SIZE_OF_GRID/2];
			}
		}
	}
	
	for (i = 0; i < SIZE_OF_GRID; i++)
		for (j = 0; j < SIZE_OF_GRID; j++)
			for (k = 0; k < SIZE_OF_GRID; k++)
				GRID[i][j][k] = F_GRID[i][j][k];
}

//Inverse Fast Fourier shift (not necessary for this code)//
void inverse_fast_fourier_shift(float GRID[SIZE_OF_GRID][SIZE_OF_GRID][SIZE_OF_GRID])
{
	int i, j, k;
	float F_GRID[SIZE_OF_GRID][SIZE_OF_GRID][SIZE_OF_GRID];

	//Fourier Shift for K//
	for (i = 0; i < SIZE_OF_GRID; i++)
	{
		for (j = 0; j < SIZE_OF_GRID; j++)
		{
			for (k = 0; k < SIZE_OF_GRID; k++)
			{	
				if (i < SIZE_OF_GRID / 2)
					F_GRID[i + SIZE_OF_GRID/2][j][k] = GRID[i][j][k];
				else
					F_GRID[i - SIZE_OF_GRID/2][j][k] = GRID[i][j][k];
			}
		}
	}

	//Fourier Shift for J//
	for (i = 0; i < SIZE_OF_GRID; i++)
	{
		for (j = 0; j < SIZE_OF_GRID; j++)
		{
			for (k = 0; k < SIZE_OF_GRID; k++)
			{	
				if (j < SIZE_OF_GRID / 2)
					GRID[i][j + SIZE_OF_GRID/2][k] = F_GRID[i][j][k];
				else
					GRID[i][j - SIZE_OF_GRID/2][k] = F_GRID[i][j][k];
			}
		}
	}

	//Fourier Shift for I//
	for (i = 0; i < SIZE_OF_GRID; i++)
	{
		for (j = 0; j < SIZE_OF_GRID; j++)
		{
			for (k = 0; k < SIZE_OF_GRID; k++)
			{	
				if (k < SIZE_OF_GRID / 2)
					F_GRID[i][j][k + SIZE_OF_GRID/2] = GRID[i][j][k];
				else
					F_GRID[i][j][k - SIZE_OF_GRID/2] = GRID[i][j][k];
			}
		}
	}

	//Saving the Grid//
	for (i = 0; i < SIZE_OF_GRID; i++)
		for (j = 0; j < SIZE_OF_GRID; j++)
			for (k = 0; k < SIZE_OF_GRID; k++)
				GRID[i][j][k] = F_GRID[i][j][k];
}

//Function for viewFileSituations
void photon_run_viewFile(float result[SIZE_OF_GRID][SIZE_OF_GRID], float dustGrid[SIZE_OF_GRID][SIZE_OF_GRID][SIZE_OF_GRID])
{
	float dir[3], norm[3];
	struct Photon phot;
	int i, j;

	for (i = 0; i < PHOTON_PACKS; i++)
	{ 
		if (i % (PHOTON_PACKS / 20) == 0)
			printf("Core %d Progress: %d %c \n", omp_get_thread_num(),  ((i*100)+100)/PHOTON_PACKS, '%');	

		//Running the route for one photon
		phot = photon_run(dustGrid);

		//Calculating if the photon hits the grid
		if (phot.pos[0] > 32)
		{
			if (phot.pos[1] > 0 && phot.pos[1] < 32 && phot.pos[2] > 0 && phot.pos[2] < 32)
			{
				photDirFloat(phot, dir); 
				norm[0] = 1; norm[1] = 0; norm[2] = 0;

				if (OBS_DOTP_LIMIT > 1 - dir[0] * norm[0])
					result[((int) phot.pos[1])][((int) phot.pos[2])] += phot.intensity;
			} 
		} 
	}
}

//Function of the route of one photon
struct Photon photon_run(float dustGrid[SIZE_OF_GRID][SIZE_OF_GRID][SIZE_OF_GRID])
{
	struct Photon phot;
	float dist, tau0, cellDen;
	int OutOfSystem = 0, test = 0;
	init_phot(&phot);

	while(OutOfSystem == 0)
	{
		//Generating tau0, the distance the photon travels before scattering
		tau0 = -log(randFloat());		

		while (tau0 > 0)
		{
			dist = distanceToNextCell(phot);
			cellDen = cellDensity(phot, dustGrid);
			
			if (cellDen < 0)
			{
				//testing for bugs
				if (test == 100000)
					printf("%f - %f\n\n",cellDen, dist);

				OutOfSystem = 1;
				break;
			}

			if (dist > tau0 / (DUST_KSCA*cellDen))
			{
				phot.pos[0] += cos(phot.dir[0]) * tau0 / (DUST_KSCA*cellDen);		
				phot.pos[1] += phot.dir[1] * tau0 / (DUST_KSCA*cellDen); 
				phot.pos[2] += sin(phot.dir[0]) * tau0 / (DUST_KSCA*cellDen); 

				tau0 = 0;
			}
			else
			{
				phot.pos[0] += cos(phot.dir[0]) * dist; 
				phot.pos[1] += phot.dir[1] * dist; 
				phot.pos[2] += sin(phot.dir[0]) * dist; 

				tau0 -= dist * DUST_KSCA * cellDen;
			}
		}

		//Testing for bugs
		if (test > 100000 && test < 100004)
		{
			print_photon(phot);
			printf("tau: %f\n", tau0);
		}

		OutOfSystem = isPhotonOutOfSystem(phot);
		if (OutOfSystem == 0) 
			scatter_photon(&phot);
		
		test++;
	}

	return phot;
}

//Photon info//
void print_photon(struct Photon phot)
{
	printf("Photon Position: %.2f, %.2f, %.2f\n", phot.pos[0], phot.pos[1], phot.pos[2]);
	printf("Photon Direction: %.2f, %.2f\n", phot.dir[0], phot.dir[1]);
	printf("Photon Intensity: %.2f\n", phot.intensity);
}

//Initialise photon//
void init_phot(struct Photon *phot)
{
	if (RANDOM_PHOTON_DIR == 0) 
	{
		(*phot).pos[0] = PHOTON_START_X; (*phot).pos[1] = PHOTON_START_Y; (*phot).pos[2] = PHOTON_START_Z; 
	}
	else
	{
		(*phot).pos[0] = PHOTON_START_X; (*phot).pos[1] = randFloat()*32; (*phot).pos[2] = randFloat()*32; 
	}

	(*phot).dir[0] = 0; (*phot).dir[1] = 0;
	(*phot).intensity = 1;
}

//Scattering//
void scatter_photon(struct Photon *phot)
{
	(*phot).dir[0] = 2 * PI * randFloat(); 
	(*phot).dir[1] = 2 *  randFloat() - 1; 
	(*phot).intensity *= DUST_ALBEDO;
}

//Function to determine distance between two points//
float distanceBetween(float vec1[3], float vec2[3])
{
	return sqrt(pow(vec2[0] - vec1[0], 2) + pow(vec1[1] - vec2[1], 2) + pow(vec1[2] - vec2[2], 2));
}

//How much distance there is to the next wall//
float distanceToNextCell(struct Photon phot)
{
	float rayPos[] = {phot.pos[0], phot.pos[1], phot.pos[2]}, distMod = 1.0;
	int wallReached = 0, dir = 1, i, rayCheck; 

	while(wallReached == 0)
	{
		rayPos[0] += dir*distMod*cos(phot.dir[0]);
		rayPos[1] += dir*distMod*phot.dir[1];
		rayPos[2] += dir*distMod*sin(phot.dir[0]);
	
		rayCheck = isRayInTheSameCell(phot, rayPos);

		if (rayCheck > 0 && dir == 1)
		{	
			dir = -1;
			distMod = distMod * 0.1;
		}
		else if (rayCheck == 0 && dir == -1)
		{
			dir = 1;
			distMod = distMod * 0.1;
		}
	        
		for (i = 0; i < 3; i++)
		{
			if (fabs(rayPos[i] - floor(rayPos[i])) < 0.001)  
				wallReached = 1; 

			if (fabs(floor(rayPos[i]) - rayPos[i]) < 0.001)  
				wallReached = 1;
		}
	}

	//Term for fixing the cases on the wrong side of 
	rayPos[0] += 0.001*cos(phot.dir[0]);
	rayPos[1] += 0.001*phot.dir[1];
	rayPos[2] += 0.001*sin(phot.dir[0]);
	
	return distanceBetween(phot.pos, rayPos);
}

//Function for checking if the ray is in the same cell than photon//
int isRayInTheSameCell(struct Photon phot, float rayPos[3])
{
	int i, check = 0;
	for (i = 0; i < 3; i++)
		if (floor(phot.pos[i]) != floor(rayPos[i]))
			check++;

	return check;
}

//Function for checking if the photon is out of the system
int isPhotonOutOfSystem(struct Photon phot)
{
	if (phot.pos[0] < 0.0 || phot.pos[0] >= 32.0)
		return 1;
	else if (phot.pos[1] < 0.0 || phot.pos[1] >= 32.0)
		return 1;
	else if (phot.pos[2] < 0.0 || phot.pos[2] >= 32.0)
		return 1;
	else
		return 0;
}

//Function for checking the cell density
float cellDensity(struct Photon phot, float dustGrid[SIZE_OF_GRID][SIZE_OF_GRID][SIZE_OF_GRID])
{
	int i, j, k;
	i = (int) floor(phot.pos[0]); j = (int) floor(phot.pos[1]); k = (int) floor(phot.pos[2]);

	if (i < 0 || i > 31 || j < 0 || j > 31 || k < 0 || k > 31)
		return -9.0;
	else
		return dustGrid[i][j][k]; 
}

float randFloat()
{
	return ((float) (rand() % 3277)) / 3277.0;
}

float dotProduct(float vec1[3], float vec2[3])
{
	return vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2];
}

void photDirFloat(struct Photon phot, float dir[3])
{
	dir[0] = cos(phot.dir[0]) * sin(PI / 2 - phot.dir[1]);
	dir[1] = sin(phot.dir[0]) * sin(PI / 2 - phot.dir[1]);
	dir[2] = cos(PI / 2 - phot.dir[1]);
}
