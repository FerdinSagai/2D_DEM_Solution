/*
Author :
DR. FERDIN SAGAI DON BOSCO
------------------------------------------------------
About Code:
This code features the DEM code in 2D

Version 0:

*/
#include<iostream>
#include<cmath>
//====================================================================================
#define PI acos(-1.0)
#define MaxBinPack 3
#define Layer 5
#define NWALLS 4

int sign(double x)
{
	int s;

	if (x == 0)
	{
		return (0);
	}
	else
	{
		s = x / fabs(x);
		return (s);
	}
}

int opposite_sign(double x)
{
	int s;

	if (x == 0)
	{
		return (0);
	}
	else
	{
		s = x / fabs(x);
		return (-s);
	}
}
double sqr(double x)
{
	return (x * x);
}

double min(double x, double y)
{
	if (x < y) return x;
	else return y;
}

double max(double x, double y)
{
	if (x > y) return x;
	else return y;
}
//====================================================================================
void initialize_variables(void);
void init_particle_position_read(void);

double Distance_between_Line_and_Point(int, double, double);

void add_particles(double, double, double, double, int, double);
void delete_particles(double, double, double, double, int);

void Nearest_Neighbour_Algorithm_Brute_Force(int);
void Nearest_Neighbour_Algorithm_L_Search(void);
void Combined_Algorithm(void);

double compute_collision_time(int, int);

void calc_force_particles(int, int);
void calculation_force_walls(int);

void calculation_force_particles(int, int);

void calculation_force_wall_left(int);
void calculation_force_wall_right(int);
void calculation_force_wall_bottom(int);

void calculation_force_wall(int, int);


void update_particle_info(int, double, double, double, double);

void save_particle_position(void);
void save_final_particle_position(void);

void DEM_void_fraction(void);
void Neighbours(int);
//=============================
const int MAX_NEIGHBOUR = 6;
const int MAXNUM = 5000;
const int IMAX = 50;
const int JMAX = 50;
/*------Determine the number of bins in the domain----------------------------------------*/
double d_max;
double BIN_S;

int Nx, Ny, Nz;				// Grid resolution for the DEM computations
int NUM;				// Number of DEM particles used in packing
double time_NNb;
double Lx, Ly;
double h, D_T, r;
int n_structures;
double CD;
double rho_g;
int nt;
int State[IMAX + 1][JMAX + 1];		// Status of each grid point. Solid or void space
double x[IMAX + 1];			// X co-ordinate of the grid point
double y[IMAX + 1];			// Y co-ordinate of the grid point
double dx[IMAX + 1];			// Cell size in the X direction
double dy[IMAX + 1];			// Cell size in the Y direction
double vfrac[IMAX + 1][JMAX + 1];		// Value of void fraction at each grid point
double sfrac[IMAX + 1][JMAX + 1];		// Value of solid fraction at each grid point (1 - void fraction))


double Surf_Xi[NWALLS], Surf_Yi[NWALLS];
double Surf_Xf[NWALLS], Surf_Yf[NWALLS];
double Surf_dx[NWALLS], Surf_dy[NWALLS];
double Surf_Angle[NWALLS];
double Surf_Normal_Angle[NWALLS];
double Surf_n[NWALLS][3];
double Surf_len[NWALLS];
//---------------------------------------------------------------------------------
// GAS FLOW COMPUTATION
double u[IMAX + 1][JMAX + 1];			// Corrected Gas velocity in the X direction at each grid point
double v[IMAX + 1][JMAX + 1];			// Corrected Gas velocity in the Y direction at each grid point
double C_x[10], C_y[10], H_block[10], L_block[10];
double gravity = 9.81;
// PROPERTIES OF PARTICLE CLASS
double rho_s;								// Density of the particle
double d[MAXNUM];							// Diameter of the particle
double mass[MAXNUM];						// Mass of the particle in kg
double MMOI[MAXNUM];						// Moment of inertia 
double xp[MAXNUM], yp[MAXNUM];				// Position for the Bed Particles
double zp[MAXNUM];
double up[MAXNUM], vp[MAXNUM];				// Velocity of Bed Particle
double acc_x[MAXNUM], acc_y[MAXNUM];		// Linear Acceleration of Bed Particle
double acc_ang[MAXNUM];						// Angular Acceleration of Bed Particle
double omega[MAXNUM];						// Angular Velocity of the bed particles
double coeff_friction;						// Coefficient of friction
double coeff_restitution;					// Coefficient of restitution
double Kn, Kt;								// Normal and Tangential spring constants of a single particle

double angle[MAXNUM][MAX_NEIGHBOUR];	// Angle between colliding pairs
int CNUM[MAXNUM][MAX_NEIGHBOUR];	// Neighbour list for a single particle

double Old_Delta_N_Particles[MAXNUM][MAXNUM];
double Old_Delta_T_Particles[MAXNUM][MAXNUM];
double Old_Delta_N_Walls[MAXNUM][NWALLS];
double Old_Delta_T_Walls[MAXNUM][NWALLS];

//double Ct;				// Damping coefficient in the Normal and Tangential direction for a single particle (Magnitude only)
double KnWall, KtWall;			// Normal and Tangential spring constants of the wall
double CnWall, CtWall;			// Damping coefficient in the Normal and Tangential direction for a wall
double mark;				// Tolerance margin for collision detection
double STATICtime;			// Standard time step
double sim_time;			// Total simulation time for DEM
int nstep;				// Number of timesteps in the DEM simulation
double rp; 				// Radius of particle.
int packing_option; 			// Packing Options
double limitH;				// Particle Generation height

double x1_t, x2_t, x3_t, x4_t;
double y1_t, y2_t, y3_t, y4_t;
double ke_system;

double Force_on_Par_X[IMAX][JMAX];
double Force_on_Par_Y[IMAX][JMAX];
double time_DEMb, time_DEMa, cpu_time_DEM;
double time_NNa;
// =============================================================================
int ini, fin;
int Particles_Added, Particles_Removed;
double generation_region;
int pack_actual, pack;
int update_count;
double insertion_min_X;
double insertion_max_X;
double insertion_min_Y;
double insertion_max_Y;

double deletion_min_X;
double deletion_max_X;
double deletion_min_Y;
double deletion_max_Y;
double discharge_rate;
double dp;
double time_step;
double cpu_time_NN;
double 	max_deformation;
double	max_gasforce;
int Particle_position_count;
int nthreads = 1;
int insertion_rate;
double phi_s;
int mode;
//=============================================================================================
#include "Inter_Particle_Collisions.h"
#include "Left_Wall.h"
#include "Right_Wall.h"
#include "Bottom_Wall.h"
#include "Wall_Particle_Collisions.h"
//=============================================================================================
int main()
{
	int i, j;
	/*------Noting the cpu start time --------------------------------------------------------*/
	time_DEMb = clock();
	/*------Variable Initializations ---------------------------------------------------------*/
	initialize_variables();

	pack_actual = NUM;
	NUM = 0;
	generation_region = (insertion_max_X - insertion_min_X) * (insertion_max_Y - insertion_min_Y);
	pack = 0.5 * generation_region / (dp * dp);
	pack = pack_actual;
	//if (pack > pack_actual) pack = pack_actual;

	nstep = (int)ceil(sim_time / STATICtime);
	nstep = 2 * nstep;
	ke_system = 0.00;
	//------Time Marching for the particles-----------------------------------------------
	for (nt = 1; nt <= nstep; nt++)
	{
		if (nt == 1)
		{
			printf("nt=%d, time=%f\n", nt, nt * STATICtime);
			cpu_time_NN = 0.00;
			max_deformation = 0.00;
			max_gasforce = 0.00;
			Particles_Added = 0;
			Particles_Removed = 0;
			ke_system = 0.00;
			
			NUM = pack_actual;
			FILE* fp1;
			fp1 = fopen("2.0.Particle_Positions_Init.inp", "r");
			for (i = 0; i < NUM; i++)
			{
				fscanf(fp1, "%lf\t%lf\t%lf\n", &xp[i], &yp[i], &d[i]);
			}
			fclose(fp1);
		}
		// Add particle at intervals
		if ((nt == 1) || (nt % 10000 == 0))
		{
			if (NUM < pack_actual)
			{
				Particles_Added = NUM;
				add_particles(insertion_min_X, insertion_max_X, insertion_min_Y, insertion_max_Y, pack, dp);
				if ((pack_actual - NUM) < pack)
				{
					pack = pack_actual - NUM;
				}
				Particles_Added = NUM - Particles_Added;
			}
		}
		
		// Delete particle at intervals
		if (nt % 10000 == 0)
		{
			Particles_Removed = NUM;
			delete_particles(deletion_min_X, deletion_max_X, deletion_min_Y, deletion_max_Y, discharge_rate);
			Particles_Removed = Particles_Removed - NUM;
		}

		// Display state at intervals
		if (nt % 10000 == 0)
		{
			printf("Update Count = %d :: nt=%d out of nstep=%d\n", update_count, nt, nstep);
			printf("\t Number of particles....................................%d \n", NUM);
			printf("\t Particles Added........................................%d \n", Particles_Added);
			printf("\t Particles Removed......................................%d \n", Particles_Removed);
			printf("\t time for DEM Process ..................................%lf seconds\n", cpu_time_NN);
			printf("\t Kinetic Energy of the system ..........................%lf Joules \n", ke_system);
			printf("\t Maximum deformation of a particle......................%lf m \n", max_deformation);
			printf("\t Maximum gas force of a particle........................%lf N \n", max_gasforce);
			save_particle_position();
		}
		if (nt % 10000 == 0)
		{
			if (nt > 100)
			{
				if (ke_system < 0.050)
				{
					//break;
				}
			}
			cpu_time_NN = 0.00;
			max_deformation = 0.00;
			max_gasforce = 0.00;
			Particles_Added = 0;
			Particles_Removed = 0;
			ke_system = 0.00;
		}
		//Nearest_Neighbour_Algorithm_Brute_Force();
		Nearest_Neighbour_Algorithm_L_Search();
		Combined_Algorithm();
	}
	//DEM_void_fraction();
	save_particle_position();
	Particle_position_count++;
	//------Saving Final position of Particles -------------------------------------------
	save_final_particle_position();
	//------Noting the cpu end time ------------------------------------------------------
	time_DEMa = clock();
	//------Noting the total cpu time for this computation--------------------------------
	cpu_time_DEM = cpu_time_DEM + ((double)(time_DEMa - time_DEMb)) / (nthreads * CLOCKS_PER_SEC);
	printf("\t Total time for the system ..........................%lf seconds \n", cpu_time_DEM);
	printf("\t FUTURE TEST PURPOSE ....... Running with Gas\n");
	return 0;
}
//==========================================================================
void initialize_variables()
{
	Lx = 0.05;
	Ly = 0.30;

	Nx = 75;
	Ny = 150;
	mark = 1.02;
	NUM = 100;
	dp = 4.00e-3;
	insertion_min_X = 0.00;
	insertion_min_Y = 0.25;
	insertion_max_X = 0.10;
	insertion_max_Y = 0.30;
	insertion_rate = 100;

	deletion_min_X = 0.00;
	deletion_min_Y = 0.00;
	deletion_max_X = 0.38;
	deletion_max_Y = 0.05;
	discharge_rate = 0;

	rho_g = 1.178;
	phi_s = 1.00;
	rho_s = 2500;
	coeff_friction = 0.3;
	coeff_restitution = 0.8;
	Kn = 1000;
	sim_time = 2.0;

	h = 5.0e-2;
	D_T = 5.00e-3;
	r = 0.00e-2;
	if ((Layer == 4) || (Layer == 5)) r = 1.00e-2;

// Left Surface
	Surf_Xi[0] = 0.00;
	Surf_Yi[0] = 0.00;
	Surf_Xf[0] = 0.00;
	Surf_Yf[0] = Ly;
// Bottom Surface
	Surf_Xi[1] = Lx;
	Surf_Yi[1] = 0.00;
	Surf_Xf[1] = 0.00;
	Surf_Yf[1] = 0.00;
// Right Surface
	Surf_Xi[2] = Lx;
	Surf_Yi[2] = Ly;
	Surf_Xf[2] = Lx;
	Surf_Yf[2] = 0.00;
	
	// Inclined plane Surface
	Surf_Xi[3] = 0.01;
	Surf_Yi[3] = 0.00;
	Surf_Xf[3] = 0.00;
	Surf_Yf[3] = 0.01;
// Surface Properties
	for (int i = 0; i < NWALLS; i++)
	{
		Surf_dx[i]		= Surf_Xf[i] - Surf_Xi[i];
		Surf_dy[i]		= Surf_Yf[i] - Surf_Yi[i];
		
		if (Surf_dx[i] == 0.00 )
		{
			if (Surf_dy[i] > 0.00 )		Surf_Angle[i]	= (PI/2); //CHeck
			else if (Surf_dy[i] < 0.00 )Surf_Angle[i]	= 3*PI/2; //CHeck
		}
		else if (Surf_dy[i] == 0.00 )
		{
			if (Surf_dx[i] > 0.00 )		Surf_Angle[i]	= 2*PI;
			else if (Surf_dx[i] < 0.00 )Surf_Angle[i]	= PI;
		}		
		else
		{
			Surf_Angle[i]	= atan(Surf_dy[i]/Surf_dx[i]);
		}
		Surf_Normal_Angle[i]= Surf_Angle[i] - (PI/2);
		if( Surf_Normal_Angle[i]<0.00) Surf_Normal_Angle[i]= Surf_Normal_Angle[i] + 2*PI;
		
		Surf_n[i][0] 		= cos(Surf_Normal_Angle[i]);
		Surf_n[i][1] 		= sin(Surf_Normal_Angle[i]);		
		Surf_n[i][2] 		= 0.00;		
		Surf_len[i] = sqrt( sqr(Surf_dx[i]) + sqr (Surf_dy[i]) );
		
		printf("%d %lf %lf \n", i, Surf_Normal_Angle[i]*180/PI, Surf_Angle[i]*180/PI);
	}

	// Initializing the key variables------------------------------------------------------------------------
	for (int i = 0; i < NUM; i++)
	{
		xp[i] = 0.00;
		yp[i] = 0.00;
		zp[i] = dp * 0.50;
		up[i] = 0.00;
		vp[i] = 0.00;
		acc_x[i] = 0.00;
		acc_y[i] = 0.00;
		acc_ang[i] = 0.00;
		omega[i] = 0.00;

		d[i] = dp;

		mass[i] = rho_s * (4.0 / 3.0) * PI * pow(d[i] / 2.0, 3.0);
		MMOI[i] = (2.0 / 5.0) * mass[i] * pow(d[i] / 2.0, 2.0);

		for (int j = 0; j < NUM; j++)
		{
			Old_Delta_N_Particles[i][j] = 0.00;
			Old_Delta_T_Particles[i][j] = 0.00;
		}
		
		for (int j = 0; j < NWALLS; j++)
		{
			Old_Delta_N_Walls[i][j] = 0.00;
			Old_Delta_T_Walls[i][j] = 0.00;
		}
	}

	double mass_min = 1000.00;
	double mass_max = 1000.0;
	d_max = 0.00;
	for (int i = 0; i < NUM; i++)
	{
		mass_min = min(mass_min, mass[i]);
		mass_max = max(mass_max, mass[i]);
		d_max = max(d_max, d[i]);
	}

	Kt = Kn;
	KnWall = Kn;
	KtWall = Kt;

	time_step = 1.0e-6;
	if (time_step > 0.1 * sqrt(mass_min / Kn))
	{
		time_step = 0.1 * sqrt(mass_min / Kn);
	}
	STATICtime = time_step;

	Particles_Added = 0;
	Particles_Removed = 0;
	ke_system = 0.00;
}
//==========================================================================
void init_particle_position_read() {
	int i, j;
	FILE* fp1;
	fp1 = fopen("1.5.Particle_Positions_Post.inp", "r");
	for (i = 0; i < NUM; i++)
	{
		fscanf(fp1, "%lf\t%lf\t%lf\n", &xp[i], &yp[i], &d[i]);
	}
	fclose(fp1);

	for (i = 0; i < NUM; i++)
	{
		d[i] = dp;
		printf("%d %lf %lf %lf\n", i, xp[i], yp[i], d[i]);
	}

}
//==========================================================================
double Distance_between_Line_and_Point(int Surf_id, double Px, double Py)
{
	double Parameter_1x  = Px - Surf_Xi[Surf_id];
	double Parameter_1y  = Py - Surf_Yi[Surf_id];
	
	double VdotPA 	= Parameter_1x * Surf_dx[Surf_id] + Parameter_1y * Surf_dy[Surf_id];
	double Vmag	  	= Surf_len[Surf_id];
	double t1		= VdotPA/(Vmag*Vmag);
	double X_int 	= Surf_dx[Surf_id] * t1 + Surf_Xi[Surf_id];
	double Y_int 	= Surf_dy[Surf_id] * t1 + Surf_Yi[Surf_id];
	double distance = sqrt( sqr(Px - X_int) + sqr(Py - Y_int) );
	
	/*
	// FOR DEBUGGING PURPOSES: DO NOT DELETE
	printf( "\t Surface :: %d \n", Surf_id);
	printf( "\t Surface :: %lf %lf \n", Surf_Xi[Surf_id], Surf_Yi[Surf_id]);
	printf( "\t Surface :: %lf %lf \n", Surf_Xf[Surf_id], Surf_Yf[Surf_id]);
	printf( "\t Surface :: %lf \n", Vmag);
	printf( "\t Direction Vec :: %lf %lf\n", Surf_dx[Surf_id], Surf_dy[Surf_id]);
	printf( "\t Particle :: %lf %lf\n", Px, Py);
	printf( "\t Numerator :: %lf\n", VdotPA);
	printf( "\t t  :: %lf\n", t1);
	printf( "\t Intersection :: %lf %lf\n", X_int, Y_int);
	printf( "Distance :: %lf\n", distance);
	*/
	return distance;
}
//==========================================================================
void add_particles(double minX, double maxX, double minY, double maxY, int N, double dia)
{
	double range_X, range_Y;
	double mass_min = 10000;
	range_X = (maxX - minX);
	range_Y = (maxY - minY);

	srand(time(0));

	for (int i = NUM; i < (NUM + N); i++)
	{
		xp[i] = minX + range_X * rand() / RAND_MAX;
		yp[i] = minY + range_Y * rand() / RAND_MAX;
		up[i] = double(rand()) / double(RAND_MAX);
		vp[i] = double(rand()) / double(RAND_MAX);
		acc_x[i] = 0.00;
		acc_y[i] = 0.00;
		acc_ang[i] = 0.00;
		omega[i] = 0.00;
		d[i] = dia;

		mass[i] = rho_s * (4.0 / 3.0) * PI * pow(d[i] / 2.0, 3.0);
		MMOI[i] = (2.0 / 5.0) * mass[i] * pow(d[i] / 2.0, 2.0);

		mass_min = min(mass_min, mass[i]);

		for (int j = 0; j < NUM; j++)
		{
			Old_Delta_N_Particles[i][j] = 0.00;
			Old_Delta_T_Particles[i][j] = 0.00;
		}
		
		for (int j = 0; j < NWALLS; j++)
		{
			Old_Delta_N_Walls[i][j] = 0.00;
			Old_Delta_T_Walls[i][j] = 0.00;
		}
		printf("Adding Particle :: %d %le %le %lf %lf %lf %lf\n", i, mass[i], MMOI[i], xp[i], yp[i], up[i], vp[i]);
	}
	NUM = NUM + N;

	if (time_step > 0.1 * sqrt(mass_min / Kn))
	{
		time_step = 0.1 * sqrt(mass_min / Kn);
	}
	STATICtime = time_step;
}
//==========================================================================
void delete_particles(double minX, double maxX, double minY, double maxY, int N)
{
	int deleted_N;
	double range_X, range_Y;

	deleted_N = 0;
	range_X = (maxX - minX);
	range_Y = (maxY - minY);

	if (N > 0)
	{
		srand(time(0));

		for (int i = 0; i < NUM; i++)
		{
			if ((xp[i] < maxX) && (xp[i] > minX) && (yp[i] < maxY) && (yp[i] > minY))
			{
				//double r = float(rand())/float(RAND_MAX);
				//if( r  > 0.50 )	
				{
					xp[i] = xp[NUM - 1];
					yp[i] = yp[NUM - 1];
					up[i] = up[NUM - 1];
					vp[i] = vp[NUM - 1];
					omega[i] = omega[NUM - 1];
					acc_x[i] = acc_x[NUM - 1];
					acc_y[i] = acc_y[NUM - 1];
					acc_ang[i] = acc_ang[NUM - 1];

					NUM = NUM - 1;
					deleted_N++;
					if (deleted_N >= N) break;
				}
			}
		}
	}
}
//==========================================================================
int SEARCH_DOMAIN[MAXNUM][5];
int	INV_P_CELL[6000][MaxBinPack];

void Nearest_Neighbour_Algorithm_L_Search()
{
	//------Noting the cpu start time --------------------------------------------------------
	time_NNb = clock();
	//------Local variables declaration--------------------------------------------------------
	int i, j, k;
	int nbin_x, nbin_y;
	int bx, by;
	double BIN_S;
	double rv_ij;
	int index;
	int N_Neigh[MAXNUM];
	int NeighCount;
	int TBINS;
	/*------Determine the number of bins in the domain----------------------------------------*/
	BIN_S = (mark * d_max);
	nbin_x = ceil(Lx / BIN_S);
	nbin_y = ceil((3.0 * Ly) / BIN_S);
	TBINS = (nbin_x + 1) * (nbin_y + 1);	
	//------Initialize array for particles in bin k------------------------------------------
	for (k = 0; k < TBINS; k++)
	{
		for (int idx = 0; idx < MaxBinPack; idx++)
		{
			INV_P_CELL[k][idx] = -1;
		}
	}
	//------Determination of search domains--------------------------------------------------
	for (i = 0; i < NUM; i++)
	{
		//------Initialize search domain for particle i----------------------------------
		SEARCH_DOMAIN[i][0] = -1;
		SEARCH_DOMAIN[i][1] = -1;
		SEARCH_DOMAIN[i][2] = -1;
		SEARCH_DOMAIN[i][3] = -1;
		SEARCH_DOMAIN[i][4] = -1;
		//------Determine the cell indices for particle i--------------------------------
		bx = ((xp[i]) / (BIN_S));
		by = ((yp[i]) / (BIN_S));
		if (bx < 0)
		{
			xp[i] = 0.5 * d[i];
			bx = ((xp[i]) / (BIN_S));
		}
		if (by < 0)
		{
			yp[i] = 0.5 * d[i];
			by = ((yp[i]) / (BIN_S));
		}
		if (bx > nbin_x)
		{
			xp[i] = Lx - (0.5 * d[i]);
			bx = ((xp[i]) / (BIN_S));
		}
		if (by > nbin_y)
		{
			yp[i] = Ly - (0.5 * d[i]);
			by = ((yp[i]) / (BIN_S));
		}
		int MC = ((by)*nbin_x) + bx;

		//------Populate the array for particles in bin k--------------------------------
		index = 0;
		for (k = 0; k < MaxBinPack; k++)
		{
			if (INV_P_CELL[MC][k] >= 0)
			{
				index = index + 1;
			}
			else
			{
				break;
			}
		}
		INV_P_CELL[MC][index] = i;

		//------Populate the search domain for particle i--------------------------------
		SEARCH_DOMAIN[i][0] = MC;
		SEARCH_DOMAIN[i][1] = (by)*nbin_x + (bx - 1); 	// WEST (SAME PLANE)
		SEARCH_DOMAIN[i][2] = (by + 1) * nbin_x + (bx - 1); 	// NORTH-WEST(SAME PLANE)
		SEARCH_DOMAIN[i][3] = (by + 1) * nbin_x + bx; 	// NORTH(SAME PLANE)
		SEARCH_DOMAIN[i][4] = (by + 1) * nbin_x + (bx + 1);  	// NORTH-EAST(SAME PLANE)	
		//------Treatment for particles near bounding walls------------------------------
		// Treatment of Western wall
		if (bx == 0)
		{
			SEARCH_DOMAIN[i][1] = -1; 			// WEST (SAME PLANE)
			SEARCH_DOMAIN[i][2] = -1; 			// NORTH-WEST(SAME PLANE)
		}

		// Treatment of Eastern wall		
		if (bx == nbin_x - 1)
		{
			SEARCH_DOMAIN[i][4] = -1; 			// NORTH-EAST(SAME PLANE)	
		}

		// Treatment of Northern wall	
		if (by == nbin_y - 1)
		{
			SEARCH_DOMAIN[i][2] = -1; 			// NORTH-WEST(SAME PLANE)
			SEARCH_DOMAIN[i][3] = -1; 			// NORTH(SAME PLANE)
			SEARCH_DOMAIN[i][4] = -1; 			// NORTH-EAST(SAME PLANE)
		}
	}
}
void Combined_Algorithm()
{
	//------Determination of neighbours for particle i by searching its search domain-------
	for (int i = 0; i < NUM; i++)
	{
		for (int N = 0; N < 5; N++)
		{
			int MC = SEARCH_DOMAIN[i][N];
			if (MC >= 0)
			{
				for (int idx = 0; idx < MaxBinPack; idx++)
				{
					int j = INV_P_CELL[MC][idx];
					if (j >= 0 && j != i)
					{
						double d_avg = 0.50 * (d[i] + d[j]);
						double R = (sqr(xp[i] - xp[j]) + sqr(yp[i] - yp[j]));
						if (R <= sqr(mark * d_avg))
						{
							calculation_force_particles(j, i);
							
							if (SEARCH_DOMAIN[i][N] != SEARCH_DOMAIN[j][N])
							{
								calculation_force_particles(i, j);	
							}
						}
					}
				}
			}
		}
		
		
		for (int k = 0; k < 4; k++)
		{
			calculation_force_wall(i,k);
		}		
		//calculation_force_wall_left(i);
		//calculation_force_wall_right(i);
		//calculation_force_wall_bottom(i);

		//acc_x[i] = acc_x[i] + body_force_x[nx][ny];
		//acc_y[i] = acc_y[i] + body_force_y[nx][ny];

		acc_y[i] = acc_y[i] - gravity;

		xp[i] = xp[i] + up[i] * STATICtime + 0.5 * acc_x[i] * pow(STATICtime, 2.0);
		yp[i] = yp[i] + vp[i] * STATICtime + 0.5 * acc_y[i] * pow(STATICtime, 2.0);
		up[i] = up[i] + acc_x[i] * STATICtime;
		vp[i] = vp[i] + acc_y[i] * STATICtime;
		omega[i] = omega[i];// +acc_ang[i] * STATICtime;

		acc_x[i] = 0.00;
		acc_y[i] = 0.00;
		double prev_ke = ke_system;
		ke_system = ke_system + 0.5 * mass[i] * (up[i] * up[i] + vp[i] * vp[i]);
	}
	//------Noting the cpu end time ------------------------------------------------------
	time_NNa = clock();
	//------Noting the total cpu time for this computation--------------------------------
	cpu_time_NN = cpu_time_NN + ((double)(time_NNa - time_NNb)) / (nthreads * CLOCKS_PER_SEC);
}
//=============================================================================================================================================
//=============================================================================================================================================
//=============================================================================================================================================
void update_particle_info(int i, double ForceX, double ForceY, double Torque, double Contact_time)
{
	acc_x[i] = acc_x[i] + (ForceX / mass[i]);
	acc_y[i] = acc_y[i] + (ForceY / mass[i]);
	acc_ang[i] = (Torque / MMOI[i]);
}
//=============================================================================================================================================
//=============================================================================================================================================
//=============================================================================================================================================
double compute_collision_time(int i, int j)
{
	double contact_time, travel_time;
	double collision_time;
	double xij, yij;
	double xfij, yfij;
	double rij, rfij;
	double vij;
	double rvij;
	double vxij, vyij;
	double det;
	double coeff_a, coeff_b, coeff_c;
	double root1, root2, discr;

	double d_avg = 0.50 * (d[i] + d[j]);
	vxij = up[i] - up[j];							// Relative Velocity between particle i and j (X direction)
	vyij = vp[i] - vp[j];							// Relative Velocity between particle i and j (Y direction)
	xij = xp[i] - xp[j];							// Distance between particle i and j (X direction)
	yij = yp[i] - yp[j];							// Distance between particle i and j (Y direction)
	xfij = xij + STATICtime * vxij;						// Predicated distance between particle i and j (X direction)
	yfij = yij + STATICtime * vyij;						// Predicated distance between particle i and j (Y direction)

	coeff_a = vxij * vxij + vyij * vyij;	 					// Square of Relative velocity  (Coefficient a) 
	coeff_b = 2 * (xij * vxij + yij * vyij); 						// Approach vector (Coefficient b) 
	coeff_c = ((xij * xij + yij * yij) - (xij * xij + yij * yij)) - (d_avg * d_avg);	// Difference between initial and contact position (Coefficient c) 

	rvij = xij * vxij + yij * vyij; 							// Approach vector
	vij = pow(vxij, 2.0) + pow(vyij, 2.0); 					// Square of Relative velocity between particles i and j
	rij = pow(xij, 2.0) + pow(yij, 2.0); 						// Square of distance between particles i and j
	rfij = pow(xfij, 2.0) + pow(yfij, 2.0);					// Expected New position of particle j wrt particle i after STATICtime
	discr = pow(coeff_b, 2.0) - 4.0 * coeff_a * coeff_c;				// Discriminant (b^2- 4ac)

	// --------------
	if (rij > pow(d_avg, 2.0))
	{
		if (rfij < pow(d_avg, 2.0))
		{	// Option A
			if (discr >= 0.0)
			{
				root1 = (-coeff_b + sqrt(discr)) / (2.0 * coeff_a);
				root2 = (-coeff_b - sqrt(discr)) / (2.0 * coeff_a);
				if ((STATICtime - root1) > 0.0)
				{
					travel_time = root1;
				}
				else
				{
					travel_time = root2;
				}
				contact_time = STATICtime - travel_time;
			}
			else if (discr < 0.0)
			{
				contact_time = 0.0;					// Anamoly condition ?
				travel_time = STATICtime;
			}
		}
		else
		{	// Option B
			travel_time = STATICtime;
			contact_time = 0.0;
		}
	}
	else
	{
		if (rfij > pow(d_avg, 2.0))
		{	// Option D		
			if (discr >= 0.0)
			{
				root1 = (-coeff_b + sqrt(discr)) / (2.0 * coeff_a);
				root2 = (-coeff_b - sqrt(discr)) / (2.0 * coeff_a);
				if ((STATICtime - root1) > 0.0)
				{
					contact_time = root1;
				}
				else
				{
					contact_time = root2;
				}
				travel_time = STATICtime - contact_time;
			}
			else if (discr < 0.0)
			{								// Anamoly condition ?
				contact_time = STATICtime;
				travel_time = 0.00;
			}
		}
		else
		{	// Option C
			contact_time = STATICtime;
			travel_time = 0.00;
		}
	}

	//Handling anamolies due to round off errors
	if (contact_time < 0.0)
	{
		contact_time = 0.0;
	}

	if (contact_time > STATICtime)
	{
		contact_time = STATICtime;
	}
	return contact_time;
}
//=============================================================================================================================================
//=============================================================================================================================================
//=============================================================================================================================================
void save_particle_position()
{
	int i, j;
	char name[80];
	double xx, yy;

	FILE* gp1;
	gp1 = fopen("3.1.Transient_Particle_Positions.dat", "a");
	if (Particle_position_count == 1)
	{
		fprintf(gp1, "TITLE = Properties at transient states\n");
		fprintf(gp1, "VARIABLES=\"X\",\"Y\",\"Dp\"\n");
	}

	fprintf(gp1, "ZONE T = \"%d\"\n", nt);
	for (i = 0; i < NUM; i++)
	{
		fprintf(gp1, "%e %e %e\n", fabs(xp[i]), fabs(yp[i]), d[i]);
	}
	fclose(gp1);
}
//=============================================================================================================================================
//=============================================================================================================================================
//=============================================================================================================================================
void save_final_particle_position()
{
	int i;
	double xc, yc, rc;
	xc = 0.02;
	yc = h;

	if ((Layer == 4) || (Layer == 5))rc = 0.50 * xc;
	if ((Layer == 3) || (Layer == 6))rc = 0.40 * xc;
	if ((Layer == 2) || (Layer == 7))rc = 0.30 * xc;
	if ((Layer == 1) || (Layer == 8))rc = 0.20 * xc;

	FILE* fp3;
	fp3 = fopen("1.1.Particle_Positions.dat", "wt");
	fprintf(fp3, "TITLE = Properties at Final State\n");
	fprintf(fp3, "VARIABLES=\"X\",\"Y\",\"Z\",\"Dp\"\n");
	for (i = 0; i < NUM; i++)
	{
		if (sqr(xp[i] - xc) + sqr(yp[i] - yc) < sqr(rc))
		{

		}
		else
		{
			fprintf(fp3, "%lf %lf %lf %lf\n", xp[i], yp[i], zp[i], d[i]);
		}
	}
	fclose(fp3);
}
//=============================================================================================================================================
//=============================================================================================================================================
//=============================================================================================================================================
void DEM_void_fraction()
{
	int i, j, k;
	double phvo, delx, dely, s;

	for (i = 0; i < IMAX - 1; i++)
	{
		for (j = 0; j < JMAX - 1; j++)
		{
			vfrac[i][j] = 1.00;
			double s = 0;
			for (int k = 0; k < NUM; k++)
			{
				if (((((x[i] <= (xp[k] + (dp / 2))) && (x[i] > (xp[k] - (dp / 2))))
					|| ((x[i + 1] < (xp[k] + (dp / 2))) && (x[i + 1] >= (xp[k] - (dp / 2)))))
					&& (((y[j] <= (yp[k] + (dp / 2))) && (y[j] > (yp[k] - (dp / 2))))
						|| ((y[j + 1] < (yp[k] + (dp / 2))) && (y[j + 1] >= (yp[k] - (dp / 2))))))
					|| (((x[i + 1] >= (xp[k] + dp / 2)) && (x[i] <= (xp[k] - dp / 2)))
						&& (((y[j + 1] >= (yp[k] + dp / 2)) && (y[j] <= (yp[k] - dp / 2))))))
				{
					delx = (min(x[i + 1], (xp[k] + (dp / 2))) - max(x[i], (xp[k] - (dp / 2))));
					dely = (min(y[j + 1], (yp[k] + (dp / 2))) - max(y[j], (yp[k] - (dp / 2))));
					phvo = (PI * delx * dely / 4);
					s = s + phvo;
				}

				vfrac[i][j] = 1 - (s / (dx[i] * dy[j]));

				if (vfrac[i][j] <= 0.3)
				{
					vfrac[i][j] = 0.3;
				}
				else if (vfrac[i][j] == 1)
				{
					vfrac[i][j] = 0.95;
				}
				else
				{
					vfrac[i][j] = vfrac[i][j];
				}
			}
		}
	}

	FILE* gp;
	gp = fopen("3.0.DEM_Void_Fraction.dat", "a");
	if (Particle_position_count == 1)
	{
		fprintf(gp, "TITLE = States of the domain\n");
		fprintf(gp, "VARIABLES=\"X (m)\",\"Y (m)\",\"State\",\"Void Fraction\"\n");
	}
	fprintf(gp, "ZONE T = \"%d\", I=%d, J=%d F=POINT\n", Particle_position_count, IMAX, JMAX);

	double yy = 0.5 * dy[1];
	for (j = 1; j < JMAX; j++)
	{
		double xx = 0.5 * dx[1];
		for (i = 1; i < IMAX; i++)
		{
			fprintf(gp, "%e %e %d %e\n", xx, yy, State[i][j], vfrac[i][j]);
			// UPDATE X POSITION
			xx = xx + 0.5 * (dx[i] + dx[i + 1]);
		}
		// UPDATE Y POSITION
		yy = yy + 0.5 * (dy[j] + dy[j + 1]);
	}
	fclose(gp);
}
//=============================================================================================================================================
//=============================================================================================================================================
//=============================================================================================================================================
void Nearest_Neighbour_Algorithm_Brute_Force(int i)
{
	printf("The particle %d has coordinates %lf %lf\n", i, xp[i], yp[i]);
	int nnb = 0;
	for (int j = 0; j < NUM; j++)
	{
		if (i != j)
		{
			double rij = sqrt(sqr(xp[i] - xp[j]) + sqr(yp[i] - yp[j]));
			double dij = 0.5 * (d[i] + d[j]);
			if (rij <= (mark * dij))
			{
				if (nnb < MAX_NEIGHBOUR)
				{
					CNUM[i][nnb] = j;
					nnb++;
				}
				else
				{
					for (int k = 0; k < nnb; k++)
					{
						int ID_1 = CNUM[i][k];
						double rij_n = sqrt(sqr(xp[i] - xp[ID_1]) + sqr(yp[i] - yp[ID_1]));
						double rij_e = sqrt(sqr(xp[i] - xp[ID_1]) + sqr(yp[i] - yp[ID_1]));
						if (rij_n <= rij_e)
						{
							CNUM[i][k] = j;
							break;
						}
					}
				}
			}
		}
	}
}