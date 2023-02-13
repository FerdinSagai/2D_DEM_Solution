#include<iostream>
//--------------------------------------------------------------------------------------------
double Lx, Ly;
int Nx, Ny;

double insertion_min_X, insertion_min_Y;
double insertion_max_X, insertion_max_Y;
double insertion_rate;

double extraction_min_X, extraction_min_Y;
double extraction_max_X, extraction_max_Y;
double discharge_rate;


// PROPERTIES OF PARTICLE CLASS
double dp;
double rho_s;								// Density of the particle
double phi_s;
double d[MAXNUM];							// Diameter of the particle
double mass[MAXNUM];						// Mass of the particle in kg
double MMOI[MAXNUM];						// Moment of inertia 
double xp[MAXNUM], yp[MAXNUM];				// Position for the Bed Particles
double zp[MAXNUM];
double up[MAXNUM], vp[MAXNUM];				// Velocity of Bed Particle
double poisson_ratio;
double coeff_friction;
double coeff_restitution;

double Kn, Kt;
double Cn, Ct;
double rho_g;

double time_DEM_Before, time_DEM_After, time_DEM_Total;
//--------------------------------------------------------------------------------------------
void input_details(void);
int main()
{
	//------Noting the cpu start time --------------------------------------------------------
	time_DEM_Before = clock();
	input_details();
	
	//------Noting the cpu finish time -------------------------------------------------------
	time_DEM_After = clock();
	//------Noting the total time ------------------------------------------------------------
	time_DEM_Total = time_DEM_After - time_DEM_Before;
	return 0;
}
void input_details()
{
	// Domain Definition
	Lx = 0.10;
	Ly = 0.30;
	
	Nx = 75;
	Ny = 150;
	
	insertion_min_X 	= 0.00;
	insertion_min_Y 	= 0.25;
	insertion_max_X 	= 0.10;
	insertion_max_Y 	= 0.30;
	insertion_rate		= 100;

	extraction_min_X	= 0.00;
	extraction_min_Y	= 0.00;
	extraction_max_X	= 0.38;
	extraction_max_Y	= 0.05;
	discharge_rate		= 0;
	
	// Gas Physical Properties
	rho_g = 1.178;
	
	// Physical Properties
	NUM					= 100;
	phi_s				= 1.00;
	rho_s				= 2500;
	dp					= 4.00e-3;
	coeff_friction		= 0.3;
	coeff_restitution	= 0.8;

	// Initializing the key variables------------------------------------------------------------------------
	for (int i = 0; i < NUM; i++)
	{
		xp[i]		= 0.00;
		yp[i]		= 0.00;
		zp[i]		= dp*0.50;
		up[i]		= 0.00;
		vp[i]		= 0.00;
		acc_x[i]	= 0.00;
		acc_y[i]	= 0.00;
		acc_ang[i]	= 0.00;
		omega[i]	= 0.00;
		d[i]		= dp;
		mass[i]			= rho_s * (4.0 / 3.0) * PI * pow(d[i] / 2.0, 3.0);
		MMOI[i]			= (2.0 / 5.0) * mass[i] * pow(d[i] / 2.0, 2.0);

		for (int j = 0; j < MAX_NEIGHBOUR+3; j++)
		{	
			CNUM[i][j] = NO_NEIGHBOUR;	
			Old_Delta_N[i][j] = 0.00;
			Old_Delta_T[i][j] = 0.00;
		}
	}
	Kn 					=	1000;
	Kt					=	Kn;
	Cn					=	2.0*(-log(coeff_restitution)/( sqrt(pow(PI,2.0) + pow(log(coeff_restitution),2.0))))*sqrt(mass[0]*Kn/2.0);;
	Ct					=	Cn;
}
