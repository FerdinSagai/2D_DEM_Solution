void calculation_force_particles(int i, int j)
{
	double delta_N, delta_T;
	double velocity_N, velocity_T;
	double sliding;
	double Contact_time, travel_time;
	double angle;

	double Cn = 2.0 * (-log(coeff_restitution) / (sqrt(pow(PI, 2.0) + pow(log(coeff_restitution), 2.0)))) * sqrt(mass[i] * Kn / 2.0);
	double Ct = Cn;
	// =======================================================================================================
	double vxij = up[i] - up[j];    				// Relative Velocity between particle i and j (X direction)
	double vyij = vp[i] - vp[j];    				// Relative Velocity between particle i and j (Y direction)
	double xij = xp[i] - xp[j];    					// Distance between particle i and j (X direction)
	double yij = yp[i] - yp[j];    					// Distance between particle i and j (Y direction)
	double rij = pow(xij, 2.0) + pow(yij, 2.0);		// Relative distance between the particles
	double d_avg = 0.5 * (d[i] + d[j]);
	if (rij < pow(mark * d_avg, 2.0))
	{
		max_deformation = 0.00;
		//  Angle computation.
		angle = (180.0 / PI) * atan2((yp[j] - yp[i]), (xp[j] - xp[i]));
		if (angle < 0.0)
		{
			angle = 360.0 + angle;
		}
		velocity_N = vxij * cos(angle * PI / 180.0) + vyij * sin(angle * PI / 180.0);
		velocity_T = vxij * sin(angle * PI / 180.0) - vyij * cos(angle * PI / 180.0);
		//  Contact time computation		
		Contact_time = compute_collision_time(i, j);
		//  Normal incremental deformation calculation
		delta_N = fabs(velocity_N * Contact_time);
		if (Contact_time == STATICtime)
		{
			delta_N = 0.5 * fabs(d_avg - sqrt(rij));
		}	
		//  Tangential incremental deformation calculation
		delta_T = (velocity_T + 0.5 * d_avg * (omega[i] - omega[j])) * Contact_time;
	
		//  Columb Criteria for Sliding. If particles slide, then only frictional forces are in effect. If sliding occurs then the sign of the frictional forces will be opposite to the deformation.															
		if (fabs(Kt * delta_T) >= coeff_friction * fabs(Kn * delta_N))			sliding = 1.00;
		//  Columb Criteria for Sliding. If particles roll, then, dashpot forces are in effect. In the absence of sliding, the tangential spring and damper are activated. Ct will always oppose motion. 
		if (fabs(Kt * delta_T) < coeff_friction * fabs(Kn * delta_N))			sliding = 0.00;		
		//  Determination of Cn and Ct in the X direction			
		double Force_Normal		= (-fabs(Kn * delta_N) + (1.0 - sliding) * Cn * fabs(velocity_N));
		double Force_Tangential	= (1.0 - sliding) * (-(Kt * delta_T) - Ct * velocity_T);
		double Force_Frictional	= sliding * -sign(velocity_T) * fabs(coeff_friction * Kn * delta_N);
		double ForceX			= Force_Normal * cos(angle * PI / 180.0) + Force_Tangential * sin(angle * PI / 180.0) + Force_Frictional * sin((angle) * PI / 180.0);
		double ForceY			= Force_Normal * sin(angle * PI / 180.0) - Force_Tangential * cos(angle * PI / 180.0) + Force_Frictional * cos((angle) * PI / 180.0);
		double Torque 			= Kt * delta_T * 0.5 * d_avg;

		//  Determination of Aggregate Displacement	
		update_particle_info(i, ForceX, ForceY, Torque, Contact_time);
	}
}