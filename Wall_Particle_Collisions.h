void calculation_force_wall(int i, int Surf_id)
{ 
	double delta_N, delta_T;
	CnWall = 2.0 * (-log(coeff_restitution) / (sqrt(pow(PI, 2.0) + pow(log(coeff_restitution), 2.0)))) * sqrt(mass[i] * Kn / 2.0); 
	CtWall = CnWall;
	
	double Wall_Angle			= Surf_Normal_Angle[Surf_id]* 180.0/PI;
	double velocity_N			= up[i] * cos(Wall_Angle * PI / 180.0) + vp[i] * sin(Wall_Angle * PI / 180.0);
	double velocity_T			= up[i] * sin(Wall_Angle * PI / 180.0) - vp[i] * cos(Wall_Angle * PI / 180.0);
	double contact_time			= 0.0;
	double travel_time			= 0.0;
	double sliding				= 0.0;
	double distance_from_wall	= Distance_between_Line_and_Point(Surf_id,xp[i],yp[i]);
	if (distance_from_wall < (0.50 * d[i]))							// In contact with the wall
	{
		travel_time = 0.00;
		contact_time = STATICtime - travel_time;
	}
	else if (distance_from_wall < mark * (0.50 * d[i]))				// Close enough to consider a hit
	{
		if (velocity_N != 0.0)
		{
			travel_time = fabs(distance_from_wall / velocity_N);
			if (travel_time < STATICtime)
			{
				contact_time = STATICtime - travel_time;
			}
			else
			{
				contact_time = 0.0;
			}
		}
		else
		{
			contact_time = 0.0;
		}
	}
	
	if (contact_time != 0.00)
	{
		//  Incremental Normal deformation calculation
		delta_N = fabs(velocity_N * contact_time);
		if (contact_time == STATICtime)
		{
			delta_N = fabs(distance_from_wall);
		}
		//  Incremental Tangential deformation calculation
		delta_T = (velocity_T + (0.5 * d[i] * omega[i])) * contact_time;
		//  Columb Criteria for Sliding. If particles slide, then, only frictional forces are in effect.
		if (fabs(KtWall * delta_T) >= coeff_friction * fabs(KnWall * delta_N))			sliding = 1.0;
		//  Columb Criteria for Sliding. If particles roll, then, frictional and tangential forces are in effect.
		if (fabs(KtWall * delta_T) < coeff_friction * fabs(KnWall * delta_N))			sliding = 0.0;
		double Force_Normal		= fabs(KnWall * delta_N) + (1.0 - sliding) * CnWall * fabs(velocity_N);
		double Force_Tangential	= (1.0 - sliding) * (-(KtWall * delta_T) - CtWall * velocity_T);
		double Force_Frictional	= sliding * -sign(velocity_T) * fabs(KnWall * delta_N * coeff_friction);
		double ForceX 			= Force_Normal * cos(Wall_Angle * PI / 180.0) + Force_Tangential * sin(Wall_Angle * PI / 180.0) + Force_Frictional * sin((Wall_Angle) * PI / 180.0);
		double ForceY 			= Force_Normal * sin(Wall_Angle * PI / 180.0) - Force_Tangential * cos(Wall_Angle * PI / 180.0) + Force_Frictional * cos((Wall_Angle) * PI / 180.0);
		double Torque 			= KtWall * delta_T * 0.50 * d[i];
		update_particle_info(i, ForceX, ForceY, Torque, contact_time);
	}
}