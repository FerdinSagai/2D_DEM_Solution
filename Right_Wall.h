/*-------------------------------------------RIGHT SIDE WALL OF VESSEL-----------------------------------------------------------*/
void calculation_force_wall_right(int i)
{ 
	int j, k;
	double Cn, Ct;
	double delta_N, delta_T;
	double velocity_N, velocity_T;
	double sliding, friction;
	double Contact_time, travel_time;
	double ForceX, ForceY, Torque;
	double Wall_Angle;
	double Particle_MinX, Particle_MinY, Particle_MaxX, Particle_MaxY;
	double Margin_MinX, Margin_MinY, Margin_MaxX, Margin_MaxY;
	double distance_from_wall;

	CnWall = 2.0 * (-log(coeff_restitution) / (sqrt(pow(PI, 2.0) + pow(log(coeff_restitution), 2.0)))) * sqrt(mass[i] * Kn / 2.0); 
	CtWall = CnWall;

	Particle_MinX = xp[i] - (0.50 * d[i]);
	Particle_MinY = yp[i] - (0.50 * d[i]);
	Particle_MaxX = xp[i] + (0.50 * d[i]);
	Particle_MaxY = yp[i] + (0.50 * d[i]);

	Margin_MinX = xp[i] - (0.50 * mark * d[i]);
	Margin_MinY = yp[i] - (0.50 * mark * d[i]);
	Margin_MaxX = xp[i] + (0.50 * mark * d[i]);
	Margin_MaxY = yp[i] + (0.50 * mark * d[i]);
	if (Margin_MaxX >= Lx)
	{
		//  Angle computation -- X not applicable for walls.. 		
		Wall_Angle = 180;
		//velocity_N = up[i] * cos(Wall_Angle * PI / 180.0) + vp[i] * sin(Wall_Angle * PI / 180.0);
		//velocity_T = up[i] * sin(Wall_Angle * PI / 180.0) - vp[i] * cos(Wall_Angle * PI / 180.0);
		velocity_N = up[i];
		velocity_T = vp[i];
		//  Contact time computation
		if (Particle_MaxX > Lx)
		{
			travel_time = 0.00;
			Contact_time = STATICtime;
		}
		else
		{
			if (velocity_N != 0.0)
			{
				distance_from_wall = Particle_MaxX;
				travel_time = fabs((Lx - Particle_MaxX) / velocity_N);
				if ((travel_time <= STATICtime) && (velocity_N > 0.0))
				{
					Contact_time = STATICtime - travel_time;
				}
				else
				{
					Contact_time = 0.0;
				}
			}
			else
			{
				Contact_time = 0.0;
			}
		}
		
		if (Contact_time == 0.00)
		{
			Old_Delta_N_Walls[i][1] = 0.00;
			Old_Delta_T_Walls[i][1] = 0.00;
		}
		//  Incremental Normal deformation calculation
		delta_N = fabs(velocity_N * Contact_time);
		if (Contact_time == STATICtime)
		{
			delta_N = fabs(Particle_MaxX - Lx);
		}
		//  Total Normal deformation calculation
		Old_Delta_N_Walls[i][1] = Old_Delta_N_Walls[i][1] + delta_N;
		
		//  Incremental Tangential deformation calculation
		delta_T = (velocity_T + (0.5 * d[i] * omega[i])) * Contact_time;
		//  Total Tangential deformation calculation
		Old_Delta_T_Walls[i][1] = Old_Delta_T_Walls[i][1] + delta_T;
		
		//  Columb Criteria for Sliding. If particles slide, then, only frictional forces are in effect.
		if (fabs(KtWall * delta_T) >= coeff_friction * fabs(KnWall * delta_N))
		{
			if (delta_T < 0.0)
			{
				friction = fabs(KnWall * delta_N * coeff_friction);
			}
			else
			{
				friction = -fabs(KnWall * delta_N * coeff_friction);
			}
			sliding = 1.0;
		}

		//  Columb Criteria for Sliding. If particles roll, then, frictional and tangential forces are in effect.
		if (fabs(KtWall * delta_T) < coeff_friction * fabs(KnWall * delta_N))
		{
			friction = -KtWall * delta_T - CtWall * velocity_T;
			sliding = 0.0;
		}
		//ForceX = -sign(velocity_N)*fabs(KnWall * delta_N) - (1.0 - sliding) * CnWall * velocity_N;
		ForceX = -fabs(KnWall * delta_N) - (1.0 - sliding) * CnWall * velocity_N;
		ForceY = friction;
		Torque = KtWall * delta_T * 0.50 * d[i];

		update_particle_info(i, ForceX, ForceY, Torque, Contact_time);
	}
}