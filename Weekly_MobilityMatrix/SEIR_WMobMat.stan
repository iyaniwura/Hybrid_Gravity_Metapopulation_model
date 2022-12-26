// applicable to weekly cases data by symptoms onset date
functions {
	real[] SEIR(real t, 
							real[] y, 
							real[] theta, 
							real[] x_r,
							int[] x_i
	) {
		int n_region = x_i[1];            // number of regions
		int N_vec[n_region];         // e vector for the populations of the regions
		
            //  
		real h1 = x_r[1];               //  rate of going from exposed E1 to  in weeks
		real h2 = x_r[2];               //  rate of going from exposed E2 to I1 in weeks
		real gamma = x_r[3];            // rate of going from E2 to I1   
		real sigma = x_r[4];          // rate of going from I2 to Recovered  
		
		real y0_vars[4*n_region];          // initial prevalence of infection as a fraction of total infected for each subgroup in E1, E2, I1, I2
		real Init_Suscept[n_region];       // initial susceptible population
		real pi_ij_Vec[30*n_region*n_region];  // a vector for the distance between the regions: a symmetric matrix put inot a vector
		real TS_Mob_Vec[30*n_region];  // time-dependent mobility vector
		
		row_vector[7*n_region] dydt;       //  an array containing the derivative of all the compartment arrange in the order S, E1, E2, I1, I2, R, and Inci 
		
		// declaring some vectors and arrays that would be use later in the program
		row_vector[n_region] Force_Infection;
		row_vector[n_region] Infect;
		row_vector[n_region] Infect_Other_Pop;
		row_vector[n_region] N_hat;
		row_vector[n_region] Pop_Vec;
		
		row_vector[n_region] beta;
		row_vector[n_region] c0;
		row_vector[n_region] theta_vec;
		row_vector[7*n_region] y0;
		row_vector[7*n_region] y0_vec;
		row_vector[30*n_region] Mobility_Vec;
		row_vector[n_region] Wkly_Mob_Vec;
		
		real N_i;
		real N_j;
		real k;
		real dt;
		int mat_len;
		real m; //used to extracting the mobility matrix
		
		row_vector[n_region] D_ij;
		// row_vector[n_region] pi_ij;
		row_vector[n_region] f_ith_row;
		row_vector[n_region] pi_ith_row;
		row_vector[30*n_region*n_region] pi_ith_rowVector;
		row_vector[n_region*n_region] pi_ij;
		
		
		// model comparments
		row_vector[n_region] S;         // susceptible pop
		row_vector[n_region] E1;        // Exposed but not infectious
		row_vector[n_region] E2;        // Pre-symptomatic
		row_vector[n_region] I1;        // Infectious (first half of the infectious period)
		row_vector[n_region] I2;        // Infectious (second half of the infectious period)
		row_vector[n_region] M;         // Recovered
		row_vector[n_region] Inci;      // Incidence (computed from those leaving E2 to I1)
		
		
		// derivatives
		row_vector[n_region] dydt1;         // susceptible pop
		row_vector[n_region] dydt2;        // Exposed but not infectious
		row_vector[n_region] dydt3;        // Pre-symptomatic
		row_vector[n_region] dydt4;        // Infectious (first half of the infectious period)
		row_vector[n_region] dydt5;        // Infectious (second half of the infectious period)
		row_vector[n_region] dydt6;         // Recovered
		row_vector[n_region] dydt7;      // Incidence (computed from those leaving E2 to I1)
		
		// Estimated parameters
		real i0	 = theta[1];    // Initial incidence
		real c0_PopMean = theta[2];
		real c0_PopVar = theta[3];
		real c0_1	 = theta[4];    // 
		real c0_2	 = theta[5];    // 
		real c0_3	 = theta[6];    // 
		real c0_4	 = theta[7];    // 
		real c0_5	 = theta[8];    // 
		real c0_6	 = theta[9];    // 
		real c0_7	 = theta[10];    // 
		real c0_8	 = theta[11];    // 
		real c0_9	 = theta[12];    // 
		real c0_10	 = theta[13];    // 
		real c0_11	 = theta[14];    // 
		real c0_12	 = theta[15];    // 
		real c0_13	 = theta[16];    // 
		real theta_1 = theta[17];    // 
		real c1 = theta[18];    // 
		
		real d_PopMean = theta[19];    // 
		real d_PopVar = theta[20];    // 
		real d1 = theta[21];    // 
		real d2 = theta[22];    // 
		real d3 = theta[23];    // 
		real d4 = theta[24];    // 
		real d5 = theta[25];    // 
		real d6 = theta[26];    // 
		real d7 = theta[27];    // 
	
		
		// real vector parameters
		N_vec = x_i[2:(n_region+1)];
		
		
		y0_vars = x_r[5:(4*n_region+4)];   	// 4*n_region elements
		Init_Suscept = x_r[(4*n_region+5):(5*n_region+4)];       // Initial population of susceptibles: n_region elements
		pi_ij_Vec = x_r[(5*n_region+5):((5*n_region+4) + (30*n_region*n_region))];   	 // 169 elements
		TS_Mob_Vec = x_r[((5*n_region+4) + (30*n_region*n_region) + 1):(30*n_region + ((5*n_region+4) + (30*n_region*n_region)))];
		
		
		
		// relating the compartments of the model to y
		S = to_row_vector(y[1:(n_region)]);
		E1 = to_row_vector(y[(n_region+1):(2*n_region)]);
		E2 = to_row_vector(y[(2*n_region+1):(3*n_region)]);
		I1 = to_row_vector(y[(3*n_region+1):(4*n_region)]);
		I2 = to_row_vector(y[(4*n_region+1):(5*n_region)]);
		M = to_row_vector(y[(5*n_region+1):(6*n_region)]);
		Inci = to_row_vector(y[(6*n_region+1):(7*n_region)]);
		
	
		Pop_Vec = to_row_vector(N_vec);
		pi_ith_rowVector = to_row_vector(pi_ij_Vec);
		mat_len = n_region*n_region;
		Mobility_Vec = to_row_vector(TS_Mob_Vec);
		
		

		if (t <= 1){
				pi_ij = pi_ith_rowVector[1:(mat_len)]; // month 1
				Wkly_Mob_Vec = Mobility_Vec[1:(n_region)]; // month 1
				
		} else if (t > 1 && t <= 2){
			pi_ij = pi_ith_rowVector[(mat_len+1):(2*mat_len)]; // month 2
			Wkly_Mob_Vec = Mobility_Vec[(n_region+1):(2*n_region)]; // month 2
			
		} else if (t > 2 && t <= 3){
			pi_ij = pi_ith_rowVector[(2*mat_len+1):(3*mat_len)]; // month 3
			Wkly_Mob_Vec = Mobility_Vec[(2*n_region+1):(3*n_region)]; // month 3
			
		} else if (t > 3 && t <= 4){
			pi_ij = pi_ith_rowVector[(3*mat_len+1):(4*mat_len)]; // month 4
			Wkly_Mob_Vec = Mobility_Vec[(3*n_region+1):(4*n_region)]; // month 4
			
		} else if (t > 4 && t <= 5){
			pi_ij = pi_ith_rowVector[(4*mat_len+1):(5*mat_len)]; // month 5
			Wkly_Mob_Vec = Mobility_Vec[(4*n_region+1):(5*n_region)]; // month 5
			
		} else if (t > 5 && t <= 6){
			pi_ij = pi_ith_rowVector[(5*mat_len+1):(6*mat_len)]; // month 5
			Wkly_Mob_Vec = Mobility_Vec[(5*n_region+1):(6*n_region)]; // month 5
			
		} else if (t > 6 && t <= 7){
			pi_ij = pi_ith_rowVector[(6*mat_len+1):(7*mat_len)]; // month 5
			Wkly_Mob_Vec = Mobility_Vec[(6*n_region+1):(7*n_region)]; // month 5
			
		} else if (t > 7 && t <= 8){
			pi_ij = pi_ith_rowVector[(7*mat_len+1):(8*mat_len)]; // month 5
			Wkly_Mob_Vec = Mobility_Vec[(7*n_region+1):(8*n_region)]; // month 5
			
		} else if (t > 8 && t <= 9){
			pi_ij = pi_ith_rowVector[(8*mat_len+1):(9*mat_len)]; // month 5
			Wkly_Mob_Vec = Mobility_Vec[(8*n_region+1):(9*n_region)]; // month 5
			
		} else if (t > 9 && t <= 10){
			pi_ij = pi_ith_rowVector[(9*mat_len+1):(10*mat_len)]; // month 6
			Wkly_Mob_Vec = Mobility_Vec[(9*n_region+1):(10*n_region)]; // month 6
			
		} else if (t > 10 && t <= 11){
			pi_ij = pi_ith_rowVector[(10*mat_len+1):(11*mat_len)]; // month 6
			Wkly_Mob_Vec = Mobility_Vec[(10*n_region+1):(11*n_region)]; // month 6
			
		} else if (t > 11 && t <= 12){
			pi_ij = pi_ith_rowVector[(11*mat_len+1):(12*mat_len)]; // month 6
			Wkly_Mob_Vec = Mobility_Vec[(11*n_region+1):(12*n_region)]; // month 6
			
		} else if (t > 12 && t <= 13){
			pi_ij = pi_ith_rowVector[(12*mat_len+1):(13*mat_len)]; // month 6
			Wkly_Mob_Vec = Mobility_Vec[(12*n_region+1):(13*n_region)]; // month 6
			
		} else if (t > 13 && t <= 14){
			pi_ij = pi_ith_rowVector[(13*mat_len+1):(14*mat_len)]; // month 6
			Wkly_Mob_Vec = Mobility_Vec[(13*n_region+1):(14*n_region)]; // month 6
			
		} else if (t > 14 && t <= 15){
			pi_ij = pi_ith_rowVector[(14*mat_len+1):(15*mat_len)]; // month 6
			Wkly_Mob_Vec = Mobility_Vec[(14*n_region+1):(15*n_region)]; // month 6
			
		} else if (t > 15 && t <= 16){
			pi_ij = pi_ith_rowVector[(15*mat_len+1):(16*mat_len)]; // month 6
			Wkly_Mob_Vec = Mobility_Vec[(15*n_region+1):(16*n_region)]; // month 6
			
		} else if (t > 16 && t <= 17){
			pi_ij = pi_ith_rowVector[(16*mat_len+1):(17*mat_len)]; // month 6
			Wkly_Mob_Vec = Mobility_Vec[(16*n_region+1):(17*n_region)]; // month 6
			
		} else if (t > 17 && t <= 18){
			pi_ij = pi_ith_rowVector[(17*mat_len+1):(18*mat_len)]; // month 6
			Wkly_Mob_Vec = Mobility_Vec[(17*n_region+1):(18*n_region)]; // month 6
			
		} else if (t > 18 && t <= 19){
			pi_ij = pi_ith_rowVector[(18*mat_len+1):(19*mat_len)]; // month 6
			Wkly_Mob_Vec = Mobility_Vec[(18*n_region+1):(19*n_region)]; // month 6
			
		} else if (t > 19 && t <= 20){
			pi_ij = pi_ith_rowVector[(19*mat_len+1):(20*mat_len)]; // month 6
			Wkly_Mob_Vec = Mobility_Vec[(19*n_region+1):(20*n_region)]; // month 6
			
		} else if (t > 20 && t <= 21){
			pi_ij = pi_ith_rowVector[(20*mat_len+1):(21*mat_len)]; // month 6
			Wkly_Mob_Vec = Mobility_Vec[(20*n_region+1):(21*n_region)]; // month 6
			
		} else if (t > 21 && t <= 22){
			pi_ij = pi_ith_rowVector[(21*mat_len+1):(22*mat_len)]; // month 6
			Wkly_Mob_Vec = Mobility_Vec[(21*n_region+1):(22*n_region)]; // month 6
			
		} else if (t > 22 && t <= 23){
			pi_ij = pi_ith_rowVector[(22*mat_len+1):(23*mat_len)]; // month 6
			Wkly_Mob_Vec = Mobility_Vec[(22*n_region+1):(23*n_region)]; // month 6
			
		} else if (t > 23 && t <= 24){
			pi_ij = pi_ith_rowVector[(23*mat_len+1):(24*mat_len)]; // month 6
			Wkly_Mob_Vec = Mobility_Vec[(23*n_region+1):(24*n_region)]; // month 6
			
		} else if (t > 24 && t <= 25){
			pi_ij = pi_ith_rowVector[(24*mat_len+1):(25*mat_len)]; // month 6
			Wkly_Mob_Vec = Mobility_Vec[(24*n_region+1):(25*n_region)]; // month 6
			
		} else if (t > 25 && t <= 26){
			pi_ij = pi_ith_rowVector[(25*mat_len+1):(26*mat_len)]; // month 6
			Wkly_Mob_Vec = Mobility_Vec[(25*n_region+1):(26*n_region)]; // month 6
			
		} else if (t > 26 && t <= 27){
			pi_ij = pi_ith_rowVector[(26*mat_len+1):(27*mat_len)]; // month 6
			Wkly_Mob_Vec = Mobility_Vec[(26*n_region+1):(27*n_region)]; // month 6
			
		} else if (t > 27 && t <= 28){
			pi_ij = pi_ith_rowVector[(27*mat_len+1):(28*mat_len)]; // month 6
			Wkly_Mob_Vec = Mobility_Vec[(27*n_region+1):(28*n_region)]; // month 6
			
		} else if (t > 28 && t <= 29){
			pi_ij = pi_ith_rowVector[(28*mat_len+1):(29*mat_len)]; // month 6
			Wkly_Mob_Vec = Mobility_Vec[(28*n_region+1):(29*n_region)]; // month 6
			
		} else if (t > 29){
			pi_ij = pi_ith_rowVector[(29*mat_len+1):(30*mat_len)]; // month 7
			Wkly_Mob_Vec = Mobility_Vec[(29*n_region+1):(30*n_region)]; // month 7
			
		}
		
	  if (t <= 4){
			dt = 0;
		} else if(t > 4 && t <= 8){
			dt = d1;
		}  else if(t > 8 && t <= 12){
			dt = d2;
		}  else if(t > 12 && t <= 16){
			dt = d3;
		}  else if(t > 16 && t <= 20){
			dt = d4;
		}  else if(t > 20 && t <= 24){
			dt = d5;
		}  else if(t > 24 && t <= 28 ){
			dt = d6;
		}	else if(t > 28 ){
			dt = d7;
		}	
		
		// compute the force of infection and N_hat
		for (i in 1:n_region){
			pi_ith_row = pi_ij[((i-1)*n_region + 1):(i*n_region)];
		 	N_hat[i] = (1 - theta_1) * Pop_Vec[i] + (theta_1 * sum(pi_ith_row .* Pop_Vec) );
			for (j in 1:n_region){
				Infect[j] = pi_ith_row[j]*(E2[j] + I1[j] + I2[j]);
			}
			Infect_Other_Pop[i] = sum(Infect)/N_hat[i];    // infections due to mixing with other population
		}
		
		// Force_Infection = 2.5*rep_row_vector(1,n_region);
		Force_Infection = (1 - theta_1) * ((E2 + I1 + I2) ./ N_hat) + (theta_1 * Infect_Other_Pop);
		c0 = [c0_1, c0_2, c0_3, c0_4, c0_5, c0_6, c0_7, c0_8, c0_9, c0_10, c0_11, c0_12, c0_13];
		beta = exp( c0 + c1* Wkly_Mob_Vec + dt);
		
		// ODE system
		dydt1 = -(beta .* S .* Force_Infection);    
		dydt2 = (beta .* S .* Force_Infection) - h1*E1;
		dydt3 = h1*E1 - h2*E2;
		dydt4 = h2*E2 - gamma*I1;
		dydt5 = gamma*I1 - sigma*I2;
		dydt6 = sigma*I2;
		dydt7 = h2*E2;
		
		// combine the derivatives into one vector
		dydt[1:(n_region)] = dydt1;
		dydt[(n_region+1):(2*n_region)] = dydt2;
		dydt[(2*n_region+1):(3*n_region)] = dydt3;
		dydt[(3*n_region+1):(4*n_region)] = dydt4;
		dydt[(4*n_region+1):(5*n_region)] = dydt5;
		dydt[(5*n_region+1):(6*n_region)] = dydt6;
		dydt[(6*n_region+1):(7*n_region)] = dydt7;
		
		return(to_array_1d(dydt)); 
		
		
	}
}



data {
	// Structure
	int n_region;
	int n_weeks;       // total number of weeks
	real t0;          //starting time
	real ts[n_weeks];  // time bins
	int doprint;
	int inference;
	real h1;
	real h2;
	real gamma;
	real sigma;
	
	
	// Data to fit
	// number of days with reported incidence
	int Rep_Reg_1[n_weeks];     
	int Rep_Reg_2[n_weeks];    
	int Rep_Reg_3[n_weeks];    
	int Rep_Reg_4[n_weeks];    
	int Rep_Reg_5[n_weeks];   
	int Rep_Reg_6[n_weeks];   
	int Rep_Reg_7[n_weeks];   
	int Rep_Reg_8[n_weeks]; 
	int Rep_Reg_9[n_weeks]; 
	int Rep_Reg_10[n_weeks]; 
	int Rep_Reg_11[n_weeks]; 
	int Rep_Reg_12[n_weeks];
	int Rep_Reg_13[n_weeks];
	
	int N_vec[n_region];  
	
		// Priors
	real p_i0[2];   // total initial prevalence
	real p_theta[2];   // total initial prevalence
	real p_c0Mean[2];   // total initial prevalence
	real p_c0Var[2];
	real p_dMean[2];   // total initial prevalence
	real p_dVar[2];
	real p_phi;

	// Fixed parameters
	  real y0_vars[4*n_region];  
		real Init_Suscept[n_region];       // initial susceptible population
		real pi_ij_Vec[30*n_region*n_region]; 
		real TS_Mob_Vec[30*n_region];

	// Fixed corrections
	real p_underreport_cases; // correction for cases reported later
}
	

transformed data {
	real x_r[((5*n_region+4) + (30*n_region*n_region) + (30*n_region)) ];
	int x_i[n_region + 1];
	
	x_i[1] = n_region;
	x_i[2:(n_region+1)] = N_vec;        // total by region
	
	x_r[1] = h1;
	x_r[2] = h2;
	x_r[3] = gamma;
	x_r[4] = sigma;

	x_r[5:(4*n_region+4)] = y0_vars;   	// 4*n_region elements
	x_r[(4*n_region+5):(5*n_region+4)] = Init_Suscept;       // Initial population of susceptibles: n_region elements
	x_r[(5*n_region+5):((5*n_region+4) + (30*n_region*n_region))] = pi_ij_Vec;   	 // 169 elements
	x_r[((5*n_region+4) + (30*n_region*n_region) + 1):(30*n_region + ((5*n_region+4) + (30*n_region*n_region)))] = TS_Mob_Vec;

}



parameters{
	real<lower=1, upper=6000> i0;      // total initial prevalence
	real<lower=0, upper=1> theta_1;
	real<lower=-2, upper=2> c0_PopMean;
	real<lower=0, upper=5> c0_PopVar;
	real<lower=-2, upper=2> c0_1;
	real<lower=-2, upper=2> c0_2;
	real<lower=-2, upper=2> c0_3;
	real<lower=-2, upper=2> c0_4;
	real<lower=-2, upper=2> c0_5;
	real<lower=-2, upper=2> c0_6;
	real<lower=-2, upper=2> c0_7;
	real<lower=-2, upper=2> c0_8;
	real<lower=-2, upper=2> c0_9;
	real<lower=-2, upper=2> c0_10;
	real<lower=-2, upper=2> c0_11;
	real<lower=-2, upper=2> c0_12;
	real<lower=-2, upper=2> c0_13;
	
	real<lower=-5, upper=5> c1;
	
	real<lower=-2, upper=2> d_PopMean;    // 
	real<lower=0, upper=5> d_PopVar;    // 
	real<lower=-2, upper=2> d1;    // 
	real<lower=-2, upper=2> d2;    // 
	real<lower=-2, upper=2> d3;    // 
	real<lower=-2, upper=2> d4;    // 
	real<lower=-2, upper=2> d5;    // 
	real<lower=-2, upper=2> d6;    // 
	real<lower=-2, upper=2> d7;    // 

	real<lower=0> phi; // 
}

transformed parameters {
	// transformed parameters
	real y[n_weeks,7*(n_region)]; // raw ODE output
	real y0[7*(n_region)];  // initial state
	row_vector[4*n_region] y0_vec;         // initial condition vector
	row_vector[7*n_region] y0_temp;         // initial condition vector
	
	// change of format for integrate_ode_rk45
	real theta[27];     // vector of parameters
	
	
	// ode outputs
	// cumulative incidence
	vector[n_weeks] Cum_inci_Reg1; // overall case incidence by day
	vector[n_weeks] Cum_inci_Reg2;   // incidence for 0-2 years
	vector[n_weeks] Cum_inci_Reg3;   // incidence for 3-4 years
	vector[n_weeks] Cum_inci_Reg4;   // incidence for 5-17 years
	vector[n_weeks] Cum_inci_Reg5;   // incidence for 18-24 years
	vector[n_weeks] Cum_inci_Reg6;   // incidence for 25-54 years
	vector[n_weeks] Cum_inci_Reg7;   // incidence for 55-64 years
	vector[n_weeks] Cum_inci_Reg8;   // incidence for 55-64 years
	vector[n_weeks] Cum_inci_Reg9;   // incidence for 55-64 years
	vector[n_weeks] Cum_inci_Reg10;   // incidence for 55-64 years
	vector[n_weeks] Cum_inci_Reg11;   // incidence for 55-64 years
	vector[n_weeks] Cum_inci_Reg12; 
	vector[n_weeks] Cum_inci_Reg13; 
	
	// incidence
	vector[n_weeks] Inci_Reg1; // overall case incidence by day
	vector[n_weeks] Inci_Reg2;   // Incidence for 0-2 years
	vector[n_weeks] Inci_Reg3;   // Incidence for 3-4 years
	vector[n_weeks] Inci_Reg4;   // Incidence for 5-17 years
	vector[n_weeks] Inci_Reg5;   // Incidence for 18-24 years
	vector[n_weeks] Inci_Reg6;   // Incidence for 25-54 years
	vector[n_weeks] Inci_Reg7;   // Incidence for 55-64 years
	vector[n_weeks] Inci_Reg8;   // Incidence for 55-64 years
	vector[n_weeks] Inci_Reg9;   // Incidence for 55-64 years
	vector[n_weeks] Inci_Reg10;   // Incidence for 55-64 years
	vector[n_weeks] Inci_Reg11;   // Incidence for 55-64 years
	vector[n_weeks] Inci_Reg12;
	vector[n_weeks] Inci_Reg13;
	
	
	// set up the initial conditions:
	y0_vec = to_row_vector(y0_vars);
	y0_temp[1:n_region] = to_row_vector(Init_Suscept);   // initial susceptibles
	y0_temp[(n_region+1):(5*n_region)] = y0_vec[1:(4*n_region)]*i0;   // E1, E2, I1, I2: i0 is initial prevalance
	y0_temp[(5*n_region+1):(6*n_region)] = rep_row_vector(0, n_region);     // R, set to zeros
	y0_temp[(6*n_region+1):(7*n_region)] =  rep_row_vector(0, n_region);     // Incidence, set to zeros
	
	y0 = to_array_1d(y0_temp);
	
	
	
	// change of format for integrate_ode_rk45
	theta[1:27] = { i0, c0_PopMean, c0_PopVar, c0_1, c0_2, c0_3, c0_4, c0_5, c0_6, c0_7, c0_8, c0_9, c0_10, c0_11, c0_12, c0_13, theta_1,
									c1, d_PopMean, d_PopVar, d1, d2, d3, d4, d5, d6, d7 };

	// // // run ODE solver
	y = integrate_ode_rk45( SEIR, // ODE function
													y0, // initial states
													t0, // t0
													ts, // evaluation dates (ts)
													theta, // parameters
													x_r, // real data
													x_i, // integer data
													1.0E-6, 1.0E-6, 1.0E3); // tolerances and maximum steps
	// y = integrate_ode_bdf(SEIR, y0, t0, ts, theta, x_r, x_i); // tolerances and maximum steps
	
	// extract and format ODE results (1.0E-9 correction to avoid negative values due to unprecise estimates of zeros as tolerance is 1.0E-10)
	
	
	// Extracting incidence from the output of the ode solver	
	Cum_inci_Reg1 = to_vector(y[,(6*n_region + 1)]) ;  
	Cum_inci_Reg2 = to_vector( y[,(6*n_region + 2)] ) ; 
	Cum_inci_Reg3 = to_vector( y[,(6*n_region + 3)] ) ; 
	Cum_inci_Reg4 = to_vector( y[,(6*n_region + 4)] ) ; 
	Cum_inci_Reg5 = to_vector( y[,(6*n_region + 5)] ) ; 
	Cum_inci_Reg6 = to_vector( y[,(6*n_region + 6)] ) ; 
	Cum_inci_Reg7 = to_vector( y[,(6*n_region + 7)] ) ; 
	Cum_inci_Reg8 = to_vector( y[,(6*n_region + 8)] ) ; 
	Cum_inci_Reg9 = to_vector( y[,(6*n_region + 9)] ) ; 
	Cum_inci_Reg10 = to_vector( y[,(6*n_region + 10)] ) ; 
	Cum_inci_Reg11 = to_vector( y[,(6*n_region + 11)] ) ; 
	Cum_inci_Reg12 = to_vector( y[,(6*n_region + 12)] ) ; 
	Cum_inci_Reg13 = to_vector( y[,(6*n_region + 13)] ) ; 
	

	Inci_Reg1[1] =  Cum_inci_Reg1[1]; 	  Inci_Reg1[2:n_weeks] =  Cum_inci_Reg1[2:n_weeks] - Cum_inci_Reg1[1:(n_weeks - 1)];
	Inci_Reg2[1] =  Cum_inci_Reg2[1];     Inci_Reg2[2:n_weeks] =  Cum_inci_Reg2[2:n_weeks] - Cum_inci_Reg2[1:(n_weeks - 1)];
	Inci_Reg3[1] =  Cum_inci_Reg3[1];     Inci_Reg3[2:n_weeks] =  Cum_inci_Reg3[2:n_weeks] - Cum_inci_Reg3[1:(n_weeks - 1)];
	Inci_Reg4[1] =  Cum_inci_Reg4[1];     Inci_Reg4[2:n_weeks] =  Cum_inci_Reg4[2:n_weeks] - Cum_inci_Reg4[1:(n_weeks - 1)];
	Inci_Reg5[1] =  Cum_inci_Reg5[1];     Inci_Reg5[2:n_weeks] =  Cum_inci_Reg5[2:n_weeks] - Cum_inci_Reg5[1:(n_weeks - 1)];
	Inci_Reg6[1] =  Cum_inci_Reg6[1];     Inci_Reg6[2:n_weeks] =  Cum_inci_Reg6[2:n_weeks] - Cum_inci_Reg6[1:(n_weeks - 1)];
	Inci_Reg7[1] =  Cum_inci_Reg7[1];     Inci_Reg7[2:n_weeks] =  Cum_inci_Reg7[2:n_weeks] - Cum_inci_Reg7[1:(n_weeks - 1)];
	Inci_Reg8[1] =  Cum_inci_Reg8[1];     Inci_Reg8[2:n_weeks] =  Cum_inci_Reg8[2:n_weeks] - Cum_inci_Reg8[1:(n_weeks - 1)];
	Inci_Reg9[1] =  Cum_inci_Reg9[1];     Inci_Reg9[2:n_weeks] =  Cum_inci_Reg9[2:n_weeks] - Cum_inci_Reg9[1:(n_weeks - 1)];
	Inci_Reg10[1] =  Cum_inci_Reg10[1];   Inci_Reg10[2:n_weeks] =  Cum_inci_Reg10[2:n_weeks] - Cum_inci_Reg10[1:(n_weeks - 1)];
	Inci_Reg11[1] =  Cum_inci_Reg11[1];   Inci_Reg11[2:n_weeks] =  Cum_inci_Reg11[2:n_weeks] - Cum_inci_Reg11[1:(n_weeks - 1)];
	Inci_Reg12[1] =  Cum_inci_Reg12[1];   Inci_Reg12[2:n_weeks] =  Cum_inci_Reg12[2:n_weeks] - Cum_inci_Reg12[1:(n_weeks - 1)];
	Inci_Reg13[1] =  Cum_inci_Reg13[1];   Inci_Reg13[2:n_weeks] =  Cum_inci_Reg13[2:n_weeks] - Cum_inci_Reg13[1:(n_weeks - 1)];
	
}

model {
	// priors
	i0 ~ lognormal(p_i0[1], p_i0[2]);
	theta_1 ~ normal(p_theta[1], p_theta[2]);
	c0_PopMean ~ normal(p_c0Mean[1], p_c0Mean[2]);
	c0_PopVar ~ normal(p_c0Var[1], p_c0Var[2]);
	c0_1 ~ normal(c0_PopMean, c0_PopVar);
	c0_2 ~ normal(c0_PopMean, c0_PopVar);
	c0_3 ~ normal(c0_PopMean, c0_PopVar);
	c0_4 ~ normal(c0_PopMean, c0_PopVar);
	c0_5 ~ normal(c0_PopMean, c0_PopVar);
	c0_6 ~ normal(c0_PopMean, c0_PopVar);
	c0_7 ~ normal(c0_PopMean, c0_PopVar);
	c0_8 ~ normal(c0_PopMean, c0_PopVar);
	c0_9 ~ normal(c0_PopMean, c0_PopVar);
	c0_10 ~ normal(c0_PopMean, c0_PopVar);
	c0_11 ~ normal(c0_PopMean, c0_PopVar);
	c0_12 ~ normal(c0_PopMean, c0_PopVar);
	c0_13 ~ normal(c0_PopMean, c0_PopVar);
	
	c1 ~ normal(0, 5);
	
	d_PopMean ~ normal(p_dMean[1], p_dMean[2]);
	d_PopVar ~ normal(p_dVar[1], p_dVar[2]);
	d1 ~ normal(d_PopMean, d_PopVar);
	d2 ~ normal(d_PopMean, d_PopVar);
	d3 ~ normal(d_PopMean, d_PopVar);
	d4 ~ normal(d_PopMean, d_PopVar);
	d5 ~ normal(d_PopMean, d_PopVar);
	d6 ~ normal(d_PopMean, d_PopVar);
	d7 ~ normal(d_PopMean, d_PopVar);

	phi ~ exponential(p_phi);

	// debug
	if(doprint == 1) {
		//print("gamma: ", gamma);
		print("y[5,]: ",y[5,]);
		
	}

		if (inference!=0){
			Rep_Reg_1 ~ neg_binomial_2((Inci_Reg1), 1/phi);
			Rep_Reg_2 ~ neg_binomial_2((Inci_Reg2), 1/phi);       // incidence for 0-2 years
			Rep_Reg_3 ~ neg_binomial_2((Inci_Reg3), 1/phi);       // incidence for 3-4 years
			Rep_Reg_4 ~ neg_binomial_2((Inci_Reg4), 1/phi);       // incidence for 5-17 years
			Rep_Reg_5 ~ neg_binomial_2((Inci_Reg5), 1/phi);       // incidence for 18-24 years
			Rep_Reg_6 ~ neg_binomial_2((Inci_Reg6), 1/phi);       // incidence for 25-54 years
			Rep_Reg_7 ~ neg_binomial_2((Inci_Reg7), 1/phi);       // incidence for 55-64 years
			Rep_Reg_8 ~ neg_binomial_2((Inci_Reg8), 1/phi);       // incidence for 65+h years
			Rep_Reg_9 ~ neg_binomial_2((Inci_Reg9), 1/phi);       // incidence for 65+c years
			Rep_Reg_10 ~ neg_binomial_2((Inci_Reg10), 1/phi);       // incidence for 65+c years
			Rep_Reg_11 ~ neg_binomial_2((Inci_Reg11), 1/phi);       // incidence for 65+c years
			Rep_Reg_12 ~ neg_binomial_2((Inci_Reg12), 1/phi);       // incidence for 65+c years
			Rep_Reg_13 ~ neg_binomial_2((Inci_Reg13), 1/phi);       // incidence for 65+c years
		}
}


generated quantities{
	
	vector[(n_weeks*n_region)] log_lik;
  int CasesData[(n_weeks*n_region)] ;
  real ModelOutput[(n_weeks*n_region)];
	
	real pred_cases_Reg_1[n_weeks];
	real pred_cases_Reg_2[n_weeks];
	real pred_cases_Reg_3[n_weeks];
	real pred_cases_Reg_4[n_weeks];
	real pred_cases_Reg_5[n_weeks];
	real pred_cases_Reg_6[n_weeks];
	real pred_cases_Reg_7[n_weeks];
	real pred_cases_Reg_8[n_weeks];
	real pred_cases_Reg_9[n_weeks];
	real pred_cases_Reg_10[n_weeks];
	real pred_cases_Reg_11[n_weeks];
	real pred_cases_Reg_12[n_weeks];
	real pred_cases_Reg_13[n_weeks];
	
	pred_cases_Reg_1 = neg_binomial_2_rng(Inci_Reg1, 1/phi);
	pred_cases_Reg_2 = neg_binomial_2_rng(Inci_Reg2, 1/phi);
	pred_cases_Reg_3 = neg_binomial_2_rng(Inci_Reg3, 1/phi);
	pred_cases_Reg_4 = neg_binomial_2_rng(Inci_Reg4, 1/phi);
	pred_cases_Reg_5 = neg_binomial_2_rng(Inci_Reg5, 1/phi);
	pred_cases_Reg_6 = neg_binomial_2_rng(Inci_Reg6, 1/phi);
	pred_cases_Reg_7 = neg_binomial_2_rng(Inci_Reg7, 1/phi);
	pred_cases_Reg_8 = neg_binomial_2_rng(Inci_Reg8, 1/phi);
	pred_cases_Reg_9 = neg_binomial_2_rng(Inci_Reg9, 1/phi);
	pred_cases_Reg_10 = neg_binomial_2_rng(Inci_Reg10, 1/phi);
	pred_cases_Reg_11 = neg_binomial_2_rng(Inci_Reg11, 1/phi);
	pred_cases_Reg_12 = neg_binomial_2_rng(Inci_Reg12, 1/phi);
	pred_cases_Reg_13 = neg_binomial_2_rng(Inci_Reg13, 1/phi);
	
	
	// Code for model comparison
	// putting the cases data into a vector
	CasesData[1:n_weeks] = to_array_1d(Rep_Reg_1);
	CasesData[(n_weeks+1):(2*n_weeks)] = to_array_1d(Rep_Reg_2);
	CasesData[(2*n_weeks+1):(3*n_weeks)] = to_array_1d(Rep_Reg_3);
	CasesData[(3*n_weeks+1):(4*n_weeks)] = to_array_1d(Rep_Reg_4);
	CasesData[(4*n_weeks+1):(5*n_weeks)] = to_array_1d(Rep_Reg_5);
	CasesData[(5*n_weeks+1):(6*n_weeks)] = to_array_1d(Rep_Reg_6);
	CasesData[(6*n_weeks+1):(7*n_weeks)] = to_array_1d(Rep_Reg_7);
	CasesData[(7*n_weeks+1):(8*n_weeks)] = to_array_1d(Rep_Reg_8);
	CasesData[(8*n_weeks+1):(9*n_weeks)] = to_array_1d(Rep_Reg_9);
	CasesData[(9*n_weeks+1):(10*n_weeks)] = to_array_1d(Rep_Reg_10);
	CasesData[(10*n_weeks+1):(11*n_weeks)] = to_array_1d(Rep_Reg_11);
	CasesData[(11*n_weeks+1):(12*n_weeks)] = to_array_1d(Rep_Reg_12);
	CasesData[(12*n_weeks+1):(13*n_weeks)] = to_array_1d(Rep_Reg_13);

	// putting the model output into a vector
	ModelOutput[1:n_weeks] = to_array_1d(Inci_Reg1);
	ModelOutput[(n_weeks+1):(2*n_weeks)] = to_array_1d(Inci_Reg2);
	ModelOutput[(2*n_weeks+1):(3*n_weeks)] = to_array_1d(Inci_Reg3);
	ModelOutput[(3*n_weeks+1):(4*n_weeks)] = to_array_1d(Inci_Reg4);
	ModelOutput[(4*n_weeks+1):(5*n_weeks)] = to_array_1d(Inci_Reg5);
	ModelOutput[(5*n_weeks+1):(6*n_weeks)] = to_array_1d(Inci_Reg6);
	ModelOutput[(6*n_weeks+1):(7*n_weeks)] = to_array_1d(Inci_Reg7);
	ModelOutput[(7*n_weeks+1):(8*n_weeks)] = to_array_1d(Inci_Reg8);
	ModelOutput[(8*n_weeks+1):(9*n_weeks)] = to_array_1d(Inci_Reg9);
	ModelOutput[(9*n_weeks+1):(10*n_weeks)] = to_array_1d(Inci_Reg10);
	ModelOutput[(10*n_weeks+1):(11*n_weeks)] = to_array_1d(Inci_Reg11);
	ModelOutput[(11*n_weeks+1):(12*n_weeks)] = to_array_1d(Inci_Reg12);
	ModelOutput[(12*n_weeks+1):(13*n_weeks)] = to_array_1d(Inci_Reg13);
	

  
  for (n in 1:(n_region*n_weeks)) {
    log_lik[n] = neg_binomial_2_lpmf( CasesData[n] | ModelOutput[n], 1/phi);
  }

	
	
}
