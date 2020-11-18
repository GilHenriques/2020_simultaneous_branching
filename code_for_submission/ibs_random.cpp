#include<vector>
#include<random>	 // random number generator
#include<ctime>		 
#include<cmath>
#include<climits>
#include<iostream>
#include<fstream>
#include<string>
#include <cstdlib>
#include <algorithm>   
#include <unistd.h> // to get process id, for seeding RNG

using namespace std;

// Parameters
const int		N = 10000;			// Population size. 
const int		tmax = 1000000000;  // Number of iterations
const int		tstep = 10000000;	// Time step for printing results. 
const double	sdv = 0.000025;    	// Standard deviation for mutation Gaussian distribution.
const double    alpha = 0.5;        // Unequal games

// Game parameters
const double    a_x = 1.0;
const double	b_x	= 1.05;
const double	c_x = 0.9;
const double	d_x = 1.65;
	// ESS x = 0.081300813; 
	
// Structure of a player
struct player {
	double xstrategy;
	double ystrategy;
	int ID;
};

// Payoff functions
double Payoff_SD_x(player A, player B){
	double sumx = A.xstrategy + B.xstrategy;
	double Benefits = a_x * (sumx - b_x * sumx * sumx);
	double Costs = c_x * (A.xstrategy - d_x * A.xstrategy * A.xstrategy);
	double P = Benefits - Costs;
	return P;
}

double Payoff_SD_y(player A, player B, double a_y, double b_y, double c_y, double d_y){
	double sumy = A.ystrategy + B.ystrategy;
	double Benefits = a_y * (sumy - b_y * sumy * sumy);
	double Costs = c_y * (A.ystrategy - d_y * A.ystrategy * A.ystrategy);
	double P = Benefits - Costs;
	return P;
}


int main()
{
	// Initiate seed with process id
	const unsigned int seed = getpid();
	
	// Create and seed an engine based on the Mersenne twister:
	mt19937_64 rng(seed);
	
	// Uniform distribution from which to pick random parameter values:
	uniform_real_distribution<> rnd_parameter(0.0, 5.0);

	// Declare parameter values
    double ab = 1.05;
    double cd = 1.485;
    double a_y;
    double c_y;
    double b_y;
    double d_y;
    double eq_y = -1.0;
    double max_y = -2.0;
    bool cs_y = false;
    bool ebp_y = false;
    
	// Sample parameter values for game 2
	while(eq_y < 0 || eq_y > max_y || cs_y == false || ebp_y == false || eq_y < 0.016)
    {
        a_y = rnd_parameter(rng);    
        c_y = rnd_parameter(rng);
        
        normal_distribution<> rnd_b{ab/a_y, 0.1*ab/a_y};
        b_y = rnd_b(rng);
        
        normal_distribution<> rnd_d{cd/c_y, 0.1*cd/c_y};
        d_y = rnd_d(rng);
        
        eq_y = (a_y - c_y)/(4*a_y*b_y - 2*c_y*d_y);
        max_y = min(1/(4*b_y), 1/(2*d_y));
        
        cs_y = -4*a_y*b_y + 2*c_y*d_y < 0.0;
        ebp_y = -2*a_y*b_y + 2*c_y*d_y > 0.0;
    }

	// Determine initial conditions
	double x_0 = 0.081300813 + 0.015*cos(2.35619);
	double y_0 = eq_y + 0.015*sin(2.35619);

	// Distributions
	uniform_int_distribution<> rnd_individual(0, N-1); // use to pick one individual from the population
	
	// Create population
	player population[N];
	
	normal_distribution<> rnd_xstart{x_0, sdv};
	normal_distribution<> rnd_ystart{y_0, sdv};
	for (int i = 0; i < N; i++){
		population[i].xstrategy = rnd_xstart(rng);
		population[i].ystrategy = rnd_ystart(rng);
		population[i].ID = i;
	}
	
	// Iteration
	int t = 0; while (t < tmax)
	{
		// Assemble first group and calculate fitness (W) of focal (f) individual
		int idF = rnd_individual(rng); int idF2 = rnd_individual(rng);
		while (idF2 == idF) idF2 = rnd_individual(rng);
		player F = population[idF]; player F2 = population[idF2];
		double fxpayoff = Payoff_SD_x(F, F2);
		double fypayoff = Payoff_SD_y(F, F2, a_y, b_y, c_y, d_y);
		double Wf = (1-alpha)*fxpayoff + alpha*fypayoff;

		// Assemble second group and calculate fitness (W) of rival (r) individual (first individual of group)
		int idR = rnd_individual(rng); int idR2 = rnd_individual(rng);
		while (idR2 == idR) idR2 = rnd_individual(rng);
		player R = population[idR]; player R2 = population[idR2];
		double rxpayoff = Payoff_SD_x(R, R2);
		double rypayoff = Payoff_SD_y(R, R2, a_y, b_y, c_y, d_y);
		double Wr = (1-alpha)*rxpayoff + alpha*rypayoff;

		// If Wr > Wf, offspring has r as a parent
		if (Wr > Wf) {
			population[F.ID].xstrategy = R.xstrategy;
			population[F.ID].ystrategy = R.ystrategy;
		}

		// Mutation
		normal_distribution<> rnd_x{population[F.ID].xstrategy, sdv};
		normal_distribution<> rnd_y{population[F.ID].ystrategy, sdv};
		double newx = rnd_x(rng);
		double newy = rnd_y(rng);
		
		// B(2x) has a maximum at x = 0.238 and C(x) has a maximum at x = 0.303. We want to avoid mutations larger than these values
        // B(2y) has a maximum at max_y. 
		bool end = 0; 	// We are also stopping the simulation when any individual reaches the edges of parameter space.
		if (newx > 0.238) {newx = 0.238; end = 1;} 
		else if (newx < 0) {newx = 0; end = 1;}
		if (newy > max_y) {newy = max_y; end = 1;}
		else if (newy < 0) {newy = 0; end = 1;}

        // cont. mutation
		population[F.ID].xstrategy = newx;
		population[F.ID].ystrategy = newy;

		// Write output
		if (t == 0) {
			cout << a_y << ", " << b_y << ", " << c_y << ", " << d_y; 
		}

		if (t % tstep == 0 || end == 1) {
			for (int i = 0; i < N; i++) cout << "\n" << t << ", " << population[i].xstrategy << ", " << population[i].ystrategy << ", NA" ;
		}

		// Iterate
		if(end == 0) t++;
		else t = tmax; // if the population approaches the edges, end the simulation
	}

	return 0;
} // close main; end
