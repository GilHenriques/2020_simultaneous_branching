#include<vector>
#include<random>	 // random number generator
#include<ctime>		 // time for seeding RNG
#include<cmath>
#include<climits>
#include<iostream>
#include<fstream>
#include<string>
#include <cstdlib>
#include <algorithm>   

using namespace std;

// Parameters
const int		N = 10000;			// Population size. 
const int		tmax = 1000000000;  // Number of iterations
const int		tstep = 10000000;	// Time step for printing results. 
const double	sdv = 0.000025;    	// Standard deviation for mutation Gaussian distribution.
const double    alpha = 0.5;        // Unequal games
double			x_0 = 0.081300813 + 0.03*cos(3.14159265359);
double			y_0 = 0.081300813 + 0.03*sin(3.14159265359);

// Game parameters
const double	b_x	= 1.05;
const double	b_y = 1.05;
const double	c_x = 0.9;
const double	c_y = 0.9;
const double	d_x = 1.65;
const double	d_y = 1.65;
	// ESS = - (c - 1)/(2*(2*b-c*d)) = 0.081300813
	
// Structure of a player
struct player {
	double xstrategy;
	double ystrategy;
	int ID;
};

// Payoff functions
double Payoff_SD_x(player A, player B){
	double sumx = A.xstrategy + B.xstrategy;
	double Benefits = sumx - b_x * sumx * sumx;
	double Costs = c_x * (A.xstrategy - d_x * A.xstrategy * A.xstrategy);
	double P = Benefits - Costs;
	return P;
}

double Payoff_SD_y(player A, player B){
	double sumy = A.ystrategy + B.ystrategy;
	double Benefits = sumy - b_y * sumy * sumy;
	double Costs = c_y * (A.ystrategy - d_y * A.ystrategy * A.ystrategy);
	double P = Benefits - Costs;
	return P;
}

int main()
{
	// Initiate seed with time and create name with seed
	const unsigned int seed = time(0);
	ofstream file;
	file.open(to_string(seed) + "_theta180_alpha0.5.txt");
	cout << "\n STARTED";
	
	// Create and seed an engine based on the Mersenne twister:
	mt19937_64 rng(seed);
	
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
		double fypayoff = Payoff_SD_y(F, F2);
		double Wf = (1-alpha)*fxpayoff + alpha*fypayoff;

		// Assemble second group and calculate fitness (W) of rival (r) individual (first individual of group)
		int idR = rnd_individual(rng); int idR2 = rnd_individual(rng);
		while (idR2 == idR) idR2 = rnd_individual(rng);
		player R = population[idR]; player R2 = population[idR2];
		double rxpayoff = Payoff_SD_x(R, R2);
		double rypayoff = Payoff_SD_y(R, R2);
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
		bool end = 0; 	// We are also stopping the simulation when any individual reaches the edges of parameter space.
		if (newx > 0.238) {newx = 0.238; end = 1;} 
		//if (abs(newx) > 0.081300813+0.03+0.005) {end = 1;} 
		else if (newx < 0) {newx = 0; end = 1;}
		if (newy > 0.238) {newy = 0.238; end = 1;}
		//if (abs(newy) > 0.081300813+0.03+0.005) {end = 1;} 
		else if (newy < 0) {newy = 0; end = 1;}

        // cont. mutation
		population[F.ID].xstrategy = newx;
		population[F.ID].ystrategy = newy;

		// Write output
		if (t % tstep == 0 | end == 1) {
			for (int i = 0; i < N; i++) file << "\n" << t << ", " << population[i].xstrategy << ", " << population[i].ystrategy;
			cout << "\n Iteration " << t << " complete!";
		}

		// Iterate
		if(end == 0) t++;
		else t = tmax; // if the population approaches the edges, end the simulation
	}

	cout << "\n FINISHED";
	file.close();   // Close the file pointer and flushing all pending IO operations to it.

	return 0;
} // close main; end
