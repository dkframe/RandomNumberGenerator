// RandomNumberGenerator.cpp : Defines the entry point for the application.
//

#include "RandomNumberGenerator.h"
#include <functional>
#include <random>
#include <chrono>
#include <cmath>
#include <vector>
#include <set>
#include <iterator>
#include <fstream>

// --------------------------------------------------------------------
const double root2 = std::sqrt(2.0);
const double pi = 3.14159265358979323846;

// --------------------------------------------------------------------
void myPause() 
{
	std::cout << std::endl << "Press any key to continue...";
	std::cin.get();
	std::cout << std::endl;
}
// --------------------------------------------------------------------
namespace randomGenerator 
{
	typedef std::set<std::pair<double, double>> cdfTable;
	typedef std::pair<double, double> Interval1D;
	typedef std::pair<Interval1D, Interval1D> Interval2D;

	double myCDF_Gaussian(double x)
	{
		// Defines the CDF for the unit Normal distribution (mean = 0; variance = 1)
		return 0.5 * (1 + std::erf(x / root2));
	}
	double myCDF_Cauchy(double x)
	{
		// Defines the CDF for the Cauchy distribution (x0 = -2, gamma = 1)
		return (std::atan(x + 2) / pi) + 0.5;
	}

	class myRNG
	{
	public:
		myRNG(double (*cdf_func)(double), int type, double xMin, double xMax, double tolerance, int max_iterations)
		{
			CDF = cdf_func;
			searchType = type;
			x_min = xMin;
			x_max = xMax;
			y_min = cdf_func(xMin);
			y_max = cdf_func(xMax);
			tol = tolerance;
			max_iters = max_iterations;

			lookupTable = initializeCDFTable();

			// Initialize Uniform Random Number Generator
			seed = static_cast<unsigned int>(std::time(NULL));
			generator.seed(seed);
			std::uniform_real_distribution<double> distribution(0, 1);

		}
		myRNG(double (*cdf_func)(double), double xMin, double xMax, double tolerance, int max_iterations)
		{
			CDF = cdf_func;
			searchType = 0;
			x_min = xMin;
			x_max = xMax;
			y_min = cdf_func(xMin);
			y_max = cdf_func(xMax);
			tol = tolerance;
			max_iters = max_iterations;

			lookupTable = initializeCDFTable();

			// Initialize Uniform Random Number Generator
			seed = static_cast<unsigned int>(std::time(NULL));
			generator.seed(seed);
			std::uniform_real_distribution<double> distribution(0, 1);

		}
		// -------------------------------------------------------------------
		void printCDFTable() 
		{
			for (auto &i : lookupTable)
			{
				std::cout << i.first << " " << i.second << std::endl;
			}
			std::cout << std::endl;
		}
		// -------------------------------------------------------------------
		double generateRandomNumber() 
		{
			// Pick a random uniform number y0
			double y0 = uniformRNG();

			// Find the interval containing y0
			Interval2D interval = lookupCDFTable(y0);

			// Search that interval for the x where the |CDF(x) - y| < tol.
			switch (searchType) {
				case 0:
					return searchBisection(y0, interval);
					break;
				case 1: 
					return searchITP(y0, interval);
					break;
				default:
					return searchBisection(y0, interval);
					break;
			}

		};
		// -------------------------------------------------------------------
	private:
		// Specified uniform random number generator
		std::default_random_engine generator;
		std::uniform_real_distribution<double> distribution;
		unsigned int seed; // seed for the uniform random number generator

		int searchType; // Type of minimization algorithm to use

		double (*CDF)(double); // Specified CDF for the generator
		cdfTable lookupTable; // Holds values for CDF at various points
		double x_min; // Left-boundary of domain of CDF
		double x_max; // Right-boundary of domain of CDF
		double y_min; // CDF(x_min): Should be almost exactly 0
		double y_max; // CDF(x_max): Should be almost exactly 1
		double tol; // Error tolerance for the minimization routines
		int max_iters; // Number of iterations required for minimization to stop prematurely

		double uniformRNG()
		{
			return distribution(generator);
		}
		cdfTable initializeCDFTable()
		{
			// Initialize a set of initial intervals for the root finding algorithm, for a specified CDF
			// Each point in the table will be separated by powers of 10, in the range between
			// x_min and x_max
			//
			// We will assume the interval [x_min, x_max] is at most either [0, 1.0E+30] or [-1.0E+30, 1.0E+30]

			cdfTable table;

			int logx_max = (int) std::log2(std::abs(x_max));
			double x;

			table.insert({ x_min, CDF(x_min) });
			table.insert({ x_max, CDF(x_max) });

			for (int i = -3; i < logx_max; i++)
			{
				x = std::pow(2, (double)i);
				if (x < x_max) table.insert({ x, CDF(x) });
				if (-x > x_min) table.insert({ -x, CDF(-x) });
			}

			return table;
		}
		Interval2D lookupCDFTable(const double y)
		{
			// Returns an interval containing the point y
			double y1, y2, x1, x2;
			bool lookupSuccess = false;

			y1 = 0;
			y2 = 0;
			x1 = 0;
			x2 = 0;

			for (auto &i : lookupTable)
			{
				y1 = y2;
				x1 = x2;
				y2 = i.second;
				x2 = i.first;
				if (y2 >= y && y1 <= y) 
				{ 
					lookupSuccess = true;
					break; 
				}
			}

			if (lookupSuccess) 
			{
				// IDEAL CASE: Function was able to find an interval containing y
				return { { x1, x2 }, {y1, y2} };
			}
			else 
			{
				// EDGE CASES: CDF(y) > y_max or CDF(y) < y_min
				// In these cases, we have to find a new interval beyond x_min or x_max
				// Hopefully, x_min and x_max are chosen so this section is never used.
				if (CDF(y) > y_max) {
					x1 = x_max;
					y1 = y_max;

					x2 = x1;
					y2 = CDF(x2);

					while (y2 < y) 
					{
						x2 *= 10.0;
						y2 = CDF(x2);
					}
					return { { x1, x2 }, {y1, y2} };
				}
				if (CDF(y) < y_min) {
					x1 = x_min;
					y1 = y_min;

					x2 = x1;
					y2 = CDF(x2);

					while (y2 < y)
					{
						x2 *= 10.0;
						y2 = CDF(x2);
					}
					return { { x1, x2 }, {y1, y2} };
				}
			}

			// Default case: should never happen
			return { {0,0}, {0,0} };
		}

		// Routines to generate a random number from a distribution with given CDF
		// We do this by generating a uniform random number in the range [ 0, 1 ),
		// then "inverting" the CDF to find a value x such that |CDF(x) - x0| < tol
		//
		// We provide 3 different minimization routines: one which uses the bisection method (type 0),
		// one which uses Brent's method (type 1), and one which uses the ITP method (type 2). 
		// 
		// If no search type is provided, the generator will default to the bisection method.
		
		double searchBisection (const double y0, const Interval2D xy_init)
		{
			
			// Set initial intervals
			double x_left = xy_init.first.first;
			double x_right = xy_init.first.second;
			double y_left = xy_init.second.first - y0;
			double y_right = xy_init.second.second - y0;
			double x, y;
			x = 0.0;
			y = 0.0;

			// Test to make sure the interval is valid
			if (!((y_left < 0) && (y_right > 0)))
			{
				std::cout << "Invalid interval: ";
				std::cout << y0 << " not in ( " << y_left << ", " << y_right << ")! " << std::endl;
				return 0.0;
			}

			// Start bisecting
			for (int i = 0; i < max_iters; i++)
			{
				x = (x_left + x_right) / 2.0; // Pick the midpoint of the interval
				y = CDF(x); // Evaluate the CDF at this point
			
				// Select which endpoint to replace
				if (y-y0 < 0) 
				{
					x_left = x;
				} else if (y-y0 > 0) 
				{
					x_right = x;
				}

				// Convergence criteria
				if (abs(y - y0) < tol) break;
			}

			// Return x that satisfies |CDF(x) - y0| < tol
			return x;

		}
		double searchITP (const double y0, const Interval2D xy_init)
		{

			// Set initial intervals
			double x1 = xy_init.first.first;
			double x2= xy_init.first.second;
			double y1 = xy_init.second.first - y0;
			double y2 = xy_init.second.second - y0;
			
			double x_mid = 0.0;
			double xf = 0.0;
			double xt = 0.0;

			double r = 0.0;
			double delta = 0.0;
			double sigma = 0.0;

			double x_ITP = 0.0;
			double y_ITP = 0.0;

			// Set ITP specific parameters
			double k1 = 0.2 / (x2 - x1);
			double k2 = 2.0;
			int n0 = 1;

			// Test to make sure the interval is valid
			if (!((y1 < 0) && (y2 > 0)))
			{
				std::cout << "Invalid interval:";
				std::cout << y0 << " not in ( " << y1 << ", " << y2 << ")! " << std::endl;
				return 0.0;
			}

			// Start ITP Method

			int n_half = (int) std::ceil(std::log2((x2 - x1) / (2 * tol)));
			int n_max = n_half + n0;
			for (int i = 0; i < n_max; i++)
			{
				// Calculate parameters
				x_mid = 0.5 * (x1 + x2);
				r = (tol * std::pow(2, n_max - i)) - (0.5 * (x2 - x1));
				delta = k1 * std::pow(x2 - x1, k2);

				// Interpolation
				xf = (y2 * x1 - y1 * x2) / (y2 - y1);

				// Truncation
				sigma = (x_mid - xf) / std::abs(x_mid - xf);
				if (delta < std::abs(x_mid - xf)) 
				{
					xt = xf + delta * sigma;
				} 
				else 
				{
					xt = x_mid;
				}

				// Projection
				if (std::abs(xt - x_mid) <= r) 
				{
					x_ITP = xt;
				}
				else
				{
					x_ITP = x_mid - sigma * r;
				}

				// Update Interval
				y_ITP = CDF(x_ITP) - y0;
				if (y_ITP > 0) 
				{
					x2 = x_ITP;
					y2 = y_ITP;
				}
				else if (y_ITP < 0)
				{
					x1 = x_ITP;
					y1 = y_ITP;
				}
				else 
				{
					x1 = x_ITP;
					x2 = x_ITP;
				}

				// Check convergence criteria
				if ((x2 - x1) < 2 * tol) break;
			}

			// Return x that satisfies |CDF(x) - y0| < tol
			return 0.5 * (x1 + x2);

		}
	};
}
// --------------------------------------------------------------------

int main()
{
	// Start timer at beginning of program
	typedef std::chrono::high_resolution_clock myclock;
	myclock::time_point t0 = myclock::now();

	// Open file for generator outputs
	std::ofstream fout1("STL_Gaussian.txt");
	std::ofstream fout2("STL_Cauchy.txt");
	std::ofstream fout3("myRNG_Gaussian_Bisection.txt");
	std::ofstream fout4("myRNG_Cauchy_Bisection.txt");
	std::ofstream fout5("myRNG_Gaussian_ITP.txt");
	std::ofstream fout6("myRNG_Cauchy_ITP.txt");
	myclock::time_point t1 = myclock::now();

	// Set parameters for our Random Number Generator
	double x_min = -1.0E+06;
	double x_max = 1.0E+06;
	double tolerance = 1.0E-07;
	int max_iterations = 10000;
	int number_of_samples = 100000;
	double rand_num = 0.0;
	bool fileIO = true;
	myclock::time_point t2 = myclock::now();

	// Initialize the STL Random Number Generators
	std::default_random_engine generatorSTL;
	std::cauchy_distribution<double> distributionSTL_cauchy(-2.0, 1.0);
	std::normal_distribution<double> distributionSTL_gaussian(0.0, 1.0);
	myclock::time_point t3 = myclock::now();

	// Test STL Generator (Gaussian)
	for (int i = 0; i < number_of_samples; i++) 
	{
		if (fileIO) {
			fout1 << distributionSTL_gaussian(generatorSTL) << std::endl;
		}
		else {
			rand_num = distributionSTL_gaussian(generatorSTL);
		}
	}
	myclock::time_point t4 = myclock::now();

	// Test STL Generator (Cauchy)
	for (int i = 0; i < number_of_samples; i++) 
	{
		if (fileIO) 
		{
			fout2 << distributionSTL_cauchy(generatorSTL) << std::endl;
		}
		else 
		{
			rand_num = distributionSTL_cauchy(generatorSTL);
		}
	}
	myclock::time_point t5 = myclock::now();

	// Initialize our Random Number Generator (Bisection search, Gaussian)
	randomGenerator::myRNG rngGaussianBisection(randomGenerator::myCDF_Gaussian, 0, x_min, x_max, tolerance, max_iterations);
	myclock::time_point t6 = myclock::now();

	// Test Generator (Bisection search, Gaussian)
	for (int i = 0; i<number_of_samples; i++) 
	{
		if (fileIO)
		{
			fout3 << rngGaussianBisection.generateRandomNumber() << std::endl;
		}
		else
		{
			rand_num = rngGaussianBisection.generateRandomNumber();
		}
	}
	myclock::time_point t7 = myclock::now();

	// Initialize our Random Number Generator (Bisection search, Cauchy)
	randomGenerator::myRNG rngCauchyBisection(randomGenerator::myCDF_Cauchy, 0, x_min, x_max, tolerance, max_iterations);
	myclock::time_point t8 = myclock::now();

	// Test Generator (Bisection search, Cauchy)
	for (int i = 0; i < number_of_samples; i++)
	{
		if (fileIO)
		{
			fout4 << rngCauchyBisection.generateRandomNumber() << std::endl;
		}
		else
		{
			rand_num = rngCauchyBisection.generateRandomNumber();
		}
	}
	myclock::time_point t9 = myclock::now();

	// Initialize our Random Number Generator (ITP search, Gaussian)
	randomGenerator::myRNG rngGaussianITP(randomGenerator::myCDF_Gaussian, 1, x_min, x_max, tolerance, max_iterations);
	myclock::time_point t10 = myclock::now();

	// Test Generator (ITP search, Gaussian)
	for (int i = 0; i < number_of_samples; i++)
	{
		if (fileIO)
		{
			fout5 << rngGaussianITP.generateRandomNumber() << std::endl;
		}
		else
		{
			rand_num = rngGaussianITP.generateRandomNumber();
		}
	}
	myclock::time_point t11 = myclock::now();

	// Initialize our Random Number Generator (ITP search, Cauchy)
	randomGenerator::myRNG rngCauchyITP(randomGenerator::myCDF_Cauchy, 1, x_min, x_max, tolerance, max_iterations);
	myclock::time_point t12 = myclock::now();

	// Test Generator (ITP search, Cauchy)
	for (int i = 0; i < number_of_samples; i++)
	{
		if (fileIO)
		{
			fout6 << rngCauchyITP.generateRandomNumber() << std::endl;
		}
		else
		{
			rand_num = rngCauchyITP.generateRandomNumber();
		}
	}
	myclock::time_point t13 = myclock::now();

	// -------------------------------------------------------------------
	//	Output results to terminal
	// -------------------------------------------------------------------

	std::chrono::duration<double, std::micro> d1 = t3 - t2;
	std::chrono::duration<double, std::micro> d2 = t4 - t3;
	std::chrono::duration<double, std::micro> d3 = t5 - t4;
	std::cout << "-------------------- STL Generator (Gaussian) -------------------------" << std::endl;
	std::cout << std::endl;
	std::cout << "Time to initialize generator: " << d1 << std::endl;
	std::cout << "Time to generate " << number_of_samples << " samples: " << d2 << std::endl;
	std::cout << "Average time per sample: " << d2 / ((double)number_of_samples) << std::endl;
	std::cout << std::endl;
	std::cout << "-------------------- STL Generator (Cauchy) ---------------------------" << std::endl;
	std::cout << std::endl;
	std::cout << "Time to initialize generator: " << d1 << std::endl;
	std::cout << "Time to generate " << number_of_samples << " samples: " << d3 << std::endl;
	std::cout << "Average time per sample: " << d3 / ((double)number_of_samples) << std::endl;
	std::cout << std::endl;

	std::chrono::duration<double, std::micro> d4 = t6 - t5;
	std::chrono::duration<double, std::micro> d5 = t7 - t6;
	std::cout << "-------------------- myRNG (Gaussian, Bisection) ----------------------" << std::endl;
	std::cout << std::endl;
	std::cout << "Time to initialize generator: " << d4 << std::endl;
	std::cout << "Time to generate " << number_of_samples << " samples: " << d5 << std::endl;
	std::cout << "Average time per sample: " << d5 / ((double)number_of_samples) << std::endl;
	std::cout << std::endl;

	std::chrono::duration<double, std::micro> d6 = t8 - t7;
	std::chrono::duration<double, std::micro> d7 = t9 - t8;
	std::cout << "-------------------- myRNG (Cauchy, Bisection) ----------------------" << std::endl;
	std::cout << std::endl;
	std::cout << "Time to initialize generator: " << d6 << std::endl;
	std::cout << "Time to generate " << number_of_samples << " samples: " << d7 << std::endl;
	std::cout << "Average time per sample: " << d7 / ((double)number_of_samples) << std::endl;
	std::cout << std::endl;

	std::chrono::duration<double, std::micro> d8 = t10 - t9;
	std::chrono::duration<double, std::micro> d9 = t11 - t10;
	std::cout << "-------------------- myRNG (Gaussian, ITP) ----------------------" << std::endl;
	std::cout << std::endl;
	std::cout << "Time to initialize generator: " << d8 << std::endl;
	std::cout << "Time to generate " << number_of_samples << " samples: " << d9 << std::endl;
	std::cout << "Average time per sample: " << d9 / ((double)number_of_samples) << std::endl;
	std::cout << std::endl;

	std::chrono::duration<double, std::micro> d10 = t12 - t11;
	std::chrono::duration<double, std::micro> d11 = t13 - t12;
	std::cout << "-------------------- myRNG (Cauchy, ITP) ----------------------" << std::endl;
	std::cout << std::endl;
	std::cout << "Time to initialize generator: " << d10 << std::endl;
	std::cout << "Time to generate " << number_of_samples << " samples: " << d11 << std::endl;
	std::cout << "Average time per sample: " << d11 / ((double)number_of_samples) << std::endl;
	std::cout << std::endl;

	// Pause before ending code
	fout1.close();
	fout2.close();
	fout3.close();
	fout4.close();
	fout5.close();
	fout6.close();

	return 0;
}
