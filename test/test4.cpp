/** test class for StopPow library to test regular dE/dx call evaluation speed
 * @author Alex Zylstra
 * @date 2014/04/08
 */


#include <stdio.h>
#include <dirent.h>

#include <iostream>
#include <vector>
#include <stdexcept>
#include <limits>
#include <ctime>
#include <random>
#include <thread>

#include "StopPow.h"
#include "StopPow_SRIM.h"
#include "StopPow_LP.h"
#include "StopPow_BetheBloch.h"
#include "StopPow_AZ.h"
#include "StopPow_Mehlhorn.h"
#include "StopPow_Zimmerman.h"
#include "StopPow_BPS.h"
#include <gsl/gsl_integration.h>

bool run_tests(std::vector<StopPow::StopPow*> models)
{
	int n = 1000;
    std::random_device rd;
    std::default_random_engine e1(rd());
    std::uniform_int_distribution<double> uniform_dist(0, 1);

	for(StopPow::StopPow* s : models)
	{
		std::clock_t start;
		double duration, randomE;
		start = std::clock();
		for(int i=0; i<n; i++)
		{
			randomE = s->get_Emin() + (s->get_Emax() - s->get_Emin()) * uniform_dist(e1);
			s->dEdx(randomE);
		}
		duration = (1e6/n)*(std::clock()-start) / (double) CLOCKS_PER_SEC;
		std::cout << "dE/dx call = " << duration << " us for model " << s->get_type() << std::endl;
	}

	return true;
}
int main(int argc, char* argv [])
{
	// Do some output
	std::cout << "========== Test Suite 4 ==========" << std::endl;
	std::cout << "   Testing dE/dx speed  " << std::endl;

	// create a variety of stopping power models to test
	std::vector<StopPow::StopPow*> models;



	// ---------------------------------------
	// 		Set up models
	// ---------------------------------------
	// create model and add to the vector:
	StopPow::StopPow* s = new StopPow::StopPow_SRIM("test0/Hydrogen in Aluminum.txt");
	models.push_back(s);
	s = new StopPow::StopPow_BetheBloch(1, 1, {26.98}, {13.0}, {6.03e22});
	models.push_back(s);
	s = new StopPow::StopPow_AZ(13);
	models.push_back(s);

	// solid density Al at 1 keV
	std::vector<double> mf {26.98};
	std::vector<double> Zf {13};
	std::vector<double> Tf {1.0};
	std::vector<double> nf {6.02e22};
	s = new StopPow::StopPow_LP(1, 1, mf, Zf, Tf, nf, 1.);
	models.push_back(s);
	std::vector<double> Zbar {7};
	s = new StopPow::StopPow_Mehlhorn(1, 1, mf, Zf, Tf, nf, Zbar, 1.);
	models.push_back(s);
	s = new StopPow::StopPow_Zimmerman(1, 1, mf, Zf, Tf, nf, Zbar, 1.);
	models.push_back(s);
	StopPow::StopPow_BPS * s2 = new StopPow::StopPow_BPS(1, 1, mf, Zf, Tf, nf, 1.);
	models.push_back(s2);

	// ---------------------------------------
	// 		Run Tests
	// ---------------------------------------
	run_tests(models);

	// Some special stuff for BPS
	int n = 1000;
    std::random_device rd;
    std::default_random_engine e1(rd());
    std::uniform_int_distribution<double> uniform_dist(0, 1);
	std::clock_t start;
	double duration, randomE;
	std::cout << "For BPS:" << std::endl;
	
	start = std::clock();
	for(int i=0; i<n; i++)
	{
		randomE = s2->get_Emin() + (s2->get_Emax() - s2->get_Emin()) * uniform_dist(e1);
		s2->dEdx_short(randomE);
	}
	duration = (1e6/n)*(std::clock()-start) / (double) CLOCKS_PER_SEC;
	std::cout << "dE/dx_short call = " << duration << " us" << std::endl;
	
	start = std::clock();
	for(int i=0; i<n; i++)
	{
		randomE = s2->get_Emin() + (s2->get_Emax() - s2->get_Emin()) * uniform_dist(e1);
		s2->dEdx_long(randomE);
	}
	duration = (1e6/n)*(std::clock()-start) / (double) CLOCKS_PER_SEC;
	std::cout << "dE/dx_long call = " << duration << " us" << std::endl;
	
	start = std::clock();
	for(int i=0; i<n; i++)
	{
		randomE = s2->get_Emin() + (s2->get_Emax() - s2->get_Emin()) * uniform_dist(e1);
		s2->dEdx_quantum(randomE);
	}
	duration = (1e6/n)*(std::clock()-start) / (double) CLOCKS_PER_SEC;
	std::cout << "dE/dx_quantum call = " << duration << " us" << std::endl;

	return 0;
}