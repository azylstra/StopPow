/** Test fitting functions
 * @author Alex Zylstra
 * @date 2014/04/11
 */

#include <iostream>
#include <vector>
#include <array>
#include <ctime>

#include "StopPow.h"
#include "StopPow_SRIM.h"
#include "StopPow_Fit.h"
#include "Fit.h"
#include "Util.h"

int main(int argc, char const *argv[])
{
	// check for verbosity flag:
	bool verbose = false;
	if( argc >= 2 )
	{
		for(int i=1; i < argc; i++)
		{
			std::string flag(argv[i]);
			if(flag == "--verbose")
			{
				verbose = true;
			}
		}
	}

	// Do some output
	std::cout << "========== Test Suite 7 ==========" << std::endl;
	std::cout << "      Test StopPow_Fit  " << std::endl;
	bool pass = true;
	bool test = true;

	std::vector<double> mf {26.98};
	std::vector<double> Zf {13};
	std::vector<double> Tf {1.0};
	std::vector<double> nf {6.02e22};
	std::vector<double> Zbar {7};
	
	StopPow::StopPow_Fit * s = new StopPow::StopPow_Fit(1, 1, mf, Zf, Tf, nf, Zbar, 1);
	// test factor:
	std::cout << "testing free-electron factor..." << std::endl;
	test &= StopPow::approx(s->dEdx(10), -0.00999, 1e-3);
	s->set_factor(2);
	test &= StopPow::approx(s->dEdx(10), -0.0164, 1e-3);
	s->set_factor(1);

	// testing be normalization
	std::cout << " testing bound-electron normalization..." << std::endl;
	test &= StopPow::approx(s->dEdx(15), -0.007302, 1e-3);
	StopPow::StopPow_SRIM * srim = new StopPow::StopPow_SRIM("SRIM/Hydrogen in Aluminum.txt");
	s->normalize_bound_e(srim, 15);
	test &= StopPow::approx(s->dEdx(15), -0.00727, 1e-3);

	// test different fe models
	std::cout << "Changing free-electron models..." << std::endl;
	s->choose_model(s->MODE_LP);
	test &= StopPow::approx(s->dEdx(15), -0.007485, 1e-3);
	s->choose_model(s->MODE_BPS);
	test &= StopPow::approx(s->dEdx(15), -0.007282, 1e-3);
	s->choose_model(s->MODE_ZIMMERMAN);
	test &= StopPow::approx(s->dEdx(15), -0.007269, 1e-3);
	std::cout << "Three tests: " << (test ? "pass" : "FAIL!") << std::endl;
	pass &= test;


	// Final result:
	if(pass)
	{
		std::cout << "PASS" << std::endl;
		return 0;
	}
	std::cout << "FAIL!" << std::endl;
	return 1;
}