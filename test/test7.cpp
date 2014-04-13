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

	// ------------ Test fitting routine ----------------
	std::cout << "testing fit routine" << std::endl;
	std::vector<double> data_x {1., 1.2, 1.4, 1.6, 1.8, 2., 2.2, 2.4, 2.6, 2.8, 3., 3.2, 3.4, 3.6, \
3.8, 4., 4.2, 4.4, 4.6, 4.8, 5., 5.2, 5.4, 5.6, 5.8, 6., 6.2, 6.4, \
6.6, 6.8, 7., 7.2, 7.4, 7.6, 7.8, 8., 8.2, 8.4, 8.6, 8.8, 9., 9.2, \
9.4, 9.6, 9.8, 10., 10.2, 10.4, 10.6, 10.8, 11., 11.2, 11.4, 11.6, \
11.8, 12., 12.2, 12.4, 12.6, 12.8, 13., 13.2, 13.4, 13.6, 13.8, 14., \
14.2, 14.4, 14.6, 14.8, 15., 15.2, 15.4, 15.6, 15.8, 16., 16.2, 16.4, \
16.6, 16.8, 17., 17.2, 17.4, 17.6, 17.8, 18., 18.2, 18.4, 18.6, 18.8, \
19., 19.2, 19.4, 19.6, 19.8, 20.};
	std::vector<double> data_y  {75259.5, 143190., 107732., -135826., 68668.8, -45779.4, -31288., \
-7906.47, -154645., -94217.2, 63672.4, -19471.6, -138668., 5475.73, \
-25031.5, 113698., 41236.6, 25076.1, -92551.9, -55113., 27626.8, \
10262.1, 164479., 152760., 52137.5, 30604.5, 45881.6, -88132.1, \
-2959.76, -55656.5, 142751., 31965.4, 249920., 157223., 228101., \
619570., 677877., 1.07328e6, 1.56837e6, 1.99958e6, 
 2.46877e6, 2.98764e6, 3.41332e6, 3.73934e6, 3.78872e6,
  3.89069e6, 3.86356e6, 3.5932e6, 3.44923e6, 2.98537e6,
  2.61213e6, 1.90998e6, 1.39943e6, 
 1.13471e6, 753923., 670349., 339406., 408169., 189216., 231251., \
19903.4, 25750.8, 49745.1, -141700., 144601., -132210., 258154., \
-186506., 95274.4, 50575.5, 3431.72, 133939., 15490., 20435., \
24945.1, -43846.3, -28191.4, -101697., -176679., 51158., -126813., \
-102399., -130652., -164509., 22325.3, 137000., -126873., 4833.01, \
-119890., -47049.6, 83393.2, 109827., -5619.18, 179639., 21938.5, \
-75947.6};
	std::vector<double> data_std {100000, 100000, 100000, 100000, 100000, 100000, 100000, 100000, \
100000, 100000, 100000, 100000, 100000, 100000, 100000, 100000, \
100000, 100000, 100000, 100000, 100000, 100000, 100000, 100000, \
100000, 100000, 100000, 100000, 100000, 100000, 100000, 100000, \
100000, 100000, 100000, 100000, 100000, 100000, 100000, 100000, \
100000, 100000, 100000, 100000, 100000, 100000, 100000, 100000, \
100000, 100000, 100000, 100000, 100000, 100000, 100000, 100000, \
100000, 100000, 100000, 100000, 100000, 100000, 100000, 100000, \
100000, 100000, 100000, 100000, 100000, 100000, 100000, 100000, \
100000, 100000, 100000, 100000, 100000, 100000, 100000, 100000, \
100000, 100000, 100000, 100000, 100000, 100000, 100000, 100000, \
100000, 100000, 100000, 100000, 100000, 100000, 100000, 100000};

	// Run:
	double dE = 0.05;
	double E0 = 14.7;
	double E0_unc = 0.05;
	double sigma0 = 0.75;
	double sigma0_unc = 0.1;
	double rhoR = 160;
	double rhoR_unc = 10;
	double chi2;
	std::vector<double> fit, fit_unc;
	forward_fit_dEdx(data_x, data_y, data_std, dE, E0, E0_unc, sigma0, sigma0_unc, rhoR, rhoR_unc, *s, chi2, fit, fit_unc, verbose);
	// print results:
	std::cout << "factor = " << fit[0] << " +/- " << fit_unc[0] << std::endl;
	std::cout << "A = " << fit[1] << " +/- " << fit_unc[1] << std::endl;

	// Final result:
	if(pass)
	{
		std::cout << "PASS" << std::endl;
		return 0;
	}
	std::cout << "FAIL!" << std::endl;
	return 1;
}