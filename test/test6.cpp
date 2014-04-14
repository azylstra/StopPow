/** Test fitting functions
 * @author Alex Zylstra
 * @date 2014/04/09
 */

#include <iostream>
#include <vector>
#include <array>
#include <ctime>

#include "StopPow.h"
#include "StopPow_SRIM.h"
#include "Util.h"
#include "Fit.h"

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
	std::cout << "========== Test Suite 6 ==========" << std::endl;
	std::cout << "      Test fitting utilities  " << std::endl;
	bool pass = true;
	bool test = true;
	
	// ------------- Test fit_Gaussian ---------------
	// Some test data:	
    std::vector<double> data_x {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
    std::vector<double> data_y {-0.000304373, 0.00335043, -0.00946906, -0.0053318, 0.00118407, 
								0.0103371, 0.0421091, 0.54167, 2.43304, 3.99037, 2.42823, 0.535774, 
								0.0466508, -0.00206453, -0.0113899, -0.00150277, -0.00766465, 
								0.00126149, 0.0118316, 0.00616385};
    std::vector<double> data_std {0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01};

    // For fit results:
    std::vector<double> fit;
    std::vector<double> fit_unc;
    double chi2;

    test = StopPow::fit_Gaussian(data_x, data_y, data_std, fit, fit_unc, chi2, verbose);
    // Check results:
	test &= StopPow::approx(fit[0], 10., 1e-2);
	test &= StopPow::approx(fit[1], 10., 1e-3);
	test &= StopPow::approx(fit[2], 1., 1e-2);

    // try fitting with larger values:
    for(int i=0; i<data_y.size(); i++)
    {
    	data_y[i] *= 1e7;
    	data_std[i] *= 1e7;
    }
    test = StopPow::fit_Gaussian(data_x, data_y, data_std, fit, fit_unc, chi2, verbose);
    // Check results:
	test &= StopPow::approx(fit[0], 1e8, 1e-2);
	test &= StopPow::approx(fit[1], 10., 1e-3);
	test &= StopPow::approx(fit[2], 1., 1e-2);

	std::cout << "fit_Gaussian: " << (test ? "pass" : "FAIL!") << std::endl;
	pass &= test;

	// ------------- Test fit_rhoR ---------------
	data_x = {1., 1.2, 1.4, 1.6, 1.8, 2., 2.2, 2.4, 2.6, 2.8, 3., 3.2, 3.4, 3.6, \
3.8, 4., 4.2, 4.4, 4.6, 4.8, 5., 5.2, 5.4, 5.6, 5.8, 6., 6.2, 6.4, \
6.6, 6.8, 7., 7.2, 7.4, 7.6, 7.8, 8., 8.2, 8.4, 8.6, 8.8, 9., 9.2, \
9.4, 9.6, 9.8, 10., 10.2, 10.4, 10.6, 10.8, 11., 11.2, 11.4, 11.6, \
11.8, 12., 12.2, 12.4, 12.6, 12.8, 13., 13.2, 13.4, 13.6, 13.8, 14., \
14.2, 14.4, 14.6, 14.8, 15., 15.2, 15.4, 15.6, 15.8, 16., 16.2, 16.4, \
16.6, 16.8, 17., 17.2, 17.4, 17.6, 17.8, 18., 18.2, 18.4, 18.6, 18.8, \
19., 19.2, 19.4, 19.6, 19.8, 20.};
	data_y =  {75259.5, 143190., 107732., -135826., 68668.8, -45779.4, -31288., \
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
	data_std = {100000, 100000, 100000, 100000, 100000, 100000, 100000, 100000, \
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

	double dE = 0.05; double E0 = 14.7; double E0_unc = 0.05;
	// stopping power model:
	StopPow::StopPow_SRIM s("SRIM/Hydrogen in Aluminum.txt");
	// For fit results
	double rhoR, rhoR_unc;

	// Do the calculation
	test = StopPow::fit_rhoR(data_x, data_y, data_std, dE, fit, fit_unc, chi2, s, E0, E0_unc, rhoR, rhoR_unc, verbose);
	// Check against pre-computed values:
	test &= StopPow::approx(rhoR, 162.548, 1e-3);
	test &= StopPow::approx(rhoR_unc, 1.513, 1e-3);
	test &= StopPow::approx(fit[0], 1e7, 2e-2);
	test &= StopPow::approx(fit[1], 10., 1e-3);
	test &= StopPow::approx(fit[2], 1., 2e-2);
	std::cout << "fit_rhoR: " << (test ? "pass" : "FAIL!") << std::endl;
	std::cout << "rhoR = " << rhoR << " +/- " << rhoR_unc << std::endl;
	pass &= test;

	// --------------- Test forward_fit_rhoR ----------
	test = StopPow::forward_fit_rhoR(data_x, data_y, data_std, dE, chi2, s, E0, E0_unc, fit, fit_unc, verbose);
	std::cout << "forward_fit_rhoR = " << fit[0] << " +/- " << fit_unc[0] << std::endl;
	test &= StopPow::approx(fit[0], 162., 2e-2);
	test &= StopPow::approx(fit[1], 1e7, 2e-2);
	test &= StopPow::approx(fit[2], 1., 2e-2);
	std::cout << "Forward fit: " << (test ? "pass" : "FAIL!") << std::endl;
	pass &= test;

	// --------------- Test deconvolve_fit_rhoR ----------
	test = StopPow::deconvolve_fit_rhoR(data_x, data_y, data_std, dE, chi2, s, E0, E0_unc, fit, fit_unc, verbose);
	std::cout << "deconvolve_fit_rhoR = " << fit[0] << " +/- " << rhoR_unc << std::endl;
	test &= StopPow::approx(fit[0], 162., 1e-2);
	test &= StopPow::approx(fit[1], 1e7, 2e-2);
	test &= StopPow::approx(fit[2], 0.75, 2e-2);
	std::cout << "Deconvolution fit: " << (test ? "pass" : "FAIL!") << std::endl;
	std::cout << fit[0] << " , " << fit[1] << " , " << fit[2] << std::endl;
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