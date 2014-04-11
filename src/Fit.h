/**
 * @brief Utilities for fitting data
 *
 * This file defines several static methods within
 * the StopPow namespace used for fitting data.
 *
 * @author Alex Zylstra
 * @date 2014/04/09
 * @copyright Alex Zylstra / MIT
 */

#ifndef FIT_H
#define FIT_H

#include <math.h>
#include <stdexcept>
#include <vector>
#include <array>
#include <iostream>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_deriv.h>
//#include <gsl/gsl_multimin.h>
#include <gsl/gsl_roots.h>

#include "StopPow.h"
#include "Spectrum.h"

namespace StopPow
{

/** Fit a Gaussian to provided data. The results are placed in variables passed by reference!
* @param data_x the independent variable (x) values
* @param data_y the dependent variable (y) values
* @param data_std error bar on y values, assumed normally distributed
* @param fit an allocated `array` to store the fit coefficients in. Will be resized to length=3. Order: [A,mu,sigma]
* @param fit_unc an allocated `array` to store the uncertainties in fit values. Will be resized to length=3. Order: [A,mu,sigma]
* @param chi2_dof the chi^2/dof for the resulting fit
* @param verbose set to true for gory details to be output to the console
* @return true if everything went OK
*/
bool fit_Gaussian(std::vector<double> & data_x, 
					std::vector<double> & data_y, 
					std::vector<double> & data_std, 
					std::vector<double> & fit,
					std::vector<double> & fit_unc,
					double & chi2_dof,
					bool verbose);

/** Use a Gaussian fit to infer rhoR from a proton spectrum. The results are placed in variables passed by reference!
* @param data_x the energy in MeV
* @param data_y the proton yield/MeV
* @param data_std the error bar on yield, assumed normally distributed
* @param dE any extra uncertainty in energy
* @param fit an allocated `array` to store the spectrum fit coefficients in. Will be resized to length=3. Order: [A,mu,sigma]
* @param fit_unc an allocated `array` to store the uncertainties in spectrum fit values. Will be resized to length=3. Order: [A,mu,sigma]
* @param chi2_dof the chi^2/dof for the resulting fit
* @param s the stopping power model to use
* @param E0 the initial (i.e. birth) proton energy [MeV]
* @param rhoR the calculated rhoR will be placed in this variable [mg/cm2]
* @param rhoR_unc the calculated uncertainty in rhoR will be placed in this variable [mg/cm2]
* @param verbose set to true for gory details to be output to the console
* @return true if everything went OK
*/
bool fit_rhoR(std::vector<double> & data_x, 
				std::vector<double> & data_y, 
				std::vector<double> & data_std,
				double & dE, 
				std::vector<double> & fit,
				std::vector<double> & fit_unc,
				double & chi2_dof,
				StopPow & s,
				double & E0,
				double & rhoR,
				double & rhoR_unc,
				bool verbose);

/** Use a Gaussian forward fit to infer rhoR from a proton spectrum. The results are placed in variables passed by reference!
* This algorithm uses a forward fit, i.e. a trial Gaussian is convolved with the rhoR downshift and compared to the data.
* @param data_x the energy in MeV
* @param data_y the proton yield/MeV
* @param data_std the error bar on yield, assumed normally distributed
* @param dE any extra uncertainty in energy
* @param chi2_dof the chi^2/dof for the resulting fit
* @param s the stopping power model to use
* @param E0 the initial (i.e. birth) proton energy [MeV]
* @param rhoR the calculated rhoR will be placed in this variable [mg/cm2]
* @param rhoR_unc the calculated uncertainty in rhoR will be placed in this variable [mg/cm2]
* @param verbose set to true for gory details to be output to the console
* @return true if everything went OK
*/
bool forward_fit_rhoR(std::vector<double> & data_x, 
						std::vector<double> & data_y, 
						std::vector<double> & data_std,
						double & dE, 
						double & chi2_dof,
						StopPow & s,
						double & E0,
						double & rhoR,
						double & rhoR_unc,
						bool verbose);

/** Use a Gaussian fit to infer rhoR from a proton spectrum. The results are placed in variables passed by reference!
* This algorithm uses a deconvolution, i.e. the observed spectrum is downshift-corrected then fit with a Gaussian.
* @param data_x the energy in MeV
* @param data_y the proton yield/MeV
* @param data_std the error bar on yield, assumed normally distributed
* @param dE any extra uncertainty in energy
* @param chi2_dof the chi^2/dof for the resulting fit
* @param s the stopping power model to use
* @param E0 the initial (i.e. birth) proton energy [MeV]
* @param rhoR the calculated rhoR will be placed in this variable [mg/cm2]
* @param rhoR_unc the calculated uncertainty in rhoR will be placed in this variable [mg/cm2]
* @param verbose set to true for gory details to be output to the console
* @return true if everything went OK
*/
bool deconvolve_fit_rhoR(std::vector<double> & data_x, 
						std::vector<double> & data_y, 
						std::vector<double> & data_std,
						double & dE, 
						double & chi2_dof,
						StopPow & s,
						double & E0,
						double & rhoR,
						double & rhoR_unc,
						bool verbose);

} // end of namespace

#endif