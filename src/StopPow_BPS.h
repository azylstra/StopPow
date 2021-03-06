// StopPow - a charged-particle stopping power library
// Copyright (C) 2014  Massachusetts Institute of Technology / Alex Zylstra

// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along
// with this program; if not, write to the Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

/**
 * @brief Calculate BPS stopping power.
 * 
 * Implement a stopping-power calculator for plasma, using
 * the BPS theory described in: L.S. Brown, D.L. Preston, and R.L. Singleton, Phys. Reports 410, 237 (2005)
 *
 * @class StopPow::StopPow_BPS
 * @author Alex Zylstra
 * @date 2014/09/17
 * @copyright Alex Zylstra / MIT
 */

#ifndef STOPPOW_BPS_H
#define STOPPOW_BPS_H

#include <math.h>

#include <vector>
#include <stdexcept>
#include <iostream>
#include <functional>
#include <thread>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_dawson.h>

#include "StopPow_Plasma.h"
#include "StopPow_Constants.h"

namespace StopPow
{

class StopPow_BPS : public StopPow_Plasma
{
public:
	/** Initialize the stopping power. Electrons should be included!
	 * @param mt the test particle mass in AMU
	 * @param Zt the test particle in charge (units of e)
	 * @param mf vector containing ordered field ion masses in AMU
	 * @param Zf vector containing ordered field ion charges in units of e
	 * @param Tf vector containing ordered field ion temperatures in units of keV
	 * @param nf vector containing ordered field ion densities in units of 1/cm3
 	 * @throws invalid_argument
	 */
	StopPow_BPS(double mt, double Zt, std::vector<double> & mf , std::vector<double> & Zf, std::vector<double> & Tf, std::vector<double> & nf) throw(std::invalid_argument);

	/** Initialize the stopping power. Electrons should be included!
	 * @param mt the test particle mass in AMU
	 * @param Zt the test particle in charge (units of e)
	 * @param field vector containing field ion info. Each row is one type of ion, then the array must contain:
	 * [mf,Zf,Tf,nf] in units of AMU, e, keV, and 1/cc
 	 * @throws invalid_argument
	 */
	StopPow_BPS(double mt, double Zt, std::vector< std::array<double,4> > & field) throw(std::invalid_argument);

	/** Initialize the stopping power. Electrons should not be included - they will be added automatically!
	 * @param mt the test particle mass in AMU
	 * @param Zt the test particle in charge (units of e)
	 * @param mf vector containing ordered field ion masses in AMU
	 * @param Zf vector containing ordered field ion charges in units of e
	 * @param Tf vector containing ordered field ion temperatures in units of keV
	 * @param nf vector containing ordered field ion densities in units of 1/cm3
	 * @param Te the electron temperature in keV
 	 * @throws invalid_argument
	 */
	StopPow_BPS(double mt, double Zt, std::vector<double> & mf , std::vector<double> & Zf, std::vector<double> & Tf, std::vector<double> & nf, double Te) throw(std::invalid_argument);

	/** Initialize the stopping power. Electrons should not be included - they will be added automatically!
	 * @param mt the test particle mass in AMU
	 * @param Zt the test particle in charge (units of e)
	 * @param field vector containing field ion info. Each row is one type of ion, then the array must contain:
	 * [mf,Zf,Tf,nf] in units of AMU, e, keV, and 1/cc
	 * @param Te the electron temperature in keV
 	 * @throws invalid_argument
	 */
	StopPow_BPS(double mt, double Zt, std::vector< std::array<double,4> > & field, double Te) throw(std::invalid_argument);

	/** Destructor */
	~StopPow_BPS();
	
	/** Calculate the total stopping power
	 * @param E the test particle energy in MeV
	 * @return stopping power in units of MeV/um
 	 * @throws invalid_argument
	 */
	double dEdx_MeV_um(double E) throw(std::invalid_argument);

	/** Calculate the total stopping power
	 * @param E the test particle energy in MeV
	 * @return stopping power in units of MeV/(mg/cm2)
 	 * @throws invalid_argument
	 */
	double dEdx_MeV_mgcm2(double E) throw(std::invalid_argument);

	/** Get stopping power due only to a specific field particle species
	* @param E the projectile energy in MeV
	* @param i the field particle index
 	* @throws invalid_argument
	*/
	double dEdx_field(double E, int i) throw(std::invalid_argument);

	/** Classical short-range stopping power (Eq. 3.3)
	* @param E the test particle energy in MeV
	* @returns the stopping power in units of MeV/um
	*/
	double dEdx_short(double E);

	/** Classical short-range stopping power (Eq. 3.3) for a single species
	* @param E the test particle energy in MeV
	* @param i the field particle index
	* @returns the stopping power in units of MeV/um
	*/
	double dEdx_short(double E, int i);

	/** Classical long-range stopping power (Eq. 3.4)
	* @param E the test particle energy in MeV
	* @returns the stopping power in units of MeV/um
	*/
	double dEdx_long(double E);

	/** Classical long-range stopping power (Eq. 3.4) for a single species
	* @param E the test particle energy in MeV
	* @param i the field particle index
	* @returns the stopping power in units of MeV/um
	*/
	double dEdx_long(double E, int i);

	/** Quantum correction to the stopping power (Eq. 3.19)
	* @param E the test particle energy in MeV
	* @returns the stopping power in units of MeV/um
	*/
	double dEdx_quantum(double E);

	/** Quantum correction to the stopping power (Eq. 3.19) for a single species
	* @param E the test particle energy in MeV
	* @param i the field particle index
	* @returns the stopping power in units of MeV/um
	*/
	double dEdx_quantum(double E, int i);

	/**
	 * Get the minimum energy that can be used for dE/dx calculations (inclusive)
	 * @return Emin in MeV
	 */
	double get_Emin();

	/**
	 * Get the maximum energy that can be used for dE/dx calculations (inclusive)
	 * @return Emax in MeV
	 */
	double get_Emax();

	/** BPS dielectric susceptibility function (Eq. 3.9), real part
	* @param u the velocity
	*/
	double Fc_real(double u);

	/** BPS dielectric susceptibility function (Eq. 3.9), imaginary part
	* @param u the velocity
	*/
	double Fc_imag(double u);

protected:
	/** Method called after field particles are changed.
	* Override if you want to do pre-calculations
	*/
	void on_field_change();

private:
	/** Initialization routine (beyond what is done by superclass constructor) */
	void init();

	// Some stored results
	std::vector<double> Mpb;
	std::vector<double> mpb;
	std::vector<double> beta_b;
	std::vector<double> kappa_b;
	double kappa_D, K;
	std::vector<double> rho_b_prefac;
	std::vector<double> eta_pb_prefac;

	/** Calculate "spectral weight" (Eq. 3.11)
	* @param i the field particle index
	* @param v the velocity in cm/s
	* @return rho_b
	*/
	double rho_b(double v, int i);

	/** Calculate total "spectral weight" (Eq. 3.10)
	* @param v the velocity in cm/s
	* @return rho_b
	*/
	double rho_tot(double v);

	/** BPS "quantum parameter" (Eq. 3.1)
	* @param vpb the velocity
	* @param i the field particle index
	* @return eta_pb
	*/
	double eta_pb(double vpb, int i);

	/** BPS dielectric susceptibility function (Eq. 3.9)
	* @param u the velocity
	*/
	gsl_complex Fc(double u);

	/** Function to integrate for long-range stopping power (Eq. 3.4)
	* @param vp the test particle velocity in cm/s
	* @param x integration parameter [=cos theta]
	* @param field particle index being evaluated
	* @return value of integrand at x
	*/
	gsl_complex dEdx_long_func(double vp, double x, int i);

	// Helper math stuff:
	
	/** Error function with a purely imaginary input: erfi = erf(i*z)/i.
	* @param z the input
	* @return erf(i*z)/i
	*/
	double erfi(double z);



	/* Minimum energy for dE/dx calculations */
	static const double Emin; 
	/* Maximum energy for dE/dx calculations */
	static const double Emax; 
};

} // end namespace StopPow

#endif