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
 * @brief Calculate Grabowski's fit to Molecular Dynamics (MD) stopping power.
 * 
 * Implement fit given in: P.E. Grabowski et al., Phys. Rev. Lett. 111, 215002 (2013)
 *
 * @class StopPow::StopPow_Grabowski
 * @author Alex Zylstra
 * @date 2014/04/01
 * @copyright Alex Zylstra / MIT
 */

#ifndef STOPPOW_GRABOWSKI_H
#define STOPPOW_GRABOWSKI_H

#include <math.h>

#include <vector>
#include <stdexcept>
 
#include <gsl/gsl_sf_erf.h>

#include "StopPow_Plasma.h"
#include "StopPow_Constants.h"

namespace StopPow
{

class StopPow_Grabowski : public StopPow_Plasma
{
public:
	/** Initialize the Li-Petrasso stopping power. Electrons should be included!
	 * @param mt the test particle mass in AMU
	 * @param Zt the test particle in charge (units of e)
	 * @param mf vector containing ordered field ion masses in AMU
	 * @param Zf vector containing ordered field ion charges in units of e
	 * @param Tf vector containing ordered field ion temperatures in units of keV
	 * @param nf vector containing ordered field ion densities in units of 1/cm3
 	 * @throws invalid_argument
	 */
	StopPow_Grabowski(double mt, double Zt, std::vector<double> & mf , std::vector<double> & Zf, std::vector<double> & Tf, std::vector<double> & nf) throw(std::invalid_argument);

	/** Initialize the stopping power. Electrons should be included!
	 * @param mt the test particle mass in AMU
	 * @param Zt the test particle in charge (units of e)
	 * @param field vector containing field ion info. Each row is one type of ion, then the array must contain:
	 * [mf,Zf,Tf,nf] in units of AMU, e, keV, and 1/cc
 	 * @throws invalid_argument
	 */
	StopPow_Grabowski(double mt, double Zt, std::vector< std::array<double,4> > & field) throw(std::invalid_argument);

	/** Initialize the stopping power.
	 * @param mt the test particle mass in AMU
	 * @param Zt the test particle in charge (units of e)
	 * @param mf vector containing ordered field particle masses in AMU
	 * @param Zf vector containing ordered field particle charges in units of e
	 * @param Tf vector containing ordered field particle temperatures in units of keV
	 * @param nf vector containing ordered field particle densities in units of 1/cm3
	 * @param Te the electron temperature in keV
 	 * @throws invalid_argument
	 */
	StopPow_Grabowski(double mt, double Zt, std::vector<double> & mf , std::vector<double> & Zf, std::vector<double> & Tf, std::vector<double> & nf, double Te) throw(std::invalid_argument);

	/** Initialize the stopping power.
	 * @param mt the test particle mass in AMU
	 * @param Zt the test particle in charge (units of e)
	 * @param field vector containing field particle info. Each row is one type of particle, then the array must contain:
	 * [mf,Zf,Tf,nf] in units of AMU, e, keV, and 1/cc
	 * @param Te the electron temperature in keV
 	 * @throws invalid_argument
	 */
	StopPow_Grabowski(double mt, double Zt, std::vector< std::array<double,4> > & field, double Te) throw(std::invalid_argument);

	/** Destructor */
	~StopPow_Grabowski();
	
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

private:
	/** Initialization routine (beyond what is done by superclass constructor) */
	void init();

	/*    Constants    */
	double gammaE {0.577216}; // Euler's constant
	double alpha {4.*exp(-2.*gammaE)};
	double G_a {1.04102e-5};
	double G_b {0.183260};
	double G_c {0.116053};
	double G_d {0.824982};
	double G_g0 {2.03301e-3};

	/*    Helper Functions from Eq 2 in paper    */
	double M1(double g, double s, double Z);
	double M2(double w, double g, double s);
	double R(double w, double g, double s, double Z);
	double G(double w);
	double H(double w);

	/* Minimum energy for dE/dx calculations */
	static const double Emin; 
	/* Maximum energy for dE/dx calculations */
	static const double Emax; 
};

} // end namespace StopPow

#endif