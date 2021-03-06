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
 * @brief Calculate Mehlhorn stopping power.
 * 
 * Implement a stopping-power calculator for partially ionized matter, using
 * the theory described in Mehlhorn 1981 J Appl Phys publication. Unlike the
 * reference, however, this uses the Li-Petrasso dE/dx for free electrons.
 *
 * @class StopPow::StopPow_Mehlhorn
 * @author Alex Zylstra
 * @date 2014/04/02
 * @copyright MIT / Alex Zylstra
 */

#ifndef STOPPOW_MEhLHORN_H
#define STOPPOW_MEhLHORN_H

#include <math.h>

#include <vector>
#include <stdexcept>

#include "StopPow_PartialIoniz.h"
#include "StopPow_Constants.h"
#include "StopPow_LP.h"
#include "AtomicData.h"

namespace StopPow
{

class StopPow_Mehlhorn : public StopPow_PartialIoniz
{
public:
	/** Initialize the stopping power. Electrons should not be included in lists - they will be added automatically!
	 * @param mt the test particle mass in AMU
	 * @param Zt the test particle in charge (units of e)
	 * @param mf vector containing ordered field ion masses in AMU
	 * @param Zf vector containing ordered field ion charges in units of e
	 * @param Tf vector containing ordered field ion temperatures in units of keV
	 * @param nf vector containing ordered field ion densities in units of 1/cm3
	 * @param Zbar a vector containing the average ionization state for each field ion. Zbar=Z corresponds to fully ionized material.
 	 * @param Te the electron temperature in keV
 	 * @throws invalid_argument
	 */
	StopPow_Mehlhorn(double mt, double Zt, std::vector<double> & mf, std::vector<double> & Zf, std::vector<double> & Tf, std::vector<double> & nf, std::vector<double> & Zbar, double Te) throw(std::invalid_argument);

	/** Initialize the stopping power. Electrons should not be included in lists - they will be added automatically!
	 * @param mt the test particle mass in AMU
	 * @param Zt the test particle in charge (units of e)
	 * @param field vector containing field ion info. Each row is one type of ion, then the array must contain:
	 * [mf,Zf,Tf,nf,Zbar] in units of AMU, e, keV, and 1/cc, e
	 * @param Te the electron temperature in keV
 	 * @throws invalid_argument
	 */
	StopPow_Mehlhorn(double mt, double Zt, std::vector< std::array<double,5> > & field, double Te) throw(std::invalid_argument);

	/** Destructor */
	~StopPow_Mehlhorn();
	
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

	/** Calculate the effective average ionization potential
	 * @param E the test particle energy in MeV
	 * @param index the field particle index
	 * @return value of Ibar in ergs
	 */
	double Ibar(double E, int index);

	/**
	 * Set the effective ionization potential for each ion
	 * @param Ibar the vector of potentials in eV
	 */
	void set_Ibar(std::vector<double> Ibar) throw(std::invalid_argument);

private:
	/** Specific initialization routines */
	void init();
	
	/** Manually specified Ibar if appropriate */
	std::vector<double> Ibar_manual;

	/** Li-Petrasso stopping power for the free electons and ions */
	StopPow * PlasmaStop;
	/** Whether to use the manual Ibar */
	bool use_manual_Ibar;

	// helper functions:

	/** Calculate the LSS low-energy stopping power
	 * @param E the test particle energy in MeV
	 * @param index the field particle index
	 * @return value of (dE/dx)_LSS in MeV/um
	 */
	double dEdx_LSS(double E, int index);

	/** Calculate the nuclear stopping power
	 * @param E the test particle energy in MeV
	 * @param index the field particle index
	 * @return value of (dE/dx)_nuc in MeV/um
	 */
	double dEdx_nuc(double E, int index);

	/** Calculate the Bethe stopping power
	 * @param E the test particle energy in MeV
	 * @param index the field particle index
	 * @return value of (dE/dx)_Bethe in MeV/um
	 */
	double dEdx_Bethe(double E, int index);

	/** Calculate the effective projectile charge
	 * @param E the test particle energy in MeV
	 * @return value of the effective charge
	 */
	double ZtEff(double E);

	/** Calculate shell correction term in log lambda for
	  * shell corrections
	  * Data are taken from Andersen and Ziegler, The Stopping and Ranges of Ions in Matter, Vol 3:
	  * Hydrogen stopping powers and ranges in all elements (1978).
	  * @param Zf the field particle atomic number
	  * @param E the test particle energy in MeV
	  * @return shell correction term
	  */
	double shell_term(double Zf, double E);


	/* Minimum energy for dE/dx calculations */
	static const double Emin; 
	/* Maximum energy for dE/dx calculations */
	static const double Emax; 
};

} // end namespace StopPow

#endif