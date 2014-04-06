/**
 * @brief Cold-matter proton stopping from fits.
 * 
 * A wrapper class for calculating proton stopping powers
 * using fits to data.
 * Data are taken from Andersen and Ziegler, The Stopping and Ranges of Ions in Matter, Vol 3:
 * Hydrogen stopping powers and ranges in all elements (1978).
 *
 * @class StopPow::StopPow_AZ
 * @author Alex Zylstra
 * @date 2013/06/06
 * @copyright MIT / Alex Zylstra
 */

#ifndef STOPPOW_AZ_H
#define STOPPOW_AZ_H

#include <math.h>

#include <cmath>
#include <string>
#include <sstream>
#include <array>
#include <stdexcept>

#include "StopPow.h"
#include "StopPow_Constants.h"
#include "AtomicData.h"

namespace StopPow
{

class StopPow_AZ : public StopPow
{
public:
	/**
	 * Constructor for Andersen-Ziegler stopping power. 
	 * Note: this constructor uses 'standard' mass densities.
	 * @param Z the atomic number (1-92) of the target element.
	 * @throws invalid_argument
	 */
	explicit StopPow_AZ(int Z) throw(std::invalid_argument);

	/**
	 * Constructor for Andersen-Ziegler stopping power. 
	 * @param Z the atomic number (1-92) of the target element.
	 * @param rho the mass density in g/cc
	 * @throws invalid_argument
	 */
	StopPow_AZ(int Z, double rho) throw(std::invalid_argument);

	/**
	 * Destructor
	 */
	~StopPow_AZ();

	/**
	 * Get stopping power from the data.
	 * @param E the particle energy in MeV
	 * @return dE/dx in MeV/um
 	 * @throws invalid_argument
	 */
	double dEdx_MeV_um(double E) throw(std::invalid_argument);

	/**
	 * Get stopping power from the data.
	 * @param E the particle energy in MeV
	 * @return dE/dx in MeV/(mg/cm2)
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

private:
	/** Atomic number */
	int Z;

	/** mass density in g/cm3 */
	double rho; 

	/** atomic number density in 1/cm3 */
	double ni; 

	/** Minimum allowable energy for calculations (MeV) */
	double Emin;

	/** Maximum allowable energy for calculations (MeV) */
	double Emax;

	static const std::array< std::array<double,12> , 92> fit_coeff;
};

} // end namespace StopPow
 #endif