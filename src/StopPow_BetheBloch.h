/**
 * @brief Calculate Bethe-Bloch stopping power.
 * 
 * Implement a stopping-power calculator for arbitrary cold matter, using
 * the simple Bethe-Bloch theory.
 *
 * @class StopPow_BetheBloch
 * @author Alex Zylstra
 * @date 2013/04/03
 * @copyright MIT / Alex Zylstra
 */

#ifndef STOPPOW_BETHEBLOCH_H
#define STOPPOW_BETHEBLOCH_H

#include <math.h>

#include <stdexcept>
#include <vector>

#include "StopPow.h"
#include "StopPow_Constants.h"

namespace StopPow
{

class StopPow_BetheBloch : public StopPow
{
public:
	/** Initialize the Bethe-Bloch calculator.
	 * @param mt the test particle mass in AMU
	 * @param Zt the test particle in charge (units of e)
	 * @param mf vector containing ordered field particle masses in AMU
	 * @param Zt vector containing ordered field particle charges in units of e
	 * @param nf vector containing ordered field particle densities in units of 1/cm3
 	 * @throws invalid_argument
	 */
	StopPow_BetheBloch(float mt, float Zt, std::vector<float> mf , std::vector<float> Zf, std::vector<float> nf);

	/** Calculate the total stopping power
	 * @param E the test particle energy in MeV
	 * @return stopping power in units of MeV/um
 	 * @throws invalid_argument
	 */
	float dEdx_MeV_um(float E);

	/** Calculate the total stopping power
	 * @param E the test particle energy in MeV
	 * @return stopping power in units of MeV/(mg/cm2)
 	 * @throws invalid_argument
	 */
	float dEdx_MeV_mgcm2(float E);

	/**
	 * Get the minimum energy that can be used for dE/dx calculations
	 * @return Emin in MeV
	 */
	float get_Emin();

	/**
	 * Get the maximum energy that can be used for dE/dx calculations
	 * @return Emax in MeV
	 */
	float get_Emax();

private:
	/** Effecive ionization potential as a function of Z.
	 * @param Zf field particle charge in units of e
	 * @return Ibar in erg
	 */
	float Ibar(float Zf);

	// data on the field particles:
	std::vector<float> mf; /** mass in atomic units */
	std::vector<float> Zf; /** charge in atomic units */
	std::vector<float> nf; /** particle density in 1/cc */
	int num; /** number of field particle species */
	float rho; /** mass density in g/cc */
	// type of test particle:
	float mt; /** mass in atomic units */
	float Zt; /** charge in atomic units */

	// ionization data
	static const float IbarData[];

	static const float Emin; /* Minimum energy for dE/dx calculations */
	static const float Emax; /* Maximum energy for dE/dx calculations */
};

} // end namespace StopPow

#endif