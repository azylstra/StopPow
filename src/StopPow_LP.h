/**
 * @brief Calculate Li-Petrasso stopping power.
 * 
 * Implement a stopping-power calculator for plasma, using
 * the Fokker-Planck theory described in Li and Petrasso, PRL 1993.
 *
 * @class StopPow::StopPow_LP
 * @author Alex Zylstra
 * @date 2013/06/04
 * @copyright MIT / Alex Zylstra
 */

#ifndef STOPPOW_LP_H
#define STOPPOW_LP_H

#include <math.h>

#include <vector>
#include <stdexcept>

#include <boost/math/special_functions.hpp>

#include "StopPow.h"
#include "StopPow_Constants.h"

namespace StopPow
{

class StopPow_LP : public StopPow
{
public:
	/** Initialize the Li-Petrasso stopping power.
	 * @param mt the test particle mass in AMU
	 * @param Zt the test particle in charge (units of e)
	 * @param mf vector containing ordered field particle masses in AMU
	 * @param Zf vector containing ordered field particle charges in units of e
	 * @param Tf vector containing ordered field particle temperatures in units of keV
	 * @param nf vector containing ordered field particle densities in units of 1/cm3
 	 * @throws invalid_argument
	 */
	StopPow_LP(float mt, float Zt, std::vector<float> mf , std::vector<float> Zf, std::vector<float> Tf, std::vector<float> nf) throw(std::invalid_argument);

	/** Destructor */
	~StopPow_LP();
	
	/** Calculate the total stopping power
	 * @param E the test particle energy in MeV
	 * @return stopping power in units of MeV/um
 	 * @throws invalid_argument
	 */
	float dEdx_MeV_um(float E) throw(std::invalid_argument);

	/** Calculate the total stopping power
	 * @param E the test particle energy in MeV
	 * @return stopping power in units of MeV/(mg/cm2)
 	 * @throws invalid_argument
	 */
	float dEdx_MeV_mgcm2(float E) throw(std::invalid_argument);

	/** Turn collective effects on or off.
	 * @param set if you want to use collective effects
	 */
	void set_collective(bool set);

	/**
	 * Get the minimum energy that can be used for dE/dx calculations (inclusive)
	 * @return Emin in MeV
	 */
	float get_Emin();

	/**
	 * Get the maximum energy that can be used for dE/dx calculations (inclusive)
	 * @return Emax in MeV
	 */
	float get_Emax();

private:
	// data on the field particles:
	/** mass in atomic units */
	std::vector<float> mf; 
	/** charge in atomic units */
	std::vector<float> Zf; 
	/** Temperature in keV */
	std::vector<float> Tf; 
	/** particle density in 1/cc */
	std::vector<float> nf; 
	/** number of field particle species */
	int num; 
	/** mass density in g/cc */
	float rho; 

	// type of test particle:
	/** mass in atomic units */
	float mt; 
	/** charge in atomic units */
	float Zt; 
	/** Use collective effects? */
	bool collective; 

	// helper functions:

	/** Calculate the Coulomb logarithm
	 * @param E the test particle energy in MeV
	 * @param index the field particle index
	 * @return value of Log(Lambda)
	 */
	float LogLambda(float E, int index);

	/** Chandrasekhar function
	 * @param E the test particle energy in MeV
	 * @param index the field particle index
	 * @return value of the Chandrasekhar function G (see L-P 1993)
	 */
	float G(float E, int index);

	/** Debye length in field plasma
	 * @return Debye length in cm
	 */
	float lDebye();

	/** Field particle thermal velocity
	 * @param index the field particle's index (for mf,Zf,Tf,nf arrays)
	 * @return thermal velocity in cm/s
	 */
	float vtf(int index);

	/** x^{t/f} parameter from Li 1993
	 * @param E test particle energy in MeV
	 * @param index the field particle's index (for mf,Zf,Tf,nf arrays)
	 * @return x^{t/f}
	 */
	float xtf(float E, int index);

	/** Relative velocity between test particle and field particle
	 * @param E test particle energy in MeV
	 * @param index the field particle's index (for mf,Zf,Tf,nf arrays)
	 * @return relative velocity in cm/s
	 */
	float u(float E, int index);

	/* Minimum energy for dE/dx calculations */
	static const float Emin; 
	/* Maximum energy for dE/dx calculations */
	static const float Emax; 
};

} // end namespace StopPow

#endif