/**
 * @brief Calculate Zimmerman stopping power.
 * 
 * Implement a stopping-power calculator for partially ionized matter, using
 * the theory described by Zimmerman in:
 * "Recent Developments in Monte Carlo Techniques", The Nuclear Explosive Code Developers' Conference, 1990.
 * UCRL-JC-105616
 *
 * @class StopPow::StopPow_Zimmerman
 * @author Alex Zylstra
 * @date 2014/04/03
 * @copyright Alex Zylstra / MIT
 */

#ifndef STOPPOW_ZIMMERMAN_H
#define STOPPOW_ZIMMERMAN_H

#include <math.h>

#include <vector>
#include <stdexcept>
#include <gsl/gsl_sf_erf.h>

#include "StopPow_PartialIoniz.h"
#include "StopPow_Constants.h"
#include "AtomicData.h"

namespace StopPow
{

class StopPow_Zimmerman : public StopPow_PartialIoniz
{
public:
	// Partially ionized constructors:

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
	StopPow_Zimmerman(float mt, float Zt, std::vector<float> & mf, std::vector<float> & Zf, std::vector<float> & Tf, std::vector<float> & nf, std::vector<float> & Zbar, float Te) throw(std::invalid_argument);

	/** Initialize the stopping power. Electrons should not be included in lists - they will be added automatically!
	 * @param mt the test particle mass in AMU
	 * @param Zt the test particle in charge (units of e)
	 * @param field vector containing field ion info. Each row is one type of ion, then the array must contain:
	 * [mf,Zf,Tf,nf,Zbar] in units of AMU, e, e, keV, and 1/cc
	 * @param Te the electron temperature in keV
 	 * @throws invalid_argument
	 */
	StopPow_Zimmerman(float mt, float Zt, std::vector< std::array<float,5> > & field, float Te) throw(std::invalid_argument);

	/** Destructor */
	~StopPow_Zimmerman();
	
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

	/** Free electron component of the stopping power
	* @param E the test particle energy in MeV
	* @return dE/dx due to electrons in MeV/um
	*/
	float dEdx_free_electron(float E);

	/** Bound electron component of the stopping power
	* @param E the test particle energy in MeV
	* @return dE/dx due to electrons in MeV/um
	*/
	float dEdx_bound_electron(float E);

	/** Ion component of the stopping power
	* @param E the test particle energy in MeV
	* @return dE/dx due to electrons in MeV/um
	*/
	float dEdx_ion(float E);

private:
	/** Specific initialization routines */
	void init();
	void init_plasma();

	// helper functions:

	/** Free electron stopping number, Eq 15
	* @param y ratio of particle/field velocity (y=Vp/Ve)
	* @param LambdaF defined in Eq 16
	* @return free electron stopping number
	*/
	float LF(float y, float LambdaF);

	/** Calculate total Debye length with all plasma components
	* @returns Debye length in cm
	*/
	float lDebye();

	/* Minimum energy for dE/dx calculations */
	static const float Emin; 
	/* Maximum energy for dE/dx calculations */
	static const float Emax; 
};

} // end namespace StopPow

#endif