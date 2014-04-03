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
	StopPow_Mehlhorn(float mt, float Zt, std::vector<float> & mf, std::vector<float> & Zf, std::vector<float> & Tf, std::vector<float> & nf, std::vector<float> & Zbar, float Te) throw(std::invalid_argument);

	/** Initialize the stopping power. Electrons should not be included in lists - they will be added automatically!
	 * @param mt the test particle mass in AMU
	 * @param Zt the test particle in charge (units of e)
	 * @param field vector containing field ion info. Each row is one type of ion, then the array must contain:
	 * [mf,Zf,Tf,nf,Zbar] in units of AMU, e, e, keV, and 1/cc
	 * @param Te the electron temperature in keV
 	 * @throws invalid_argument
	 */
	StopPow_Mehlhorn(float mt, float Zt, std::vector< std::array<float,5> > & field, float Te) throw(std::invalid_argument);

	/** Destructor */
	~StopPow_Mehlhorn();
	
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

	/** Calculate the effective average ionization potential
	 * @param E the test particle energy in MeV
	 * @param index the field particle index
	 * @return value of Ibar in ergs
	 */
	float Ibar(float E, int index);

	/**
	 * Set the effective ionization potential for each ion
	 * @param Ibar the vector of potentials in eV
	 */
	void set_Ibar(std::vector<float> Ibar) throw(std::invalid_argument);

private:
	/** Specific initialization routines */
	void init();
	
	/** Manually specified Ibar if appropriate */
	std::vector<float> Ibar_manual;

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
	float dEdx_LSS(float E, int index);

	/** Calculate the nuclear stopping power
	 * @param E the test particle energy in MeV
	 * @param index the field particle index
	 * @return value of (dE/dx)_nuc in MeV/um
	 */
	float dEdx_nuc(float E, int index);

	/** Calculate the Bethe stopping power
	 * @param E the test particle energy in MeV
	 * @param index the field particle index
	 * @return value of (dE/dx)_Bethe in MeV/um
	 */
	float dEdx_Bethe(float E, int index);

	/** Calculate the effective projectile charge
	 * @param E the test particle energy in MeV
	 * @return value of the effective charge
	 */
	float ZtEff(float E);

	/** Calculate shell correction term in log lambda for
	  * shell corrections
	  * Data are taken from Andersen and Ziegler, The Stopping and Ranges of Ions in Matter, Vol 3:
	  * Hydrogen stopping powers and ranges in all elements (1978).
	  * @param Zf the field particle atomic number
	  * @param E the test particle energy in MeV
	  * @return shell correction term
	  */
	float shell_term(float Zf, float E);


	/* Minimum energy for dE/dx calculations */
	static const float Emin; 
	/* Maximum energy for dE/dx calculations */
	static const float Emax; 
};

} // end namespace StopPow

#endif