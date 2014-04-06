/**
* @brief Calculate Bethe-Bloch stopping power.
* 
* Implement a stopping-power calculator for arbitrary cold matter, using
* the simple Bethe-Bloch theory.
*
* @class StopPow::StopPow_BetheBloch
* @author Alex Zylstra
* @date 2014/03/06
* @copyright MIT / Alex Zylstra
*/

#ifndef STOPPOW_BETHEBLOCH_H
#define STOPPOW_BETHEBLOCH_H

#include <math.h>
#include <cmath>

#include <stdexcept>
#include <vector>
#include <array>

#include "StopPow.h"
#include "StopPow_Constants.h"
#include "AtomicData.h"

namespace StopPow
{

class StopPow_BetheBloch : public StopPow
{
public:
	/** Initialize the Bethe-Bloch calculator.
	 * @param mt the test particle mass in AMU
	 * @param Zt the test particle in charge (units of e)
	 * @param mf vector containing ordered field particle masses in AMU
	 * @param Zf vector containing ordered field particle charges in units of e
	 * @param nf vector containing ordered field particle densities in units of 1/cm3
 	 * @throws invalid_argument
	 */
	StopPow_BetheBloch(double mt, double Zt, std::vector<double> mf , std::vector<double> Zf, std::vector<double> nf) throw(std::invalid_argument);

	/** Destructor */
	~StopPow_BetheBloch();
	
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

	/**
	  * Turn shell corrections on or off in the model.
	  * @param enabled set to true to use shell corrections
	  */
	void use_shell_correction(bool enabled);

	/**
	  * Get whether the model is currently using shell corrections
	  * @return true if shell corrections are enabled
	  */
	bool using_shell_correction();

	/**
	 * Set the effective ionization potential manually.
	 * @param Ibar the value to use in eV for each field particle
	 */
	void set_Ibar(std::vector<double> Ibar) throw(std::invalid_argument);

	/** Effecive ionization potential as a function of Z.
	 * @param Zf field particle charge in units of e
	 * @return Ibar in erg
	 */
	double Ibar(double Zf);

private:
	/** Calculate shell correction term in log lambda for
	  * shell corrections
	  * Data are taken from Andersen and Ziegler, The Stopping and Ranges of Ions in Matter, Vol 3:
	  * Hydrogen stopping powers and ranges in all elements (1978).
	  * @param Zf the field particle atomic number
	  * @param E the test particle energy in MeV
	  * @return shell correction term
	  */
	double shell_term(double Zf, double E);

	// data on the field particles:
	/** mass in atomic units */
	std::vector<double> mf; 
	/** charge in atomic units */
	std::vector<double> Zf; 
	/** particle density in 1/cc */
	std::vector<double> nf; 
	/** For manual setting of Ibar in eV */
	std::vector<double> Ibar_manual;
	/** number of field particle species */
	int num; 
	/** mass density in g/cc */
	double rho; 

	// type of test particle:
	/** mass in atomic units */
	double mt; 
	/** charge in atomic units */
	double Zt; 

	/* Minimum energy for dE/dx calculations */
	double Emin; 
	/* Maximum energy for dE/dx calculations */
	double Emax; 

	/** Represent the type of model described by this class */
	static const std::string TYPE;
	/** Some information about the model, stored as string */
	std::string info;

	/** If shell corrections should be used */
	bool use_shell_corr;
	/** If manual Ibar should be used */
	bool use_manual_Ibar;
};

} // end namespace StopPow

#endif