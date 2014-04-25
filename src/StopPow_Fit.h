/**
 * @brief Adjustable and combination model for fitting experimental data
 * 
 * A wrapper class for calculating stopping powers
 * using a variety of models for WDM conditions. Also
 * allows adjustment of a single parameter for use in
 * fitting routines. 
 *
 * The constructors take standard plasma conditions and initialize
 * a default model. The model adjustment is then accomplished
 * by calling the appropriate class methods as needed.
 *
 * The total stopping power is written:
 * dE/dx = (dE/dx)_fe + (dE/dx)_be + (dE/dx)_i
 * for free electrons, bound electrons, and ions respectively.
 *
 * Bound electrons and ions both use the Zimmerman model. The bound electron
 * dE/dx may be normalized to a reference case (e.g. SRIM).
 *
 * For free electrons, the following models may be used: 
 * Zimmerman (default), Li-Petrasso, BPS, Grabowski, Grabowski w/ quantum BPS
 * The entire free-electron dE/dx is scaled.
 *
 * @class StopPow::StopPow_Fit
 * @author Alex Zylstra
 * @date 2014/04/21
 * @copyright Alex Zylstra / MIT
 */

#ifndef STOPPOW_FIT_H
#define STOPPOW_FIT_H

#include <vector>
#include <stdexcept>
#include <math.h>

#include "StopPow.h"
#include "StopPow_PartialIoniz.h"
#include "StopPow_Zimmerman.h"
#include "StopPow_LP.h"
#include "StopPow_BPS.h"
#include "StopPow_Grabowski.h"

namespace StopPow
{

class StopPow_Fit : public StopPow_PartialIoniz
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
	StopPow_Fit(double mt, double Zt, std::vector<double> & mf, std::vector<double> & Zf, std::vector<double> & Tf, std::vector<double> & nf, std::vector<double> & Zbar, double Te) throw(std::invalid_argument);

	/** Initialize the stopping power. Electrons should not be included in lists - they will be added automatically!
	 * @param mt the test particle mass in AMU
	 * @param Zt the test particle in charge (units of e)
	 * @param field vector containing field ion info. Each row is one type of ion, then the array must contain:
	 * [mf,Zf,Tf,nf,Zbar] in units of AMU, e, e, keV, and 1/cc
	 * @param Te the electron temperature in keV
 	 * @throws invalid_argument
	 */
	StopPow_Fit(double mt, double Zt, std::vector< std::array<double,5> > & field, double Te) throw(std::invalid_argument);

	/**
	 * Destructor
	 */
	~StopPow_Fit();

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

	/** Normalize the bound-electron stopping to a reference case at a given proton energy.
	 * @param ref The cold-matter reference dE/dx to use
	 * @param Ep normalize the stopping at this energy
	 */
	void normalize_bound_e(StopPow * ref, double Ep);

	/** Free electron models that may be used */
	static const int MODE_ZIMMERMAN;
	static const int MODE_LP;
	static const int MODE_BPS;
	static const int MODE_GRABOWSKI;
	static const int MODE_QUANTUM_GRABOWSKI;
	/** Choose the free-electron stopping model.
	* @param new_model the new model to be used, must pass one of `MODE_ZIMMERMAN`, `MODE_LP`, `MODE_BPS`, `MODE_GRABOWSKI`, or `MODE_QUANTUM_GRABOWSKI`.
	*/
	void choose_model(int new_model) throw(std::invalid_argument);

	/** Adjust the stopping, to be used for fitting
	* @param factor set the free-electron adjustment factor
	*/
	void set_factor(double factor);

	/** Get the current adjustment factor
	@return current free-electron adjustment factor
	*/
	double get_factor();

private:
	/** Initialization with default parameters */
	void init();

	/** Scaling factors */
	double be_factor {1};
	double fe_factor {1};

	/** Models to be used */
	StopPow_Zimmerman * z {NULL};
	StopPow * fe {NULL};
	StopPow * fe2 {NULL};
	int fe_model {MODE_ZIMMERMAN};

};

} // end namespace StopPow
 #endif