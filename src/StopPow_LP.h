/**
 * @brief Calculate Li-Petrasso stopping power.
 * 
 * Implement a stopping-power calculator for plasma, using
 * the Fokker-Planck theory described in: C.K. Li and R.D. Petrasso, Phys. Rev. Lett. 70, 3059 (1993).
 *
 * @class StopPow::StopPow_LP
 * @author Alex Zylstra
 * @date 2014/04/01
 * @copyright Alex Zylstra / MIT
 */

#ifndef STOPPOW_LP_H
#define STOPPOW_LP_H

#include <math.h>

#include <vector>
#include <stdexcept>
#include <iostream>

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_fermi_dirac.h>

#include "StopPow_Plasma.h"
#include "StopPow_Constants.h"

namespace StopPow
{

class StopPow_LP : public StopPow_Plasma
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
	StopPow_LP(float mt, float Zt, std::vector<float> & mf , std::vector<float> & Zf, std::vector<float> & Tf, std::vector<float> & nf) throw(std::invalid_argument);

	/** Initialize the stopping power. Electrons should be included!
	 * @param mt the test particle mass in AMU
	 * @param Zt the test particle in charge (units of e)
	 * @param field vector containing field ion info. Each row is one type of ion, then the array must contain:
	 * [mf,Zf,Tf,nf] in units of AMU, e, keV, and 1/cc
 	 * @throws invalid_argument
	 */
	StopPow_LP(float mt, float Zt, std::vector< std::array<float,4> > & field) throw(std::invalid_argument);

	/** Initialize the Li-Petrasso stopping power. Electrons should not be included - they will be added automatically!
	 * @param mt the test particle mass in AMU
	 * @param Zt the test particle in charge (units of e)
	 * @param mf vector containing ordered field ion masses in AMU
	 * @param Zf vector containing ordered field ion charges in units of e
	 * @param Tf vector containing ordered field ion temperatures in units of keV
	 * @param nf vector containing ordered field ion densities in units of 1/cm3
	 * @param Te the electron temperature in keV
 	 * @throws invalid_argument
	 */
	StopPow_LP(float mt, float Zt, std::vector<float> & mf , std::vector<float> & Zf, std::vector<float> & Tf, std::vector<float> & nf, float Te) throw(std::invalid_argument);

	/** Initialize the stopping power. Electrons should not be included - they will be added automatically!
	 * @param mt the test particle mass in AMU
	 * @param Zt the test particle in charge (units of e)
	 * @param field vector containing field ion info. Each row is one type of ion, then the array must contain:
	 * [mf,Zf,Tf,nf] in units of AMU, e, keV, and 1/cc
	 * @param Te the electron temperature in keV
 	 * @throws invalid_argument
	 */
	StopPow_LP(float mt, float Zt, std::vector< std::array<float,4> > & field, float Te) throw(std::invalid_argument);

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

	/** Turn quantum effective temperature correction on or off.
	 * @param set if you want to use quantum correction
	 */
	void set_quantum(bool set);

	/** Set the factor for for the thermal velocity in x^{t/f}, i.e. vf = sqrt(a*kB*Tf/mf)
	 * WARNING: Changing this not recommended unless you know what you're doing
	 * @param a the constant for the field particle velocity in x^{t/f} for binary collisions
	 */
	void set_xtf_factor(float a);

	/** Set the factor for for the thermal velocity for collective effects in x^{t/f}, i.e. vf = sqrt(a*kB*Tf/mf)
	 * WARNING: Changing this not recommended unless you know what you're doing
	 * @param a the constant for the field particle velocity in x^{t/f} for collective effects
	 */
	void set_xtf_collective_factor(float a);

	/** Set the factor for for the thermal velocity for calculating relative velocity, i.e. vf = sqrt(a*kB*Tf/mf)
	 * WARNING: Changing this not recommended unless you know what you're doing
	 * @param a the constant for the field particle velocity in x^{t/f} for calculating relative velocity
	 */
	void set_u_factor(float a);

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
	/** Initialization routine (beyond what is done by superclass constructor) */
	void init();

	// Some configurable options
	/** Use collective effects? */
	bool collective {true}; 
	/** Use quantum temperature correction? */
	bool quantumT {true};
	/** Factor for thermal velocity in binary collision x^{t/f} */
	float xtf_factor {2.0};
	/** Factor for thermal velocity in collective effects x^{t/f} */
	float xtf_collective_factor {1.0};
	/** Factor for thermal velocity in collective effects x^{t/f} */
	float u_factor {8./M_PI};

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

	/** Field particle thermal velocity. This one defaults to using 2kT/m.
	 * @param index the field particle's index (for mf,Zf,Tf,nf arrays)
	 * @return thermal velocity in cm/s
	 */
	float vtf(int index);

	/** Field particle thermal velocity with specified constant. EG:
	 *		For problems that use kT/m (eg lDebye*wpe), use constant=1
	 *		For most probable speed, use constant=2
	 *		For Maxwell-averaged, use constant=8/Pi
	 *		For RMS velocity, use constant=3
	 * @param index the field particle's index (for mf,Zf,Tf,nf arrays)
	 * @param constant the multiplicitive factor. vtf = c*sqrt(constant*k*T/m)
	 * @return thermal velocity in cm/s
	 */
	float vtf(int index, float constant);

	/** x^{t/f} parameter from Li 1993. For Chandrasekhar function.
	 * @param E test particle energy in MeV
	 * @param index the field particle's index (for mf,Zf,Tf,nf arrays)
	 * @return x^{t/f}
	 */
	float xtf(float E, int index);

	/** x^{t/f} parameter from Li 1993. For collective effects.
	 * @param E test particle energy in MeV
	 * @param index the field particle's index (for mf,Zf,Tf,nf arrays)
	 * @return x^{t/f}
	 */
	float xtf_collective(float E, int index);

	/** Relative velocity between test particle and field particle
	 * @param E test particle energy in MeV
	 * @param index the field particle's index (for mf,Zf,Tf,nf arrays)
	 * @return relative velocity in cm/s
	 */
	float u(float E, int index);

	/** Temperature to use for calculations, taking into account quantum correction (or not)
	* @param index the field particle's index
	* @return effective temperature in keV
	*/
	float Tq(int index);

	/* Minimum energy for dE/dx calculations */
	static const float Emin; 
	/* Maximum energy for dE/dx calculations */
	static const float Emax; 
};

} // end namespace StopPow

#endif