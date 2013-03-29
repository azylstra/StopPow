/** 
 * @brief Generic class for stopping power calculators. 
 * 
 * In addition to setting the abstract template for stopping power calculators, this also includes several generic methods.
 * The stopping power utilities here can be called as functions of linear distance or areal density. To specify which,
 * the mode must be set correctly.
 * 
 * @class StopPow
 * @author Alex Zylstra
 * @date 2013/03/29
 */

#ifndef STOPPOW_H
#define STOPPOW_H

#include <math.h>
#include <stdexcept>
#include <sstream>

using namespace std;

class StopPow
{
public:
	/** Simple constructor for the generic class */
	StopPow();
	/**
	 * Construct a new StopPow object given a starting mode
	 * @param set_mode the mode you want to use (defined using class constants)
	 */
	StopPow(int set_mode);

	/**
	 * Calculate stopping power. Return units depend on mode.
	 * @param E the particle energy in MeV
	 * @return dE/dx in MeV/um [MeV/(mg/cm2)]
 	 * @throws invalid_argument
	 */
	float dEdx(float E);

	/* Extending classes must implement these two dEdx functions: */
	virtual float dEdx_MeV_um(float E) = 0;
	virtual float dEdx_MeV_mgcm2(float E) = 0;

 	/**
 	 * Get energy downshift for a particle.
 	 * @param E the particle energy in MeV
 	 * @param x thickness of material in um [mg/cm2]
 	 * @return final particle energy in MeV
 	 * @throws invalid_argument
 	 */
	float Eout(float E, float x);

 	/**
 	 * Get incident energy for a particle.
 	 * @param E the particle energy in MeV
 	 * @param x thickness of material in um [mg/cm2]
 	 * @return initial particle energy in MeV
 	 * @throws invalid_argument
 	 */
	float Ein(float E, float x);

 	/**
 	 * Get thickness of material traversed.
 	 * @param E1 the initial particle energy in MeV
 	 * @param E2 the final particle energy in MeV
 	 * @return material thickness in um [mg/cm2]
 	 * @throws invalid_argument
 	 */
	float Thickness(float E1, float E2);

	/** Get the current step sized being used for calculations.
	 * @return dx the step size in um [mg/cm2]
	 */
	float get_dx();
	/** Set the step size for calculations
	 * @param new_dx the new step size to use, in um [mg/cm2]
 	 * @throws invalid_argument
	 */
	void set_dx(float new_dx);

	/** Get the current mode being used for calculations.
	 * @return mode Either StopPow.MODE_LENGTH or StopPow.MODE_RHOR
	 */
	int get_mode();
	/** Set the mode for calculations
	 * @param new_mode Either StopPow.MODE_LENGTH or StopPow.MODE_RHOR
 	 * @throws invalid_argument
	 */
	void set_mode(int new_mode);

	static const float DEFAULT_DX; /* default step size for length-based calculations */
	static const float DEFAULT_DRHOR; /* default step size for areal-density calculations */
	static const int MODE_LENGTH; /* perform calculations as functions of length (um) */
	static const int MODE_RHOR; /* perform calculations as functions of rhoR (mg/cm2) */

protected:
	float dx; // step size
	int mode; // current mode for calculations
};

#endif