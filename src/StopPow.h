/** 
 * @brief Generic class for stopping power calculators.
 *
 * In addition to setting the abstract template for stopping power calculators, this also includes several generic methods.
 * The stopping power utilities here can be called as functions of linear distance or areal density. To specify which,
 * the mode must be set correctly.
 *
 * @class StopPow::StopPow
 * @author Alex Zylstra
 * @date 2013/06/04
 * @copyright MIT / Alex Zylstra
 */

#ifndef STOPPOW_H
#define STOPPOW_H

#include <math.h>

#include <stdexcept>
#include <sstream>
#include <limits>

/** @namespace StopPow */
namespace StopPow
{

class StopPow
{
public:
	/** Simple constructor for the generic class */
	StopPow();

	/**
	 * Virtual destructor is necessary so that deriving classes
	 * can do garbage collection when a StopPow pointer is deleted.
	 */
	virtual ~StopPow(){};

	/**
	 * Construct a new StopPow object given a starting mode
	 * @param set_mode the mode you want to use (defined using class constants)
	 */
	explicit StopPow(int set_mode);

	/**
	 * Calculate stopping power. Return units depend on mode.
	 * @param E the particle energy in MeV
	 * @return dE/dx in MeV/um [MeV/(mg/cm2)]
 	 * @throws invalid_argument
	 */
	double dEdx(double E) throw(std::invalid_argument);

	/* Extending classes must implement these two dEdx functions: */
	virtual double dEdx_MeV_um(double E) = 0;
	virtual double dEdx_MeV_mgcm2(double E) = 0;
	/* Extending classes must also implement defined (inclusive) energy limits: */
	virtual double get_Emin() = 0;
	virtual double get_Emax() = 0;

	/**
	  * Get the type of stopping power model described by this class
	  * @return a std::string type descriptor
	  */
	std::string get_type();

	/**
	  * Get some information about the model
	  * @return a std::string containing useful info (in this case, the file name of SRIM data used)
	  */
	std::string get_info();

 	/**
	 * Get energy downshift for a particle. If the particle energy
	 * drops below the model's minimum energy, the particle is considered
	 * ranged out and this method returns 0.
	 * @param E the particle energy in MeV
	 * @param x thickness of material in um [mg/cm2]
	 * @return final particle energy in MeV
	 */
	double Eout(double E, double x) throw(std::invalid_argument);

 	/**
	 * Get incident energy for a particle. If the particle energy
	 * goes above the model's maximum energy, then quiet not-a-number is
	 * returned, i.e. {@code std::numeric_limits<double>::quiet_NaN()}.
	 * @param E the particle energy in MeV
	 * @param x thickness of material in um [mg/cm2]
	 * @return initial particle energy in MeV
	 */
	double Ein(double E, double x) throw(std::invalid_argument);

 	/**
	 * Get thickness of material traversed.
	 * @param E1 the initial particle energy in MeV
	 * @param E2 the final particle energy in MeV
	 * @throws std::invalid_argument
	 * @return material thickness in um [mg/cm2]
	 */
	double Thickness(double E1, double E2) throw(std::invalid_argument);

	/**
	 * Get the range of a particle with given energy
	 * @param E the particle energy in MeV
	 * @return range in um [mg/cm2]
	 * @throws invalid_argument
	*/
	double Range(double E) throw(std::invalid_argument);

	/** Get the current step sized being used for calculations.
	 * @return dx the step size in um [mg/cm2]
	 */
	double get_dx();
	/** Set the step size for calculations
	 * @param new_dx the new step size to use, in um [mg/cm2]
 	 * @throws invalid_argument
	 */
	void set_dx(double new_dx) throw(std::invalid_argument);

	/** Get the current mode being used for calculations.
	 * @return mode Either StopPow.MODE_LENGTH or StopPow.MODE_RHOR
	 */
	int get_mode();
	/** Set the mode for calculations
	 * @param new_mode Either StopPow.MODE_LENGTH or StopPow.MODE_RHOR
 	 * @throws invalid_argument
	 */
	void set_mode(int new_mode);

	/** default step size for length-based calculations */
	static const double DEFAULT_DX;
	/** default step size for areal-density calculations */
	static const double DEFAULT_DRHOR;
	/** perform calculations as functions of length (um) */
	static const int MODE_LENGTH;
	/** perform calculations as functions of rhoR (mg/cm2) */
	static const int MODE_RHOR;

protected:
	/** step size in um [mg/cm2] */
	double dx;
	/** current mode for calculations */
	int mode;

	/** Represent the type of model described by this class */
	std::string model_type;
	/** Some information about the model, stored as string */
	std::string info;
};

} // end namespace StopPow

#endif
