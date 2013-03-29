/* Generic class for stopping power calculators. Also includes several generic methods.
 * @author Alex Zylstra
 * @date 2013/03/29
 */

#include "StopPow.h"

const float StopPow::DEFAULT_DX = 0.1; /* default step size for length-based calculations */
const float StopPow::DEFAULT_DRHOR = 0.1; /* default step size for areal-density calculations */
const int StopPow::MODE_LENGTH = 0; /* perform calculations as functions of length (um) */
const int StopPow::MODE_RHOR = 1; /* perform calculations as functions of rhoR (mg/cm2) */

/* Basic constructor, which simply sets dx and uses length as default mode*/
StopPow::StopPow()
{
	dx = DEFAULT_DX; // set step size
	mode = MODE_LENGTH;
}
/* Constructor which takes an initial mode */
StopPow::StopPow(int set_mode)
{
	dx = DEFAULT_DX; // set step size
	mode = set_mode;
}

/**
 * Calculate stopping power. Return units depend on mode.
 * @param E the particle energy in MeV
 * @return dE/dx in MeV/um [MeV/(mg/cm2)]
 * @throws invalid_argument
 */
float StopPow::dEdx(float E)
{
	if( mode == MODE_LENGTH )
		return dEdx_MeV_um(E);
	if( mode == MODE_RHOR )
		return dEdx_MeV_mgcm2(E);
	return 0;
}

/**
 * Get energy downshift for a particle.
 * @param E the particle energy in MeV
 * @param x thickness of material in um [mg/cm2]
 * @return final particle energy in MeV
 */
float StopPow::Eout(float E, float x)
{
	// sanity checking:
	if( E < 0 || x < 0)
	{
		std::stringstream msg;
		msg << "Energies passed to StopPow::Eout are bad: " << E << "," << x;
		throw new std::invalid_argument(msg.str());
	}
	
	float ret = E; // return value

	// iterate through total thickness.
	// if the energy is too low, stop looping
	for( float i = 0; i < x && ret > 0; i+=dx )
	{
		ret += dx*dEdx(ret);
	}

	// make sure we do not return a negative energy:
	return fmax( ret , 0.0 );
}




/**
 * Get incident energy for a particle.
 * @param E the particle energy in MeV
 * @param x thickness of material in um [mg/cm2]
 * @return initial particle energy in MeV
 */
float StopPow::Ein(float E, float x)
{
	// sanity checking:
	if( E < 0 || x < 0)
	{
		std::stringstream msg;
		msg << "Args passed to StopPow::Ein are bad: " << E << "," << x;
		throw new std::invalid_argument(msg.str());
	}

	float ret = E; // return value

	// iterate through total thickness.
	for( float i = 0; i < x; i+=dx )
	{
		ret -= dx*dEdx(ret);
	}

	return ret;
}

/**
 * Get thickness of material traversed.
 * @param E1 the initial particle energy in MeV
 * @param E2 the final particle energy in MeV
 * @throws std::invalid_argument
 * @return material thickness in um [mg/cm2]
 */
float StopPow::Thickness(float E1, float E2)
{
	// sanity checking:
	if (E1 <= 0 || E2 <= 0 || E2 > E1)
	{
		std::stringstream msg;
		msg << "Energies passed to StopPow::Thickness are bad: " << E1 << "," << E2;
		throw new std::invalid_argument(msg.str());
	}

	float ret = 0; // calculated thickness
	float E = E1; // working energy variable

	// iterate, increasing thickness, until
	// E is <= the "final particle energy" E2
	while (E > E2)
	{
		E += dx*dEdx(E);
		ret += dx;
	}

	return ret;
}

/** Get the current step sized being used for calculations.
* @return dx the step size in um [mg/cm2]
*/
float StopPow::get_dx()
{
	return dx;
}
/** Set the step size for calculations
* @param new_dx the new step size to use, in um [mg/cm2]
 * @throws std::invalid_argument
*/
void StopPow::set_dx(float new_dx)
{
	// sanity checking:
	if (new_dx <= 0)
	{
		std::stringstream msg;
		msg << "Non-positive step size passed to StopPow::set_dx: " << new_dx;
		throw new std::invalid_argument(msg.str());
	}

	dx = new_dx;
}

/** Get the current mode being used for calculations.
 * @return mode Either StopPow.MODE_LENGTH or StopPow.MODE_RHOR
 */
int StopPow::get_mode()
{
	return mode;
}
/** Set the mode for calculations
 * @param new_mode Either StopPow.MODE_LENGTH or StopPow.MODE_RHOR
 * @throws invalid_argument
 */
void StopPow::set_mode(int new_mode)
{
	if( new_mode == MODE_LENGTH )
	{
		mode = MODE_LENGTH;
		dx = DEFAULT_DX;
	}
	else if ( new_mode == MODE_RHOR )
	{
		mode = MODE_RHOR;
		dx = DEFAULT_DRHOR;
	}
	// if we get inside this else, new_mode was invalid:
	else
	{
		std::stringstream msg;
		msg << "Invalid mode passed to StopPow::set_mode: " << new_mode;
		throw new std::invalid_argument(msg.str());
	}
}