#include "StopPow.h"

namespace StopPow
{

const float StopPow::DEFAULT_DX = 0.1; /* default step size for length-based calculations */
const float StopPow::DEFAULT_DRHOR = 0.1; /* default step size for areal-density calculations */
const int StopPow::MODE_LENGTH = 0; /* perform calculations as functions of length (um) */
const int StopPow::MODE_RHOR = 1; /* perform calculations as functions of rhoR (mg/cm2) */

/* Basic constructor, which simply sets dx and uses length as default mode*/
StopPow::StopPow()
{
	dx = DEFAULT_DX; // set step size
	mode = MODE_LENGTH;
	range_Emin = 0.02; // default 20 keV
}
/* Constructor which takes an initial mode */
StopPow::StopPow(int set_mode)
{
	dx = DEFAULT_DX; // set step size
	mode = set_mode;
	range_Emin = 0.02; // default 20 keV
}

/**
 * Calculate stopping power. Return units depend on mode.
 * @param E the particle energy in MeV
 * @return dE/dx in MeV/um [MeV/(mg/cm2)]
 * @throws invalid_argument
 */
float StopPow::dEdx(float E) throw(std::invalid_argument)
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
float StopPow::Eout(float E, float x) throw(std::invalid_argument)
{
	// sanity checking:
	if( E < get_Emin() || E > get_Emax() || x < 0 )
	{
		std::stringstream msg;
		msg << "Energies passed to StopPow::Eout are bad: " << E << "," << x;
		throw std::invalid_argument(msg.str());
	}
	
	float ret = E; // return value

	// iterate through total thickness.
	// if the energy is too low, stop looping
	for( float i = 0; i < x && ret >= get_Emin(); i+=dx )
	{
		ret += dx*dEdx(ret);
	}

	// Account for remainder:
	ret += fmod(x,dx)*dEdx(ret);

	// make sure we do not return a negative energy:
	return fmax( ret , 0.0 );
}

/**
 * Get incident energy for a particle.
 * @param E the particle energy in MeV
 * @param x thickness of material in um [mg/cm2]
 * @return initial particle energy in MeV
 */
float StopPow::Ein(float E, float x) throw(std::invalid_argument)
{
	// sanity checking:
	if( E < get_Emin() || E > get_Emax() || x < 0 )
	{
		std::stringstream msg;
		msg << "Args passed to StopPow::Ein are bad: " << E << "," << x;
		throw std::invalid_argument(msg.str());
	}

	float ret = E; // return value

	// iterate through total thickness.
	for( float i = 0; i < x && ret <= get_Emax(); i+=dx )
	{
		ret -= dx*dEdx(ret);
	}

	// Account for remainder:
	ret -= fmod(x,dx)*dEdx(ret);

	return ret;
}

/**
 * Get thickness of material traversed.
 * @param E1 the initial particle energy in MeV
 * @param E2 the final particle energy in MeV
 * @throws std::invalid_argument
 * @return material thickness in um [mg/cm2]
 */
float StopPow::Thickness(float E1, float E2) throw(std::invalid_argument)
{
	// sanity checking:
	if (E1 < get_Emin() || E1 > get_Emax() || 
		E2 < get_Emin() || E2 > get_Emax()
		 || E2 > E1)
	{
		std::stringstream msg;
		msg << "Energies passed to StopPow::Thickness are bad: " << E1 << "," << E2;
		throw std::invalid_argument(msg.str());
	}

	float ret = 0; // calculated thickness
	float E = E1; // working energy variable

	// iterate, increasing thickness, until
	// E is <= the "final particle energy" E2
	// also require that dEdx < 0 to prevent infinite loops
	while (E > E2 && dEdx(E) < 0)
	{
		E += dx*dEdx(E);
		ret += dx;
	}

	// Account for remainder:
	ret += (E2-E) / dEdx(E);

	return ret;
}

/**
 * Get the range of a particle with given energy
 * @param E the particle energy in MeV
 * @return range in um [mg/cm2]
 * @throws invalid_argument
*/
float StopPow::Range(float E) throw(std::invalid_argument)
{
	// sanity checking:
	if ( E < get_Emin() || E > get_Emax() )
	{
		std::stringstream msg;
		msg << "Energy passed to StopPow::Range is bad: " << E;
		throw std::invalid_argument(msg.str());
	}

	// use either 0 or Emin for the lower cutoff:
	float E2 = fmax( get_Emin() , range_Emin );
	
	// use the Thickness method:
	return Thickness(E,E2);
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
void StopPow::set_dx(float new_dx) throw(std::invalid_argument)
{
	// sanity checking:
	if (new_dx <= 0)
	{
		std::stringstream msg;
		msg << "Non-positive step size passed to StopPow::set_dx: " << new_dx;
		throw std::invalid_argument(msg.str());
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
		throw std::invalid_argument(msg.str());
	}
}

} // end namespace StopPow