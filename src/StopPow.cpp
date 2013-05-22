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
}
/* Constructor which takes an initial mode */
StopPow::StopPow(int set_mode)
{
	dx = DEFAULT_DX; // set step size
	mode = set_mode;
}

// Calculate stopping power:
float StopPow::dEdx(float E) throw(std::invalid_argument)
{
	// return depending on mode:
	if( mode == MODE_LENGTH )
		return dEdx_MeV_um(E);
	if( mode == MODE_RHOR )
		return dEdx_MeV_mgcm2(E);
	// return NAN if mode selection fails
	return std::numeric_limits<float>::quiet_NaN();
}

// Calculate energy downshift:
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
	for( float i = 0; i < x && ret >= get_Emin() && dEdx(ret) < 0; i+=dx )
	{
		ret += dx*dEdx(ret);
	}

	// Account for remainder if possible:
	if( ret > get_Emin() )
		ret += fmod(x,dx)*dEdx(ret);
	else // if we've dropped below min energy, call it zero:
		ret = 0;

	// make sure we do not return a negative energy:
	return fmax( ret , 0.0 );
}

// Calculate energy upshift
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
	for( float i = 0; i < x && ret <= get_Emax() && dEdx(ret) < 0; i+=dx )
	{
		ret -= dx*dEdx(ret);
	}

	// Account for remainder if possible:
	if( ret < get_Emax() )
		ret -= fmod(x,dx)*dEdx(ret);
	else // energy got too big
		ret = std::numeric_limits<float>::quiet_NaN();

	return ret;
}

// Calculate thickness of material traversed for a given energy shift
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
	// and that E stay within E limits
	while (E > E2 && 
		E > get_Emin() && E < get_Emax()
		&& dEdx(E) < 0)
	{
		E += dx*dEdx(E);
		ret += dx;
	}

	// Account for remainder if possible:
	//if( E > get_Emin() && E < get_Emax() )
	//	ret += (E2-E) / dEdx(E);

	return ret;
}

// Calculate the range of a particle with given energy
float StopPow::Range(float E) throw(std::invalid_argument)
{
	// sanity checking:
	if ( E < get_Emin() || E > get_Emax() )
	{
		std::stringstream msg;
		msg << "Energy passed to StopPow::Range is bad: " << E;
		throw std::invalid_argument(msg.str());
	}

	// use either self-defined range cutoff or Emin for the lower cutoff:
	float E2 = fmax( get_Emin() , 0 );
	
	// sanity check:
	if( E <= E2 )
		return 0;

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