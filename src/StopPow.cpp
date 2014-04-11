#include "StopPow.h"

namespace StopPow
{

const int StopPow::MODE_LENGTH = 0; /* perform calculations as functions of length (um) */
const int StopPow::MODE_RHOR = 1; /* perform calculations as functions of rhoR (mg/cm2) */

/* Basic constructor, which simply sets dx and uses length as default mode*/
/*StopPow::StopPow()
 : StopPow(MODE_LENGTH)
{}*/ // commented out to allow compilation with gcc 4.6.x
 StopPow::StopPow()
{
	mode = MODE_LENGTH;

	// set the default type and info strings to empty:
	model_type = "";
	info = "";
}
/* Constructor which takes an initial mode */
StopPow::StopPow(int set_mode)
{
	mode = set_mode;

	// set the default type and info strings to empty:
	model_type = "";
	info = "";
}

// Calculate stopping power:
double StopPow::dEdx(double E) throw(std::invalid_argument)
{
	// return depending on mode:
	if( mode == MODE_LENGTH )
		return dEdx_MeV_um(E);
	if( mode == MODE_RHOR )
		return dEdx_MeV_mgcm2(E);
	// return NAN if mode selection fails
	return std::numeric_limits<double>::quiet_NaN();
}

// set up function for GSL, used in calculating meta stuff (Eout, Thickness, Range)
auto Eout_func = [] (double t, const double y[], double dydt[], void * params)
{
	StopPow::StopPow * s = (StopPow::StopPow *)params;
	dydt[0] = s->dEdx(y[0]);
	return (int)GSL_SUCCESS;
};

// set up function for GSL, used in calculating meta stuff (Ein)
auto Ein_func = [] (double t, const double y[], double dydt[], void * params)
{
	StopPow::StopPow * s = (StopPow::StopPow *)params;
	dydt[0] = -1.*(s->dEdx(y[0]));
	return (int)GSL_SUCCESS;
};

// Calculate energy downshift:
double StopPow::Eout(double E, double x) throw(std::invalid_argument, std::runtime_error)
{
	// sanity checking:
	if( E < get_Emin() || E > get_Emax() || x < 0 )
	{
		std::stringstream msg;
		msg << "Energies passed to StopPow::Eout are bad: " << E << "," << x;
		throw std::invalid_argument(msg.str());
	}
	
	// set up GSL ODE solver
	gsl_odeiv2_system sys = {Eout_func, NULL, 1, this};
	gsl_odeiv2_driver * d = 
		gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk4,
				  1e-6, 1e-6, 0.0);
	
	double x0 = 0.; double x1 = x; // range of ODE integration
	double y[1] = { E };
	int status;
	try
	{
		status = gsl_odeiv2_driver_apply(d, &x0, x1, y);
	}
	catch(std::invalid_argument e)
	{
		return 0; // got to the end of the range
	}

	// check for errors:
	if( status != GSL_SUCCESS )
	{
		throw std::runtime_error::runtime_error("GSL RK4 ODE integration failed in StopPow::Eout!");
	}

	gsl_odeiv2_driver_free (d);

	// make sure we do not return a negative energy:
	return fmax( y[0] , 0.0 );
}

// Calculate energy upshift
double StopPow::Ein(double E, double x) throw(std::invalid_argument, std::runtime_error)
{
	// sanity checking:
	if( E < get_Emin() || E > get_Emax() || x < 0 )
	{
		std::stringstream msg;
		msg << "Args passed to StopPow::Ein are bad: " << E << "," << x;
		throw std::invalid_argument(msg.str());
	}
	
	// set up GSL ODE solver
	gsl_odeiv2_system sys = {Ein_func, NULL, 1, this};
	gsl_odeiv2_driver * d = 
		gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk4,
				  1e-6, 1e-6, 0.0);
	
	double x0 = 0.; double x1 = x; // range of ODE integration
	double y[1] = { E };
	int status;
	// run the solver:
	try
	{
		status = gsl_odeiv2_driver_apply(d, &x0, x1, y);
	}
	catch(std::invalid_argument e)
	{
		return get_Emax(); // got to Emax
	}

	// check for errors:
	if( status != GSL_SUCCESS )
	{
		throw std::runtime_error::runtime_error("GSL RK4 ODE integration failed in StopPow::Ein!");
	}

	gsl_odeiv2_driver_free (d);

	return y[0];
}

// Calculate thickness of material traversed for a given energy shift
double StopPow::Thickness(double E1, double E2) throw(std::invalid_argument)
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

	// ODE system to solve:
	gsl_odeiv2_system sys = {Eout_func, NULL, 1, this};

	// set up GSL ODE solver: stepping done manually for thickness
	const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rk4;
	gsl_odeiv2_step * step = gsl_odeiv2_step_alloc (T, 1);
	gsl_odeiv2_control * c = gsl_odeiv2_control_y_new (1e-6, 0.0);
	gsl_odeiv2_evolve * e = gsl_odeiv2_evolve_alloc (1);

	int status; double x = 0;
	// step size for thickness iteration, corresponds to 50 keV change
	double dx = -0.05 / dEdx(E1);
	// step size for RK ODE solver, set at 1/100 of previous
	double h = dx / 100.; 
	double y[1] = { E1 };
	double y_last = E1;

	// Loop until we overshoot, i.e. energy calculated becomes lower than E2
	do
	{
		dx = -0.05 / dEdx(y_last);
		y_last = y[0];
		try
		{
			status = gsl_odeiv2_evolve_apply (e, c, step, &sys, &x, x+dx, &h, y);
		}
		catch(std::invalid_argument e)
		{
			break; // if we get to the end of the particle's range, this is reached
		}

		// check for errors:
		if( status != GSL_SUCCESS )
			throw std::runtime_error::runtime_error("GSL RK4 ODE integration failed in StopPow::Ein!");
	} while( y[0] > E2 );

	// Do a linear interpolation between current point and previous point to get
	// the most accurate value of thickness:
	double slope = dEdx(y[0]);
	double thick = x + (E2-y[0])/slope;

	return thick;
}

// Calculate the range of a particle with given energy
double StopPow::Range(double E) throw(std::invalid_argument)
{
	// sanity checking:
	if ( E < get_Emin() || E > get_Emax() )
	{
		std::stringstream msg;
		msg << "Energy passed to StopPow::Range is bad: " << E;
		throw std::invalid_argument(msg.str());
	}

	// use either self-defined range cutoff or Emin for the lower cutoff:
	double E2 = fmax( get_Emin() , 0 );

	// sanity check:
	if( E <= E2 )
		return 0;

	// use the Thickness method:
	return Thickness(E,E2);
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
	}
	else if ( new_mode == MODE_RHOR )
	{
		mode = MODE_RHOR;
	}
	// if we get inside this else, new_mode was invalid:
	else
	{
		std::stringstream msg;
		msg << "Invalid mode passed to StopPow::set_mode: " << new_mode;
		throw std::invalid_argument(msg.str());
	}
}


// Get the type of model
std::string StopPow::get_type()
{
	return model_type;
}

// Get an info string for this model
std::string StopPow::get_info()
{
	return info;
}

} // end namespace StopPow
