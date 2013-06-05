
#include "PlotGen.h"

namespace StopPow
{

// Generate a plot of dEdx vs E for a given model. Results stored in data.
bool get_dEdx_vs_E( 	StopPow & model ,
						std::vector< std::vector<float> > & data )
{
	// call overloaded function using default step size:
	return get_dEdx_vs_E(	model,
							PLOT_DEFAULT_NUM_POINTS,
							data );
}
// Generate a plot of dEdx vs E for a given model. Results stored in data.
bool get_dEdx_vs_E( 	StopPow & model ,
						int num_points ,
						std::vector< std::vector<float> > & data )
{
	// call overloaded function using model's energy limits:
	return get_dEdx_vs_E(	model,
							model.get_Emin() ,
							model.get_Emax() ,
							num_points ,
							data );
}
// Generate a plot of dEdx vs E for a given model. Results stored in data.
bool get_dEdx_vs_E( 	StopPow & model ,
						float Emin ,
						float Emax ,
						int num_points ,
						std::vector< std::vector<float> > & data )
{
	// clear data
	data.clear();

	// set up 2 vectors for E and dE/dx
	std::vector<float> E_vector; std::vector<float> dEdx_vector;
	data.push_back( E_vector );
	data.push_back( dEdx_vector );

	// step size:
	float dE = (Emax-Emin) / ((float)num_points);

	// now generate the plot values:
	for(float E = Emin; E <= Emax; E+=dE )
	{
		try
		{
			data[0].push_back(E);
			data[1].push_back( model.dEdx(E) );
		}
		catch( std::exception e )
		{
			// if an error occurs, abort and return false:
			return false;
		}
	}
	
	return true;
}

// generate a plot of range vs E for a given model. Results stored in data
bool get_Range_vs_E( 	StopPow & model ,
						std::vector< std::vector<float> > & data )
{
	return get_Range_vs_E(	model ,
							PLOT_DEFAULT_NUM_POINTS ,
							data );
}

// generate a plot of range vs E for a given model. Results stored in data
bool get_Range_vs_E( 	StopPow & model ,
						int num_points ,
						std::vector< std::vector<float> > & data )
{
	return get_Range_vs_E(	model,
							model.get_Emin(),
							model.get_Emax(),
							num_points,
							data );
}

// generate a plot of range vs E for a given model. Results stored in data
bool get_Range_vs_E( 	StopPow & model ,
						float Emin ,
						float Emax ,
						int num_points ,
						std::vector< std::vector<float> > & data )
{
	// clear data
	data.clear();

	// set up 2 vectors for E and dE/dx
	std::vector<float> E_vector; std::vector<float> Range_vector;
	data.push_back( E_vector );
	data.push_back( Range_vector );

	// step size:
	float dE = (Emax-Emin) / ((float)num_points);

	// now generate the plot values:
	for(float E = Emin; E <= Emax; E+=dE )
	{
		try
		{
			// calculate first element directly:
			if( data[0].size() == 0 )
			{
				data[0].push_back(E);
				data[1].push_back( model.Range(E) );
			}
			else
			{
				float Range_last = data[1].back();

				data[0].push_back(E);
				// instead of calculating whole range, use previous result
				// to speed up calculation
				float range = Range_last + model.Thickness(E, fmax(E-dE,Emin) );
				data[1].push_back( range );
			}
		}
		catch( std::exception e )
		{
			// if an error occurs, abort and return false:
			return false;
		}
	}

	return true;
}

// Eout vs Ein given thickness
bool get_Eout_vs_Ein(	StopPow & model ,
						float Thickness ,
						std::vector< std::vector<float> > & data )
{	
	return get_Eout_vs_Ein(	model,
							PLOT_DEFAULT_NUM_POINTS,
							Thickness,
							data );
}

// Eout vs Ein given thickness
bool get_Eout_vs_Ein(	StopPow & model ,
						int num_points ,
						float Thickness ,
						std::vector< std::vector<float> > & data )
{	
	// calculate limits. Lower cutoff based on particle
	// that just makes it through
	float Emin, Emax;
	try
	{
		Emin = model.Ein(model.get_Emin(),Thickness);
		Emax = model.get_Emax();
	}
	catch( std::exception e)
	{
		return false;
	}
	return get_Eout_vs_Ein(	model,
							Emin,
							Emax,
							num_points,
							Thickness,
							data );
}

// Eout vs Ein given thickness
bool get_Eout_vs_Ein(	StopPow & model ,
						float Emin ,
						float Emax ,
						int num_points ,
						float Thickness ,
						std::vector< std::vector<float> > & data )
{	
	// clear data
	data.clear();

	// set up 2 vectors for E and dE/dx
	std::vector<float> E_vector; std::vector<float> Eout_vector;
	data.push_back( E_vector );
	data.push_back( Eout_vector );

	// step size:
	float dE = (Emax-Emin) / ((float)num_points);

	// now generate the plot values:
	for(float Ein = Emin; Ein <= Emax; Ein+=dE )
	{
		try
		{
			data[0].push_back(Ein);
			data[1].push_back( model.Eout(Ein,Thickness) );
		}
		catch( std::exception e )
		{
			// if an error occurs, abort and return false:
			return false;
		}
	}

	return true;
}

// Calculate Eout vs thickness for given input energy
bool get_Eout_vs_Thickness(	StopPow & model ,
							float Ein ,
							std::vector< std::vector<float> > & data )
{
	return get_Eout_vs_Thickness(	model,
									PLOT_DEFAULT_NUM_POINTS,
									Ein,
									data );
}

// Calculate Eout vs thickness for given input energy
bool get_Eout_vs_Thickness(	StopPow & model ,
							int num_points ,
							float Ein ,
							std::vector< std::vector<float> > & data )
{
	// calculate limits. Lower cutoff based on particle
	// that just makes it through
	float Tmin, Tmax;
	try
	{
		Tmin = 0;
		Tmax = model.Range(Ein);
	}
	catch( std::exception e)
	{
		return false;
	}
	return get_Eout_vs_Thickness(	model,
									Tmin,
									Tmax,
									num_points,
									Ein,
									data );
}

// Calculate Eout vs thickness for given input energy
bool get_Eout_vs_Thickness(	StopPow & model ,
							float Tmin ,
							float Tmax ,
							int num_points ,
							float Ein ,
							std::vector< std::vector<float> > & data )
{
	// clear data
	data.clear();

	// set up 2 vectors for E and dE/dx
	std::vector<float> E_vector; std::vector<float> Thickness_vector;
	data.push_back( E_vector );
	data.push_back( Thickness_vector );

	// step size:
	float dT = (Tmax-Tmin) / ((float)num_points);

	// now generate the plot values:
	for(float T = Tmin; T <= Tmax; T+=dT )
	{
		try
		{
			// first data point:
			if(data[0].size() == 0)
			{
				data[0].push_back(T);
				data[1].push_back( model.Eout(Ein,T) );
			}
			else
			{
				float Elast = data[1].back();
				// do an incremental calculation, which is faster
				data[0].push_back(T);
				data[1].push_back( model.Eout(Elast,dT) );
			}
		}
		catch( std::exception e )
		{
			// if an error occurs, abort and return false:
			return false;
		}
	}

	return true;
}

// Get dataset for Ein vs Eout
bool get_Ein_vs_Eout(	StopPow & model ,
						float Thickness ,
						std::vector< std::vector<float> > & data )
{
	return get_Ein_vs_Eout(	model,
							PLOT_DEFAULT_NUM_POINTS,
							Thickness,
							data );
}

// Get dataset for Ein vs Eout
bool get_Ein_vs_Eout(	StopPow & model ,
						int num_points ,
						float Thickness ,
						std::vector< std::vector<float> > & data )
{
	// calculate energy limits
	// lower limit = particle almost ranged out
	// upper limit = Emax shifted through Thickness
	float Emin, Emax;
	try
	{
		Emin = model.get_Emin();
		Emax = model.Eout( model.get_Emax() , Thickness );
	}
	catch( std::exception e )
	{
		return false;
	}
	return get_Ein_vs_Eout(	model,
							Emin,
							Emax,
							num_points,
							Thickness,
							data );

}

// Get dataset for Ein vs Eout
bool get_Ein_vs_Eout(	StopPow & model ,
						float Emin ,
						float Emax ,
						int num_points ,
						float Thickness ,
						std::vector< std::vector<float> > & data )
{	
	// clear data
	data.clear();

	// set up 2 vectors
	std::vector<float> E_vector; std::vector<float> Ein_vector;
	data.push_back( E_vector );
	data.push_back( Ein_vector );

	// step size:
	float dE = (Emax-Emin) / ((float)num_points);

	// now generate the plot values:
	for(float Eout = Emin; Eout <= Emax; Eout+=dE )
	{
		try
		{
			data[0].push_back(Eout);
			data[1].push_back( model.Ein( fmin(Eout,Emax),Thickness) );
		}
		catch( std::exception e )
		{
			// if an error occurs, abort and return false:
			return false;
		}
	}

	return true;

}

// Get dataset for Ein vs thickness
bool get_Ein_vs_Thickness(	StopPow & model ,
							float Eout ,
							std::vector< std::vector<float> > & data )
{
	return get_Ein_vs_Thickness(model,
								PLOT_DEFAULT_NUM_POINTS,
								Eout,
								data );
}

// Get dataset for Ein vs thickness
bool get_Ein_vs_Thickness(	StopPow & model ,
							int num_points ,
							float Eout ,
							std::vector< std::vector<float> > & data )
{
	// calculate limits. Lower cutoff is zero
	// Upper cutoff corresponds to thickness that ranges
	// the max model energy to the specified Eout
	float Tmin, Tmax;
	try
	{
		Tmin = 0;
		Tmax = model.Thickness(model.get_Emax(),Eout);
	}
	catch( std::exception e)
	{
		return false;
	}
	return get_Ein_vs_Thickness(	model,
									Tmin,
									Tmax,
									num_points,
									Eout,
									data );
}

// Get dataset for Ein vs thickness
bool get_Ein_vs_Thickness(	StopPow & model ,
							float Tmin ,
							float Tmax ,
							int num_points ,
							float Eout ,
							std::vector< std::vector<float> > & data )
{
	// clear data
	data.clear();

	// set up 2 vectors for E and dE/dx
	std::vector<float> E_vector; std::vector<float> Thickness_vector;
	data.push_back( E_vector );
	data.push_back( Thickness_vector );

	// step size:
	float dT = (Tmax-Tmin) / ((float)num_points);

	// now generate the plot values:
	for(float T = Tmin; T <= Tmax; T+=dT )
	{
		try
		{
			if( data[0].size() == 0 )
			{
				data[0].push_back(T);
				data[1].push_back( model.Ein(Eout,T) );
			}
			else
			{
				float Ein_last = data[1].back();
				float Ein = model.Ein(Ein_last,dT);
				data[0].push_back(T);
				data[1].push_back(Ein);
			}
		}
		catch( std::exception e )
		{
			// if an error occurs, abort and return false:
			return false;
		}
	}

	return true;
}

// Get data for thickess and a function of energy out given energy in
bool get_Thickness_vs_Eout(	StopPow & model ,
							float Ein ,
							std::vector< std::vector<float> > & data )
{
	return get_Thickness_vs_Eout(model,
								PLOT_DEFAULT_NUM_POINTS,
								Ein,
								data);
}

// Get data for thickess and a function of energy out given energy in
bool get_Thickness_vs_Eout(	StopPow & model ,
							int num_points ,
							float Ein ,
							std::vector< std::vector<float> > & data )
{
	// energy limits
	// lower limit is 0 (i.e. ranged out)
	// upper limit is Ein (i.e. thickness = 0)
	return get_Thickness_vs_Eout(model,
								model.get_Emin(),
								Ein,
								num_points,
								Ein,
								data );
}

// Get data for thickess and a function of energy out given energy in
bool get_Thickness_vs_Eout(	StopPow & model ,
							float Emin ,
							float Emax ,
							int num_points ,
							float Ein ,
							std::vector< std::vector<float> > & data )
{
	// clear data
	data.clear();

	// set up 2 vectors for E and dE/dx
	std::vector<float> E_vector; std::vector<float> Thickness_vector;
	data.push_back( Thickness_vector );
	data.push_back( E_vector );

	// step size:
	float dE = (Emax-Emin) / ((float)num_points);

	// now generate the plot values:
	for(float E = Emin; E <= Emax; E+=dE )
	{
		try
		{
			// for first point, do direct calculation:
			if( data[0].size() == 0)
			{
				data[0].push_back(E);
				data[1].push_back( model.Thickness(Ein,E) );
			}
			else
			{
				float Tlast = data[1].back();
				// do an incremental calculation based on last point:
				data[0].push_back(E);
				data[1].push_back( Tlast - model.Thickness(E, fmax(Emin,E-dE)) );
			}
		}
		catch( std::invalid_argument e )
		{
			// if an error occurs, abort and return false:
			return false;
		}
	}

	return true;
}

// Get dataset for thickness as a function of energy in given energy out
bool get_Thickness_vs_Ein(	StopPow & model ,
							float Eout ,
							std::vector< std::vector<float> > & data )
{
	return get_Thickness_vs_Ein(model,
								PLOT_DEFAULT_NUM_POINTS,
								Eout,
								data );
}

// Get dataset for thickness as a function of energy in given energy out
bool get_Thickness_vs_Ein(	StopPow & model ,
							int num_points ,
							float Eout ,
							std::vector< std::vector<float> > & data )
{
	// set energy limits
	// minimum energy is Eout (thickness -> 0)
	// maximum energy is model's Emax
	return get_Thickness_vs_Ein(model,
								Eout,
								model.get_Emax(),
								num_points,
								Eout,
								data);
}

// Get dataset for thickness as a function of energy in given energy out
bool get_Thickness_vs_Ein(	StopPow & model ,
							float Emin ,
							float Emax ,
							int num_points ,
							float Eout ,
							std::vector< std::vector<float> > & data )
{
	// clear data
	data.clear();

	// set up 2 vectors for E and dE/dx
	std::vector<float> E_vector; std::vector<float> Thickness_vector;
	data.push_back( Thickness_vector );
	data.push_back( E_vector );

	// step size:
	float dE = (Emax-Emin) / ((float)num_points);

	// now generate the plot values:
	for(float E = Emin; E <= Emax; E+=dE )
	{
		try
		{
			// for first point, do direct calculation:
			if( data[0].size() == 0 )
			{
				data[0].push_back(E);
				data[1].push_back( model.Thickness(E,Eout) );
			}
			else
			{
				float Tlast = data[1].back();
				// otherwise do incremental based off of last point
				data[0].push_back(E);
				data[1].push_back( Tlast + model.Thickness(E, fmax(E-dE,Emin)) );
			}
		}
		catch( std::exception e )
		{
			// if an error occurs, abort and return false:
			return false;
		}
	}

	return true;
}

} // end of namespace
