
#include "Spectrum.h"

namespace StopPow
{

void shift(StopPow & model, float thickness, std::vector<float> & data_E, std::vector<float> & data_Y) throw(std::invalid_argument)
{
	// Use the other function with dummy error bars:
	std::vector<float> data_err;
	data_err.reserve(data_E.size());
	for(int i=0; i<data_E.size(); i++)
		data_err[i] = 0.;
	shift(model, thickness, data_E, data_Y, data_err);
}

void shift(StopPow & model, float thickness, std::vector<float> & data_E, std::vector<float> & data_Y, std::vector<float> & data_err) throw(std::invalid_argument)
{
	/*
	* Various sanity checks
	*/
	// Must be same dimensions:
	if( data_Y.size() != data_E.size() || data_err.size() != data_E.size())
	{
		throw std::invalid_argument("StopPow::shift - data vectors of different sizes");
	}
	// Spectrum must be regularly spaced:
	float dE = data_E[1] - data_E[0];
	for(int i=0; i<data_E.size()-1; i++)
	{
		if( !approx(data_E[i]-data_E[i+0], dE, 1e-4f) )
		{
			throw std::invalid_argument("StopPow::shift - Energy bins invalid.");
		}
	}

	/*
	* Perform actual calculation
	*/
	// Temporary values:
	std::vector< std::array<float,3> > ret;
	ret.reserve(data_E.size());
	for(int i=0; i<data_E.size(); i++)
	{
		ret[i][0] = data_E[i];
		ret[i][1] = 0.;
		ret[i][2] = 0.;
	}

	// Now do the calculation for each bin in the original dataset
	for(int i=0; i<data_E.size(); i++)
	{
		// energy limits of bin:
		float Emin, E, Emax;
		Emin = data_E[i]-dE/2;
		E = data_E[i];
		Emax = data_E[i]+dE/2;

		// New yield/MeV and error:
		float Y = data_Y[i];
		float err = data_err[i];

		// split initial bin into 50 pieces and shift, rebin each:
		float n = 50;
		float dE2 = dE/n;
		float Eshift; int index;
		for(float E2=Emin+dE2/2; E2<Emax; E2+=dE2)
		{
			if(thickness<0)
				Eshift = model.Ein(E2, -1.*thickness);
			else if(thickness>0)
				Eshift = model.Eout(E2, thickness);
			index = (int)((Eshift-data_E[0])/dE);
			if(index>0 && index<data_E.size())
			{
				ret[index][1] += Y/n;
				ret[index][2] += err/n;
			}
		}
	}

	// Copy values back from ret:
	for(int i=0; i<data_E.size(); i++)
	{
		data_E[i] = ret[i][0];
		data_Y[i] = ret[i][1];
		data_err[i] = ret[i][2];
	}
}

} // end of namespace StopPow