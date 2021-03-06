// StopPow - a charged-particle stopping power library
// Copyright (C) 2014  Massachusetts Institute of Technology / Alex Zylstra

// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along
// with this program; if not, write to the Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

#include "Spectrum.h"

namespace StopPow
{

void shift(StopPow & model, double thickness, std::vector<double> & data_E, std::vector<double> & data_Y) throw(std::invalid_argument)
{
	// Use the other function with dummy error bars:
	std::vector<double> data_err;
	for(int i=0; i<data_E.size(); i++)
		data_err.push_back(0);
	shift(model, thickness, data_E, data_Y, data_err);
}

void shift(StopPow & model, double thickness, std::vector<double> & data_E, std::vector<double> & data_Y, std::vector<double> & data_err) throw(std::invalid_argument)
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
	double dE = data_E[1] - data_E[0];
	for(int i=0; i<data_E.size()-1; i++)
	{
		if( !approx(data_E[i+1]-data_E[i], dE, 1e-4) )
		{
			throw std::invalid_argument("StopPow::shift - Energy bins invalid.");
		}
	}

	/*
	* Perform actual calculation
	*/
	// Temporary values:
	std::vector< std::array<double,3> > ret;
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
		double Emin, E, Emax;
		Emin = data_E[i]-dE/2;
		E = data_E[i];
		Emax = data_E[i]+dE/2;

		// New yield/MeV and error:
		double Y = data_Y[i];
		double err = data_err[i];

		// split initial bin into 50 pieces and shift, rebin each:
		double n = 50;
		double dE2 = dE/n;
		double Eshift; int index;
		for(double E2=Emin+dE2/2.; E2<Emax; E2+=dE2)
		{
			if(thickness<0)
				Eshift = model.Ein(E2, -1.*thickness);
			else if(thickness>0)
				Eshift = model.Eout(E2, thickness);
			index = floor( (Eshift-(data_E[0]-dE/2.)) / dE );
			if(index>=0 && index<data_E.size())
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