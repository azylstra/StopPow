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

#include "StopPow_Fit.h"

namespace StopPow
{

const int StopPow_Fit::MODE_ZIMMERMAN = 0; 
const int StopPow_Fit::MODE_LP = 1; 
const int StopPow_Fit::MODE_BPS = 2;
const int StopPow_Fit::MODE_GRABOWSKI = 3;
const int StopPow_Fit::MODE_QUANTUM_GRABOWSKI = 4;
const int StopPow_Fit::MODE_LP_PUB = 5;

// Constructors use those in PartialIoniz
StopPow_Fit::StopPow_Fit(double mt_in, double Zt_in, std::vector<double> & mf_in, std::vector<double> & Zf_in, std::vector<double> & Tf_in, std::vector<double> & nf_in, std::vector<double> & Zbar_in, double Te_in) throw(std::invalid_argument)
	: StopPow_PartialIoniz(mt_in, Zt_in, mf_in, Zf_in, Tf_in, nf_in, Zbar_in, Te_in) 
{
	init();
}

StopPow_Fit::StopPow_Fit(double mt_in, double Zt_in, std::vector< std::array<double,5> > & field_in, double Te_in) throw(std::invalid_argument)
	: StopPow_PartialIoniz(mt_in, Zt_in, field_in, Te_in) 
{
	init();
}

// Destructor
StopPow_Fit::~StopPow_Fit()
{
	if( fe != z && fe != NULL )
		delete fe;
	if( fe2 != NULL )
		delete fe2;
	delete z;
}

void StopPow_Fit::init()
{
	// Set up the Zimmerman model:
	z = new StopPow_Zimmerman(mt, Zt, mf, Zf, Tf, nf, Zbar, Te);
	fe = z;
}

// Calculate stopping power
double StopPow_Fit::dEdx_MeV_um(double E) throw(std::invalid_argument)
{
	// three components:
	double dEdx_i, dEdx_be, dEdx_fe;

	dEdx_i = z->dEdx_ion(E);
	dEdx_be = be_factor * z->dEdx_bound_electron(E);
	// Zimmerman model uses a single StopPow object, necessitates calling
	// the free electron function
	if( fe == z )
	{
		dEdx_fe = fe_factor * ((StopPow_Zimmerman*)fe)->dEdx_free_electron(E);
	}
	else // LP or BPS or Grabowski
	{
		dEdx_fe = fe_factor * ((StopPow_Plasma*)fe)->dEdx_plasma_electrons(E);
	}

	// Quantum Grabowski needs quantum correction:
	if( fe_model == MODE_QUANTUM_GRABOWSKI )
	{
		dEdx_fe += fe_factor * ((StopPow_BPS*)fe2)->dEdx_quantum(E,0);
	}

	return (dEdx_i + dEdx_be + dEdx_fe);
}

// Calculate stoppign power
double StopPow_Fit::dEdx_MeV_mgcm2(double E) throw(std::invalid_argument)
{
	return (dEdx_MeV_um(E)*1e4) / (rho*1e3);
}

// Get energy limits, set by most restrictive of two sub-models:
double StopPow_Fit::get_Emin()
{
	return fmax(z->get_Emin(), fe->get_Emin());
}
double StopPow_Fit::get_Emax()
{
	return fmin(z->get_Emax(), fe->get_Emax());
}

// Normalize the bound-electron stopping to a reference case at a given proton energy.
void StopPow_Fit::normalize_bound_e(StopPow * ref, double Ep)
{
	// Need to create a dummy zero-ionization version of the Zimmerman model being used:
	std::vector<double> Zbar_0(Zbar);
	for(int i=0; i<Zbar_0.size(); i++)
		Zbar_0[i] = 0.0; // no ionization
	StopPow_Zimmerman * z2 = new StopPow_Zimmerman(mt, Zt, mf, Zf, Tf, nf, Zbar_0, Te);
	be_factor = ref->dEdx(Ep) / z2->dEdx(Ep);
	delete z2;
}

// Choose the free-electron stopping model.
void StopPow_Fit::choose_model(int new_model) throw(std::invalid_argument)
{
	// Plasma conditions for the free electrons:
	std::vector<double> mf_fe {me/amu};
	std::vector<double> Zf_fe {-1};
	std::vector<double> Tf_fe {Te};
	std::vector<double> nf_fe {ne};

	// Add in ions:
	for(int i=0; i<mf.size(); i++)
	{
		mf_fe.push_back(mf[i]);
		Zf_fe.push_back(Zbar[i]);
		Tf_fe.push_back(Tf[i]);
		nf_fe.push_back(nf[i]);
	}

	// Configure the new model based on the given flag:
	switch(new_model)
	{
		case MODE_ZIMMERMAN:
			if( fe != z )
				delete fe;
			fe = z;
			break;
		case MODE_LP:
			if( fe != z )
				delete fe;
			fe = new StopPow_LP(mt, Zt, mf_fe, Zf_fe, Tf_fe, nf_fe);
			break;
		case MODE_BPS:
			if( fe != z )
				delete fe;
			fe = new StopPow_BPS(mt, Zt, mf_fe, Zf_fe, Tf_fe, nf_fe);
			break;
		case MODE_GRABOWSKI:
			if( fe != z )
				delete fe;
			fe = new StopPow_Grabowski(mt, Zt, mf_fe, Zf_fe, Tf_fe, nf_fe);
			break;
		case MODE_QUANTUM_GRABOWSKI:
			if( fe != z )
				delete fe;
			if( fe2 != NULL )
				delete fe2;
			fe = new StopPow_Grabowski(mt, Zt, mf_fe, Zf_fe, Tf_fe, nf_fe);
			fe2 = new StopPow_BPS(mt, Zt, mf_fe, Zf_fe, Tf_fe, nf_fe);
			break;
		case MODE_LP_PUB:
			if( fe != z )
				delete fe;
			fe = new StopPow_LP(mt, Zt, mf_fe, Zf_fe, Tf_fe, nf_fe);
			((StopPow_LP*)fe)->set_xtf_factor(2.);
			((StopPow_LP*)fe)->set_u_factor(2.);
			((StopPow_LP*)fe)->use_published_collective(true);
			break;
		default:
			throw std::invalid_argument("Model choice passed to StopPow_Fit::choose_model is invalid");
	}
	fe_model = new_model;
}

// Adjust the stopping, to be used for fitting
void StopPow_Fit::set_factor(double factor)
{
	fe_factor = factor;
}

// Get the current adjustment factor
double StopPow_Fit::get_factor()
{
	return fe_factor;
}

} // end of namespace