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

#include "StopPow_Mehlhorn.h"

namespace StopPow
{

const double StopPow_Mehlhorn::Emin = 0.1; /* Minimum energy for dE/dx calculations */
const double StopPow_Mehlhorn::Emax = 30; /* Maximum energy for dE/dx calculations */


// Initialization routines specific to this model
void StopPow_Mehlhorn::init()
{
	if(ne > 0)
	{
		// set up Li-Petrasso for the free electrons and ions:
		std::vector<double> plasma_mf {me/mp};
		std::vector<double> plasma_Zf {-1.};
		std::vector<double> plasma_Tf {Te};
		std::vector<double> plasma_nf {ne};
		// add the ions:
		for(int i=0; i<num; i++)
		{
			if(Zbar[i] > 0)
			{
				plasma_mf.push_back( mf[i] );
				plasma_Zf.push_back( Zbar[i] );
				plasma_Tf.push_back( Tf[i] );
				plasma_nf.push_back( nf[i] );
			}
		}
		PlasmaStop = new StopPow_LP(mt, Zt, plasma_mf, plasma_Zf, plasma_Tf, plasma_nf);
	}
	else
	{
		PlasmaStop = NULL;
	}

	// set the info string:
	model_type = "Mehlhorn";
	info = "";

	use_manual_Ibar = false;
}
// constructors
StopPow_Mehlhorn::StopPow_Mehlhorn(double mt_in, double Zt_in, std::vector<double> & mf_in, std::vector<double> & Zf_in, std::vector<double> & Tf_in, std::vector<double> & nf_in, std::vector<double> & Zbar_in, double Te_in) throw(std::invalid_argument)
	: StopPow_PartialIoniz(mt_in, Zt_in, mf_in, Zf_in, Tf_in, nf_in, Zbar_in, Te_in)
{
	init();
}
StopPow_Mehlhorn::StopPow_Mehlhorn(double mt, double Zt, std::vector< std::array<double,5> > & field, double Te_in) throw(std::invalid_argument)
	: StopPow_PartialIoniz(mt, Zt, field, Te_in)
{
	init();
}

// Destructor
StopPow_Mehlhorn::~StopPow_Mehlhorn()
{
	delete PlasmaStop;
}

/** Calculate the total stopping power
 * @param E the test particle energy in MeV
 * @return stopping power in units of MeV/um
 * @throws invalid_argument
 */
double StopPow_Mehlhorn::dEdx_MeV_um(double E) throw(std::invalid_argument)
{
	// sanity check:
	if( E < Emin || E > Emax )
	{
		std::stringstream msg;
		msg << "Energy passed to StopPow_Mehlhorn::dEdx is bad: " << E;
		throw std::invalid_argument(msg.str());
	}

	double cold = 0; // return value for cold component
	// components:
	double Bethe, LSS, nuc;

	//iterate over all field particles:
	for(int i=0; i < num; i++)
	{
		// If fully ionized, only return plasma!
		if(Zbar[i] < Zf[i])
		{
			// call private helper functions:
			Bethe = dEdx_Bethe(E, i);
			LSS = dEdx_LSS(E, i);
			nuc = dEdx_nuc(E, i);
			cold += fmax(Bethe, LSS) + nuc;
		}
	}

	// calculate the free electron contribution:
	double hot = 0;
	if(PlasmaStop != NULL)
		hot = PlasmaStop->dEdx_MeV_um(E);

	return (cold+hot);
}

/** Calculate the total stopping power
 * @param E the test particle energy in MeV
 * @return stopping power in units of MeV/(mg/cm2)
 * @throws invalid_argument
 */
double StopPow_Mehlhorn::dEdx_MeV_mgcm2(double E) throw(std::invalid_argument)
{
	return (dEdx_MeV_um(E)*1e4) / (rho*1e3);
}

/**
 * Get the minimum energy that can be used for dE/dx calculations
 * @return Emin in MeV
 */
double StopPow_Mehlhorn::get_Emin()
{
	return Emin;
}

/**
 * Get the maximum energy that can be used for dE/dx calculations
 * @return Emax in MeV
 */
double StopPow_Mehlhorn::get_Emax()
{
	return Emax;
}

// Stopping for low energy ions, from LSS theory
double StopPow_Mehlhorn::dEdx_LSS(double E, int index)
{
	// See Eq 3 of T. Mehlhorn, C. Appl. Phys. 52, 6522 (1981)

	double A = mf[index] / mt;
	double K = 0.0793 * pow(ZtEff(E), 2/3) * sqrt(Zf[index]) * pow(1+A, 1.5)
				/ ( pow(pow(ZtEff(E), 2/3) + pow(Zf[index], 2/3), 0.75) * sqrt(mf[index]) );
	double a = 4.683e-9 / sqrt( pow(ZtEff(E), 2/3) + pow(Zf[index], 2/3) ); // cm
	double El = (1 + A)*Zf[index]*ZtEff(E)*e*e/(A*a); // erg
	double RL = pow(1+A, 2) / (4*M_PI*A*nf[index]*a*a);
	double C_LSS = K * sqrt(El/1.602e-9) / (RL*1e4); // KeV^1/2 / um

	// convert. Sign is for convention in this program that dE/dx < 0
	C_LSS = -1 * C_LSS * sqrt(1e3); // MeV^1.2 / um
	// and return the stopping power:
	return C_LSS * sqrt(E);
}

// Nuclear stopping power
double StopPow_Mehlhorn::dEdx_nuc(double E, int index)
{
	// See Eq 4 of T. Mehlhorn, C. Appl. Phys. 52, 6522 (1981)
	double C = E / mt; // MeV/amu
	double Cn = 4.14e6 * pow(mt/(mt+mf[index]), 1.5) * sqrt(ZtEff(E)*Zf[index]/mf[index])
				/ pow( pow(ZtEff(E), 2/3) + pow(Zf[index], 2/3), 0.75);
	double Cnprime = mf[index]*mt/(mf[index]+mt) * (1./(ZtEff(E)*Zf[index])) / sqrt( pow(ZtEff(E), 2/3) + pow(Zf[index], 2/3) );
	// given in paper in terms of areal density:
	double dEdr = Cn * sqrt(C) * exp(-45.2*pow(Cnprime*C, 0.277));
	// return with conversion to MeV / um:
	return (dEdr *1e4) / (rho*1e3);
}

// Bethe stopping power, with Mehlhorn's adjustments
double StopPow_Mehlhorn::dEdx_Bethe(double E, int index)
{
	double Ekev = E * 1e3; // energy in keV for convenience
	double ret = 0; // return value

	// see Eq 1 of T. Mehlhorn, C. Appl. Phys. 52, 6522 (1981)
	double rho_i, LogLamda, vt, beta, gamma, prefac;
	rho_i = nf[index] * mf[index] / Na; // mass density in g/cm3
	LogLamda = 0.0; // initialize

	vt = c*sqrt(2.0*Ekev/(mt*mpc2)); // test particle velocity
	beta = vt/c; // normalized to c
	gamma = 1.0/sqrt(1-pow(beta,2)); // relativistic gamma factor

	prefac = 4.0*M_PI*Na*rho_i*pow(ZtEff(E)*e*e,2)*(Zf[index]-Zbar[index]) / (me*c*c*beta*beta*mf[index]);
	LogLamda += log(2.0*me*pow(c*beta*gamma,2)/Ibar(E, index));
	LogLamda -= pow(beta,2);
	LogLamda -= shell_term(Zf[index],E);
	// no polarization effects included for now
	ret -= prefac*LogLamda*(1e-13)/(1.602e-19); // MeV/cm

	return ret*1e-4; // MeV/um
}

// effective ionization potential of partially ionized matter 
// for use in Bethe stopping power
double StopPow_Mehlhorn::Ibar(double E, int index)
{
	if( use_manual_Ibar )
	{
		return Ibar_manual[index]*1.602e-12;
	}
	// otherwise default to using Mehlhorn model
	// See Eq 7 of T. Mehlhorn, C. Appl. Phys. 52, 6522 (1981)

	// need to calculate Ibar(Z-Zbar), which is generally nonintegral.
	// calculate bounding indices:
	int i1 = floor( Zf[index] - Zbar[index] );
	int i2 =  ceil( Zf[index] - Zbar[index] );
	double di = Zf[index] - Zbar[index] - (double)i1;

	// sanity check:
	if( i1<0 || i2>=AtomicData::n )
	{
		std::stringstream msg;
		msg << "Out of range in StopPow_Mehlhorn::Ibar. Got i1, i2 = ";
		msg << i1 << "," << i2;
		throw std::invalid_argument(msg.str());
	}

	// check if we have an integer value:
	double Ibar;
	if(i1==i2)
	{
		Ibar = AtomicData::get_mean_ionization(i1);
	}
	else
	{
		// use linear interpolation to get the right value
		double Ibar1 = AtomicData::get_mean_ionization(i1);
		double Ibar2 = AtomicData::get_mean_ionization(i2);

		// there's no data in the atomic table for 0, so:
		if( i1 == 0 )
			Ibar1 = 0.;

		Ibar = Ibar1 + di*(Ibar2-Ibar1)/( (double)(i2-i1) );
	}
	double ret_eV = pow(Zf[index], 2) * Ibar / pow( Zf[index] - Zbar[index] , 2);
	return ret_eV * 1.602e-12;
}

// Set manual Ibars:
void StopPow_Mehlhorn::set_Ibar(std::vector<double> Ibar) throw(std::invalid_argument)
{
	if(Ibar.size() == Zf.size())
	{
		Ibar_manual = std::vector<double>(Ibar);
		use_manual_Ibar = true;
	}
	else
	{
		throw std::invalid_argument("StopPow_Mehlhorn::set_Ibar got wrong number of elements passed to it");
	}
}

// effective projectile charge
double StopPow_Mehlhorn::ZtEff(double E)
{
	// See Eq 6 of T. Mehlhorn, C. Appl. Phys. 52, 6522 (1981)
	// test particle velocity
	double beta = sqrt(2e3*E/(mt*mpc2)); // normalized to c

	return Zt * (1 - 1.034*exp(-137.04*beta/pow(Zt, 0.69)));
}

// Calculate shell correction term in log lambda 
double StopPow_Mehlhorn::shell_term(double Zf, double E)
{
	int Z = (int)Zf;

	// sanity check:
	if( Z < 1 || Z > AtomicData::n || E < Emin || E > Emax)
		return 0;

	// get coefficients:
	std::array<double,5> coeff = AtomicData::get_shell_coeff(Z);

	double shell = 0;
	// for each coefficient:
	for(int i=0; i < coeff.size(); i++)
		shell += coeff[i] * pow( log(1e3*E/mt) , i); // convert E to keV

	return shell;
}

} // end namespace StopPow
