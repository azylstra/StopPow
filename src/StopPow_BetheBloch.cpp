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

#include "StopPow_BetheBloch.h"

namespace StopPow
{

 /** Initialize the Bethe-Bloch calculator.
 * @param mt the test particle mass in AMU
 * @param Zt the test particle in charge (units of e)
 * @param mf vector containing ordered field particle masses in AMU
 * @param Zt vector containing ordered field particle charges in units of e
 * @param nf vector containing ordered field particle densities in units of 1/cm3
 * @throws invalid_argument
*/
StopPow_BetheBloch::StopPow_BetheBloch(double mt_in, double Zt_in, std::vector<double> mf_in, std::vector<double> Zf_in, std::vector<double> nf_in) throw(std::invalid_argument)
{
	// default mode for B-B:
	set_mode(MODE_LENGTH);

	// infer size of the field particle arrays:
	num = mf_in.size();

	// sanity checking. 
	bool args_ok = true;
	// Make sure mt and Zt are positive,
	// and that all field particle arrays have same size
	if( mt_in <= 0 || Zt_in <= 0
		|| Zf_in.size() != num
		|| nf_in.size() != num )
	{
		args_ok = false;
	}

	// now do sanity checking on the field particle values:
	for(int i=0; i<num; i++)
	{
		args_ok = args_ok && mf_in[i] > 0;
		args_ok = args_ok && nf_in[i] > 0;
	}

	// throw an exception if necessary:
	if( !args_ok )
	{
		std::stringstream msg;
		// start constructing message, add info on mt and Zt:
		msg << "Values passed to StopPow_LP constructor are bad: " 
		 << mt_in << "," << Zt_in << "," << std::endl;

		std::vector<double>::iterator it; // to iterate over field particles

		// add each element in mf:
		msg << "mf = ";
		for(it=mf.begin(); it<mf.end(); it++)
		 	msg << (*it) << ",";

		// add each element in Zf:
		msg << std::endl << "Zf = ";
		for(it=Zf.begin(); it<Zf.end(); it++)
		 	msg << (*it) << ",";

		// add each element in nf:
		msg << std::endl << "nf = ";
		for(it=nf.begin(); it<nf.end(); it++)
		 	msg << (*it) << ",";

		 // throw the exception:
		throw std::invalid_argument(msg.str());
	}

	// set class variables:
	mt = mt_in;
	Zt = Zt_in;
	mf = mf_in;
	Zf = Zf_in;
	nf = nf_in;

	// calculate the field particle mass density:
	rho = 0; // g/cm3
	// iterate over field particles:
	for(int i=0; i<num; i++)
	{
		rho += mf[i] * mp * nf[i];
	}

	// default state is to have shell corrections on:
	use_shell_corr = true;

	// Minimum energy for dE/dx calculations
	Emin = 0.6 * mt; // see Andersen and Ziegler for justification of this limit  
	// Maximum energy for dE/dx calculations
	Emax = 30; 

	// set the info strings:
	model_type = "Bethe-Bloch";
	info = "";

	use_manual_Ibar = false;
}

// Destructor
StopPow_BetheBloch::~StopPow_BetheBloch()
{
	// nothing to do
}

/** Calculate the total stopping power
 * @param E the test particle energy in MeV
 * @return stopping power in units of MeV/um
 * @throws invalid_argument
*/
double StopPow_BetheBloch::dEdx_MeV_um(double E) throw(std::invalid_argument)
{
	// sanity check:
	if( E < Emin || E > Emax || std::isnan(E) )
	{
		std::stringstream msg;
		msg << "Energy passed to StopPow_BetheBloch::dEdx is bad: " << E;
		throw std::invalid_argument(msg.str());
	}

	double Ekev = E * 1e3; // energy in keV for convenience

	double ret = 0;

	// iterate over field particles:
	double rho, LogLamda, vt, beta, gamma, prefac;
	for(int i=0; i < num; i++)
	{
		rho = nf[i] * mf[i] / Na; // mass density in g/cm3
		LogLamda = 0.0; // initialize

		vt = c*sqrt(2.0*Ekev/(mt*mpc2)); // test particle velocity
		beta = vt/c; // normalized to c
		gamma = 1.0/sqrt(1-pow(beta,2)); // relativistic gamma factor

		prefac = 4.0*M_PI*Na*rho*pow(Zt*e*e,2)*Zf[i] / (me*c*c*beta*beta*mf[i]);
		LogLamda += log(2.0*me*pow(c*beta*gamma,2)/Ibar(Zf[i]));
		LogLamda -= pow(beta,2);
		if(use_shell_corr)
			LogLamda -= shell_term(Zf[i],E);
		// TODO: polarization effects and shell corrections for now
		ret -= prefac*LogLamda*(1e-13)/(1.602e-19); // MeV/cm
	}

	return ret*1e-4;
}

/** Calculate the total stopping power
 * @param E the test particle energy in MeV
 * @return stopping power in units of MeV/(mg/cm2)
 * @throws invalid_argument
 */
double StopPow_BetheBloch::dEdx_MeV_mgcm2(double E) throw(std::invalid_argument)
{
	return (dEdx_MeV_um(E)*1e4) / (rho*1e3);
}

/**
 * Get the minimum energy that can be used for dE/dx calculations
 * @return Emin in MeV
 */
double StopPow_BetheBloch::get_Emin()
{
	return Emin;
}

/**
 * Get the maximum energy that can be used for dE/dx calculations
 * @return Emax in MeV
 */
double StopPow_BetheBloch::get_Emax()
{
	return Emax;
}

/** Effecive ionization potential as a function of Z.
 * @param Z field particle charge in units of e
 * @return Ibar in erg
 */
double StopPow_BetheBloch::Ibar(double Z)
{
	// use manual value if it has been set:
	if( use_manual_Ibar )
	{
		int i = 0;
		// find the right one for Zf:
		for(int j=0; j<Zf.size(); j++)
		{
			if(Zf[i] == Z)
				return Ibar_manual[i]*1.602e-12;  // return in erg
		}
	}
	// default is to use built-in atomic data:
	int iZ = (int) Z;
	if( iZ > 0 && iZ < AtomicData::n  )
		return AtomicData::get_mean_ionization(iZ)*1.602e-12;
	return 0;
}

// Set the ionization potential manually
void StopPow_BetheBloch::set_Ibar(std::vector<double> Ibar) throw(std::invalid_argument)
{
	if(Ibar.size() == Zf.size())
	{
		Ibar_manual = std::vector<double>(Ibar);
		use_manual_Ibar = true;
	}
	else
	{
		throw std::invalid_argument("StopPow_BetheBloch::set_Ibar got wrong number of elements passed to it");
	}
}

// Calculate shell correction term in log lambda 
double StopPow_BetheBloch::shell_term(double Zf, double E)
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

// Turn shell corrections on or off
void StopPow_BetheBloch::use_shell_correction(bool enabled)
{
	use_shell_corr = enabled;
}

// get current state of shell corrections
bool StopPow_BetheBloch::using_shell_correction()
{
	return use_shell_corr;
}
} // end namespace StopPow