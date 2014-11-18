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

#include "StopPow_Zimmerman.h"

namespace StopPow
{

const double StopPow_Zimmerman::Emin = 0.01; /* Minimum energy for dE/dx calculations */
const double StopPow_Zimmerman::Emax = 30; /* Maximum energy for dE/dx calculations */


// Initialization routines specific to this model
void StopPow_Zimmerman::init()
{
	// set the info string:
	model_type = "Zimmerman";
	info = "";
}
// constructors for partial ionized material:
StopPow_Zimmerman::StopPow_Zimmerman(double mt_in, double Zt_in, std::vector<double> & mf_in, std::vector<double> & Zf_in, std::vector<double> & Tf_in, std::vector<double> & nf_in, std::vector<double> & Zbar_in, double Te_in) throw(std::invalid_argument)
	: StopPow_PartialIoniz(mt_in, Zt_in, mf_in, Zf_in, Tf_in, nf_in, Zbar_in, Te_in)
{
	init();
}
StopPow_Zimmerman::StopPow_Zimmerman(double mt, double Zt, std::vector< std::array<double,5> > & field, double Te_in) throw(std::invalid_argument)
	: StopPow_PartialIoniz(mt, Zt, field, Te_in)
{
	init();
}

// Destructor
StopPow_Zimmerman::~StopPow_Zimmerman()
{
	// nothing to do
}

// Calculate stopping power
double StopPow_Zimmerman::dEdx_MeV_um(double E) throw(std::invalid_argument)
{
	// sanity check:
	if( E < Emin || E > Emax )
	{
		std::stringstream msg;
		msg << "Energy passed to StopPow_Zimmerman::dEdx is bad: " << E;
		throw std::invalid_argument(msg.str());
	}

	return dEdx_free_electron(E) + dEdx_ion(E) + dEdx_bound_electron(E);
}

// total stopping in rhoR units
double StopPow_Zimmerman::dEdx_MeV_mgcm2(double E) throw(std::invalid_argument)
{
	return (dEdx_MeV_um(E)*1e4) / (rho*1e3);
}

// For solving for mu:
// See Atzeni p 329-330
// Code uses GSL Newton's method; therefore requires several functions to work
struct mu_params
{
double kT, lth, ne;
};

double mu_f (double x, void * params) { 
	struct mu_params *p = (struct mu_params *) params;
	gsl_sf_result result;
	int code = gsl_sf_fermi_dirac_half_e(x/p->kT, &result);
	// gamma functions fix normalization diff between Atzeni and GSL
	if( code == 0 )
		return result.val * gsl_sf_gamma(1.5) / gsl_sf_gamma(0.5) - pow(p->lth,3.)*p->ne/2.;
	return -pow(p->lth,3.)*p->ne/2.;
};
double mu_df (double x, void * params) { 
	//struct mu_params *p = (struct mu_params *) params;
	double result, abserr;
	gsl_function F;
	F.function = &mu_f;
	F.params = params;
	gsl_deriv_central (&F, x, 1e-12, &result, &abserr);
	return result;
};
void mu_fdf(double x, void * params, 
               double *y, double *dy) {
	*y = mu_f(x,params);
	*dy = mu_df(x,params);
}

// Free electron stopping power
double StopPow_Zimmerman::dEdx_free_electron(double E)
{
	// Sanity check, if there are no electrons, dE/dx=0
	if(ne == 0)
		return 0.;
	
	// test particle velocity
	double vt = c*sqrt(2e3*E/(mt*mpc2));
	// y parameter just ratio of test / thermal velocity
	// for thermal velocity, depends on if we are using quantum effects:
	double vth = sqrt(2.*kB*Te*keVtoK/me); // standard nondegenerate (Eq 19)
	if(quantum)
	{
		// solve for mu:
		// set up finder using Newton's method:
		const gsl_root_fdfsolver_type *Tsolver;
		gsl_root_fdfsolver *s;
		gsl_function_fdf Fmu;
		Fmu.f = &mu_f;
		Fmu.df = &mu_df;
		Fmu.fdf = &mu_fdf;
		double lth = sqrt(2*M_PI*hbar*hbar/(me*kB*Te*keVtoK));
		struct mu_params params = {kB*Te*keVtoK, lth, ne};
		Fmu.params = &params;
		Tsolver = gsl_root_fdfsolver_newton;
		s = gsl_root_fdfsolver_alloc (Tsolver);
		gsl_root_fdfsolver_set (s, &Fmu, 0);
		double mu = 0.; int status;
		do
		{
		  status = gsl_root_fdfsolver_iterate (s);
		  mu = (gsl_root_fdfsolver_root (s));
		  status = (mu_f(mu, &params) <= 1e-12);
		}
		while (status == GSL_CONTINUE);

		// Zimmerman Eq 18 gives a quantum expression for vth, but it is only really applicable
		// if greater than the usual thermal velocity, thus taking the max below:
		vth = fmax(vth, (h/(2.*sqrt(M_PI)*me)) * pow( 4*ne*(1 + exp(-mu/(kB*Te*keVtoK))) , 1./3 ));

		gsl_root_fdfsolver_free(s);
	}
	double y = vt/vth;
	double omega_pe = sqrt(4*M_PI*esu*esu*ne/me);
	// Eq 16:
	double LambdaF = (4*M_PI*me*pow(vth,2)/(h*omega_pe)) * (0.321+0.259*pow(y,2)+0.0707*pow(y,4)+0.05*pow(y,6))/(1.+0.130*pow(y,2)+0.05*pow(y,4));
	double dEdx_F = 4*M_PI*(1./me)*pow((double)esu,4)*pow(Zt/vt,2)*ne*LF(y,LambdaF);
	return -1.*dEdx_F * 624150.934 * 1e-4; // MeV/um
}

// Bound electron stopping power
double StopPow_Zimmerman::dEdx_bound_electron(double E)
{
	// test particle velocity
	double vt = c*sqrt(2e3*E/(mt*mpc2));

	double dEdx_BE = 0.;
	double Ibar, ZiB, LiB;
	double prefac = pow((double)e,4)*(4.*M_PI*pow(Zt,2) / (me*pow(vt,2)));
	// have to loop over all ions
	for(int i=0; i<num; i++)
	{
		ZiB = (Zf[i] - Zbar[i]); // number of bound electrons
		if(ZiB > 0.)
		{
			// Eq 20:
			Ibar = Zf[i] * (0.024 - 0.013 * ZiB/Zf[i]) / sqrt(ZiB/Zf[i]); // in keV
			Ibar = Ibar * 1e3 * 1.60217e-12; // in erg
			LiB = log(2. * me * pow(vt,2) / Ibar);
			dEdx_BE +=  prefac * nf[i] * ZiB * LiB;
		}
	}
	return -1. * dEdx_BE * 624150.934 * 1e-4; // MeV/um
}

// Ion stopping power
double StopPow_Zimmerman::dEdx_ion(double E)
{
	// test particle velocity
	double vt = c*sqrt(2e3*E/(mt*mpc2));

	// need to loop over field ions
	double dEdx_I = 0.;
	double mr, bi, Li;
	double prefac = (4*M_PI*pow(esu,4)*pow(Zt/vt,2)) / amu;
	for(int i=0; i<num; i++)
	{
		mr = amu*mf[i]*mt/(mf[i]+mt);
		// Eq 14, effective minimum impact param for ion stopping
		bi = sqrt( pow(h/(4*M_PI*mr*vt),2) + pow(esu*esu*Zf[i]*Zt/(mr*vt*vt),2) ); 
		Li = log(lDebye()/bi); // ion stopping number
		// Eq 12:
		dEdx_I += prefac * (nf[i]*Zf[i]*Zf[i]*Li/(mf[i]));
	}
	return -1.*dEdx_I * 624150.934 * 1e-4; // MeV/um
}

// whether to use quantum correction
void StopPow_Zimmerman::set_quantum(bool set)
{
	quantum = set;
}

// Minimum energy limit
double StopPow_Zimmerman::get_Emin()
{
	return Emin;
}

// Maximum energy limit
double StopPow_Zimmerman::get_Emax()
{
	return Emax;
}

// Electron stopping number
double StopPow_Zimmerman::LF(double y, double LambdaF)
{
	return 0.5 * log(1.+pow(LambdaF,2)) * (gsl_sf_erf(y) - (2./sqrt(M_PI))*y*exp(-y*y));
}

// Debye length in field plasma
double StopPow_Zimmerman::lDebye()
{
	double ret = 0; // temporary return value
	//iterate over all field ions:
	for(int i=0; i < num; i++)
	{
		ret += (4.*M_PI*nf[i]*pow(Zf[i]*e,2.0) / (kB*Tf[i]*keVtoK) );
	}
	// electrons:
	ret += (4.*M_PI*ne*pow(e,2.0) / (kB*Te*keVtoK) );
	return 1.0/sqrt(ret);
}

} // end namespace StopPow
