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

#include "StopPow_BPS.h"

namespace StopPow
{

const double StopPow_BPS::Emin = 0.01; /* Minimum energy/A for dE/dx calculations */
const double StopPow_BPS::Emax = 50; /* Maximum energy/A for dE/dx calculations */

// L-P specific initialization stuff:
void StopPow_BPS::init()
{
	// set the info string:
	model_type = "BPS";
	info = "";
	// call helper method which precomputs some stuff:
	on_field_change();
}

// Li-Petrasso constructor primarily relies on superclass constructor
StopPow_BPS::StopPow_BPS(double mt_in, double Zt_in, std::vector<double> & mf_in, std::vector<double> & Zf_in, std::vector<double> & Tf_in, std::vector<double> & nf_in) throw(std::invalid_argument)
	: StopPow_Plasma(mt_in, Zt_in, mf_in, Zf_in, Tf_in, nf_in)
{
	init();
}
StopPow_BPS::StopPow_BPS(double mt, double Zt, std::vector< std::array<double,4> > & field) throw(std::invalid_argument)
	: StopPow_Plasma(mt, Zt, field)
{
	init();
}
StopPow_BPS::StopPow_BPS(double mt_in, double Zt_in, std::vector<double> & mf_in, std::vector<double> & Zf_in, std::vector<double> & Tf_in, std::vector<double> & nf_in, double Te) throw(std::invalid_argument)
	: StopPow_Plasma(mt_in, Zt_in, mf_in, Zf_in, Tf_in, nf_in, Te)
{
	init();
}
StopPow_BPS::StopPow_BPS(double mt, double Zt, std::vector< std::array<double,4> > & field, double Te) throw(std::invalid_argument)
	: StopPow_Plasma(mt, Zt, field, Te)
{
	init();
}

// Destructor
StopPow_BPS::~StopPow_BPS()
{
	// nothing to do
}

// Calculate the total stopping power
double StopPow_BPS::dEdx_MeV_um(double E) throw(std::invalid_argument)
{
	// sanity check:
	if( E < Emin || E > Emax )
	{
		std::stringstream msg;
		msg << "Energy passed to StopPow_BPS::dEdx is bad: " << E;
		throw std::invalid_argument(msg.str());
	}

	double ret_short, ret_long, ret_quantum; // return values for the terms

	// Use threading to speed up:
	auto shortF = [this,&E,&ret_short] () {ret_short = this->dEdx_short(E);};
	std::thread t1(shortF);
	auto longF = [this,&E,&ret_long] () {ret_long = this->dEdx_long(E);};
	std::thread t2(longF);
	auto quantumF = [this,&E,&ret_quantum] () {ret_quantum = this->dEdx_quantum(E);};
	std::thread t3(quantumF);

	t1.join();
	t2.join();
	t3.join();

	return ret_short + ret_long + ret_quantum;
}

// Calculate the total stopping power
double StopPow_BPS::dEdx_MeV_mgcm2(double E) throw(std::invalid_argument)
{
	return (dEdx_MeV_um(E)*1e4) / (rho*1e3);
}

// Get stopping power due only to a specific field particle species
double StopPow_BPS::dEdx_field(double E, int i) throw(std::invalid_argument)
{
	double ret_short, ret_long, ret_quantum; // return values for the terms

	// Use threading to speed up:
	auto shortF = [this,&E,&ret_short,&i] () {ret_short = this->dEdx_short(E,i);};
	std::thread t1(shortF);
	auto longF = [this,&E,&ret_long,&i] () {ret_long = this->dEdx_long(E,i);};
	std::thread t2(longF);
	auto quantumF = [this,&E,&ret_quantum,&i] () {ret_quantum = this->dEdx_quantum(E,i);};
	std::thread t3(quantumF);

	t1.join();
	t2.join();
	t3.join();

	return ret_short + ret_long + ret_quantum;
}

// For use in GSL integrations
struct int_params
{
	double beta_b, vp, ep, eb, kappa_b, K, mp, mb, mpb, Mpb;
};

// Function to integrate for short-range stopping power, i.e. integrand of Eq 3.3
double dEdxcs_func (double u, void * params) {
	struct int_params *p = (struct int_params *) params; // for GSL-style params, have to cast from void*
	double term1 = sqrt(u) * exp(-0.5 * p->beta_b * p->mb * pow(p->vp,2) * u);
	double term2 = ( (-log(p->beta_b * (fabs(p->ep * p->eb) * p->K/(4. * M_PI)) * (p->mb/p->mpb) * u/(1. - u)) + 2. - 2. * 0.5772)
	 	 * (p->beta_b * p->Mpb * pow(p->vp,2) - 1./u) + 2./u );
	return term1 * term2;
}

// Classical short-range stopping power (Eq. 3.3) for one species
double StopPow_BPS::dEdx_short(double E, int i)
{
	double vp = c*sqrt(2e3*E/(mt*mpc2)); // test particle velocity
	double ret = 0.;
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (100);
	double prefac, result, err;

	// Set up and integrate the integrand using GSL numerical integration library
	struct int_params params = {beta_b[i], vp, Zt*e_LH, Zf[i]*e_LH, kappa_b[i], K, mt*amu, mf[i]*amu, mpb[i], Mpb[i]};
	prefac = (pow(Zt*e_LH,2)/(4.*M_PI)) * (pow(kappa_b[i],2)/(mt*amu*vp)) * sqrt(mf[i]*amu/(2.*M_PI*beta_b[i]));
	gsl_function Fa;
	Fa.function = dEdxcs_func;
	Fa.params = &params;
	gsl_integration_qag(&Fa, 0, 1, 1e-8, 1e-4, 100, 1, w, &result, &err);

	ret += prefac * result;

	gsl_integration_workspace_free (w);
	// Convert from erg/cm to MeV/um and also flip sign for consistency
	return -1 * ret * 624150.934 * 1e-4; // MeV/um
}

// Classical short-range stopping power (Eq. 3.3)
double StopPow_BPS::dEdx_short(double E)
{
	double ret = 0;
	// Loop over field particles:
	for(int i=0; i<num; i++)
	{
		ret += dEdx_short(E,i);
	}
	return ret;
}

// Helper function for calculating long-range stopping power, i.e. integrand of Eq 3.4
gsl_complex StopPow_BPS::dEdx_long_func(double vp, double x, int i)
{
	// variable substitution: x = cos(theta)

	// Evaluate F:
	gsl_complex F1 = Fc(vp*x); // F(vp*costheta)

	// prefactor:
	double part1 = (pow(Zt*e_LH,2)/(8*M_PI*M_PI)) * x * (rho_b(vp*x, i)/rho_tot(vp*x)); // this part is nice and real

	// need to use complex numbers for the result:
	gsl_complex ret = gsl_complex_rect(0, part1); // includes factor of i in equation!
	ret = gsl_complex_mul(ret, F1);
	gsl_complex logval = gsl_complex_div(F1, gsl_complex_rect(pow(K,2), 0.));
	ret = gsl_complex_mul(ret, gsl_complex_log(logval));

	return ret;
}

// Classical long-range stopping power (Eq 3.4) for a single species
double StopPow_BPS::dEdx_long(double E, int i)
{
	double vp = c*sqrt(2e3*E/(mt*mpc2)); // test particle velocity
	double du = 0.025; // step size in numerical integration

	double ret = 0.; // return value

	// original version
	// Evaluation of the first term of the equation (part with integral):
	gsl_complex dEdx_cR_1 = gsl_complex_rect(0,0);
	// do integration manually
	for(double u=-1+du/2; u<1.; u+=du) {
		gsl_complex temp = gsl_complex_mul(dEdx_long_func(vp, u, i) , gsl_complex_rect(du,0));
		if(!isnan(GSL_REAL(temp)) && !isinf(GSL_REAL(temp))) {
			dEdx_cR_1 = gsl_complex_add(temp, dEdx_cR_1); // add value
		}
	}
	dEdx_cR_1 = gsl_complex_mul(dEdx_cR_1, gsl_complex_rect(624150.934 * 1e-4,0)); // Convert to MeV/um

	// Evaluation of the second term of the equation
	gsl_complex prefac2 = gsl_complex_rect(0, (pow(Zt*e_LH,2)/(8*M_PI*M_PI)) * (1./(beta_b[i]*mt*amu*pow(vp,2))) 
		* ( rho_b(vp, i) / rho_tot(vp) ));
	// Calculate F and F*
	gsl_complex Fv = Fc(vp);
	gsl_complex Fvc = gsl_complex_conjugate(Fv); // == Fc(-1*vp)
	// calculate two complex terms inside the square brackets of Eq 3.4
	gsl_complex term1 = gsl_complex_mul( Fv, gsl_complex_log( gsl_complex_div(Fv,gsl_complex_rect(pow(K,2),0)) ) );
	gsl_complex term2 = gsl_complex_mul( Fvc, gsl_complex_log( gsl_complex_div(Fvc, gsl_complex_rect(pow(K,2),0)) ) );
	gsl_complex term = gsl_complex_sub(term1,term2);
	term = gsl_complex_mul(term, prefac2);
	// final result:
	gsl_complex dEdx_cR_2 = gsl_complex_mul(term, gsl_complex_rect(624150.934*1e-4,0)); // Convert to MeV/um

	// combine two terms
	if(!isnan( GSL_REAL(dEdx_cR_1) ))
		ret += GSL_REAL(dEdx_cR_1);
	if(!isnan( GSL_REAL(dEdx_cR_2) ))
		ret -= GSL_REAL(dEdx_cR_2);

	// sign flip for consistency
	return -1 * ret;
}

// Classical long-range stopping power (Eq 3.4)
double StopPow_BPS::dEdx_long(double E)
{
	double ret = 0;
	// Loop over all field species:
	for(int i=0; i<num; i++)
	{
		ret += dEdx_long(E,i);
	}
	return ret;
}

// For use in GSL quantum integrations
struct quantum_params
{
	double beta_b, vp, ep, eb, kappa_b, K, mp, mb, mpb, Mpb;
	std::function<double(double)> eta_pb;
};

// Function to integrate for quantum correction (integrand of Eq 3.19)
double dEdxQ_func(double vpb, void * params) 
{
	struct quantum_params *p = (struct quantum_params *) params;
	double eta_pb = p->eta_pb(vpb);

	// Calculate first part of the integrand
	double theta = atan2(eta_pb,1); // need to correct for GSL normalization weirdness
	// gsl_sf_psi_1piy is psi(1+i*y)
	double term1 = 2.*gsl_sf_psi_1piy(eta_pb)*cos(theta) - log(pow(eta_pb,2));
	// terms in curly braces:
	double term2a = (1. + ((p->Mpb * p->vp)/(p->mb * vpb)) * (1./(p->beta_b * p->mb * p->vp * vpb) - 1.) ) 
		* exp(-0.5 * p->beta_b * p->mb * pow(p->vp - vpb, 2));
	double term2b = (1. + ((p->Mpb * p->vp)/(p->mb * vpb)) * (1./(p->beta_b * p->mb * p->vp * vpb) + 1.) ) 
		* exp(-0.5 * p->beta_b * p->mb * pow(p->vp + vpb, 2));
	return -1 * term1 * (term2a - term2b);
};

// Evaluate BPS quantum correction (Eq. 3.19) for a single species
double StopPow_BPS::dEdx_quantum(double E, int i)
{
	// test particle velocity
	double vp = c*sqrt(2e3*E/(mt*mpc2));

	double ret = 0.;
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (100);
	double prefac, result, err;
	double vb, v_min, v_max;

	// prefactor:
	prefac = (pow(Zt*e_LH,2)/(4*M_PI)) * (pow(kappa_b[i],2)/(2*beta_b[i]*mt*amu*pow(vp,2))) * sqrt(beta_b[i]*mf[i]*amu/(2*M_PI));

	// field particle thermal velocity
	vb = sqrt(3*kB*Tf[i]*keVtoK / (mf[i]*amu));
	// limits for numerical evaluation:
	v_min = fmin(vb, vp)/5.;
	v_max = fmax(vb, vp)*5.;

	// use a lambda to wrap eta_pb evaluation for field particle i
	auto etaFunc = [this,&i] (double vpb) {return this->eta_pb(vpb, i);};
	// parameters for numerical integration
	struct quantum_params params = {beta_b[i], vp, Zt*e_LH, Zf[i]*e_LH, kappa_b[i], K, mt*amu, mf[i]*amu, mpb[i], Mpb[i], etaFunc};

	// Use GSL library for evaluation:
	gsl_function Fc;
	Fc.function = dEdxQ_func;
	Fc.params = &params;
	gsl_integration_qag(&Fc, v_min, v_max, 1e-12, 1e-4, 100, 2, w, &result, &err);

	// with conversion from erg/cm to meV/um:
	ret += result * prefac * 624150.934 * 1e-4; // MeV/um

	gsl_integration_workspace_free (w);
	return ret;
}

// Evaluate BPS quantum correction (Eq. 3.19)
double StopPow_BPS::dEdx_quantum(double E)
{
	double ret = 0;
	// Loop over field particles:
	for(int i=0; i<num; i++)
	{
		ret += dEdx_quantum(E,i);
	}
	return ret;
}

// Get the minimum energy that can be used for dE/dx calculations
double StopPow_BPS::get_Emin()
{
	return Emin * mt;
}

// Get the maximum energy that can be used for dE/dx calculations
double StopPow_BPS::get_Emax()
{
	return Emax * mt;
}

// Function to evaluate rho_b (Eq 3.11 in paper)
double StopPow_BPS::rho_b(double v, int i) {
	return rho_b_prefac[i] * v * exp(-0.5 * beta_b[i] * mf[i]*amu * pow(v,2));
}

// Evaluate rho_total (Eq 3.10)
double StopPow_BPS::rho_tot(double v) {
	double ret = 0.;
	// sum over all field plasma species:
	for(int i=0; i<num; i++)
		ret += rho_b(v, i);
	return ret;
}

// Quantum parameter (Eq 3.1)
double StopPow_BPS::eta_pb(double vpb, int i) {
	return eta_pb_prefac[i] / vpb;
}

// Error function with a purely imaginary input: erfi = erf(i*z)/i
double StopPow_BPS::erfi(double z) {
	// Evaluation using GSL Dawson's Function
	// erfi(z) = Dawson(z) * 2/sqrt(pi) * exp(z^2)
	return gsl_sf_dawson(z) * 2/sqrt(M_PI) * exp(pow(z,2));
}

// Calculate the dielectric susceptibility integral (Eq 3.9):
gsl_complex StopPow_BPS::Fc(double u) {
	// Need to sum over plasma species. Simplification was made:
	// rho_b(v) = rho * v * exp(-a*v^2)
	// convenient form for evaluation, because then a Cauchy analysis (via Mathematica) is:
	// LaTeX syntax:
	// \int_{-\infty }^{\infty } \frac{\rho  v \exp \left(-a v^2\right)}{-i \eta +u-v} \, dv
	// \rho  \left(-\frac{\sqrt{\pi }}{\sqrt{a}}+(u-i \eta ) e^{-a (u-i \eta )^2} \left(\pi  \text{erfi}\left(\sqrt{a} (u-i \eta )\right)+\log (-u+i \eta )-\log (u-i \eta )\right)\right),\Im(u)\neq \Re(\eta )\land \Re(a)>0
	// Mathematica input:
	// Integrate[(\[Rho]*v*Exp[-a*v^2])/(u - v - I*\[Eta]), {v, -Infinity, Infinity}]
	// \[Rho] (-(Sqrt[\[Pi]]/Sqrt[a]) + E^(-a (u - I \[Eta])^2) (u - I \[Eta]) (\[Pi] Erfi[Sqrt[a] (u - I \[Eta])] - Log[u - I \[Eta]] + Log[-u + I \[Eta]])), Im[u] != Re[\[Eta]] && Re[a] > 0
	
	// Return value:
	gsl_complex ret = gsl_complex_rect(0,0);

	// loop over species:
	for(int i=0; i<num; i++)
	{
		double rho = pow(kappa_b[i],2)*sqrt(beta_b[i]*mf[i]*amu/(2*M_PI));
		double a = 0.5*beta_b[i]*mf[i]*amu;

		// for numerics
		double eta = 1e-6;

		// need to calculate with complex numbers:
		gsl_complex rhoc = gsl_complex_rect(rho, 0);
		gsl_complex acm = gsl_complex_rect(-1.*a, 0);
		gsl_complex pi_a = gsl_complex_rect(sqrt(M_PI/a), 0); // -sqrt(pi/a)
		gsl_complex uc = gsl_complex_rect(u, -1*eta);
		gsl_complex ucm = gsl_complex_rect(-1.*u, eta);
		gsl_complex uc2 = gsl_complex_pow(gsl_complex_rect(u, -1*eta), gsl_complex_rect(2,0));
		gsl_complex erf_term = gsl_complex_rect( M_PI*erfi(sqrt(a)*u), 0);
		gsl_complex exp_term = gsl_complex_exp(gsl_complex_mul(acm,uc2));

		// start building up return:
		gsl_complex temp = gsl_complex_sub(gsl_complex_log(ucm), gsl_complex_log(uc));
		temp = gsl_complex_add(temp, erf_term);
		temp = gsl_complex_mul(temp, uc);
		temp = gsl_complex_mul(temp, exp_term);
		temp = gsl_complex_sub(temp, pi_a);
		temp = gsl_complex_mul(rhoc, temp);
		temp = gsl_complex_mul(temp, gsl_complex_rect(-1,0));

		// can have trouble with ions...so add some sanity checks
		if(!isnan(GSL_REAL(temp)))
			ret = gsl_complex_add_real(ret, GSL_REAL(temp));
		if(!isnan(GSL_IMAG(temp)))
			ret = gsl_complex_add_imag(ret, GSL_IMAG(temp));
	}

	// need to flip sign:
	gsl_complex ret2 = gsl_complex_rect( GSL_REAL(ret) , -1*GSL_IMAG(ret) );

	return ret2;
}

double StopPow_BPS::Fc_real(double u)
{
	gsl_complex temp = Fc(u);
	return GSL_REAL(temp);
}

double StopPow_BPS::Fc_imag(double u)
{
	gsl_complex temp = Fc(u);
	return GSL_IMAG(temp);
}

// Do some precalculation of parameters
void StopPow_BPS::on_field_change()
{
	// allocate:
	Mpb.resize(num);
	mpb.resize(num);
	beta_b.resize(num);
	kappa_b.resize(num);
	rho_b_prefac.resize(num);
	eta_pb_prefac.resize(num);

	// reset Debye wavenumber
	kappa_D = 0.;
	// find electron index:
	int ie = 0;

	// loop over field species:
	for(int i=0; i<num; i++)
	{
		// Calculate total mass (Eq. 3.8)
		Mpb[i] = amu * (mt + mf[i]);
		// Calculate reduced mass (Eq. 3.7)
		mpb[i] = 1. / ( 1./(amu*mt) + 1./(amu*mf[i]) );
		// Calculate inverse temperature in energy units
		beta_b[i] = 1. / (kB * Tf[i] * keVtoK);
		// Calculate Debye wave number for a species (Eq. 3.5)
		kappa_b[i] = sqrt( beta_b[i] * pow(Zf[i]*e_LH,2) * nf[i] );

		// Calculate total Debye wave number (Eq. 3.6)
		kappa_D += pow(kappa_b[i],2);

		ie = (mf[i] < 0.9) ? i : ie;

		// prefactor for rho_b
		rho_b_prefac[i] = pow(kappa_b[i],2) * sqrt(beta_b[i] * mf[i]*amu/(2.*M_PI));

		// prefactor for eta_pb
		eta_pb_prefac[i] = Zf[i]*e_LH * Zt*e_LH / (4.*M_PI * hbar);
	}

	kappa_D = sqrt(kappa_D);
	// An arbitrary choice:
	K = 1.0*kappa_b[ie];
}

} // end namespace StopPow
