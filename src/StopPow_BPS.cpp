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
}

// Li-Petrasso constructor primarily relies on superclass constructor
StopPow_BPS::StopPow_BPS(double mt_in, double Zt_in, std::vector<double> & mf_in, std::vector<double> & Zf_in, std::vector<double> & Tf_in, std::vector<double> & nf_in) throw(std::invalid_argument)
	: StopPow_Plasma::StopPow_Plasma(mt_in, Zt_in, mf_in, Zf_in, Tf_in, nf_in)
{
	init();
}
StopPow_BPS::StopPow_BPS(double mt, double Zt, std::vector< std::array<double,4> > & field) throw(std::invalid_argument)
	: StopPow_Plasma::StopPow_Plasma(mt, Zt, field)
{
	init();
}
StopPow_BPS::StopPow_BPS(double mt_in, double Zt_in, std::vector<double> & mf_in, std::vector<double> & Zf_in, std::vector<double> & Tf_in, std::vector<double> & nf_in, double Te) throw(std::invalid_argument)
	: StopPow_Plasma::StopPow_Plasma(mt_in, Zt_in, mf_in, Zf_in, Tf_in, nf_in, Te)
{
	init();
}
StopPow_BPS::StopPow_BPS(double mt, double Zt, std::vector< std::array<double,4> > & field, double Te) throw(std::invalid_argument)
	: StopPow_Plasma::StopPow_Plasma(mt, Zt, field, Te)
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

	return dEdx_short(E) + dEdx_long(E) + dEdx_quantum(E); // MeV/um
}

// Calculate the total stopping power
double StopPow_BPS::dEdx_MeV_mgcm2(double E) throw(std::invalid_argument)
{
	return (dEdx_MeV_um(E)*1e4) / (rho*1e3);
}

// For use in GSL integrations
struct int_params
{
	double beta_b, vp, ep, eb, kappa_b, K, mp, mb, mpb, Mpb;
};

// Function to integrate for short-range stopping power
double dEdxcs_func (double u, void * params) {
	struct int_params *p = (struct int_params *) params;
	double term1 = sqrt(u) * exp(-0.5 * p->beta_b * p->mb * pow(p->vp,2) * u);
	double term2 = ( (-log(p->beta_b * (fabs(p->ep * p->eb) * p->K/(4. * M_PI)) * (p->mb/p->mpb) * u/(1. - u)) + 2. - 2. * 0.5772)
	 	 * (p->beta_b * p->Mpb * pow(p->vp,2) - 1./u) + 2./u );
	//std::cout << term1 << " , " << term2 << " , " << p->beta_b * (fabs(p->ep * p->eb) * p->K/(4 * M_PI)) * (p->mb/p->mpb) * u/(1. - u) << std::endl;
	return term1 * term2;
}
// Classical short-range stopping power (Eq. 3.3)
double StopPow_BPS::dEdx_short(double E)
{
	double vp = c*sqrt(2e3*E/(mt*mpc2));
	double ret = 0.;
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (100);
	double prefac, result, err;

	// Loop over field particles:
	for(int i=0; i<num; i++)
	{
		struct int_params params = {beta_b(i), vp, Zt*e_LH, Zf[i]*e_LH, kappa_b(i), K(i), mt*amu, mf[i]*amu, mpb(i), Mpb(i)};
		prefac = (pow(Zt*e_LH,2)/(4.*M_PI)) * (pow(kappa_b(i),2)/(mt*amu*vp)) * sqrt(mf[i]*amu/(2.*M_PI*beta_b(i)));
		gsl_function Fa;
		Fa.function = dEdxcs_func;
		Fa.params = &params;
		gsl_integration_qag(&Fa, 0, 1, 1e-8, 1e-6, 100, 1, w, &result, &err);

		ret += prefac * result;
	}

	gsl_integration_workspace_free (w);
	return ret * 624150.934 * 1e-4; // MeV/um
}

// Helper function for calculating long-range stopping power
gsl_complex StopPow_BPS::dEdx_long_func(double E, double x, int i)
{
	double vp = c*sqrt(2e3*E/(mt*mpc2));

	double r_b = rho_b(vp*x, i);
	double r_tot = rho_tot(vp*x);

	// Evaluate F:
	gsl_complex F1 = Fc(vp*x); // F(vp*costheta)

	// prefactor:
	double part1 = (pow(Zt*e_LH,2)/(8*M_PI*M_PI)) * x * (r_b/r_tot); // this part is nice and real

	// need to use complex for the second part:
	gsl_complex ret = gsl_complex_rect(0, part1); // includes factor of i in equation!
	ret = gsl_complex_mul(ret, F1);
	gsl_complex logval = gsl_complex_div(F1, gsl_complex_rect(pow(K(i),2), 0.));
	ret = gsl_complex_mul(ret, gsl_complex_log(logval));

	return ret;
}

// Classical long-range stopping power (Eq 3.4)
double StopPow_BPS::dEdx_long(double E)
{
	double vp = c*sqrt(2e3*E/(mt*mpc2)); // test particle velocity
	double du = 0.01; // step size in numerical integration

	double dEdx_cR = 0.; // return value

	// Loop over all field species:
	for(int i=0; i<num; i++)
	{
		// Evaluation of the first term of the equation (part with integral):
		gsl_complex dEdx_cR_1 = gsl_complex_rect(0,0); // first term (one with integral)
		for(double u=-1+du/2; u<1.; u+=du) {
			gsl_complex temp = gsl_complex_mul(dEdx_long_func(E, u, i) , gsl_complex_rect(du,0));
			if(!isnan(GSL_REAL(temp)) and !isinf(GSL_REAL(temp))) {
				dEdx_cR_1 = gsl_complex_add(temp, dEdx_cR_1); // add value
			}
		}
		dEdx_cR_1 = gsl_complex_mul(dEdx_cR_1, gsl_complex_rect(624150.934 * 1e-4,0)); // Convert to MeV/um

		// Evaluation of the second term of the equation
		gsl_complex prefac2 = gsl_complex_rect(0, (pow(Zt*e_LH,2)/(8*M_PI*M_PI)) * (1./(beta_b(i)*mt*amu*pow(vp,2))) 
			* ( rho_b(vp, i) / rho_tot(vp) ));
		// Calculate F and F*
		gsl_complex Fv = Fc(vp);
		gsl_complex Fvc = gsl_complex_conjugate(Fv); //Fc(-1*vp);
		// calculate two complex terms inside the square brackets of Eq 3.4
		gsl_complex term1 = gsl_complex_mul( Fv, gsl_complex_log( gsl_complex_div(Fv,gsl_complex_rect(pow(K(i),2),0)) ) );
		gsl_complex term2 = gsl_complex_mul( Fvc, gsl_complex_log( gsl_complex_div(Fvc, gsl_complex_rect(pow(K(i),2),0)) ) );
		gsl_complex term = gsl_complex_sub(term1,term2);
		term = gsl_complex_mul(term, prefac2);
		// std::cout << std::endl << "Fv = " << GSL_REAL(Fv) << " , " << GSL_IMAG(Fv) << std::endl;
		// std::cout << "prefac = " << GSL_REAL(prefac2) << " , " << GSL_IMAG(prefac2) << std::endl;
		// std::cout << vp << std::endl;
		//cout << GSL_REAL(term) << " + " << GSL_IMAG(term) << " i"  << endl;
		gsl_complex dEdx_cR_2 = gsl_complex_mul(term, gsl_complex_rect(624150.934*1e-4,0)); // Convert to MeV/um

		dEdx_cR = GSL_REAL(dEdx_cR_1) - GSL_REAL(dEdx_cR_2);
		//std::cout << std::endl << "long term: " << GSL_REAL(dEdx_cR_1) << " , " << GSL_REAL(dEdx_cR_2) << std::endl;
	}

	return dEdx_cR;
}

// For use in GSL quantum integrations
struct quantum_params
{
	double beta_b, vp, ep, eb, kappa_b, K, mp, mb, mpb, Mpb;
	std::function<double(double)> eta_pb;
};

// Function to integrate for quantum correction
double dEdxQ_func(double vpb, void * params) 
{
	struct quantum_params *p = (struct quantum_params *) params;
	double eta_pb = p->eta_pb(vpb);

	// Calculate first part of the integrand
	double theta = atan2(eta_pb,1); // need to correct for GSL normalization weirdness
	double term1 = 2.*gsl_sf_psi_1piy(eta_pb)*cos(theta) - log(pow(eta_pb,2));
	double term2a = (1. + ((p->Mpb * p->vp)/(p->mb * vpb)) * (1./(p->beta_b * p->mb * p->vp * vpb) - 1.) ) 
		* exp(-0.5 * p->beta_b * p->mb * pow(p->vp - vpb, 2));
	double term2b = (1. + ((p->Mpb * p->vp)/(p->mb * vpb)) * (1./(p->beta_b * p->mb * p->vp * vpb) + 1.) ) 
		* exp(-0.5 * p->beta_b * p->mb * pow(p->vp + vpb, 2));
	return term1 * (term2a - term2b);
};

double StopPow_BPS::dEdx_quantum(double E)
{
	// test particle velocity
	double vp = c*sqrt(2e3*E/(mt*mpc2));

	double ret = 0.;
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (100);
	double prefac, result, err;

	// Loop over field particles:
	for(int i=0; i<num; i++)
	{
		// prefactor:
		prefac = (pow(Zt*e_LH,2)/(4*M_PI)) * (pow(kappa_b(i),2)/(2*beta_b(i)*mt*amu*pow(vp,2))) * sqrt(beta_b(i)*mf[i]*amu/(2*M_PI));

		// field particle thermal velocity
		double vb = sqrt(3*kB*Tf[i]*keVtoK / (mf[i]*amu));
		// limits for numerical evaluation:
		double v_min = fmin(vb, vp)/5.;
		double v_max = fmax(vb, vp)*5.;

		auto etaFunc = [this,&i] (double vpb) {return this->eta_pb(vpb, i);};
		struct quantum_params params = {beta_b(i), vp, Zt*e_LH, Zf[i]*e_LH, kappa_b(i), K(i), mt*amu, mf[i]*amu, mpb(i), Mpb(i), etaFunc};

		gsl_function Fc;
		Fc.function = dEdxQ_func;
		Fc.params = &params;
		gsl_integration_qag(&Fc, v_min, v_max, 1e-12, 1e-6, 100, 2, w, &result, &err);
		ret += result * prefac * 624150.934 * 1e-4; // MeV/um
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

// Calculate total mass (Eq. 3.8)
double StopPow_BPS::Mpb(int i)
{
	return amu * (mt + mf[i]);
}

// Calculate reduced mass (Eq. 3.7)
double StopPow_BPS::mpb(int i)
{
	return 1. / ( 1./(amu*mt) + 1./(amu*mf[i]) );
}

// Calculate normalized temperature
double StopPow_BPS::beta_b(int i)
{
	return 1. / (kB * Tf[i] * keVtoK);

}

// Calculate Debye wave number for a species (Eq. 3.5)
double StopPow_BPS::kappa_b(int i)
{
	return sqrt( beta_b(i) * pow(Zf[i]*e_LH,2) * nf[i] );

}

// Calculate total Debye wave number (Eq. 3.6)
double StopPow_BPS::kappa_D()
{
	double ret = 0.;
	for(int i=0; i<num; i++)
		ret += pow(kappa_b(i),2);
	return sqrt(ret);
}

// Calculate wave number to use
double StopPow_BPS::K(int i)
{
	// find electron index:
	int ie = 0;
	for(int j=0; j<num; j++)
		ie = (mf[j] < 0.9) ? j : ie;
	return 1.0*kappa_b(ie);
}

// Function to evaluate rho_b (Eq 3.11 in paper)
double StopPow_BPS::rho_b(double v, int i) {
	return pow(kappa_b(i),2) * sqrt(beta_b(i) * mf[i]*amu/(2.*M_PI)) * v * exp(-0.5 * beta_b(i) * mf[i]*amu * pow(v,2));
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
	return Zf[i]*e_LH * Zt*e_LH / (4.*M_PI * hbar * vpb);
}

// Error function with a purely imaginary input: erfi = erf(i*z)/i
double StopPow_BPS::erfi(double z) {
	// Evaluation using GSL Dawson's Function
	// erfi(z) = Dawson(z) * 2/sqrt(pi) * exp(z^2)
	return gsl_sf_dawson(z) * 2/sqrt(M_PI) * exp(pow(z,2));
}

// Imaginary error function with a complex argument
gsl_complex StopPow_BPS::erfi(gsl_complex z) {
	// double real_part = erfi(GSL_REAL(z));
	// double zi = GSL_IMAG(z);
	// double imag_part = (1./sqrt(M_PI)) * (2.*zi - (2./3.)*pow(zi,3) + (1./5.)*pow(zi,5) - (1./21.)*pow(zi,7) + (1./108.)*pow(zi,9) - (1./660.)*pow(zi,11) + (1./4680.)*pow(zi,13) - (1./37800.)*pow(zi,15));
	// return gsl_complex_rect(real_part, imag_part);

	// z = x + i*y        x,y real
	double x = GSL_REAL(z);
	double y = GSL_IMAG(z);
	double real_part = ( x*(2 - 2*pow(y,2) + pow(y,4)) + pow(x,3)*( 2./3. - 2.*pow(y,2) + (5./3.)*pow(y,4)) ) / sqrt(M_PI);
	double imag_part = ( y*(2 - 2*pow(y,3)/3.) + pow(x,2)*(2*y - 2*pow(y,3)) + pow(x,4)*(y - (5./3.)*pow(y,3)) ) / sqrt(M_PI);
	return gsl_complex_rect(real_part, imag_part);
}

// Calculate the dielectric susceptibility integral:
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
		double rho = pow(kappa_b(i),2)*sqrt(beta_b(i)*mf[i]*amu/(2*M_PI));
		double a = 0.5*beta_b(i)*mf[i]*amu;

		// for numerics
		double eta = 1e-6;

		// need to calculate with complex numbers:
		gsl_complex rhoc = gsl_complex_rect(rho, 0);
		gsl_complex acm = gsl_complex_rect(-1.*a, 0);
		gsl_complex pi_a = gsl_complex_rect(sqrt(M_PI/a), 0); // -sqrt(pi/a)
		gsl_complex uc = gsl_complex_rect(u, -1*eta);
		gsl_complex ucm = gsl_complex_rect(-1.*u, eta);
		gsl_complex uc2 = gsl_complex_pow(gsl_complex_rect(u, -1*eta), gsl_complex_rect(2,0));
		//gsl_complex erf_arg = gsl_complex_mul_real(uc, sqrt(a));
		//gsl_complex erf_term = gsl_complex_mul_real(erfi(erf_arg), M_PI);
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

		// can have trouble with ions...
		if(!isnan(GSL_REAL(temp)))
			ret = gsl_complex_add_real(ret, GSL_REAL(temp));
		if(!isnan(GSL_IMAG(temp)))
			ret = gsl_complex_add_imag(ret, GSL_IMAG(temp));
		//ret = gsl_complex_add(ret, temp);
	}

	gsl_complex ret2 = gsl_complex_rect( GSL_REAL(ret) , -1*GSL_IMAG(ret) );

	return ret2;
}

} // end namespace StopPow
