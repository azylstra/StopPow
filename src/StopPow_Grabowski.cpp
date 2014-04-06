#include "StopPow_Grabowski.h"

namespace StopPow
{

const double StopPow_Grabowski::Emin = 0.1; /* Minimum energy/A for dE/dx calculations */
const double StopPow_Grabowski::Emax = 30; /* Maximum energy/A for dE/dx calculations */

// specific initialization stuff:
void StopPow_Grabowski::init()
{
	// set the info string:
	model_type = "Grabowski";
	info = "";
}

// constructors primarily rely on superclass constructors
StopPow_Grabowski::StopPow_Grabowski(double mt_in, double Zt_in, std::vector<double> & mf_in, std::vector<double> & Zf_in, std::vector<double> & Tf_in, std::vector<double> & nf_in) throw(std::invalid_argument)
	: StopPow_Plasma::StopPow_Plasma(mt_in, Zt_in, mf_in, Zf_in, Tf_in, nf_in)
{
	init();
}
StopPow_Grabowski::StopPow_Grabowski(double mt, double Zt, std::vector< std::array<double,4> > & field) throw(std::invalid_argument)
	: StopPow_Plasma::StopPow_Plasma(mt, Zt, field)
{
	init();
}
StopPow_Grabowski::StopPow_Grabowski(double mt_in, double Zt_in, std::vector<double> & mf_in, std::vector<double> & Zf_in, std::vector<double> & Tf_in, std::vector<double> & nf_in, double Te) throw(std::invalid_argument)
	: StopPow_Plasma::StopPow_Plasma(mt_in, Zt_in, mf_in, Zf_in, Tf_in, nf_in, Te)
{
	init();
}
StopPow_Grabowski::StopPow_Grabowski(double mt, double Zt, std::vector< std::array<double,4> > & field, double Te) throw(std::invalid_argument)
	: StopPow_Plasma::StopPow_Plasma(mt, Zt, field, Te)
{
	init();
}

// Destructor
StopPow_Grabowski::~StopPow_Grabowski()
{
	// nothing to do
}

// Calculate the total stopping power
double StopPow_Grabowski::dEdx_MeV_um(double E) throw(std::invalid_argument)
{
	// sanity check:
	if( E < Emin || E > Emax )
	{
		std::stringstream msg;
		msg << "Energy passed to StopPow_Grabowski::dEdx is bad: " << E;
		throw std::invalid_argument(msg.str());
	}

	double ret = 0; // return value

	//iterate over all field particles:
	for(int i=0; i < num; i++)
	{
		// for convenience:
		double Tf_K = Tf[i] * keVtoK;

		// test particle velocity:
		double v = c*sqrt(2000.*E/(mt*mpc2));

		// Wigner-Seitz radius
		double G_r0 = pow(4*M_PI*nf[i]/3., -1/3.);

		// thermal velocity:
		double vth = sqrt(kB*Tf_K/(amu*mf[i]));

		// intratarget coupling parameter Gamma = qe^2 / (G_r0*Te)
		double Gamma = pow(Zf[i]*esu,2) / (G_r0*kB*Tf_K);
		// top left of pg 2
		double g = sqrt(3)*fabs(Zt)*pow(Gamma,1.5);
		double s = G_d*pow(1.+G_c*g, 1/3);
		double w = v/(vth*s);

		// normalization factor for stopping power
		// (Z^2 qe^2 / lD^2) / (1+g)^2/3
		double lD = sqrt(kB*Tf_K / (4*M_PI*nf[i]*esu*esu));
		double norm = pow(Zt*Zf[i]*esu/lD, 2) / pow(1+g, 2/3);

		ret += (-R(w,g,s,Zt) * (G(w)*log(pow(M_E,.5) + (alpha+w*w)/G_g0) + H(w))) * norm * 624150.934; // MeV/cm
	}

	return ret*1e-4; // MeV/um
}

// Calculate the total stopping power
double StopPow_Grabowski::dEdx_MeV_mgcm2(double E) throw(std::invalid_argument)
{
	return (dEdx_MeV_um(E)*1e4) / (rho*1e3);
}

// Get the minimum energy that can be used for dE/dx calculations
double StopPow_Grabowski::get_Emin()
{
	return Emin * mt;
}

// Get the maximum energy that can be used for dE/dx calculations
double StopPow_Grabowski::get_Emax()
{
	return Emax * mt;
}

// Eq 2 in paper
double StopPow_Grabowski::M1(double g, double s, double Z) {
	return s*log(1.+alpha*pow(M_E,-0.5)/(g*(1+G_a*Z*Z*g))) / log(1. + alpha*pow(M_E,-0.5)/G_g0);
}

// Eq 2 in paper
double StopPow_Grabowski::M2(double w, double g, double s) {
	return (1./pow(s,2)) * log(1 + pow(s*w,3)/g) / log(1+pow(w,3)/G_g0);
}

// Eq 2 in paper
double StopPow_Grabowski::R(double w, double g, double s, double Z) {
	return (M1(g,s,Z) + G_b*M2(w,g,s)*pow(w,2))*pow(1+g,2/3) / (w*w*(1.+G_b*w*w));
}

// Eq 2 in paper
double StopPow_Grabowski::G(double w) {
	return gsl_sf_erf(w/sqrt(2)) - sqrt(2/M_PI)*w*exp(-w*w/2.);
}

// Eq 2 in paper
double StopPow_Grabowski::H(double w) {
	return pow(w,4.)*log(w)/(12.+pow(w,4)) - pow(w,3)*exp(-w*w/2.)/(3*sqrt(2*M_PI));
}

} // end namespace StopPow
