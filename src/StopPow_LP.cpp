#include "StopPow_LP.h"

namespace StopPow
{

const float StopPow_LP::Emin = 0.1; /* Minimum energy/A for dE/dx calculations */
const float StopPow_LP::Emax = 30; /* Maximum energy/A for dE/dx calculations */

// L-P specific initialization stuff:
void StopPow_LP::init()
{
	// set the info string:
	model_type = "Li-Petrasso";
	info = "";
}

// Li-Petrasso constructor primarily relies on superclass constructor
StopPow_LP::StopPow_LP(float mt_in, float Zt_in, std::vector<float> & mf_in, std::vector<float> & Zf_in, std::vector<float> & Tf_in, std::vector<float> & nf_in) throw(std::invalid_argument)
	: StopPow_Plasma::StopPow_Plasma(mt_in, Zt_in, mf_in, Zf_in, Tf_in, nf_in)
{
	init();
}
StopPow_LP::StopPow_LP(float mt, float Zt, std::vector< std::array<float,4> > & field) throw(std::invalid_argument)
	: StopPow_Plasma::StopPow_Plasma(mt, Zt, field)
{
	init();
}
StopPow_LP::StopPow_LP(float mt_in, float Zt_in, std::vector<float> & mf_in, std::vector<float> & Zf_in, std::vector<float> & Tf_in, std::vector<float> & nf_in, float Te) throw(std::invalid_argument)
	: StopPow_Plasma::StopPow_Plasma(mt_in, Zt_in, mf_in, Zf_in, Tf_in, nf_in, Te)
{
	init();
}
StopPow_LP::StopPow_LP(float mt, float Zt, std::vector< std::array<float,4> > & field, float Te) throw(std::invalid_argument)
	: StopPow_Plasma::StopPow_Plasma(mt, Zt, field, Te)
{
	init();
}

// Destructor
StopPow_LP::~StopPow_LP()
{
	// nothing to do
}

// Calculate the total stopping power
float StopPow_LP::dEdx_MeV_um(float E) throw(std::invalid_argument)
{
	// sanity check:
	if( E < Emin || E > Emax )
	{
		std::stringstream msg;
		msg << "Energy passed to StopPow_LP::dEdx is bad: " << E;
		throw std::invalid_argument(msg.str());
	}

	float ret = 0; // return value

	//iterate over all field particles:
	for(int i=0; i < num; i++)
	{
		float dEdx_single; // stopping power due to this species
		dEdx_single = LogLambda(E,i)*G(E,i); // standard term
		// collective effects:
		if(collective)
		{
			//dEdx_single += 0.5*log(1.261*xtf(E,i));
			float xInvSqrt = 1./sqrt(xtf_collective(E,i));
			float LogLambdaC = gsl_sf_bessel_K0(xInvSqrt) 
							* gsl_sf_bessel_K1(xInvSqrt) * xInvSqrt;
			dEdx_single += LogLambdaC;
		}

		// calculate prefactor for the term:
		float vt = c*sqrt(2*E*1e3/(mt*mpc2)); // test particle velocity
		float tmp = pow(Zt*e/vt,2.0);
		float wpf = sqrt(4*M_PI*nf[i]*pow(Zf[i]*e,2.0)/(mf[i]*mp));// plasma frequency
		dEdx_single = -tmp*wpf*wpf*dEdx_single; // erg/cm
		dEdx_single = dEdx_single*(1e-13)/(1.602e-19); // MeV/cm
		ret += dEdx_single;
	}

	return ret*1e-4; // MeV/um
}

// Calculate the total stopping power
float StopPow_LP::dEdx_MeV_mgcm2(float E) throw(std::invalid_argument)
{
	return (dEdx_MeV_um(E)*1e4) / (rho*1e3);
}

// Turn collective effects on or off.
void StopPow_LP::set_collective(bool set)
{
	collective = set;
}

// Turn quantum effects on or off.
void StopPow_LP::set_quantum(bool set)
{
	quantumT = set;
}

// Set factor for calculating binary collision xtf
void StopPow_LP::set_xtf_factor(float a)
{
	xtf_factor = a;
}

// Set factor for calculating collective effects xtf
void StopPow_LP::set_xtf_collective_factor(float a)
{
	xtf_factor = a;
}

// Set factor for calculating u
void StopPow_LP::set_u_factor(float a)
{
	u_factor = a;
}

// Get the minimum energy that can be used for dE/dx calculations
float StopPow_LP::get_Emin()
{
	return Emin * mt;
}

// Get the maximum energy that can be used for dE/dx calculations
float StopPow_LP::get_Emax()
{
	return Emax * mt;
}

// Calculate the Coulomb logarithm
float StopPow_LP::LogLambda(float E, int index)
{
	// reduced mass:
	float mr = mp*mt*mf[index]/(mt+mf[index]);
	// relative velocity:
	float u1 = u(E,index);
	// classical b90:
	float pperp = Zf[index]*e*Zt*e / (mr*u1*u1);
	// L-P style quantum b:
	float pmin = sqrt( pow(pperp,2.0) + pow(hbar/(2*mr*u1),2.0) );

	// calculate LogLambda:
	float LogLambda = 0.5*log(1 + pow(lDebye()/pmin,2.0) );
	// sanity. LogLambda cannot be negative:
	if( LogLambda > 0. )
		return LogLambda;
	return 0;
}

// Chandrasekhar function
float StopPow_LP::G(float E, int index)
{
	float rat = mf[index] / mt; // mass ratio
	float mu = 1.12838*sqrt(xtf(E,index))*exp(-xtf(E,index));
	float erfunc = erf( sqrt(xtf(E,index)) );
	return (erfunc  - mu) - rat*(mu - erfunc/LogLambda(E,index)) ;
}

// Debye length in field plasma
float StopPow_LP::lDebye()
{
	float ret = 0; // temporary return value
	//iterate over all field particles:
	for(int i=0; i < num; i++)
	{
		ret += (4.*M_PI*nf[i]*pow(Zf[i]*e,2.0) / (kB*Tq(i)*keVtoK) );
	}
	return 1.0/sqrt(ret);
}

/* Field particle thermal velocity */
float StopPow_LP::vtf(int index)
{
	return vtf(index, 2.);
}

/* Field particle thermal velocity with specified constant */
float StopPow_LP::vtf(int index, float constant)
{
	return c*sqrt(constant*Tq(index)/(mpc2*mf[index]));
}

/* x^{t/f} parameter from Li 1993 */
float StopPow_LP::xtf(float E, int index)
{
	// test particle velocity:
	float vt = c*sqrt(2*E*1e3/(mt*mpc2));
	float vf = vtf(index, xtf_factor); // sqrt(2kT/m) by default
	return pow( vt/vf ,2);
}

/* x^{t/f} parameter from Li 1993 for the collective effects term */
float StopPow_LP::xtf_collective(float E, int index)
{
	// test particle velocity:
	float vt = c*sqrt(2*E*1e3/(mt*mpc2));
	float vf = vtf(index, xtf_collective_factor); // sqrt(kT/m) by default
	return pow( vt/vf ,2);
}

/* Relative velocity between test particle and field particle */
float StopPow_LP::u(float E, int index)
{
	float vt = c*sqrt(2.*E*1e3/(mt*mpc2));
	float vf = vtf(index, u_factor); // sqrt(8kT/pi*m) by default
	
	// simple model:
	//return sqrt( (pow(vt,2) + pow(vf,2)) );
	
	// more exact:
	float u = (vf/2.)*exp(-4.*vt*vt/(M_PI*vf*vf))
		+ vt*(1.+M_PI*vf*vf/(8.*vt*vt))
		* erf(sqrt(4.*vt*vt/(M_PI*vf*vf)));

	return u;
}

/* Quantum-corrected temperature */
float StopPow_LP::Tq(int index)
{
	if(quantumT)
	{
		// degeneracy parameter: 
		double TF = (1/kB)*((double)hbar*(double)hbar/(2.*mf[index]*amu))*pow(3*M_PI*M_PI*nf[index], 2./3.);
		float theta = Tf[index]*keVtoK / TF;
		// chemical potential over kT: Drake Eq 3.20
		float mukT = -1.5*log(theta) + log(4./(3.*sqrt(M_PI))) + (0.25054*pow(theta,-1.858) + 0.072*pow(theta,-1.858/2))/(1+0.25054*pow(theta,-0.858));
		// Effective temperature from Drake Eq 3.22
		return Tf[index] * gsl_sf_fermi_dirac_3half(mukT)/gsl_sf_fermi_dirac_half(mukT);
	}
	return Tf[index];
}

} // end namespace StopPow
