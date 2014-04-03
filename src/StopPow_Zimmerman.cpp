#include "StopPow_Zimmerman.h"

namespace StopPow
{

const float StopPow_Zimmerman::Emin = 0.1; /* Minimum energy for dE/dx calculations */
const float StopPow_Zimmerman::Emax = 30; /* Maximum energy for dE/dx calculations */


// Initialization routines specific to this model
void StopPow_Zimmerman::init()
{
	// set the info string:
	model_type = "Zimmerman";
	info = "";
}
// constructors for partial ionized material:
StopPow_Zimmerman::StopPow_Zimmerman(float mt_in, float Zt_in, std::vector<float> & mf_in, std::vector<float> & Zf_in, std::vector<float> & Tf_in, std::vector<float> & nf_in, std::vector<float> & Zbar_in, float Te_in) throw(std::invalid_argument)
	: StopPow_PartialIoniz::StopPow_PartialIoniz(mt_in, Zt_in, mf_in, Zf_in, Tf_in, mf_in, Zbar_in, Te)
{
	init();
}
StopPow_Zimmerman::StopPow_Zimmerman(float mt, float Zt, std::vector< std::array<float,5> > & field, float Te) throw(std::invalid_argument)
	: StopPow_PartialIoniz::StopPow_PartialIoniz(mt, Zt, field, Te)
{
	init();
}

// Destructor
StopPow_Zimmerman::~StopPow_Zimmerman()
{
	// nothing to do
}

// Calculate stopping power
float StopPow_Zimmerman::dEdx_MeV_um(float E) throw(std::invalid_argument)
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
float StopPow_Zimmerman::dEdx_MeV_mgcm2(float E) throw(std::invalid_argument)
{
	return (dEdx_MeV_um(E)*1e4) / (rho*1e3);
}

// Free electron stopping power
float StopPow_Zimmerman::dEdx_free_electron(float E)
{
	// test particle velocity
	float vt = c*sqrt(2e3*E/(mt*mpc2));
	// y parameter just ratio of test / thermal velocity
	float vth = sqrt(2.*kB*Te*keVtoK/me); // nondegenerate! He also gives Eq 18 for degenerate
	float y = vt/vth;
	float omega_pe = sqrt(4*M_PI*esu*esu*ne/me);
	// Eq 16:
	float LambdaF = (4*M_PI*me*pow(vth,2)/(h*omega_pe)) * (0.321+0.259*pow(y,2)+0.0707*pow(y,4)+0.05*pow(y,6))/(1.+0.130*pow(y,2)+0.05*pow(y,4));
	float dEdx_F = 4*M_PI*(1./me)*pow(esu,4)*pow(Zt/vt,2)*ne*LF(y,LambdaF);
	return -1.*dEdx_F * 624150.934 * 1e-4; // MeV/um
}

// Bound electron stopping power
float StopPow_Zimmerman::dEdx_bound_electron(float E)
{
	// test particle velocity
	float vt = c*sqrt(2e3*E/(mt*mpc2));

	float dEdx_BE = 0.;
	float Ibar, ZiB, LiB;
	double prefac = (4.*M_PI*pow(e,4)*pow(Zt,2) / (me*pow(vt,2)));
	// have to loop over all ions
	for(int i=0; i<num; i++)
	{
		ZiB = (Zf[i] - Zbar[i]); // number of bound electrons
		if(ZiB > 0.)
		{
			// Eq 20:
			Ibar = Zf[i] * (0.024 - 0.013 * ZiB/Zf[i]) / sqrt(ZiB/Zf[i]);
			LiB = log(2. * me * pow(vt,2) / Ibar);
			dEdx_BE +=  prefac * nf[i] * ZiB * LiB;
		}
	}
	return -1.*dEdx_BE * 624150.934 * 1e-4; // MeV/um
}

// Ion stopping power
float StopPow_Zimmerman::dEdx_ion(float E)
{
	// test particle velocity
	float vt = c*sqrt(2e3*E/(mt*mpc2));

	// need to loop over field ions
	float dEdx_I = 0.;
	float mr, bi, Li;
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

// Minimum energy limit
float StopPow_Zimmerman::get_Emin()
{
	return Emin;
}

// Maximum energy limit
float StopPow_Zimmerman::get_Emax()
{
	return Emax;
}

// Electron stopping number
float StopPow_Zimmerman::LF(float y, float LambdaF)
{
	return 0.5 * log(1.+pow(LambdaF,2)) * (gsl_sf_erf(y) - (2./sqrt(M_PI))*y*exp(-y*y));
}

// Debye length in field plasma
float StopPow_Zimmerman::lDebye()
{
	float ret = 0; // temporary return value
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
