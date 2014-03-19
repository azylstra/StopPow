#include "StopPow_LP.h"

namespace StopPow
{

const float StopPow_LP::Emin = 0.1; /* Minimum energy for dE/dx calculations */
const float StopPow_LP::Emax = 30; /* Maximum energy for dE/dx calculations */

/** Initialize the Li-Petrasso stopping power.
 * @param mt the test particle mass in AMU
 * @param Zt the test particle in charge (units of e)
 * @param mf vector containing ordered field particle masses in AMU
 * @param Zt vector containing ordered field particle charges in units of e
 * @param Tf vector containing ordered field particle temperatures in units of keV
 * @param nf vector containing ordered field particle densities in units of 1/cm3
 * @throws invalid_argument
 */
StopPow_LP::StopPow_LP(float mt_in, float Zt_in, std::vector<float> mf_in, std::vector<float> Zf_in, std::vector<float> Tf_in, std::vector<float> nf_in) throw(std::invalid_argument)
{
	// default mode for LP:
	set_mode(MODE_LENGTH);
	
	// infer size of the field particle arrays:
	num = mf_in.size();

	// sanity checking. 
	bool args_ok = true;
	// Make sure mt and Zt are positive,
	// and that all field particle arrays have same size
	if( mt_in <= 0 || Zt_in <= 0
		|| Zf_in.size() != num 
		|| Tf_in.size() != num
		|| nf_in.size() != num )
	{
		args_ok = false;
	}

	// now do sanity checking on the field particle values:
	for(int i=0; i<num; i++)
	{
		args_ok = args_ok && mf_in[i] > 0;
		args_ok = args_ok && Tf_in[i] > 0;
		args_ok = args_ok && nf_in[i] > 0;
	}

	// throw an exception if necessary:
	if( !args_ok )
	{
		std::stringstream msg;
		// start constructing message, add info on mt and Zt:
		msg << "Values passed to StopPow_LP constructor are bad: " 
		 << mt_in << "," << Zt_in << "," << std::endl;

		std::vector<float>::iterator it; // to iterate over field particles

		// add each element in mf:
		msg << "mf = ";
		for(it=mf.begin(); it<mf.end(); it++)
		 	msg << (*it) << ",";

		// add each element in Zf:
		msg << std::endl << "Zf = ";
		for(it=Zf.begin(); it<Zf.end(); it++)
		 	msg << (*it) << ",";

		// add each element in Tf:
		msg << std::endl << "Tf = ";
		for(it=Tf.begin(); it<Tf.end(); it++)
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
	Tf = Tf_in;
	nf = nf_in;
	collective = true;

	// calculate the field particle mass density:
	rho = 0; // g/cm3
	// iterate over field particles:
	for(int i=0; i<num; i++)
	{
		rho += mf[i] * mp * nf[i];
	}

	// set the info string:
	model_type = "Li-Petrasso";
	info = "";
}

// Destructor
StopPow_LP::~StopPow_LP()
{
	// nothing to do
}

/** Calculate the total stopping power
 * @param E the test particle energy in MeV
 * @return stopping power in units of MeV/um
 * @throws invalid_argument
 */
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
		float dEdx_single = 0; // stopping power due to this species
		dEdx_single += LogLambda(E,i)*G(E,i); // standard term
		// collective effects:
		if(collective)
		{
			//dEdx_single += 0.5*log(1.261*xtf(E,i));
			float xInvSqrt = 1/sqrt(xtf_collective(E,i));
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

/** Calculate the total stopping power
 * @param E the test particle energy in MeV
 * @return stopping power in units of MeV/(mg/cm2)
 * @throws invalid_argument
 */
float StopPow_LP::dEdx_MeV_mgcm2(float E) throw(std::invalid_argument)
{
	return (dEdx_MeV_um(E)*1e4) / (rho*1e3);
}

/** Turn collective effects on or off.
* @param set if you want to use collective effects
*/
void StopPow_LP::set_collective(bool set)
{
	collective = set;
}


/**
 * Get the minimum energy that can be used for dE/dx calculations
 * @return Emin in MeV
 */
float StopPow_LP::get_Emin()
{
	return Emin;
}

/**
 * Get the maximum energy that can be used for dE/dx calculations
 * @return Emax in MeV
 */
float StopPow_LP::get_Emax()
{
	return Emax;
}


/** Calculate the Coulomb logarithm
 * @param E the test particle energy in MeV
 * @param index the field particle index
 * @return value of Log(Lambda)
 */
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

/** Chandrasekhar function
 * @param E the test particle energy in MeV
 * @param index the field particle index
 * @return value of the Chandrasekhar function G (see L-P 1993)
 */
float StopPow_LP::G(float E, int index)
{
	float rat = mf[index] / mt; // mass ratio
	float mu = 1.12838*sqrt(xtf(E,index))*exp(-xtf(E,index));
	float erfunc = erf( sqrt(xtf(E,index)) );
	return (erfunc  - mu) - rat*(mu - erfunc/LogLambda(E,index)) ;
}

/** Debye length in field plasma
 * @return Debye length in cm
 */
float StopPow_LP::lDebye()
{
	float ret = 0; // temporary return value
	//iterate over all field particles:
	for(int i=0; i < num; i++)
	{
		ret += (4.*M_PI*nf[i]*pow(Zf[i]*e,2.0) / (kB*Tf[i]*keVtoK) );
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
	return c*sqrt(constant*Tf[index]/(mpc2*mf[index]));
}

/** x^{t/f} parameter from Li 1993
 * @param E test particle energy in MeV
 * @param index the field particle's index (for mf,Zf,Tf,nf arrays)
 * @return x^{t/f}
 */
float StopPow_LP::xtf(float E, int index)
{
	// test particle velocity:
	float vt = c*sqrt(2*E*1e3/(mt*mpc2));
	float vf = vtf(index, 2.); // sqrt(2kT/m)
	return pow( vt/vf ,2);
}

/* x^{t/f} parameter from Li 1993 for the collective effects term */
float StopPow_LP::xtf_collective(float E, int index)
{
	// test particle velocity:
	float vt = c*sqrt(2*E*1e3/(mt*mpc2));
	float vf = vtf(index, 1.); // sqrt(kT/m)
	return pow( vt/vf ,2);
}

/** Relative velocity between test particle and field particle
 * @param E test particle energy in MeV
 * @param index the field particle's index (for mf,Zf,Tf,nf arrays)
 * @return relative velocity in cm/s
 */
float StopPow_LP::u(float E, int index)
{
	float vt = c*sqrt(2.*E*1e3/(mt*mpc2));
	float vf = vtf(index, 8./M_PI); // sqrt(8kT/pi*m)
	
	// simple model:
	//return sqrt( (pow(vt,2) + pow(vf,2)) );
	
	// more exact:
	float u = (vf/2.)*exp(-4.*vt*vt/(M_PI*vf*vf))
		+ vt*(1.+M_PI*vf*vf/(8.*vt*vt))
		* erf(sqrt(4*vt*vt/(M_PI*vf*vf)));

	return u;
}

} // end namespace StopPow
