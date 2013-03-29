/**
 * @brief Calculate Bethe-Bloch stopping power.
 * 
 * Implement a stopping-power calculator for arbitrary cold matter, using
 * the simple Bethe-Bloch theory.
 *
 * @class StopPow_BetheBloch
 * @author Alex Zylstra
 * @date 2013/03/25
 */


 #include "StopPow_BetheBloch.h"

const float StopPow_BetheBloch::IbarData[] = {19.0f,21.0f,16.0f,15.0f,15.0f,13.0f};

 /** Initialize the Bethe-Bloch calculator.
 * @param mt the test particle mass in AMU
 * @param Zt the test particle in charge (units of e)
 * @param mf vector containing ordered field particle masses in AMU
 * @param Zt vector containing ordered field particle charges in units of e
 * @param nf vector containing ordered field particle densities in units of 1/cm3
 * @throws invalid_argument
*/
StopPow_BetheBloch::StopPow_BetheBloch(float mt_in, float Zt_in, vector<float> mf_in, vector<float> Zf_in, vector<float> nf_in)
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
		stringstream msg;
		msg << "Values passed to StopPow_LP constructor are bad: " 
		 << mt_in << "," << Zt_in << ","; // << mf_in << "," << Zf_in << "," << "," << nf_in;
		throw new invalid_argument(msg.str());
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
}

/** Calculate the total stopping power
 * @param E the test particle energy in MeV
 * @return stopping power in units of MeV/um
 * @throws invalid_argument
*/
float StopPow_BetheBloch::dEdx_MeV_um(float E)
{
	float Ekev = E * 1e3; // energy in keV for convenience

	float ret = 0;

	// iterate over field particles:
	float rho, LogLamda, vt, beta, gamma, prefac;
	for(int i=0; i < num; i++)
	{
		rho = nf[i] * mf[i] / Na; // mass density in g/cm3
		LogLamda = 0.0; // initialize

		vt = c*sqrt(2.0*Ekev/mpc2); // test particle velocity
		beta = vt/c; // normalized to c
		gamma = 1.0/sqrt(1-pow(beta,2)); // relativistic gamma factor

		prefac = 4.0*M_PI*Na*rho*pow(Zt*e*e,2)*Zf[i] / (me*c*c*beta*beta*mf[i]);
		LogLamda += log(2.0*me*pow(c*beta*gamma,2)/Ibar(Zf[i]));
		LogLamda -= pow(beta,2);
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
float StopPow_BetheBloch::dEdx_MeV_mgcm2(float E)
{
	return (dEdx_MeV_um(E)*1e4) / (rho*1e3);
}

/** Effecive ionization potential as a function of Z.
 * @param Zf field particle charge in units of e
 * @return Ibar in erg
 */
float StopPow_BetheBloch::Ibar(float Zf)
{
	if( Zf <= sizeof(IbarData)/sizeof(float) )
		return (IbarData[(int)Zf - 1])*Zf*1.602e-12;
	return 0;
}