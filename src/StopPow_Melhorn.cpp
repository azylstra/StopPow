#include "StopPow_Melhorn.h"

namespace StopPow
{

const float StopPow_Melhorn::Emin = 0.1; /* Minimum energy for dE/dx calculations */
const float StopPow_Melhorn::Emax = 30; /* Maximum energy for dE/dx calculations */

// constructor
StopPow_Melhorn::StopPow_Melhorn(float mt_in, float Zt_in, std::vector<float> mf_in, std::vector<float> Zf_in, std::vector<float> Zbar_in, std::vector<float> nf_in, float Te_in) throw(std::invalid_argument)
{
	// default mode:
	set_mode(MODE_LENGTH);
	
	// infer size of the field particle arrays:
	num = mf_in.size();

	// sanity checking. 
	bool args_ok = true;
	// Make sure mt and Zt are positive,
	// and that all field particle arrays have same size
	if( mt_in <= 0 || Zt_in <= 0
		|| Zf_in.size() != num 
		|| Zbar_in.size() != num
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
		msg << "Values passed to StopPow_Melhorn constructor are bad: " 
		 << mt_in << "," << Zt_in << "," << std::endl;

		std::vector<float>::iterator it; // to iterate over field particles

		// add each element in mf:
		msg << "mf = ";
		for(it=mf_in.begin(); it<mf_in.end(); it++)
		 	msg << (*it) << ",";

		// add each element in Zf:
		msg << std::endl << "Zf = ";
		for(it=Zf_in.begin(); it<Zf_in.end(); it++)
		 	msg << (*it) << ",";

		// add each element in Zbar:
		msg << std::endl << "Zbar = ";
		for(it=Zbar_in.begin(); it<Zbar_in.end(); it++)
		 	msg << (*it) << ",";

		// add each element in nf:
		msg << std::endl << "nf = ";
		for(it=nf_in.begin(); it<nf_in.end(); it++)
		 	msg << (*it) << ",";

		 // throw the exception:
		throw std::invalid_argument(msg.str());
	}

	// set class variables:
	mt = mt_in;
	Zt = Zt_in;
	mf = mf_in;
	Zf = Zf_in;
	Zbar = Zbar_in;
	nf = nf_in;

	// and free electron stuff:
	Te = Te_in;
	// electron density:
	ne = 0;
	for(int i=0; i<num; i++)
	{
		ne += nf[i] * Zbar[i];
	}

	// calculate the field particle mass density:
	rho = 0; // g/cm3
	// iterate over field particles:
	for(int i=0; i<num; i++)
	{
		rho += mf[i] * mp * nf[i];
	}
	// set up Li-Petrasso for the free electrons and ions:
	std::vector<float> plasma_mf {{me/mp}};
	std::vector<float> plasma_Zf {{-1.}};
	std::vector<float> plasma_Tf {{Te}};
	std::vector<float> plasma_nf {{ne}};
	// add the ions:
	for(int i=0; i<num; i++)
	{
		plasma_mf.push_back( mf[i] );
		plasma_Zf.push_back( Zbar[i] );
		plasma_Tf.push_back( Te );
		plasma_nf.push_back( nf[i] );
	}
	PlasmaStop = new StopPow_LP(mt, Zt, plasma_mf, plasma_Zf, plasma_Tf, plasma_nf);

	// set the info string:
	model_type = "Melhorn";
	info = "";
}

// Destructor
StopPow_Melhorn::~StopPow_Melhorn()
{
	delete PlasmaStop;
}

/** Calculate the total stopping power
 * @param E the test particle energy in MeV
 * @return stopping power in units of MeV/um
 * @throws invalid_argument
 */
float StopPow_Melhorn::dEdx_MeV_um(float E) throw(std::invalid_argument)
{
	// sanity check:
	if( E < Emin || E > Emax )
	{
		std::stringstream msg;
		msg << "Energy passed to StopPow_Melhorn::dEdx is bad: " << E;
		throw std::invalid_argument(msg.str());
	}

	float cold = 0; // return value for cold component
	// components:
	float Bethe, LSS, nuc;

	//iterate over all field particles:
	for(int i=0; i < num; i++)
	{
		// call private helper functions:
		Bethe = dEdx_Bethe(E, i);
		LSS = dEdx_LSS(E, i);
		nuc = dEdx_nuc(E, i);
		cold += fmax(Bethe, LSS) + nuc;
	}

	// calculate the free electron contribution:
	float hot = PlasmaStop->dEdx_MeV_um(E);

	return (cold+hot);
}

/** Calculate the total stopping power
 * @param E the test particle energy in MeV
 * @return stopping power in units of MeV/(mg/cm2)
 * @throws invalid_argument
 */
float StopPow_Melhorn::dEdx_MeV_mgcm2(float E) throw(std::invalid_argument)
{
	return (dEdx_MeV_um(E)*1e4) / (rho*1e3);
}

/**
 * Get the minimum energy that can be used for dE/dx calculations
 * @return Emin in MeV
 */
float StopPow_Melhorn::get_Emin()
{
	return Emin;
}

/**
 * Get the maximum energy that can be used for dE/dx calculations
 * @return Emax in MeV
 */
float StopPow_Melhorn::get_Emax()
{
	return Emax;
}

// Stopping for low energy ions, from LSS theory
float StopPow_Melhorn::dEdx_LSS(float E, int index)
{
	// See Eq 3 of T. Melhorn, C. Appl. Phys. 52, 6522 (1981)

	float A = mf[index] / mt;
	float K = 0.0793 * pow(ZtEff(E), 2/3) * sqrt(Zf[index]) * pow(1+A, 1.5)
				/ ( pow(pow(ZtEff(E), 2/3) + pow(Zf[index], 2/3), 0.75) * sqrt(mf[index]) );
	float a = 4.683e-9 / sqrt( pow(ZtEff(E), 2/3) + pow(Zf[index], 2/3) ); // cm
	float El = (1 + A)*Zf[index]*ZtEff(E)*e*e/(A*a); // erg
	float RL = pow(1+A, 2) / (4*M_PI*A*nf[index]*a*a);
	float C_LSS = K * sqrt(El/1.602e-9) / (RL*1e4); // KeV^1/2 / um

	// convert. Sign is for convention in this program that dE/dx < 0
	C_LSS = -1 * C_LSS * sqrt(1e3); // MeV^1.2 / um
	// and return the stopping power:
	return C_LSS * sqrt(E);
}

// Nuclear stopping power
float StopPow_Melhorn::dEdx_nuc(float E, int index)
{
	// See Eq 4 of T. Melhorn, C. Appl. Phys. 52, 6522 (1981)
	float C = E / mt; // MeV/amu
	float Cn = 4.14e6 * pow(mt/(mt+mf[index]), 1.5) * sqrt(ZtEff(E)*Zf[index]/mf[index])
				/ pow( pow(ZtEff(E), 2/3) + pow(Zf[index], 2/3), 0.75);
	float Cnprime = mf[index]*mt/(mf[index]+mt) * (1./(ZtEff(E)*Zf[index])) / sqrt( pow(ZtEff(E), 2/3) + pow(Zf[index], 2/3) );
	// given in paper in terms of areal density:
	float dEdr = Cn * sqrt(C) * exp(-45.2*pow(Cnprime*C, 0.277));
	// return with conversion to MeV / um:
	return (dEdr *1e4) / (rho*1e3);
}

// Bethe stopping power, with Melhorn's adjustments
float StopPow_Melhorn::dEdx_Bethe(float E, int index)
{
	float Ekev = E * 1e3; // energy in keV for convenience
	float ret = 0; // return value

	// see Eq 1 of T. Melhorn, C. Appl. Phys. 52, 6522 (1981)
	float rho_i, LogLamda, vt, beta, gamma, prefac;
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
float StopPow_Melhorn::Ibar(float E, int index)
{
	// See Eq 7 of T. Melhorn, C. Appl. Phys. 52, 6522 (1981)

	// need to calculate Ibar(Z-Zbar), which is generally nonintegral.
	// calculate bounding indices:
	int i1 = floor( Zf[index] - Zbar[index] );
	int i2 =  ceil( Zf[index] - Zbar[index] );
	float di = Zf[index] - Zbar[index] - (float)i1;

	// sanity check:
	if( i1<0 || i2>=AtomicData::n )
	{
		std::stringstream msg;
		msg << "Out of range in StopPow_Melhorn::Ibar. Got i1, i2 = ";
		msg << i1 << "," << i2;
		throw std::invalid_argument(msg.str());
	}

	// check if we have an integer value:
	float Ibar;
	if(i1==i2)
	{
		Ibar = AtomicData::get_mean_ionization(i1);
	}
	else
	{
		// use linear interpolation to get the right value
		float Ibar1 = AtomicData::get_mean_ionization(i1);
		float Ibar2 = AtomicData::get_mean_ionization(i2);

		// there's no data in the atomic table for 0, so:
		if( i1 == 0 )
			Ibar1 = 0.;

		Ibar = Ibar1 + di*(Ibar2-Ibar1)/( (float)(i2-i1) );
	}
	float ret_eV = pow(Zf[index], 2) * Ibar / pow( Zf[index] - Zbar[index] , 2);
	return ret_eV * 1.602e-12;
}

// effective projectile charge
float StopPow_Melhorn::ZtEff(float E)
{
	// See Eq 6 of T. Melhorn, C. Appl. Phys. 52, 6522 (1981)
	// test particle velocity
	float beta = sqrt(2e3*E/(mt*mpc2)); // normalized to c

	return Zt * (1 - 1.034*exp(-137.04*beta/pow(Zt, 0.69)));
}

// Calculate shell correction term in log lambda 
float StopPow_Melhorn::shell_term(float Zf, float E)
{
	int Z = (int)Zf;

	// sanity check:
	if( Z < 1 || Z > AtomicData::n || E < Emin || E > Emax)
		return 0;

	// get coefficients:
	std::array<float,5> coeff = AtomicData::get_shell_coeff(Z);

	float shell = 0;
	// for each coefficient:
	for(int i=0; i < coeff.size(); i++)
		shell += coeff[i] * pow( log(1e3*E/mt) , i); // convert E to keV

	return shell;
}

} // end namespace StopPow
