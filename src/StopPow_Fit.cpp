#include "StopPow_Fit.h"

namespace StopPow
{

const int StopPow_Fit::MODE_ZIMMERMAN = 0; 
const int StopPow_Fit::MODE_LP = 1; 
const int StopPow_Fit::MODE_BPS = 2;

// Constructors use those in PartialIoniz
StopPow_Fit::StopPow_Fit(double mt_in, double Zt_in, std::vector<double> & mf_in, std::vector<double> & Zf_in, std::vector<double> & Tf_in, std::vector<double> & nf_in, std::vector<double> & Zbar_in, double Te_in) throw(std::invalid_argument)
	: StopPow_PartialIoniz::StopPow_PartialIoniz(mt_in, Zt_in, mf_in, Zf_in, Tf_in, nf_in, Zbar_in, Te_in) 
{
	init();
}

StopPow_Fit::StopPow_Fit(double mt_in, double Zt_in, std::vector< std::array<double,5> > & field_in, double Te_in) throw(std::invalid_argument)
	: StopPow_PartialIoniz::StopPow_PartialIoniz(mt_in, Zt_in, field_in, Te_in) 
{
	init();
}

// Destructor
StopPow_Fit::~StopPow_Fit()
{
	delete z;
	delete fe;
}

void StopPow_Fit::init()
{
	// Set up the Zimmerman model:
	z = new StopPow_Zimmerman(mt, Zt, mf, Zf, Tf, nf, Zbar, Te);
	fe = z;
}

// Calculate stopping power
double StopPow_Fit::dEdx_MeV_um(double E) throw(std::invalid_argument)
{
	// three components:
	double dEdx_i, dEdx_be, dEdx_fe;

	dEdx_i = z->dEdx_ion(E);
	dEdx_be = be_factor * z->dEdx_bound_electron(E);
	// Zimmerman model uses a single StopPow object, necessitates calling
	// the free electron function
	if( fe == z )
		dEdx_fe = fe_factor * ((StopPow_Zimmerman*)fe)->dEdx_free_electron(E);
	else if( fe_model == MODE_BPS )
	{
		// Scale only the quantum correction:
		double dEdx_short, dEdx_long, dEdx_quantum;
		dEdx_short = ((StopPow_BPS*)fe)->dEdx_short(E);
		dEdx_long = ((StopPow_BPS*)fe)->dEdx_long(E);
		dEdx_quantum = ((StopPow_BPS*)fe)->dEdx_quantum(E);
		dEdx_fe = dEdx_short + dEdx_long + fe_factor * dEdx_quantum;
	}
	else // assume LP
		dEdx_fe = fe_factor * fe->dEdx(E);

	return (dEdx_i + dEdx_be + dEdx_fe);
}

// Calculate stoppign power
double StopPow_Fit::dEdx_MeV_mgcm2(double E) throw(std::invalid_argument)
{
	return (dEdx_MeV_um(E)*1e4) / (rho*1e3);
}

// Get energy limits, set by most restrictive of two sub-models:
double StopPow_Fit::get_Emin()
{
	return fmax(z->get_Emin(), fe->get_Emin());
}
double StopPow_Fit::get_Emax()
{
	return fmin(z->get_Emax(), fe->get_Emax());
}

// Normalize the bound-electron stopping to a reference case at a given proton energy.
void StopPow_Fit::normalize_bound_e(StopPow * ref, double Ep)
{
	// Need to create a dummy zero-ionization version of the Zimmerman model being used:
	std::vector<double> Zbar_0(Zbar);
	for(int i=0; i<Zbar_0.size(); i++)
		Zbar_0[i] = 0.0; // no ionization
	StopPow_Zimmerman * z2 = new StopPow_Zimmerman(mt, Zt, mf, Zf, Tf, nf, Zbar_0, Te);
	be_factor = ref->dEdx(Ep) / z2->dEdx(Ep);
	delete z2;
}

// Choose the free-electron stopping model.
void StopPow_Fit::choose_model(int new_model) throw(std::invalid_argument)
{
	std::vector<double> mf_fe {me/amu};
	std::vector<double> Zf_fe {-1};
	std::vector<double> Tf_fe {Te};
	std::vector<double> nf_fe {ne};

	switch(new_model)
	{
		case MODE_ZIMMERMAN:
			if( fe != z )
				delete fe;
			fe = z;
			break;
		case MODE_LP:
			if( fe != z )
				delete fe;
			fe = new StopPow_LP(mt, Zt, mf_fe, Zf_fe, Tf_fe, nf_fe);
			break;
		case MODE_BPS:
			if( fe != z )
				delete fe;
			fe = new StopPow_BPS(mt, Zt, mf_fe, Zf_fe, Tf_fe, nf_fe);
			break;
		default:
			throw std::invalid_argument("Model choice passed to StopPow_Fit::choose_model is invalid");
	}
	fe_model = new_model;
}

// Adjust the stopping, to be used for fitting
void StopPow_Fit::set_factor(double factor)
{
	fe_factor = factor;
}

// Get the current adjustment factor
double StopPow_Fit::get_factor()
{
	return fe_factor;
}

} // end of namespace