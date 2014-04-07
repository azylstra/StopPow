#include "StopPow_PartialIoniz.h"

namespace StopPow
{

// Default plasma model constructor:
StopPow_PartialIoniz::StopPow_PartialIoniz(double mt_in, double Zt_in, std::vector<double> & mf_in, std::vector<double> & Zf_in, std::vector<double> & Tf_in, std::vector<double> & nf_in, std::vector<double> & Zbar_in, double Te_in) throw(std::invalid_argument)
{
	// default mode is length:
	set_mode(MODE_LENGTH);
	
	// call helper methods:
	set_particle(mt_in, Zt_in);
	set_field(mf_in, Zf_in, Tf_in, nf_in, Zbar_in, Te_in);
}

// constructor using vector of arrays
StopPow_PartialIoniz::StopPow_PartialIoniz(double mt_in, double Zt_in, std::vector< std::array<double,5> > & field, double Te_in) throw(std::invalid_argument)
{
	// default mode is length:
	set_mode(MODE_LENGTH);
	
	// call helper method for test particle:
	set_particle(mt_in, Zt_in);
	set_field(field, Te_in);
}

// Destructor
StopPow_PartialIoniz::~StopPow_PartialIoniz()
{
	// nothing to do
}

// Method to set test particle info
void StopPow_PartialIoniz::set_particle(double mt_in, double Zt_in)
{
	// sanity check
	if( mt_in <= 0 || isnan(mt_in)
		|| Zt_in <= 0 || isnan(Zt_in) )
	{
		std::stringstream msg;
		// start constructing message, add info on mt and Zt:
		msg << "Values passed to StopPow_PartialIoniz::set_particle are bad: " 
		 << mt_in << "," << Zt_in << "," << std::endl;

		 // throw the exception:
		throw std::invalid_argument(msg.str());
	}

	// set class variables:
	mt = mt_in;
	Zt = Zt_in;
}

// Method to set field particle info
void StopPow_PartialIoniz::set_field(std::vector< std::array<double,5> > & field, double Te) throw(std::invalid_argument)
{
	// Convert field particle info:
	std::vector<double> mf, Zf, Tf, nf, Zbar;
	for( std::array<double,5> row : field)
	{
		mf.push_back(row[0]);
		Zf.push_back(row[1]);
		Tf.push_back(row[2]);
		nf.push_back(row[3]);
		Zbar.push_back(row[4]);
	}
	set_field(mf, Zf, Tf, nf, Zbar, Te);
}

// Method to set field particle info
void StopPow_PartialIoniz::set_field(std::vector<double> & mf_in, std::vector<double> & Zf_in, std::vector<double> & Tf_in, std::vector<double> & nf_in, std::vector<double> & Zbar_in, double Te_in)
{
	// infer size of the field particle arrays:
	num = mf_in.size();

	// sanity checking. 
	bool args_ok = true;
	// Make sure that all field particle arrays have same size
	if( Zf_in.size() != num 
		|| Tf_in.size() != num
		|| nf_in.size() != num 
		|| Zbar_in.size() != num )
	{
		args_ok = false;
	}

	// now do sanity checking on the field particle values:
	for(int i=0; i<num; i++)
	{
		args_ok = args_ok && mf_in[i] > 0;
		args_ok = args_ok && Tf_in[i] > 0;
		args_ok = args_ok && nf_in[i] > 0;
		args_ok = args_ok && Zbar_in[i] >= 0;
		args_ok = args_ok && Zbar_in[i] <= Zf_in[i];
	}

	// throw an exception if necessary:
	if( !args_ok )
	{
		std::stringstream msg;
		// start constructing message, add info on mt and Zt:
		msg << "Values passed to StopPow_PartialIoniz constructor are bad: " << std::endl;

		std::vector<double>::iterator it; // to iterate over field particles

		// add each element in mf:
		msg << "mf = ";
		for(it=mf_in.begin(); it<mf_in.end(); it++)
		 	msg << (*it) << ",";

		// add each element in Zf:
		msg << std::endl << "Zf = ";
		for(it=Zf_in.begin(); it<Zf_in.end(); it++)
		 	msg << (*it) << ",";

		// add each element in Tf:
		msg << std::endl << "Tf = ";
		for(it=Tf_in.begin(); it<Tf_in.end(); it++)
		 	msg << (*it) << ",";

		// add each element in nf:
		msg << std::endl << "nf = ";
		for(it=nf_in.begin(); it<nf_in.end(); it++)
		 	msg << (*it) << ",";

		// add each element in Zbar:
		msg << std::endl << "Zbar = ";
		for(it=Zbar_in.begin(); it<Zbar_in.end(); it++)
		 	msg << (*it) << ",";

		 // throw the exception:
		throw std::invalid_argument(msg.str());
	}
	// invoke copy constructor for vectors:
	mf = std::vector<double>(mf_in);
	Zf = std::vector<double>(Zf_in);
	Tf = std::vector<double>(Tf_in);
	nf = std::vector<double>(nf_in);
	Zbar = std::vector<double>(Zbar_in);

	// calculate the field particle mass density:
	rho = 0; // g/cm3
	// iterate over field particles:
	for(int i=0; i<num; i++)
	{
		rho += mf[i] * mp * nf[i];
	}

	// Calculate electron number density:
	ne = 0.;
	for(int i=0; i < num; i++)
		ne += Zbar[i] * nf[i]; // including ionization state
	Te = Te_in;
}

} // end of namespace