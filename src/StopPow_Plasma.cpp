#include "StopPow_Plasma.h"

namespace StopPow
{

// Constructor with individual vectors for field particles, electrons specified manually
StopPow_Plasma::StopPow_Plasma(float mt_in, float Zt_in, std::vector<float> & mf_in, std::vector<float> & Zf_in, std::vector<float> & Tf_in, std::vector<float> & nf_in) throw(std::invalid_argument)
{
	// default mode is length:
	set_mode(MODE_LENGTH);
	
	// call helper methods:
	set_particle(mt_in, Zt_in);
	set_field(mf_in, Zf_in, Tf_in, nf_in);
}

// Constructor with 2-D field info, electrons specified manually
StopPow_Plasma::StopPow_Plasma(float mt_in, float Zt_in, std::vector< std::array<float,4> > & field) throw(std::invalid_argument)
{
	// default mode is length:
	set_mode(MODE_LENGTH);
	
	// call helper methods:
	set_particle(mt_in, Zt_in);
	set_field(field);
}

// Constructor with individual vectors for field particles, electrons specified quasi-automatically
StopPow_Plasma::StopPow_Plasma(float mt_in, float Zt_in, std::vector<float> & mf_in, std::vector<float> & Zf_in, std::vector<float> & Tf_in, std::vector<float> & nf_in, float Te) throw(std::invalid_argument)
{
	// default mode is length:
	set_mode(MODE_LENGTH);
	
	// call helper methods:
	set_particle(mt_in, Zt_in);
	set_field(mf_in, Zf_in, Tf_in, nf_in, Te);
}

// Constructor with 2-D field info, electrons specified quasi-automatically
StopPow_Plasma::StopPow_Plasma(float mt_in, float Zt_in, std::vector< std::array<float,4> > & field, float Te) throw(std::invalid_argument)
{
	// default mode is length:
	set_mode(MODE_LENGTH);
	
	// call helper methods:
	set_particle(mt_in, Zt_in);
	set_field(field, Te);
}

// Destructor
StopPow_Plasma::~StopPow_Plasma()
{
	// nothing to do
}

// Method to set test particle info
void StopPow_Plasma::set_particle(float mt_in, float Zt_in) throw(std::invalid_argument)
{
	// sanity check
	if( mt_in <= 0 || isnan(mt_in)
		|| Zt_in <= 0 || isnan(Zt_in) )
	{
		std::stringstream msg;
		// start constructing message, add info on mt and Zt:
		msg << "Values passed to StopPow_Plasma::set_particle are bad: " 
		 << mt_in << "," << Zt_in << "," << std::endl;

		 // throw the exception:
		throw std::invalid_argument(msg.str());
	}

	// set class variables:
	mt = mt_in;
	Zt = Zt_in;
}

// Method to set field particle info
void StopPow_Plasma::set_field(std::vector<float> & mf_in, std::vector<float> & Zf_in, std::vector<float> & Tf_in, std::vector<float> & nf_in) throw(std::invalid_argument)
{
	// infer size of the field particle arrays:
	num = mf_in.size();

	// sanity checking. 
	bool args_ok = true;
	// Make sure that all field particle arrays have same size
	if( Zf_in.size() != num 
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
		msg << "Values passed to StopPow_Plasma constructor are bad: " << std::endl;

		std::vector<float>::iterator it; // to iterate over field particles

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

		 // throw the exception:
		throw std::invalid_argument(msg.str());
	}
	// invoke copy constructor for vectors:
	mf = std::vector<float>(mf_in);
	Zf = std::vector<float>(Zf_in);
	Tf = std::vector<float>(Tf_in);
	nf = std::vector<float>(nf_in);

	// calculate the field particle mass density:
	rho = 0; // g/cm3
	// iterate over field particles:
	for(int i=0; i<num; i++)
	{
		rho += mf[i] * mp * nf[i];
	}
}

// Method to set field particle info with quasi-automatic electrons
void StopPow_Plasma::set_field(std::vector<float> & mf_in, std::vector<float> & Zf_in, std::vector<float> & Tf_in, std::vector<float> & nf_in, float Te) throw(std::invalid_argument)
{
	// call other version first:
	set_field(mf_in, Zf_in, Tf_in, nf_in);

	// Calculate electron number density:
	float ne = 0.;
	for(int i=0; i < num; i++)
		ne += Zf[i] * nf[i]; // fully ionized
	// Add to field particle vectors:
	mf.push_back(me/amu);
	Zf.push_back(-1.);
	Tf.push_back(Te);
	nf.push_back(ne);
}

void StopPow_Plasma::set_field(std::vector< std::array<float,4> > & field) throw(std::invalid_argument)
{
	// Convert field particle info:
	std::vector<float> mf, Zf, Tf, nf;
	for( std::array<float,4> row : field)
	{
		mf.push_back(row[0]);
		Zf.push_back(row[1]);
		Tf.push_back(row[2]);
		nf.push_back(row[3]);
	}
	set_field(mf, Zf, Tf, nf);
}

void StopPow_Plasma::set_field(std::vector< std::array<float,4> > & field, float Te) throw(std::invalid_argument)
{
	// Convert field particle info:
	std::vector<float> mf, Zf, Tf, nf;
	for( std::array<float,4> row : field)
	{
		mf.push_back(row[0]);
		Zf.push_back(row[1]);
		Tf.push_back(row[2]);
		nf.push_back(row[3]);
	}
	set_field(mf, Zf, Tf, nf, Te);
}

} // end of namespace