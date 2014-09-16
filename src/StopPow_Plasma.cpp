// StopPow - a charged-particle stopping power library
// Copyright (C) 2014  Massachusetts Institute of Technology / Alex Zylstra

// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along
// with this program; if not, write to the Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

#include "StopPow_Plasma.h"

namespace StopPow
{

// Constructor with individual vectors for field particles, electrons specified manually
StopPow_Plasma::StopPow_Plasma(double mt_in, double Zt_in, std::vector<double> & mf_in, std::vector<double> & Zf_in, std::vector<double> & Tf_in, std::vector<double> & nf_in) throw(std::invalid_argument)
{
	// default mode is length:
	set_mode(MODE_LENGTH);
	
	// call helper methods:
	set_particle(mt_in, Zt_in);
	set_field(mf_in, Zf_in, Tf_in, nf_in);
}

// Constructor with 2-D field info, electrons specified manually
StopPow_Plasma::StopPow_Plasma(double mt_in, double Zt_in, std::vector< std::array<double,4> > & field) throw(std::invalid_argument)
{
	// default mode is length:
	set_mode(MODE_LENGTH);
	
	// call helper methods:
	set_particle(mt_in, Zt_in);
	set_field(field);
}

// Constructor with individual vectors for field particles, electrons specified quasi-automatically
StopPow_Plasma::StopPow_Plasma(double mt_in, double Zt_in, std::vector<double> & mf_in, std::vector<double> & Zf_in, std::vector<double> & Tf_in, std::vector<double> & nf_in, double Te_in) throw(std::invalid_argument)
{
	// default mode is length:
	set_mode(MODE_LENGTH);
	
	// call helper methods:
	set_particle(mt_in, Zt_in);
	set_field(mf_in, Zf_in, Tf_in, nf_in, Te_in);
}

// Constructor with 2-D field info, electrons specified quasi-automatically
StopPow_Plasma::StopPow_Plasma(double mt_in, double Zt_in, std::vector< std::array<double,4> > & field, double Te_in) throw(std::invalid_argument)
{
	// default mode is length:
	set_mode(MODE_LENGTH);
	
	// call helper methods:
	set_particle(mt_in, Zt_in);
	set_field(field, Te_in);
}

// Destructor
StopPow_Plasma::~StopPow_Plasma()
{
	// nothing to do
}

// Get stopping power due only to electrons
double StopPow_Plasma::dEdx_plasma_electrons(double E) throw(std::invalid_argument)
{
	// loop over all field particles
	for(int i=0; i < num; i++)
	{
		// electrons identified by mass:
		if( approx(mf[i], me/mp, 1e-2) )
			return dEdx_field(E,i);
	}
	return 0; // no electrons!
}

// Get stopping power due only to ions
double StopPow_Plasma::dEdx_plasma_ions(double E) throw(std::invalid_argument)
{
	double ret = 0;
	// loop over all field particles
	for(int i=0; i < num; i++)
	{
		// electrons excluded by mass:
		if( !approx(mf[i], me/mp, 1e-2) )
			ret += dEdx_field(E,i);
	}
	return ret;
}

// Method to set test particle info
void StopPow_Plasma::set_particle(double mt_in, double Zt_in) throw(std::invalid_argument)
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
void StopPow_Plasma::set_field(std::vector<double> & mf_in, std::vector<double> & Zf_in, std::vector<double> & Tf_in, std::vector<double> & nf_in) throw(std::invalid_argument)
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

		 // throw the exception:
		throw std::invalid_argument(msg.str());
	}
	// invoke copy constructor for vectors:
	mf = std::vector<double>(mf_in);
	Zf = std::vector<double>(Zf_in);
	Tf = std::vector<double>(Tf_in);
	nf = std::vector<double>(nf_in);

	// calculate the field particle mass density:
	rho = 0; // g/cm3
	// iterate over field particles:
	for(int i=0; i<num; i++)
	{
		rho += mf[i] * mp * nf[i];
	}

	on_field_change();
}

// Method to set field particle info with quasi-automatic electrons
void StopPow_Plasma::set_field(std::vector<double> & mf_in, std::vector<double> & Zf_in, std::vector<double> & Tf_in, std::vector<double> & nf_in, double Te) throw(std::invalid_argument)
{
	// call other version first:
	set_field(mf_in, Zf_in, Tf_in, nf_in);

	// Calculate electron number density:
	double ne = 0.;
	for(int i=0; i < num; i++)
		ne += Zf_in[i] * nf_in[i]; // fully ionized
	// Add to field particle vectors:
	mf.push_back(me/amu);
	Zf.push_back(-1.);
	Tf.push_back(Te);
	nf.push_back(ne);
	num++;

	on_field_change();
}

void StopPow_Plasma::set_field(std::vector< std::array<double,4> > & field) throw(std::invalid_argument)
{
	// Convert field particle info:
	std::vector<double> mf, Zf, Tf, nf;
	for( std::array<double,4> row : field)
	{
		mf.push_back(row[0]);
		Zf.push_back(row[1]);
		Tf.push_back(row[2]);
		nf.push_back(row[3]);
	}
	set_field(mf, Zf, Tf, nf);
}

void StopPow_Plasma::set_field(std::vector< std::array<double,4> > & field, double Te) throw(std::invalid_argument)
{
	// Convert field particle info:
	std::vector<double> mf, Zf, Tf, nf;
	for( std::array<double,4> row : field)
	{
		mf.push_back(row[0]);
		Zf.push_back(row[1]);
		Tf.push_back(row[2]);
		nf.push_back(row[3]);
	}
	set_field(mf, Zf, Tf, nf, Te);
}

void StopPow_Plasma::on_field_change(){}


} // end of namespace