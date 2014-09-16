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

/** test class for StopPow library
 * @author Alex Zylstra
 * @date 2014/04/07
 */


#include <stdio.h>
#include <dirent.h>

#include <iostream>
#include <vector>
#include <array>
#include <stdexcept>
#include <limits>

#include "StopPow.h"
#include "StopPow_Plasma.h"
#include "StopPow_SRIM.h"
#include "StopPow_LP.h"
#include "StopPow_AZ.h"
#include "StopPow_BetheBloch.h"
#include "StopPow_Grabowski.h"
#include "StopPow_Zimmerman.h"
#include "StopPow_BPS.h"
#include "StopPow_Constants.h"

#include <gsl/gsl_sf_fermi_dirac.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_sf_dawson.h>

// Read in plasma conditions and test dE/dx data from a file
void read_partial_ioniz_file(std::string fname, double & mt, double & Zt, std::vector< std::array<double,5> > & field_data, double & Te, std::vector< std::array<double,2> > & test_data)
{
	// Clear input vectors:
	field_data.clear();
	test_data.clear();

	std::ifstream file( fname.c_str() );
	
	if( file.is_open() )
	{
		std::string line;
		while( !file.eof() )
		{
			getline(file,line);
			// ignore name and comment lines:
			if( line.find("#") == std::string::npos && line.size() > 0)
			{
				// use stringstream for line:
				std::stringstream ss(line);
				std::string val;

				// Case depends on flag
				if(line[0] == 'f') // field particle info on this line
				{
					// discard header
					getline(ss,val,',');
					// read data:
					int i=0;
					std::array<double,5> new_line;
					while(getline(ss,val,',')) {
						new_line[i] = atof(val.c_str());
						i++;
					}
					field_data.push_back(new_line);
				}
				else if(line[0] == 't') // test particle info
				{
					// discard header
					getline(ss,val,',');

					// read mt:
					getline(ss,val,',');
					mt = atof(val.c_str());
					
					// read Zt:
					getline(ss,val,',');
					Zt = atof(val.c_str());
				}
				else if(line[0] == 'T' && line[1] == 'e') // electron temp
				{
					// discard header
					getline(ss,val,',');

					// read Te:
					getline(ss,val,',');
					Te = atof(val.c_str());
				}
				else // data on this line
				{
					std::array<double,2> new_line;

					// read E:
					getline(ss,val,',');
					new_line[0] = atof(val.c_str());
					
					// read dE/dx:
					getline(ss,val,',');
					new_line[1] = atof(val.c_str());

					test_data.push_back(new_line);
				}
			}
		}
	}
}


// Read in plasma conditions and test dE/dx data from a file
void read_plasma_file(std::string fname, double & mt, double & Zt, std::vector< std::array<double,4> > & field_data, std::vector< std::array<double,2> > & test_data)
{
	// Clear input vectors:
	field_data.clear();
	test_data.clear();

	std::ifstream file( fname.c_str() );
	
	if( file.is_open() )
	{
		std::string line;
		while( !file.eof() )
		{
			getline(file,line);
			// ignore name and comment lines:
			if( line.find("#") == std::string::npos && line.size() > 0)
			{
				// use stringstream for line:
				std::stringstream ss(line);
				std::string val;

				// Case depends on flag
				if(line[0] == 'f') // field particle info on this line
				{
					// discard header
					getline(ss,val,',');
					// read data:
					int i=0;
					std::array<double,4> new_line;
					while(getline(ss,val,',')) {
						new_line[i] = atof(val.c_str());
						i++;
					}
					field_data.push_back(new_line);
				}
				else if(line[0] == 't') // test particle info
				{
					// discard header
					getline(ss,val,',');

					// read mt:
					getline(ss,val,',');
					mt = atof(val.c_str());
					
					// read Zt:
					getline(ss,val,',');
					Zt = atof(val.c_str());
				}
				else // data on this line
				{
					std::array<double,2> new_line;

					// read E:
					getline(ss,val,',');
					new_line[0] = atof(val.c_str());
					
					// read dE/dx:
					getline(ss,val,',');
					new_line[1] = atof(val.c_str());

					test_data.push_back(new_line);
				}
			}
		}
	}
}

bool run_test(StopPow::StopPow* model, std::vector< std::array<double,2> > data, double tol, bool verbose)
{
	if(verbose)
	{
		std::cout << "E (MeV) , dE/dx (ref) , dE/dx (calc)" << std::endl;
		verbose = true;
	}

	// ---------------------------------------
	// 		test dE/dx
	// ---------------------------------------
	double result, delta;
	bool pass = true;
	for( std::array<double,2> row : data )
	{
		result = model->dEdx_MeV_um(row[0]);
		delta = fabs(result - row[1]) / fabs(row[1]);
		pass = pass && (delta <= tol);
		if(verbose)
			std::cout << row[0] << " , " << row[1] << " , " << result << " -> " << delta << " , " << tol << " , " << pass << std::endl;
	}

	return pass;
}

// Helper function for testing a model, loads all tests from files and runs them
// Templated to allow for creation of all models with some splitting (below method) for plasma vs partial ioniz models
template<class T>
bool test_plasma_model(std::string dir_name, std::string name, float tol, bool verbose)
{
	std::cout << "Testing " + name + " model..." << std::endl;
	bool pass = true;
	int n = 0;
	// Load all BPS test data directory
	DIR *model_dir = opendir(dir_name.c_str());
	if(model_dir) // dir is open
	{
		// loop over all files:
		dirent* result = readdir(model_dir);
		while( (result=readdir(model_dir)) != NULL )
		{
			// try to load if an actual file:
			if(std::string(result->d_name).find("csv") != std::string::npos)
			{
				// construct relative file path/name:
				std::string fname(dir_name);
				fname.append("/");
				fname.append(result->d_name);
				
				// test data points:
				std::vector< std::array<double,2> > test_data;
				// storage for data from file:
				double mt, Zt;
				
				// necessary dimension for field data depends on model type:
				// Call appropriate helper function to read data
				std::vector< std::array<double,4> > field_data;
				read_plasma_file(fname, mt, Zt, field_data, test_data);
				StopPow::StopPow * s = new T(mt, Zt, field_data);

				// create model and test:
				bool test = run_test(s, test_data, tol, verbose);
				std::cout << fname << ": " << (test ? "pass" : "FAIL") << std::endl;
				pass = pass && test;
				n++;

				// free memory:
				delete s;
			}
		}
	}
	std::cout << n << " " + name + " model(s) tested: " << (pass ? "pass" : "FAIL") << std::endl << std::endl;

	return pass;
}

// Helper function for testing a partially ionized model, loads all tests from files and runs them
template<class T>
bool test_partial_ioniz_model(std::string dir_name, std::string name, float tol, bool verbose)
{
	std::cout << "Testing " + name + " model..." << std::endl;
	bool pass = true;
	int n = 0;
	// Load all BPS test data directory
	DIR *model_dir = opendir(dir_name.c_str());
	if(model_dir) // dir is open
	{
		// loop over all files:
		dirent* result = readdir(model_dir);
		while( (result=readdir(model_dir)) != NULL )
		{
			// try to load if an actual file:
			if(std::string(result->d_name).find("csv") != std::string::npos)
			{
				// construct relative file path/name:
				std::string fname(dir_name);
				fname.append("/");
				fname.append(result->d_name);
				
				// test data points:
				std::vector< std::array<double,2> > test_data;
				// storage for data from file:
				double mt, Zt, Te;
				
				// necessary dimension for field data depends on model type:
				// Call appropriate helper function to read data
				std::vector< std::array<double,5> > field_data;
				read_partial_ioniz_file(fname, mt, Zt, field_data, Te, test_data);
				StopPow::StopPow * s = new T(mt, Zt, field_data, Te);

				// create model and test:
				bool test = run_test(s, test_data, tol, verbose);
				std::cout << fname << ": " << (test ? "pass" : "FAIL") << std::endl;
				pass = pass && test;
				n++;

				// free memory:
				delete s;
			}
		}
	}
	std::cout << n << " " + name + " model(s) tested: " << (pass ? "pass" : "FAIL") << std::endl << std::endl;

	return pass;
}

// Helper function for testing a cold models against SRIM
template<class T>
bool test_cold_model(std::string srim_name, std::string name, double mt, double Zt, std::vector<double> mf, std::vector<double> Zf, std::vector<double> nf, float tol, bool verbose)
{
	StopPow::StopPow * SRIM = new StopPow::StopPow_SRIM(srim_name);
	StopPow::StopPow * s = new T(mt, Zt, mf, Zf, nf);
	bool pass = true;
	for(double E=1.0; E<=20.; E+=0.05)
	{
		double dEdx = s->dEdx_MeV_um(E);
		double dEdx_SRIM = SRIM->dEdx_MeV_um(E);
		double delta = fabs(dEdx - dEdx_SRIM) / fabs(dEdx_SRIM);
		pass = pass && (delta < 3e-2);
		if(verbose)
			std::cout << E << " , " << dEdx_SRIM << " , " << dEdx << " -> " << pass << std::endl;
	}
	std::cout << name + " model tested: " << (pass ? "pass" : "FAIL") << std::endl;

	return pass;
}

// Andersen-Ziegler model has unique constructor, requires its own helper function
bool test_AZ_model(std::string srim_name, std::string name, double Z, float tol, bool verbose)
{
	StopPow::StopPow * SRIM = new StopPow::StopPow_SRIM(srim_name);
	StopPow::StopPow * s = new StopPow::StopPow_AZ(Z);
	bool pass = true;
	for(double E=1.0; E<=20.; E+=0.05)
	{
		double dEdx = s->dEdx_MeV_um(E);
		double dEdx_SRIM = SRIM->dEdx_MeV_um(E);
		double delta = fabs(dEdx - dEdx_SRIM) / fabs(dEdx_SRIM);
		pass = pass && (delta < 3e-2);
		if(verbose)
			std::cout << E << " , " << dEdx_SRIM << " , " << dEdx << " -> " << pass << std::endl;
	}
	std::cout << name + " model tested: " << (pass ? "pass" : "FAIL") << std::endl;

	return pass;
}

int main(int argc, char* argv [])
{
	// check for verbosity flag:
	bool verbose = false;
	if( argc >= 2 )
	{
		for(int i=1; i < argc; i++)
		{
			std::string flag(argv[i]);
			if(flag == "--verbose")
			{
				verbose = true;
			}
		}
	}

	// Do some output
	std::cout << "========== Test Suite 0 ==========" << std::endl;
	std::cout << "  Testing computational fidelity  " << std::endl << std::endl;
	bool all_pass = true;


	// ---------------------------------------
	// 	   Test cold models against SRIM
	// ---------------------------------------
	std::cout << "Testing Bethe-Bloch model against SRIM..." << std::endl;
	all_pass &= test_cold_model<StopPow::StopPow_BetheBloch>("SRIM/Hydrogen in Beryllium.txt", "Bethe-Bloch (Be)", 1, 1, {9.012}, {4.0}, {1.235e23}, 3e-2, verbose);
	all_pass &= test_cold_model<StopPow::StopPow_BetheBloch>("SRIM/Hydrogen in Aluminum.txt", "Bethe-Bloch (Al)", 1, 1, {26.98}, {13.0}, {6.03e22}, 3e-2, verbose);
	all_pass &= test_cold_model<StopPow::StopPow_BetheBloch>("SRIM/Hydrogen in Tantalum.txt", "Bethe-Bloch (Ta)", 1, 1, {180.95}, {73.0}, {5.525e22}, 3e-2, verbose);
	std::cout << std::endl;
	
	std::cout << "Testing Andersen-Ziegler model against SRIM..." << std::endl;
	all_pass &= test_AZ_model("SRIM/Hydrogen in Beryllium.txt", "A-Z (Be)", 4., 3e-2, verbose);
	all_pass &= test_AZ_model("SRIM/Hydrogen in Aluminum.txt", "A-Z (Al)", 13., 3e-2, verbose);
	all_pass &= test_AZ_model("SRIM/Hydrogen in Tantalum.txt", "A-Z (Ta)", 73., 3e-2, verbose);
	std::cout << std::endl;

	// ---------------------------------------
	//       Set up and test various models
	// ---------------------------------------
	all_pass &= test_plasma_model<StopPow::StopPow_LP>("test0/Li-Petrasso", "Li-Petrasso", 2e-2, verbose);
	all_pass &= test_plasma_model<StopPow::StopPow_Grabowski>("test0/Grabowski", "Grabowski", 1e-2, verbose);
	all_pass &= test_partial_ioniz_model<StopPow::StopPow_Zimmerman>("test0/Zimmerman", "Zimmerman", 3e-2, verbose);
	all_pass &= test_plasma_model<StopPow::StopPow_BPS>("test0/BPS", "BPS", 7e-2, verbose); // BPS limits looser because DataThief is not particularly accurate

	std::cout << "RESULT: " << (all_pass ? "PASS" : "FAIL") << std::endl;

	return (all_pass? 0 : 1);
	
}