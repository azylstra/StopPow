/** test class for StopPow library
 * @author Alex Zylstra
 * @date 2014/04/02
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
#include "StopPow_Constants.h"

#include <gsl/gsl_sf_fermi_dirac.h>


// Read in plasma conditions and test dE/dx data from a file
void read_partial_ioniz_file(std::string fname, float & mt, float & Zt, std::vector< std::array<float,5> > & field_data, float & Te, std::vector< std::array<float,2> > & test_data)
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
					std::array<float,5> new_line;
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
					std::array<float,2> new_line;

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
void read_plasma_file(std::string fname, float & mt, float & Zt, std::vector< std::array<float,4> > & field_data, std::vector< std::array<float,2> > & test_data)
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
					std::array<float,4> new_line;
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
					std::array<float,2> new_line;

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
bool run_test(StopPow::StopPow* model, std::vector< std::array<float,2> > data, float tol, bool verbose)
{
	if(verbose)
	{
		std::cout << "E (MeV) , dE/dx (ref) , dE/dx (calc)" << std::endl;
		verbose = true;
	}

	// ---------------------------------------
	// 		test dE/dx
	// ---------------------------------------
	float result, delta;
	bool pass = true;
	for( std::array<float,2> row : data )
	{
		result = model->dEdx_MeV_um(row[0]);
		delta = fabs(result - row[1]) / fabs(row[1]);
		pass = pass && (delta <= tol);
		if(verbose)
			std::cout << row[0] << " , " << row[1] << " , " << result << " -> " << delta << " , " << tol << " , " << pass << std::endl;
	}

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

	// create a variety of stopping power models to test
	std::vector<StopPow::StopPow*> models;

	// ---------------------------------------
	// 		Set up and test Li-Petrasso
	// ---------------------------------------
	std::cout << "Testing Li-Petrasso model..." << std::endl;
	bool pass = true;
	int n = 0;
	// Load all Li-Petrasso test data directory
	std::string dir_name("test0/Li-Petrasso");
	DIR *LP_dir = opendir(dir_name.c_str());
	if(LP_dir) // dir is open
	{
		// loop over all files:
		dirent* result = readdir(LP_dir);
		while( (result=readdir(LP_dir)) != NULL )
		{
			// try to load if an actual file:
			if(std::string(result->d_name) != "..")
			{
				// construct relative file path/name:
				std::string fname(dir_name);
				fname.append("/");
				fname.append(result->d_name);

				// storage for data:
				float mt, Zt;
				std::vector< std::array<float,4> > field_data;
				std::vector< std::array<float,2> > test_data;

				// Call helper function to read data
				read_plasma_file(fname, mt, Zt, field_data, test_data);

				// create model and test:
				StopPow::StopPow* s = new StopPow::StopPow_LP(mt, Zt, field_data);
				bool test = run_test(s, test_data, 2e-2, verbose);
				std::cout << fname << ": " << (test ? "pass" : "FAIL") << std::endl;
				pass = pass && test;
				n++;

				// free memory:
				delete s;
			}
		}
	}
	std::cout << n << " Li-Petrasso model(s) tested: " << (pass ? "pass" : "FAIL") << std::endl << std::endl;

	// ---------------------------------------
	// 		Set up and test Grabowski
	// ---------------------------------------
	std::cout << "Testing Grabowski model..." << std::endl;
	pass = true;
	n = 0;
	// Load all Grabowski test data directory
	dir_name = "test0/Grabowski";
	DIR *Grabowski_dir = opendir(dir_name.c_str());
	if(Grabowski_dir) // dir is open
	{
		// loop over all files:
		dirent* result = readdir(Grabowski_dir);
		while( (result=readdir(Grabowski_dir)) != NULL )
		{
			// try to load if an actual file:
			if(std::string(result->d_name) != ".." && result->d_name[0] != '.')
			{
				// construct relative file path/name:
				std::string fname(dir_name);
				fname.append("/");
				fname.append(result->d_name);

				// storage for data:
				float mt, Zt;
				std::vector< std::array<float,4> > field_data;
				std::vector< std::array<float,2> > test_data;

				// Call helper function to read data
				read_plasma_file(fname, mt, Zt, field_data, test_data);

				// create model and test:
				StopPow::StopPow* s = new StopPow::StopPow_Grabowski(mt, Zt, field_data);
				bool test = run_test(s, test_data, 1e-2, verbose);
				std::cout << fname << ": " << (test ? "pass" : "FAIL") << std::endl;
				pass = pass && test;
				n++;

				// free memory:
				delete s;
			}
		}
	}
	std::cout << n << " Grabowski model(s) tested: " << (pass ? "pass" : "FAIL") << std::endl << std::endl;
	

	// ---------------------------------------
	// 	    Test Bethe-Bloch against SRIM
	// ---------------------------------------
	std::cout << "Testing Bethe-Bloch model..." << std::endl;
	StopPow::StopPow * SRIM = new StopPow::StopPow_SRIM("test0/Hydrogen in Aluminum.txt");
	StopPow::StopPow * BB = new StopPow::StopPow_BetheBloch(1, 1, {26.98}, {13.0}, {6.026e22});
	pass = true;
	for(float E=1.0; E<=20.; E+=0.05)
	{
		float dEdx_BB = BB->dEdx_MeV_um(E);
		float dEdx_SRIM = SRIM->dEdx_MeV_um(E);
		float delta = fabs(dEdx_BB - dEdx_SRIM) / fabs(dEdx_SRIM);
		pass = pass && (delta < 3e-2);
		if(verbose)
			std::cout << E << " , " << dEdx_SRIM << " , " << dEdx_BB << " -> " << pass << std::endl;
	}
	std::cout << " Bethe-Bloch model tested: " << (pass ? "pass" : "FAIL") << std::endl << std::endl;
	

	// ---------------------------------------
	// 	    Test A-Z against SRIM
	// ---------------------------------------
	std::cout << "Testing Andersen-Ziegler model..." << std::endl;
	StopPow::StopPow * AZ = new StopPow::StopPow_AZ(13);
	pass = true;
	for(float E=1.0; E<=20.; E+=0.05)
	{
		float dEdx_AZ = AZ->dEdx_MeV_um(E);
		float dEdx_SRIM = SRIM->dEdx_MeV_um(E);
		float delta = fabs(dEdx_AZ - dEdx_SRIM) / fabs(dEdx_SRIM);
		pass = pass && (delta < 3e-2);
		if(verbose)
			std::cout << E << " , " << dEdx_SRIM << " , " << dEdx_AZ << " -> " << pass << std::endl;
	}
	std::cout << " Andersen-Ziegler model tested: " << (pass ? "pass" : "FAIL") << std::endl << std::endl;

	// ---------------------------------------
	// 		Set up and test Zimmerman
	// ---------------------------------------
	std::cout << "Testing Zimmerman model..." << std::endl;
	pass = true;
	n = 0;
	// Load all Zimmerman test data directory
	dir_name = "test0/Zimmerman";
	DIR *Zimmerman_dir = opendir(dir_name.c_str());
	if(Zimmerman_dir) // dir is open
	{
		// loop over all files:
		dirent* result = readdir(Zimmerman_dir);
		while( (result=readdir(Zimmerman_dir)) != NULL )
		{
			// try to load if an actual file:
			if(std::string(result->d_name) != ".." && result->d_name[0] != '.')
			{
				// construct relative file path/name:
				std::string fname(dir_name);
				fname.append("/");
				fname.append(result->d_name);

				// storage for data:
				float mt, Zt, Te;
				std::vector< std::array<float,5> > field_data;
				std::vector< std::array<float,2> > test_data;

				// Call helper function to read data
				read_partial_ioniz_file(fname, mt, Zt, field_data, Te, test_data);

				// create model and test:
				StopPow::StopPow* s = new StopPow::StopPow_Zimmerman(mt, Zt, field_data, Te);
				bool test = run_test(s, test_data, 3e-2, verbose);
				std::cout << fname << ": " << (test ? "pass" : "FAIL") << std::endl;
				pass = pass && test;
				n++;

				// free memory:
				delete s;
			}
		}
	}
	std::cout << n << " Zimmerman model(s) tested: " << (pass ? "pass" : "FAIL") << std::endl << std::endl;

	return 0;
	
}