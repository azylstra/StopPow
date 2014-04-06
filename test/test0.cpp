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
				double mt, Zt;
				std::vector< std::array<double,4> > field_data;
				std::vector< std::array<double,2> > test_data;

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
				double mt, Zt;
				std::vector< std::array<double,4> > field_data;
				std::vector< std::array<double,2> > test_data;

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
	for(double E=1.0; E<=20.; E+=0.05)
	{
		double dEdx_BB = BB->dEdx_MeV_um(E);
		double dEdx_SRIM = SRIM->dEdx_MeV_um(E);
		double delta = fabs(dEdx_BB - dEdx_SRIM) / fabs(dEdx_SRIM);
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
	for(double E=1.0; E<=20.; E+=0.05)
	{
		double dEdx_AZ = AZ->dEdx_MeV_um(E);
		double dEdx_SRIM = SRIM->dEdx_MeV_um(E);
		double delta = fabs(dEdx_AZ - dEdx_SRIM) / fabs(dEdx_SRIM);
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
				double mt, Zt, Te;
				std::vector< std::array<double,5> > field_data;
				std::vector< std::array<double,2> > test_data;

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


	// Basic BPS test
	std::cout << "Fig 6" << std::endl;
	std::vector<double> mf = {1.008, StopPow::me/StopPow::mp};
	std::vector<double> Zf = {1.0, -1};
	std::vector<double> Tf = {1.0, 1};
	std::vector<double> nf = {5e25, 5e25};
	double vb = sqrt(3*StopPow::kB*Tf[0]*StopPow::keVtoK/(StopPow::me));
	StopPow::StopPow_BPS * s = new StopPow::StopPow_BPS(1., 1., mf, Zf, Tf, nf);
	for(double E=0.1; E<=30.5; E+=1.0)
	{
		double vt = StopPow::c * sqrt(2e3*E/StopPow::mpc2);
		std::cout << vt/vb << " , " << s->dEdx_short(E) << " , " << s->dEdx_long(E) << " , " << s->dEdx_quantum(E) << std::endl;
	}
	// std::cout << "F(1) = " << GSL_REAL(s->Fc(1)) << " + " << GSL_IMAG(s->Fc(1)) << "i" << std::endl;
	// std::cout << "F(-1) = " << GSL_REAL(s->Fc(-1)) << " + " << GSL_IMAG(s->Fc(-1)) << "i" << std::endl;
	// std::cout << "F(7.5e9) = " << GSL_REAL(s->Fc(7.5e9)) << " + " << GSL_IMAG(s->Fc(7.5e9)) << "i" << std::endl;
	// std::cout << "F(-7.5e9) = " << GSL_REAL(s->Fc(-7.5e9)) << " + " << GSL_IMAG(s->Fc(-7.5e9)) << "i" << std::endl;

	std::cout << "---------------------" << std::endl;
	std::cout << "Fig 8" << std::endl;
	std::vector<double> mf2 = {StopPow::me/StopPow::mp};
	std::vector<double> Zf2 = {-1};
	std::vector<double> Tf2 = {3};
	std::vector<double> nf2 = {1e25};
	for(double E=0.1; E<=3.5; E+=0.2)
	{
		StopPow::StopPow_BPS * s = new StopPow::StopPow_BPS(4., 2., mf2, Zf2, Tf2, nf2);
		std::cout << E << " , " << s->dEdx(E) << std::endl;
	}

	std::cout << "---------------------" << std::endl;
	std::cout << "Fig 11" << std::endl;
	std::vector<double> mf3 = {3.016,2.014,StopPow::me/StopPow::mp};
	std::vector<double> Zf3 = {1,1,-1};
	std::vector<double> Tf3 = {30,30,30};
	std::vector<double> nf3 = {5e26,5e26,1e27};
	for(double E=0.05; E<=3.5; E+=0.05)
	{
		StopPow::StopPow_BPS * s = new StopPow::StopPow_BPS(4., 2., mf3, Zf3, Tf3, nf3);
		std::cout << E << " , " << s->dEdx(E) << std::endl;
	}

	std::cout << "---------------------" << std::endl;
	std::cout << "Fig 1" << std::endl;
	std::vector<double> mf4 = {StopPow::me/StopPow::mp};
	std::vector<double> Zf4 = {-1};
	std::vector<double> Tf4 = {(double)1.6e5/StopPow::keVtoK};
	std::vector<double> nf4 = {1.1e20};
	double Zt = 5;
	double At = 100;
	vb = sqrt(3*StopPow::kB*Tf4[0]*StopPow::keVtoK/(StopPow::me));
	for(double E=0.1; E<=40; E+=0.5)
	{
		double vt = StopPow::c * sqrt(2e3*E/(At*StopPow::mpc2));
		StopPow::StopPow_BPS * s = new StopPow::StopPow_BPS(At, Zt, mf4, Zf4, Tf4, nf4);
		std::cout << vt/vb << " , " << s->dEdx_short(E)+s->dEdx_long(E) << std::endl;
	}

	std::cout << "---------------------" << std::endl;
	std::cout << "for Johan" << std::endl;
	std::vector<double> mf5 = {StopPow::me/StopPow::mp};
	std::vector<double> Zf5 = {-1};
	std::vector<double> Tf5 = {0.5};
	std::vector<double> nf5 = {5e22};
	StopPow::StopPow_BPS * s5 = new StopPow::StopPow_BPS(1, 1, mf5, Zf5, Tf5, nf5);
	for(double E=0.1; E<=15; E+=0.1)
	{
		std::cout << E << " , " << s5->dEdx(E)*1e4 << std::endl;
	}

	// some BPS testing stuff
	// std::vector<double> mfZ = {StopPow::me/StopPow::mp};
	// std::vector<double> ZfZ = {-1.0};
	// std::vector<double> TfZ = {0.5};
	// std::vector<double> nfZ = {5e25};
	// StopPow::StopPow_BPS* s = new StopPow::StopPow_BPS(1, 1, mfZ, ZfZ, TfZ, nfZ);
	// std::cout << s->dEdx_short(1)*1e4 << "," << s->dEdx_long(1)*1e4 << "," << s->dEdx_quantum(1)*1e4 << std::endl;
	// std::cout << "F(1) = " << GSL_REAL(s->Fc(1)) << " + " << GSL_IMAG(s->Fc(1)) << "i" << std::endl;
	// std::cout << "F(-1) = " << GSL_REAL(s->Fc(-1)) << " + " << GSL_IMAG(s->Fc(-1)) << "i" << std::endl;
	// std::cout << "F(1e9) = " << GSL_REAL(s->Fc(1e9)) << " + " << GSL_IMAG(s->Fc(1e9)) << "i" << std::endl;
	// std::cout << "F(-1e9) = " << GSL_REAL(s->Fc(-1e9)) << " + " << GSL_IMAG(s->Fc(-1e9)) << "i" << std::endl;
	// gsl_complex x = gsl_complex_rect(1,1);
	// std::cout << "erfi(1+i) = " << GSL_REAL(s->erfi(x)) << " + " << GSL_IMAG(s->erfi(x)) << "i" << std::endl;
	// std::cout << "dawson(1) = " << gsl_sf_dawson(1.) << std::endl;
	// std::cout << "dawson(1)*(2/sqrt(pi))*e^1 = " << gsl_sf_dawson(1.)*(2./sqrt(M_PI))*M_E << std::endl;
	return 0;
	
}