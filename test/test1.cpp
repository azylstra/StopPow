/** test class for StopPow library
 * @author Alex Zylstra
 * @date 2013/05/24
 */


#include <stdio.h>
#include <dirent.h>

#include <iostream>
#include <vector>
#include <stdexcept>
#include <limits>

#include "StopPow.h"
#include "StopPow_SRIM.h"
#include "StopPow_LP.h"
#include "StopPow_BetheBloch.h"

std::vector< std::vector<double> > read_model_file(std::string fname)
{
	std::vector< std::vector<double> > file_data;
	std::ifstream file( fname.c_str() );
	if( file.is_open() )
	{
		std::string line;
		while( !file.eof() )
		{
			getline(file,line);
			// ignore name and comment lines:
			if( line.find("Model") == std::string::npos
				&& line.find("model") == std::string::npos
				&& line.find("#") == std::string::npos)
			{
				std::vector<double> new_line;
				// use stringstream for line:
				std::stringstream ss(line);
				std::string val;
				while( getline(ss,val,',') )
					new_line.push_back( atof(val.c_str()) );

				file_data.push_back(new_line);
			}
		}
	}

	return file_data;
}

bool run_tests(std::vector<StopPow::StopPow*> models, int argc, char* argv [])
{
	// check for verbosity flag:
	bool verbose = false;
	if( argc >= 2 )
	{
		for(int i=1; i < argc; i++)
		{
			std::string flag(argv[i]);
			if(flag == "--verbose")
				verbose = true;
		}
	}

	// ---------------------------------------
	// 		test dE/dx
	// ---------------------------------------
	std::cout << "Testing dE/dx functions under normal conditions" << std::flush;
	for( StopPow::StopPow * s : models )
	{
		// do 20 steps:
		double dE = (s->get_Emax()-s->get_Emin())/200.;
		try
		{
			for(double E=s->get_Emin(); E<s->get_Emax(); E+=dE)
			{
				if(verbose)
					std::cout << "dEdx_MeV_um(" << E << ") = " << std::flush;
				double temp = s->dEdx_MeV_um(E);
				if(verbose)
					std::cout << temp << std::endl;

				if(verbose)
					std::cout << "dEdx_MeV_mgcm2(" << E << ") = " << std::flush;
				temp = s->dEdx_MeV_mgcm2(E);
				if(verbose)
					std::cout << temp << std::endl;
			}
		}
		// if any errors are thrown here, it is a problem:
		catch(...) {return false;}
		std::cout << "." << std::flush;
	}
	std::cout << "done" << std::endl;

	std::cout << "Testing dE/dx functions under abnormal conditions" << std::flush;
	for( StopPow::StopPow * s : models )
	{
		double dE = 0.1;
		double temp;

		// call at a variety of abnormal conditions and make sure it behaves OK
		// should throw std::invalid_argument exceptions, so wrap in try/catch
		// call below min E:
		if(verbose)
			std::cout << "Calling dEdx below Emin" << std::endl;
		try {temp = s->dEdx_MeV_um(s->get_Emin()-dE);} 
		catch(std::invalid_argument e){} // this is what we want
		catch(...){return false;} // any other exceptions are failure
		try {temp = s->dEdx_MeV_mgcm2(s->get_Emin()-dE);} 
		catch(std::invalid_argument e){} // this is what we want
		catch(...){return false;} // any other exceptions are failure

		// call above max E:
		if(verbose)
			std::cout << "Calling dEdx above Emax" << std::endl;
		try {temp = s->dEdx_MeV_um(s->get_Emax()+dE);} 
		catch(std::invalid_argument e){} // this is what we want
		catch(...){return false;} // any other exceptions are failure
		try {temp = s->dEdx_MeV_mgcm2(s->get_Emax()+dE);} 
		catch(std::invalid_argument e){} // this is what we want
		catch(...){return false;} // any other exceptions are failure

		// call with 0 energy:
		if(verbose)
			std::cout << "Calling dEdx at E=0" << std::endl;
		try {temp = s->dEdx_MeV_um(0);} 
		catch(std::invalid_argument e){} // this is what we want
		catch(...){return false;} // any other exceptions are failure
		try {temp = s->dEdx_MeV_mgcm2(0);} 
		catch(std::invalid_argument e){} // this is what we want
		catch(...){return false;} // any other exceptions are failure

		// call with NAN energy:
		if(verbose)
			std::cout << "Calling dEdx at NaN" << std::endl;
		double nan = std::numeric_limits<double>::quiet_NaN();
		try {temp = s->dEdx_MeV_um(nan);} 
		catch(std::invalid_argument e){} // this is what we want
		catch(...){return false;} // any other exceptions are failure
		try {temp = s->dEdx_MeV_mgcm2(nan);} 
		catch(std::invalid_argument e){} // this is what we want
		catch(...){return false;} // any other exceptions are failure

		// call with inf energy:
		if(verbose)
			std::cout << "Calling dEdx at inf" << std::endl;
		double inf = std::numeric_limits<double>::infinity();
		try {temp = s->dEdx_MeV_um(inf);} 
		catch(std::invalid_argument e){} // this is what we want
		catch(...){return false;} // any other exceptions are failure
		try {temp = s->dEdx_MeV_mgcm2(inf);} 
		catch(std::invalid_argument e){} // this is what we want
		catch(...){return false;} // any other exceptions are failure

		std::cout << "." << std::flush;
	}
	std::cout << "done" << std::endl;

	// ---------------------------------------
	// 		test Eout
	// ---------------------------------------
	std::cout << "Testing Eout functions" << std::flush;
	for( StopPow::StopPow * s : models )
	{
		// do 20 steps:
		double dE = (s->get_Emax()-s->get_Emin())/20.;

		//try
		//{
			// loop over energies:
			for(double E=s->get_Emin(); E<s->get_Emax(); E+=dE)
			{
				// get range for this energy:
				double range = s->Range(E);
				double dr = range / 20.;
				// try a variety of thicknesses between 0 and range:
				for(double x=0; x<range; x+=dr)
				{
					if(verbose)
						std::cout << "Eout(" << E << "," << x << ") = " << std::flush;
					double temp = s->Eout(E,x);
					if(verbose)
						std::cout << temp << std::endl;
				}

				// also test outside limits for range:
				if(verbose)
					std::cout << "testing Eout outside E range" << std::endl;
				double temp = s->Eout(E,0);
				temp = s->Eout(E,range*2);
			}
		//}
		//catch(...) {return false;} // failed

		std::cout << "." << std::flush;
	}
	std::cout << "done" << std::endl;

	// ---------------------------------------
	// 		test Ein
	// ---------------------------------------
	std::cout << "Testing Ein functions" << std::flush;
	for( StopPow::StopPow * s : models )
	{
		// do 20 steps:
		double dE = (s->get_Emax()-s->get_Emin())/20.;

		try
		{
			// loop over energies:
			for(double E=s->get_Emin(); E<s->get_Emax(); E+=dE)
			{
				// get range for this energy:
				double range = s->Range(s->get_Emax());
				double dr = range / 20.;
				// try a variety of thicknesses between 0 and range:
				for(double x=0; x<range; x+=dr)
				{
					if(verbose)
						std::cout << "Ein(" << E << "," << x << ") = " << std::flush;
					double temp = s->Ein(E,x);
					if(verbose)
						std::cout << temp << std::endl;
				}

				// also test outside limits for range:
				if(verbose)
					std::cout << "testing Ein outside E range" << std::endl;
				double temp = s->Ein(E,0);
				temp = s->Ein(E,range*2);
			}
		}
		catch(...) {return false;} // failed

		std::cout << "." << std::flush;
	}
	std::cout << "done" << std::endl;

	// ---------------------------------------
	// 		test Thickness
	// ---------------------------------------
	std::cout << "Testing Thickness functions" << std::flush;
	for( StopPow::StopPow * s : models )
	{
		// do 20 steps:
		double dE = (s->get_Emax()-s->get_Emin())/20.;

		try
		{
			// loop over output energies:
			for(double E2=s->get_Emin(); E2<s->get_Emax(); E2+=dE)
			{
				// also loop over possible input energies
				for(double E1=E2; E1<s->get_Emax(); E1+=dE)
				{
					if(verbose)
						std::cout << "Thickness(" << E1 << "," << E2 << ") = " << std::flush;
					double temp = s->Thickness(E1,E2);
					if(verbose)
						std::cout << temp << std::endl;
				}
			}
		}
		catch(...) {return false;} // failed

		std::cout << "." << std::flush;
	}
	std::cout << "done" << std::endl;

	// ---------------------------------------
	// 		test Range
	// ---------------------------------------
	std::cout << "Testing Range functions" << std::flush;
	for( StopPow::StopPow * s : models )
	{
		// do 20 steps:
		double dE = (s->get_Emax()-s->get_Emin())/20.;

		try
		{
			// loop over energies:
			for(double E=s->get_Emin(); E<s->get_Emax(); E+=dE)
			{
				if(verbose)
					std::cout << "Range(" << E << ") = " << std::flush;
				double temp = s->Range(E);
				if(verbose)
					std::cout << temp << std::endl;
			}
		}
		catch(...){return false;} //failed

		std::cout << "." << std::flush;
	}
	std::cout << "done" << std::endl;


	// ---------------------------------------
	// 		Output some test cases
	// ---------------------------------------
	if( verbose )
	{
		std::vector<StopPow::StopPow*>::iterator it;
		for(it=models.begin(); it<models.end(); it++)
		{
			std::cout << "Emin = " << (*it)->get_Emin() << " MeV" << std::endl;
			std::cout << "Emax = " << (*it)->get_Emax() << " MeV" << std::endl;

			if( (*it)->get_Emin() < 1 && (*it)->get_Emax() > 10)
			{
				(*it)->set_mode((*it)->MODE_LENGTH);
				std::cout << "dEdx(10 MeV) = " << (*it)->dEdx(10) << " MeV/um" << std::endl;
				(*it)->set_mode((*it)->MODE_RHOR);
				std::cout << "dEdx(10 MeV) = " << (*it)->dEdx(10) << " MeV/(mg/cm2)" << std::endl;

				(*it)->set_mode((*it)->MODE_LENGTH);

				std::cout << "Eout(10 MeV, 100um) = " << (*it)->Eout(10,100) << std::endl;
				std::cout << "Ein(10 MeV, 100um) = " << (*it)->Ein(10,100) << std::endl;
				std::cout << "Thickness(10 MeV, 9 MeV) = " << (*it)->Thickness(10,9) << std::endl;
				std::cout << "Thickness(10 MeV, 1 MeV) = " << (*it)->Thickness(10,1) << std::endl;
				std::cout << "Range(10 MeV) = " << (*it)->Range(10) << std::endl;
			}
			std::cout << "-----------" << std::endl;
		}
	}

	return true;
}
int main(int argc, char* argv [])
{
	// Do some output
	std::cout << "========== Test Suite 1 ==========" << std::endl;
	std::cout << "   Testing computational aspects  " << std::endl;

	// create a variety of stopping power models to test
	std::vector<StopPow::StopPow*> models;



	// ---------------------------------------
	// 		Set up SRIM models
	// ---------------------------------------
	// Load all SRIM models from SRIM directory
	std::string dir_name("SRIM");
	DIR *SRIM_dir = opendir(dir_name.c_str());
	if(SRIM_dir) // dir is open
	{
		// loop over all files:
		dirent* result = readdir(SRIM_dir);
		while( (result=readdir(SRIM_dir)) != NULL )
		{
			// try to load SRIM:
			try
			{
				// construct relative file path/name:
				std::string fname(dir_name);
				fname.append("/");
				fname.append(result->d_name);

				// create model and add to the vector:
				StopPow::StopPow* s = new StopPow::StopPow_SRIM(fname);
				models.push_back(s);
			}
			// catch and ignore all exceptions
			catch(...){}
		}
	}
	std::cout << models.size() << " SRIM model(s) loaded" << std::endl;
	// test the SRIM models:
	bool SRIM_passed = run_tests(models,argc,argv);
	if(SRIM_passed)
		std::cout << "Passed!" << std::endl;
	else
	{
		std::cout << "FAILED!" << std::endl;
		return 1;
	}
	// clear the models:
	for( StopPow::StopPow* s : models )
		delete s;
	models.clear();


	// ---------------------------------------
	// 		Set up Li-Petrasso models
	// ---------------------------------------
	// plasma parameters read from file:
	std::string fname("test1/LiPetrasso.csv");
	std::vector< std::vector<double> > file_data = read_model_file(fname);
	
	// construct the models:
	double mt = 1;
	double Zt = 1;
	for(int i=0; i<file_data.size()-3; i+=4)
	{
		StopPow::StopPow * s = new StopPow::StopPow_LP(mt,Zt,
			file_data[i],file_data[i+1],
			file_data[i+2],file_data[i+3] );
		models.push_back(s);
	}
	std::cout << models.size() << " Li-Petrasso model(s) loaded" << std::endl;
	// test the L-P models:
	bool LP_passed = run_tests(models,argc,argv);
	if(LP_passed)
		std::cout << "Passed!" << std::endl;
	else
	{
		std::cout << "FAILED!" << std::endl;
		return 1;
	}
	// clear the models:
	for( StopPow::StopPow* s : models )
		delete s;
	models.clear();

	// ---------------------------------------
	// 		Set up Li-Petrasso models
	// ---------------------------------------
	// plasma parameters read from file:
	fname = "test1/BetheBloch.csv";
	file_data = read_model_file(fname);
	
	// construct the models:
	mt = 1;
	Zt = 1;
	for(int i=0; i<file_data.size()-2; i+=4)
	{
		StopPow::StopPow * s = new StopPow::StopPow_BetheBloch(mt,Zt,
			file_data[i],file_data[i+1],
			file_data[i+2] );
		models.push_back(s);
	}
	std::cout << models.size() << " Bethe-Bloch model(s) loaded" << std::endl;
	// test the L-P models:
	bool BB_passed = run_tests(models,argc,argv);
	if(BB_passed)
		std::cout << "Passed!" << std::endl;
	else
	{
		std::cout << "FAILED!" << std::endl;
		return 1;
	}

	std::cout << "All tests passed" << std::endl;
	return 0;
}