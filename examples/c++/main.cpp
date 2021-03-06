/** test class for StopPow library
 * @author Alex Zylstra
 * @date 2013/05/07
 */

/*
 * Building this example:
 * First, compile the library (lib directory in the top-level directory)
 * Copy the library binary to this directory
 * Use the included makefile to compile the example
 */

#include <stdio.h>

#include <iostream>
#include <vector>
#include <stdexcept>

#include "StopPow.h"
#include "StopPow_SRIM.h"
#include "StopPow_LP.h"
#include "StopPow_BetheBloch.h"

int main(int argc, char* argv [])
{
	// create a variety of stopping power models to test
	std::vector<StopPow::StopPow*> models;

	// cold matter:
	// proton in solid aluminum evaluated with SRIM
	StopPow::StopPow* s = new StopPow::StopPow_SRIM("data/Hydrogen in Aluminum.txt");
	models.push_back(s);
	
	// Li-Petrasso
	// proton in hydrogen plasma at 1e24 ions/cc and 1keV temp
	std::vector<double> mf(2);
	mf[0] = 1.0;
	mf[1] = 1/1800.;
	std::vector<double> Zf(2);
	Zf[0] = 1.0;
	Zf[1] = -1.;
	std::vector<double> Tf(2);
	Tf[0] = 1.0;
	Tf[1] = 1.0;
	std::vector<double> nf(2);
	nf[0] = 1e24;
	nf[1] = 1e24;
	StopPow::StopPow* s2 = new StopPow::StopPow_LP(1,1,mf,Zf,Tf,nf);
	models.push_back(s2);

	// Bethe-Bloch
	// protons in cold diamond
	std::vector<double> mf2;
	mf2.push_back(12.);
	std::vector<double> Zf2;
	Zf2.push_back(6);
	std::vector<double> nf2;
	nf2.push_back(1.76e23);
	StopPow::StopPow* s3 = new StopPow::StopPow_BetheBloch(1,1,mf2,Zf2,nf2);
	models.push_back(s3);

	// test models
	std::vector<StopPow::StopPow*>::iterator it;
	for(it=models.begin(); it<models.end(); it++)
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
		std::cout << "-----------" << std::endl;
	}
	return 0;
}