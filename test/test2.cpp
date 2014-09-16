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
 * this one tests plot generation utilities
 * @author Alex Zylstra
 * @date 2013/05/24
 */


#include <stdio.h>
#include <dirent.h>

#include <iostream>
#include <vector>
#include <stdexcept>
#include <limits>
#include <ctime>

#include "StopPow.h"
#include "StopPow_SRIM.h"
#include "StopPow_LP.h"
#include "PlotGen.h"

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
				verbose = true;
		}
	}

	bool pass = true;

	// Do some output
	std::cout << "========== Test Suite 2 ==========" << std::endl;
	std::cout << "   Testing plot generators  " << std::endl;

	// create a model to use for the test:
	std::string fname("SRIM/Hydrogen in Aluminum.txt");
	//StopPow::StopPow_SRIM s(fname);
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
	nf[0] = 1e24; nf[1] = 1e24;
	StopPow::StopPow_LP s(1,1,mf,Zf,Tf,nf);


	std::vector< std::vector<double> > dEdx_plot;
	bool ret = StopPow::get_dEdx_vs_E( s , dEdx_plot );
	if(ret)
		std::cout << "dE/dx plot generated successfully" << std::endl;
	else
		std::cout << "ERROR: could not generate dE/dx plot" << std::endl;
	pass &= ret;

	std::vector< std::vector<double> > Range_plot;
	ret = StopPow::get_Range_vs_E( s , Range_plot );
	if(ret)
		std::cout << "Range plot generated successfully" << std::endl;
	else
		std::cout << "ERROR: could not generate Range plot" << std::endl;
	pass &= ret;

	double thickness = 100; // um
	std::vector< std::vector<double> > Eout_plot_1;
	ret = StopPow::get_Eout_vs_Ein( s , thickness , Eout_plot_1 );
	if(ret)
		std::cout << "Eout vs Ein plot generated successfully" << std::endl;
	else
		std::cout << "ERROR: could not generate Eout vs Ein plot" << std::endl;
	pass &= ret;

	double Ein = 15; // MeV
	std::vector< std::vector<double> > Eout_plot_2;
	ret = StopPow::get_Eout_vs_Thickness( s , Ein , Eout_plot_2 );
	if(ret)
		std::cout << "Eout vs Thickness plot generated successfully" << std::endl;
	else
		std::cout << "ERROR: could not generate Eout vs Thickness plot" << std::endl;
	pass &= ret;

	thickness = 100; // um
	std::vector< std::vector<double> > Ein_plot_1;
	ret = StopPow::get_Ein_vs_Eout( s , thickness , Ein_plot_1 );
	if(ret)
		std::cout << "Ein vs Eout plot generated successfully" << std::endl;
	else
		std::cout << "ERROR: could not generate Ein vs Eout plot" << std::endl;
	pass &= ret;

	double Eout = 15; // MeV
	std::vector< std::vector<double> > Ein_plot_2;
	ret = StopPow::get_Ein_vs_Thickness( s , Eout , Ein_plot_2 );
	if(ret)
		std::cout << "Ein vs Thickness plot generated successfully" << std::endl;
	else
		std::cout << "ERROR: could not generate Ein vs Thickness plot" << std::endl;
	pass &= ret;

	Ein = 15;
	std::vector< std::vector<double> > Thickness_plot_1;
	ret = StopPow::get_Thickness_vs_Eout( s , Ein , Thickness_plot_1 );
	if(ret)
		std::cout << "Thickness vs Eout plot generated successfully" << std::endl;
	else
		std::cout << "ERROR: could not generate Thickness vs Eout plot" << std::endl;
	pass &= ret;

	Eout = 5;
	std::vector< std::vector<double> > Thickness_plot_2;
	ret = StopPow::get_Thickness_vs_Ein( s , Eout , Thickness_plot_2 );
	if(ret)
		std::cout << "Thickness vs Ein plot generated successfully" << std::endl;
	else
		std::cout << "ERROR: could not generate Thickness vs Ein plot" << std::endl;
	pass &= ret;

	// print if requested
	if( verbose) 
	{
		std::cout << "E (MeV) , dE/dx" << std::endl;
		for(int j=0; j<dEdx_plot[0].size(); j++)
		{
			std::cout << dEdx_plot[0][j] << "," << dEdx_plot[1][j] << std::endl;
		}
		std::cout << "--------------------------------" << std::endl;


		std::cout << "E (MeV) , Range" << std::endl;
		for(int j=0; j<Range_plot[0].size(); j++)
		{
			std::cout << Range_plot[0][j] << "," << Range_plot[1][j] << std::endl;
		}
		std::cout << "--------------------------------" << std::endl;

		std::cout << "Ein (MeV) , Eout (MeV)" << std::endl;
		for(int j=0; j<Eout_plot_1[0].size(); j++)
		{
			std::cout << Eout_plot_1[0][j] << "," << Eout_plot_1[1][j] << std::endl;
		}
		std::cout << "--------------------------------" << std::endl;

		std::cout << "Thickness , Eout (MeV)" << std::endl;
		for(int j=0; j<Eout_plot_2[0].size(); j++)
		{
			std::cout << Eout_plot_2[0][j] << "," << Eout_plot_2[1][j] << std::endl;
		}
		std::cout << "--------------------------------" << std::endl;

		std::cout << "Eout (MeV) , Ein (MeV)" << std::endl;
		for(int j=0; j<Ein_plot_1[0].size(); j++)
		{
			std::cout << Ein_plot_1[0][j] << "," << Ein_plot_1[1][j] << std::endl;
		}
		std::cout << "--------------------------------" << std::endl;

		std::cout << "Thickness , Ein (MeV)" << std::endl;
		for(int j=0; j<Ein_plot_2[0].size(); j++)
		{
			std::cout << Ein_plot_2[0][j] << "," << Ein_plot_2[1][j] << std::endl;
		}
		std::cout << "--------------------------------" << std::endl;

		std::cout << "Eout (MeV) , Thickness" << std::endl;
		for(int j=0; j<Thickness_plot_1[0].size(); j++)
		{
			std::cout << Thickness_plot_1[0][j] << "," << Thickness_plot_1[1][j] << std::endl;
		}
		std::cout << "--------------------------------" << std::endl;

		std::cout << "Ein (MeV) , Thickness" << std::endl;
		for(int j=0; j<Thickness_plot_2[0].size(); j++)
		{
			std::cout << Thickness_plot_2[0][j] << "," << Thickness_plot_2[1][j] << std::endl;
		}
		std::cout << "--------------------------------" << std::endl;

	}

	// ---------------------------------------
	//				Speed tests
	// ---------------------------------------
	std::cout << "Speed tests (ms / generation):" << std::endl;
	int n = 10;

	std::clock_t start;
	double duration;
	start = std::clock();
	for(int i=0; i<n; i++)
		ret = StopPow::get_dEdx_vs_E( s , dEdx_plot );
	// duration per call in ms:
	duration = (1000./n)*(std::clock()-start) / (double) CLOCKS_PER_SEC;
	std::cout << "dE/dx vs E = " << duration << " ms" << std::endl;

	start = std::clock();
	for(int i=0; i<n; i++)
		ret = StopPow::get_Range_vs_E( s , Range_plot );
	// duration per call in ms:
	duration = (1000./n)*(std::clock()-start) / (double) CLOCKS_PER_SEC;
	std::cout << "Range vs E = " << duration << " ms" << std::endl;

	start = std::clock();
	for(int i=0; i<n; i++)
		ret = StopPow::get_Eout_vs_Ein( s , thickness , Eout_plot_1 );
	// duration per call in ms:
	duration = (1000./n)*(std::clock()-start) / (double) CLOCKS_PER_SEC;
	std::cout << "Eout vs Ein = " << duration << " ms" << std::endl;

	start = std::clock();
	for(int i=0; i<n; i++)
		ret = StopPow::get_Eout_vs_Thickness( s , Ein , Eout_plot_2 );
	// duration per call in ms:
	duration = (1000./n)*(std::clock()-start) / (double) CLOCKS_PER_SEC;
	std::cout << "Eout vs Thickness = " << duration << " ms" << std::endl;

	start = std::clock();
	for(int i=0; i<n; i++)
		ret = StopPow::get_Ein_vs_Eout( s , thickness , Ein_plot_1 );
	// duration per call in ms:
	duration = (1000./n)*(std::clock()-start) / (double) CLOCKS_PER_SEC;
	std::cout << "Ein vs Eout = " << duration << " ms" << std::endl;

	start = std::clock();
	for(int i=0; i<n; i++)
		ret = StopPow::get_Ein_vs_Thickness( s , Eout , Ein_plot_2 );
	// duration per call in ms:
	duration = (1000./n)*(std::clock()-start) / (double) CLOCKS_PER_SEC;
	std::cout << "Ein vs Thickness = " << duration << " ms" << std::endl;

	start = std::clock();
	for(int i=0; i<n; i++)
		ret = StopPow::get_Thickness_vs_Eout( s , Ein , Thickness_plot_1 );
	// duration per call in ms:
	duration = (1000./n)*(std::clock()-start) / (double) CLOCKS_PER_SEC;
	std::cout << "Thickness vs Eout = " << duration << " ms" << std::endl;

	start = std::clock();
	for(int i=0; i<n; i++)
		ret = StopPow::get_Thickness_vs_Ein( s , Eout , Thickness_plot_2 );
	// duration per call in ms:
	duration = (1000./n)*(std::clock()-start) / (double) CLOCKS_PER_SEC;
	std::cout << "Thickness vs Ein = " << duration << " ms" << std::endl;
	
	if(pass)
	{
		std::cout << "PASS" << std::endl;
		return 0;
	}
	std::cout << "FAIL!" << std::endl;
	return 1;
}