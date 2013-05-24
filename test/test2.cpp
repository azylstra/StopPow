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

	// Do some output
	std::cout << "========== Test Suite 2 ==========" << std::endl;
	std::cout << "   Testing plot generators  " << std::endl;

	// create a model to use for the test:
	std::string fname("test2/Hydrogen in Aluminum.txt");
	//StopPow::StopPow_SRIM s(fname);
	std::vector<float> mf(2);
	mf[0] = 1.0;
	mf[1] = 1/1800.;
	std::vector<float> Zf(2);
	Zf[0] = 1.0;
	Zf[1] = -1.;
	std::vector<float> Tf(2);
	Tf[0] = 1.0;
	Tf[1] = 1.0;
	std::vector<float> nf(2);
	nf[0] = 1e24;
	nf[1] = 1e24;
	nf[0] = 1e24; nf[1] = 1e24;
	StopPow::StopPow_LP s(1,1,mf,Zf,Tf,nf);


	std::vector< std::vector<float> > dEdx_plot;
	bool ret = StopPow::get_dEdx_vs_E( s , dEdx_plot );
	if(ret)
		std::cout << "dE/dx plot generated successfully" << std::endl;
	else
		std::cout << "ERROR: could not generate dE/dx plot" << std::endl;

	std::vector< std::vector<float> > Range_plot;
	ret = StopPow::get_Range_vs_E( s , Range_plot );
	if(ret)
		std::cout << "Range plot generated successfully" << std::endl;
	else
		std::cout << "ERROR: could not generate Range plot" << std::endl;

	float thickness = 100; // um
	std::vector< std::vector<float> > Eout_plot_1;
	ret = StopPow::get_Eout_vs_Ein( s , thickness , Eout_plot_1 );
	if(ret)
		std::cout << "Eout vs Ein plot generated successfully" << std::endl;
	else
		std::cout << "ERROR: could not generate Eout vs Ein plot" << std::endl;

	float Ein = 15; // MeV
	std::vector< std::vector<float> > Eout_plot_2;
	ret = StopPow::get_Eout_vs_Thickness( s , Ein , Eout_plot_2 );
	if(ret)
		std::cout << "Eout vs Thickness plot generated successfully" << std::endl;
	else
		std::cout << "ERROR: could not generate Eout vs Thickness plot" << std::endl;

	thickness = 100; // um
	std::vector< std::vector<float> > Ein_plot_1;
	ret = StopPow::get_Ein_vs_Eout( s , thickness , Ein_plot_1 );
	if(ret)
		std::cout << "Ein vs Eout plot generated successfully" << std::endl;
	else
		std::cout << "ERROR: could not generate Ein vs Eout plot" << std::endl;

	float Eout = 15; // MeV
	std::vector< std::vector<float> > Ein_plot_2;
	ret = StopPow::get_Ein_vs_Thickness( s , Eout , Ein_plot_2 );
	if(ret)
		std::cout << "Ein vs Thickness plot generated successfully" << std::endl;
	else
		std::cout << "ERROR: could not generate Ein vs Thickness plot" << std::endl;

	Ein = 15;
	std::vector< std::vector<float> > Thickness_plot_1;
	ret = StopPow::get_Thickness_vs_Eout( s , Ein , Thickness_plot_1 );
	if(ret)
		std::cout << "Thickness vs Eout plot generated successfully" << std::endl;
	else
		std::cout << "ERROR: could not generate Thickness vs Eout plot" << std::endl;

	Eout = 5;
	std::vector< std::vector<float> > Thickness_plot_2;
	ret = StopPow::get_Thickness_vs_Ein( s , Eout , Thickness_plot_2 );
	if(ret)
		std::cout << "Thickness vs Ein plot generated successfully" << std::endl;
	else
		std::cout << "ERROR: could not generate Thickness vs Ein plot" << std::endl;

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

	std::clock_t start;
	float duration;
	start = std::clock();
	for(int i=0; i<100; i++)
		ret = StopPow::get_dEdx_vs_E( s , dEdx_plot );
	// duration per call in ms:
	duration = (1000./100.)*(std::clock()-start) / (float) CLOCKS_PER_SEC;
	std::cout << "dE/dx vs E = " << duration << " ms" << std::endl;

	start = std::clock();
	for(int i=0; i<100; i++)
		ret = StopPow::get_Range_vs_E( s , Range_plot );
	// duration per call in ms:
	duration = (1000./100.)*(std::clock()-start) / (float) CLOCKS_PER_SEC;
	std::cout << "Range vs E = " << duration << " ms" << std::endl;

	start = std::clock();
	for(int i=0; i<100; i++)
		ret = StopPow::get_Eout_vs_Ein( s , thickness , Eout_plot_1 );
	// duration per call in ms:
	duration = (1000./100.)*(std::clock()-start) / (float) CLOCKS_PER_SEC;
	std::cout << "Eout vs Ein = " << duration << " ms" << std::endl;

	start = std::clock();
	for(int i=0; i<100; i++)
		ret = StopPow::get_Eout_vs_Thickness( s , Ein , Eout_plot_2 );
	// duration per call in ms:
	duration = (1000./100.)*(std::clock()-start) / (float) CLOCKS_PER_SEC;
	std::cout << "Eout vs Thickness = " << duration << " ms" << std::endl;

	start = std::clock();
	for(int i=0; i<100; i++)
		ret = StopPow::get_Ein_vs_Eout( s , thickness , Ein_plot_1 );
	// duration per call in ms:
	duration = (1000./100.)*(std::clock()-start) / (float) CLOCKS_PER_SEC;
	std::cout << "Ein vs Eout = " << duration << " ms" << std::endl;

	start = std::clock();
	for(int i=0; i<100; i++)
		ret = StopPow::get_Ein_vs_Thickness( s , Eout , Ein_plot_2 );
	// duration per call in ms:
	duration = (1000./100.)*(std::clock()-start) / (float) CLOCKS_PER_SEC;
	std::cout << "Ein vs Thickness = " << duration << " ms" << std::endl;

	start = std::clock();
	for(int i=0; i<100; i++)
		ret = StopPow::get_Thickness_vs_Eout( s , Ein , Thickness_plot_1 );
	// duration per call in ms:
	duration = (1000./100.)*(std::clock()-start) / (float) CLOCKS_PER_SEC;
	std::cout << "Thickness vs Eout = " << duration << " ms" << std::endl;

	start = std::clock();
	for(int i=0; i<100; i++)
		ret = StopPow::get_Thickness_vs_Ein( s , Eout , Thickness_plot_2 );
	// duration per call in ms:
	duration = (1000./100.)*(std::clock()-start) / (float) CLOCKS_PER_SEC;
	std::cout << "Thickness vs Ein = " << duration << " ms" << std::endl;
	
	return 0;
}