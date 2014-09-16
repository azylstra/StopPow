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

/** test class for StopPow library to test regular dE/dx call evaluation speed
 * @author Alex Zylstra
 * @date 2014/04/09
 */


#include <stdio.h>

#include <iostream>
#include <vector>
#include <ctime>

#include "StopPow.h"
#include "StopPow_SRIM.h"
#include "Util.h"

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
	std::cout << "========== Test Suite 5 ==========" << std::endl;
	std::cout << "      Numerical accuracy of  " << std::endl;
	std::cout << "     basic StopPow functions " << std::endl;
	bool pass = true;
	bool test;
	
	StopPow::StopPow * s = new StopPow::StopPow_SRIM("SRIM/Hydrogen in Aluminum.txt");

	// Test Eout method against some precomputed test cases:
	double E = 14.7, E1;
	double tol = 1e-4;
	bool Eout_pass = true;
	std::vector<double> thicknesses {1, 10, 100, 500};
	std::vector<double> expected {14.6818, 14.6304, 14.0049, 10.9047};
	for(int i=0; i<thicknesses.size(); i++)
	{
		E1 = s->Eout(E, thicknesses[i]);
		test = StopPow::approx(E1, expected[i], tol);
		if(verbose || !test)
		{
			std::cout << "Eout test: " << E << "->" << E1 << " thru " << thicknesses[i] << ", expected: " << expected[i];
			if(test)
				std::cout << " pass";
			else
				std::cout << " FAIL!";
			std::cout << std::endl;
		}

		Eout_pass &= test;
	}
	std::cout << "Eout tests: " << (Eout_pass ? "pass" : "FAIL!") << std::endl;
	pass &= Eout_pass;

	// Test Ein method against precomputed cases:
	bool Ein_pass = true;
	expected = {14.7082, 14.7693, 15.3701, 17.8445};
	for(int i=0; i<thicknesses.size(); i++)
	{
		E1 = s->Ein(E, thicknesses[i]);
		test = StopPow::approx(E1, expected[i], tol);
		if(verbose || !test)
		{
			std::cout << "Ein test: " << E << "->" << E1 << " thru " << thicknesses[i] << ", expected: " << expected[i];
			if(test)
				std::cout << " pass";
			else
				std::cout << " FAIL!";
			std::cout << std::endl;
		}

		Ein_pass &= test;
	}
	std::cout << "Ein tests: " << (Ein_pass ? "pass" : "FAIL!") << std::endl;
	pass &= Ein_pass;

	// Test thickness method against precomputed cases:
	bool T_pass = true;
	std::vector<double> E2 = {14.5, 13.0, 12.0, 11.0};
	expected = {29.2148, 238.305, 367.704, 488.965}; 
	double Ttest;
	for(int i=0; i<E2.size(); i++)
	{
		Ttest = s->Thickness(E, E2[i]);
		test = StopPow::approx(Ttest, expected[i], tol);
		if(verbose || !test)
		{
			std::cout << "Thickness test: " << E << "->" << E2[i] << " thru " << Ttest << ", expected: " << expected[i];
			if(test)
				std::cout << " pass";
			else
				std::cout << " FAIL!";
			std::cout << std::endl;
		}

		T_pass &= test;
	}
	std::cout << "Thickness tests: " << (T_pass ? "pass" : "FAIL!") << std::endl;
	pass &= T_pass;

	// Test range method against precomputed cases:
	bool R_pass = true;
	E2 = {5, 10, 15};
	expected = {189.9, 625, 1271.4}; 
	double Rtest;
	for(int i=0; i<E2.size(); i++)
	{
		Rtest = s->Range(E2[i]);
		test = StopPow::approx(Rtest, expected[i], tol);
		if(verbose || !test)
		{
			std::cout << "Range test: " << E2[i] << "->" << 0 << " thru " << Rtest << ", expected: " << expected[i];
			if(test)
				std::cout << " pass";
			else
				std::cout << " FAIL!";
			std::cout << std::endl;
		}

		R_pass &= test;
	}
	std::cout << "Range tests: " << (R_pass ? "pass" : "FAIL!") << std::endl;
	pass &= R_pass;

	if(pass)
	{
		std::cout << "PASS" << std::endl;
		return 0;
	}
	std::cout << "FAIL" << std::endl;
	return 1;
}