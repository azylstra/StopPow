/** test class for StopPow library
 * this one tests atomic data utilities
 * @author Alex Zylstra
 * @date 2013/06/05
 */


#include <stdio.h>
#include <iostream>
#include <math.h>
#include <cmath>
#include <string>
#include <limits>
#include <vector>
#include "AtomicData.h"

/** Test if two numbers are approximately (1%) equal
  * @param a the first number
  * @param b the second number
  * @return true if a and b are within 1% of each other */
bool approx(double a, double b)
{
	double diff = abs(a-b);
	double avg = (a+b)/2.0;
	return (diff/avg < 0.01);
}

/**
  * Test data of a single element
  * @param Z the element number
  * @param AMU the expected atomic weight in AMU
  * @param rho the expected mass density in g/cc
  * @param symbol the expected atomic symbol
  * @param name the expected common name
  * @param ioniz the expected average ionization potential
  * @return true if the retreived data matches the expectation
  */
bool test_elem(int Z, double AMU, double rho, std::string symbol, std::string name, double ioniz)
{
	bool test =		approx(StopPow::AtomicData::get_AMU(Z) , AMU)
				&&	approx(StopPow::AtomicData::get_rho(Z) , rho )
				&&	( StopPow::AtomicData::get_symbol(Z) == symbol )
				&&	( StopPow::AtomicData::get_name(Z) == name )
				&&	approx(StopPow::AtomicData::get_mean_ionization(Z) , ioniz );

	if( !test )
	{
		std::cout << "Failed " + name << std::endl;
		std::cout << "   AMU: expected " << AMU << ", got " << StopPow::AtomicData::get_AMU(Z) << std::endl;
		std::cout << "   rho: expected " << rho << ", got " << StopPow::AtomicData::get_rho(Z) << std::endl;
		std::cout << "   Symbol: expected " << symbol << ", got " << StopPow::AtomicData::get_symbol(Z) << std::endl;
		std::cout << "   Name: expected " << name << ", got " << StopPow::AtomicData::get_name(Z) << std::endl;
		std::cout << "   Ionization: expected " << ioniz << ", got " << StopPow::AtomicData::get_mean_ionization(Z) << std::endl;
	}

	return test;
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
				verbose = true;
		}
	}

	// if the tests passes
	bool pass = true;

	// Do some output
	std::cout << "========== Test Suite 3 ==========" << std::endl;
	std::cout << "   Testing atomic data  " << std::endl;

	// test that we can retreive data for all elements:
	std::cout << "Testing data retrieval..." << std::flush;
	for(int i=1; i<=StopPow::AtomicData::n; i++)
	{
		try
		{
			StopPow::AtomicData::get_AMU(i);
			StopPow::AtomicData::get_rho(i);
			StopPow::AtomicData::get_symbol(i);
			StopPow::AtomicData::get_name(i);
			StopPow::AtomicData::get_mean_ionization(i);
		}
		catch(...)
		{
			pass = false;
			std::cout << std::endl;
			std::cout << "ERROR: could not retreive element " + std::to_string(i) << std::endl;
		}
	}
	// test data retreival for mischevous input
	std::vector<int> v {-1,0,(int)1e5,std::numeric_limits<int>::quiet_NaN(),std::numeric_limits<int>::infinity(),-std::numeric_limits<int>::infinity()};
	for(int val : v)
	{
		if( !std::isnan(StopPow::AtomicData::get_AMU(val)) )
		{
			std::cout << std::endl << "StopPow::AtomicData::get_AMU did not properly return NaN \
				when called with " << val << std::endl;
			pass = false;
		}
		if( !std::isnan(StopPow::AtomicData::get_rho(val)) )
		{
			std::cout << std::endl << "StopPow::AtomicData::get_rho did not properly return NaN \
				when called with " << val << std::endl;
			pass = false;
		}
		if( StopPow::AtomicData::get_symbol(val) != "" )
		{
			std::cout << std::endl << "StopPow::AtomicData::get_symbol did not properly return "" \
				when called with " << val << std::endl;
			pass = false;
		}
		if( StopPow::AtomicData::get_name(val) != "" )
		{
			std::cout << std::endl << "StopPow::AtomicData::get_name did not properly return "" \
				when called with " << val << std::endl;
			pass = false;
		}
		if( !std::isnan(StopPow::AtomicData::get_mean_ionization(val)) )
		{
			std::cout << std::endl << "StopPow::AtomicData::get_mean_ionization did not properly return NaN \
				when called with " << val << std::endl;
			pass = false;
		}
	}
	std::cout << "done!" << std::endl;

	// Run some hard-coded test cases:
	std::cout << "Running test cases..." << std::flush;

	// hydrogen:
	bool test = test_elem(1,1.008,8.99e-5,"H","Hydrogen",18.8);
	pass = pass && test;
	if( !test )
		std::cout << std::endl << "Failed hydrogen!" << std::endl;

	// aluminum:
	test = test_elem(13,26.98,2.7,"Al","Aluminum",162.);
	pass = pass && test;
	if( !test )
		std::cout << std::endl << "Failed aluminum!" << std::endl;

	// tantalum
	test = test_elem(73,180.94,16.69,"Ta","Tantalum",684);
	pass = pass && test;
	if( !test )
		std::cout << std::endl << "Failed tantalum!" << std::endl;

	std::cout << "done!" << std::endl;


	// end result:
	if( pass )
		std::cout << "Passed test!" << std::endl;
	else
		std::cout << "FAILED" << std::endl;
	
	return (pass ? 0 : 1);
}