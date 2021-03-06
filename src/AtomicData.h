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

/*
 * Physical constants for atomic data. 
 * 
 * @author Alex Zylstra
 * @date 2013/06/05
 * @copyright MIT / Alex Zylstra
 */

#ifndef ATOMICDATA_H
#define ATOMICDATA_H
 
#include <string>
#include <array>
#include <limits>

namespace StopPow
{

/**
  * Provide access to atomic data.
  * Various parameters are available, retreived by atomic number (1-92).
  * 
  * @author Alex Zylstra
  * @date 2013/06/05
  * @copyright MIT / Alex Zylstra
  * @class AtomicData
  */
class AtomicData
{
public:
	/**
	  * Get the atomic mass
	  * @param Z the atomic number (1-92)
	  * @return the atomic mass for element Z in AMU
	  * if an invalid Z is given, returns NaN
	  */
	static double get_AMU(int Z);

	/**
	  * Get the mass density at standard conditions
	  * @param Z the atomic number (1-92)
	  * @return density in g/cm^3
	  * if an invalid Z is given, returns NaN
	  */
	static double get_rho(int Z);

	/**
	  * Get the symbol (e.g. "H" for Z=1 hydrogen)
	  * @param Z the atomic number (1-92)
	  * @return a string containing the symbol for element Z
	  * if an invalid Z is given, returns ""
	  */
	static std::string get_symbol(int Z);

	/**
	  * Look up an element by symbol
	  * @param symbol the common symbol (e.g. "Al")
	  * @return the atomic number (e.g. 13)
	  * returns -1 if it is not found
	  */
	static int get_num_from_symbol(std::string symbol);

	/**
	  * Get the name (e.g. "Hydrogen" for Z=1)
	  * @param Z the atomic number (1-92)
	  * @return a string containing the common element name
	  * if an invalid Z is given, returns ""
	  */
	static std::string get_name(int Z);

	/**
	  * Look up an element by name
	  * @param name the common name (e.g. "Aluminum")
	  * @return the atomic number (e.g. 13)
	  * returns -1 if it is not found
	  */
	static int get_num_from_name(std::string name);

	/**
	  * Get the mean ionization potential. This is for use in
	  * Bethe-Bloch style calculations. The data are not normalized
	  * to Z.
	  * Data are taken from Andersen and Ziegler, The Stopping and Ranges of Ions in Matter, Vol 3:
	  * Hydrogen stopping powers and ranges in all elements (1978).
	  * @param Z the atomic number (1-92)
	  * @return the mean ionization potential in eV
	  * if an invalid Z is given, returns NaN
	  */
	static double get_mean_ionization(int Z);

	/**
	  * Get the shell correction coefficients for Bethe-Bloch
	  * style calculations.
	  * Data are taken from Andersen and Ziegler, The Stopping and Ranges of Ions in Matter, Vol 3:
	  * Hydrogen stopping powers and ranges in all elements (1978).
	  * @param Z the atomic number (1-92)
	  * @return the five shell coefficients {A0,A1,A2,A3,A4}
	  * if an invalid Z is given, returns all values as NaN
	  */
	static std::array<double,5> get_shell_coeff(int Z);
	
	/** Number of elements defined */
	static const int n = 92;

	/* data */
private:
	/** Atomic mass data */
	static const std::array<double,n> AMU;

	/** Mass density data */
	static const std::array<double,n> rho;

	/** Symbols */
	static const std::array<std::string,n> symbol;

	/** Names */
	static const std::array<std::string,n> name;

	/** Mean ionization potential data */
	static const std::array<double,n> ioniz;

	/** Shell correction factors, for Bethe-Bloch style problems */
	static const std::array< std::array<double,5> , n> shell;
};

} // end namespace StopPow

#endif