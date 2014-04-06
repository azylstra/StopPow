/**
 * @brief Cold-matter tabulated stopping.
 * 
 * A wrapper class for calculating stopping powers
 * using tabulated SRIM data (stored in csv files)
 * Linear interpolation is performed between data points.
 *
 * @class StopPow::StopPow_SRIM
 * @author Alex Zylstra
 * @date 2013/06/04
 * @copyright MIT / Alex Zylstra
 */

#ifndef STOPPOW_SRIM_H
#define STOPPOW_SRIM_H

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <ios>
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <math.h>
#include <algorithm>

#include "StopPow.h"

namespace StopPow
{

class StopPow_SRIM : public StopPow
{
public:
	/**
	 * Constructor for SRIM object. Data file must be in standard SRIM format.
	 * @param fname file name (or relative path) for the data
	 * @throws ios_base::failure
	 */
	explicit StopPow_SRIM(std::string fname) throw(std::ios_base::failure);

	/**
	 * Destructor
	 */
	~StopPow_SRIM();

	/**
	 * Get stopping power from the data.
	 * @param E the particle energy in MeV
	 * @return dE/dx in MeV/um
 	 * @throws invalid_argument
	 */
	double dEdx_MeV_um(double E) throw(std::invalid_argument);

	/**
	 * Get stopping power from the data.
	 * @param E the particle energy in MeV
	 * @return dE/dx in MeV/(mg/cm2)
 	 * @throws invalid_argument
	 */
	double dEdx_MeV_mgcm2(double E) throw(std::invalid_argument);

	/**
	 * Get the minimum energy that can be used for dE/dx calculations (inclusive)
	 * @return Emin in MeV
	 */
	double get_Emin();

	/**
	 * Get the maximum energy that can be used for dE/dx calculations (inclusive)
	 * @return Emax in MeV
	 */
	double get_Emax();

private:
	/**
	 * Compare two vectors by first element.
	 */
	static bool vector_compare(const std::vector<double>& v1, const std::vector<double>& v2);
	/**
	* Function to find a data point based on energy.
	*/
	static bool find_compare(const std::vector<double>& v, const double& E);
	/**
	 * Parse utility for the SRIM file's header
	 */
	void parse_header(std::stringstream& header);
	/**
	 * Parse utility for the SRIM file's body
	 */
	void parse_body(std::stringstream& body);
	/**
	 * Parse utility for the SRIM file's footer
	 */
	void parse_footer(std::stringstream& footer);

	/** The stopping power data from SRIM */
	std::vector< std::vector<double> > data;
	/** mass density in g/cm3 */
	double rho; 
	/** atomic number density in 1/cm3 */
	double ni; 
	/** Scale factor to convert data to keV/um */
	double scale_keV_um; 
	/** Scale factor to convert data to MeV/(mg/cm2) */
	double scale_Mev_mgcm2; 

	// some consts to define various things:
	/** Whitespace character used in SRIM file */
	static const char WHITESPACE;
	/** Represent a separator for the SRIM file header and main body */
	static const std::string header_sep;
	/** Represent a separator for the SRIM file body and footer */
	static const std::string footer_sep;
	/** Represent a key to look for to identify density in SRIM file */
	static const std::string KEY_DENSITY;
};

} // end namespace StopPow
 #endif