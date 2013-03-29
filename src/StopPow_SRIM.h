/**
 * @class StopPow_SRIM
 * @brief Cold-matter tabulated stopping.
 * 
 * A wrapper class for calculating stopping powers
 * using tabulated SRIM data (stored in csv files)
 * Linear interpolation is performed between data points.
 *
 * @author Alex Zylstra
 * @date 2013/03/29
 */

#ifndef STOPPOW_SRIM_H
#define STOPPOW_SRIM_H

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include "StopPow.h"

using namespace std;

class StopPow_SRIM : public StopPow
{
public:
	/**
	 * Constructor for SRIM object. Data file must be in standard SRIM format.
	 * @param fname file name (or relative path) for the data
	 * @throws ios_base::failure
	 */
	StopPow_SRIM(string fname);

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
	float dEdx_MeV_um(float E);

	/**
	 * Get stopping power from the data.
	 * @param E the particle energy in MeV
	 * @return dE/dx in MeV/(mg/cm2)
 	 * @throws invalid_argument
	 */
	float dEdx_MeV_mgcm2(float E);

private:
	/**
	 * Compare two vectors by first element.
	 */
	static bool vector_compare(const vector<float>& v1, const vector<float>& v2);
	/**
	* Function to find a data point based on energy.
	*/
	static bool find_compare(const vector<float>& v, const float& E);
	/**
	 * Parse utility for the SRIM file's header
	 */
	void parse_header(stringstream& header);
	/**
	 * Parse utility for the SRIM file's body
	 */
	void parse_body(stringstream& body);
	/**
	 * Parse utility for the SRIM file's footer
	 */
	void parse_footer(stringstream& footer);

	/**
	 * The stopping power data from SRIM
	 */
	vector< vector<float> > data;
	float rho; /** mass density in g/cm3 */
	float ni; /** atomic number density in 1/cm3 */
	float scale_keV_um; /** Scale factor to convert data to keV/um */
	float scale_Mev_mgcm2; /** Scale factor to convert data to MeV/(mg/cm2) */

	// some consts to define various things:
	static const char WHITESPACE;
	static const string header_sep;
	static const string footer_sep;
	static const string KEY_DENSITY;
};

 #endif