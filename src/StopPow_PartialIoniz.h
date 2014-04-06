/**
 * @brief Wrapper class for partially-ionized plasma stopping power theories.
 * 
 * This does not actually implement a theory, but defines an interface
 * for all theories which can handle partial ionization.
 *
 * @class StopPow::StopPow_PartialIoniz
 * @author Alex Zylstra
 * @date 2014/04/02
 * @copyright Alex Zylstra / MIT
 */

#ifndef STOPPOW_PartialIoniz_H
#define STOPPOW_PartialIoniz_H

#include <iostream>
#include <vector>
#include <array>
#include <stdexcept>

#include "StopPow.h"
#include "StopPow_Constants.h"

namespace StopPow
{

class StopPow_PartialIoniz : public StopPow
{
public:
	/** Initialize the stopping power. Electrons should not be included in lists - they will be added automatically!
	 * @param mt the test particle mass in AMU
	 * @param Zt the test particle in charge (units of e)
	 * @param mf vector containing ordered field ion masses in AMU
	 * @param Zf vector containing ordered field ion charges in units of e
	 * @param Tf vector containing ordered field ion temperatures in units of keV
	 * @param nf vector containing ordered field ion densities in units of 1/cm3
	 * @param Zbar a vector containing the average ionization state for each field ion. Zbar=Z corresponds to fully ionized material.
 	 * @param Te the electron temperature in keV
 	 * @throws invalid_argument
	 */
	StopPow_PartialIoniz(double mt, double Zt, std::vector<double> & mf, std::vector<double> & Zf, std::vector<double> & Tf, std::vector<double> & nf, std::vector<double> & Zbar, double Te) throw(std::invalid_argument);

	/** Initialize the stopping power. Electrons should not be included in lists - they will be added automatically!
	 * @param mt the test particle mass in AMU
	 * @param Zt the test particle in charge (units of e)
	 * @param field vector containing field ion info. Each row is one type of ion, then the array must contain:
	 * [mf,Zf,Tf,nf,Zbar] in units of AMU, e, keV, and 1/cc, e
	 * @param Te the electron temperature in keV
 	 * @throws invalid_argument
	 */
	StopPow_PartialIoniz(double mt, double Zt, std::vector< std::array<double,5> > & field, double Te) throw(std::invalid_argument);

	/** Destructor */
	~StopPow_PartialIoniz();
	
	/** Modify the test particle used in the theory
	 * @param mt the test particle mass in AMU
	 * @param Zt the test particle in charge (units of e)
	 */
	void set_particle(double mt, double Zt);
	
	/** Modify the field particles used in the theory. Electrons should not be included in lists - they will be added automatically!
	 * @param mf vector containing ordered field ions masses in AMU
	 * @param Zf vector containing ordered field ions charges in units of e
	 * @param Tf vector containing ordered field ions temperatures in units of keV
	 * @param nf vector containing ordered field ions densities in units of 1/cm3
	 * @param Zbar a vector containing the average ionization state for each field ion. Zbar=Z corresponds to fully ionized material.
 	 * @param Te the electron temperature in keV
	 */
	void set_field(std::vector<double> & mf, std::vector<double> & Zf, std::vector<double> & Tf, std::vector<double> & nf, std::vector<double> & Zbar, double Te);
	
	/** Modify the field ionss used in the theory. Electrons should not be included in lists - they will be added automatically!
	 * @param field vector containing field ions info. Each row is one type of ions, then the array must contain:
	 * [mf,Zf,Tf,nf,Zbar] in units of AMU, e, keV, and 1/cc, e
 	 * @param Te the electron temperature in keV
 	 * @throws invalid_argument
	 */
	void set_field(std::vector< std::array<double,5> > & field, double Te) throw(std::invalid_argument);

protected:
	// data on the field ions:
	/** mass in atomic units */
	std::vector<double> mf; 
	/** charge in atomic units */
	std::vector<double> Zf; 
	/** Temperature in keV */
	std::vector<double> Tf; 
	/** particle density in 1/cc */
	std::vector<double> nf; 
	/** particle ionization state */
	std::vector<double> Zbar;
	/** Free electron number density */
	double ne;
	/** Free electron temperature */
	double Te;
	/** number of field particle species */
	int num; 
	/** mass density in g/cc */
	double rho; 

	// type of test particle:
	/** mass in atomic units */
	double mt; 
	/** charge in atomic units */
	double Zt; 
};

} // end namespace StopPow

#endif