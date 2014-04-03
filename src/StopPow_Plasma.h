/**
 * @brief Wrapper class for pure plasma stopping power theories.
 * 
 * This does not actually implement a theory, but defines an interface
 * for all plasma-based theories.
 *
 * @class StopPow::StopPow_Plasma
 * @author Alex Zylstra
 * @date 2014/04/02
 * @copyright Alex Zylstra / MIT
 */

#ifndef STOPPOW_PLASMA_H
#define STOPPOW_PLASMA_H

#include <vector>
#include <array>
#include <stdexcept>

#include "StopPow.h"
#include "StopPow_Constants.h"

namespace StopPow
{

class StopPow_Plasma : public StopPow
{
public:
	/** Initialize the stopping power. Electrons should be included!
	 * @param mt the test particle mass in AMU
	 * @param Zt the test particle in charge (units of e)
	 * @param mf vector containing ordered field ion masses in AMU
	 * @param Zf vector containing ordered field ion charges in units of e
	 * @param Tf vector containing ordered field ion temperatures in units of keV
	 * @param nf vector containing ordered field ion densities in units of 1/cm3
 	 * @throws invalid_argument
	 */
	StopPow_Plasma(float mt, float Zt, std::vector<float> & mf, std::vector<float> & Zf, std::vector<float> & Tf, std::vector<float> & nf) throw(std::invalid_argument);

	/** Initialize the stopping power. Electrons should be included!
	 * @param mt the test particle mass in AMU
	 * @param Zt the test particle in charge (units of e)
	 * @param field vector containing field ion info. Each row is one type of ion, then the array must contain:
	 * [mf,Zf,Tf,nf] in units of AMU, e, keV, and 1/cc
 	 * @throws invalid_argument
	 */
	StopPow_Plasma(float mt, float Zt, std::vector< std::array<float,4> > & field) throw(std::invalid_argument);

	/** Initialize the stopping power. Electrons should not be included - they will be added automatically!
	 * @param mt the test particle mass in AMU
	 * @param Zt the test particle in charge (units of e)
	 * @param mf vector containing ordered field ion masses in AMU
	 * @param Zf vector containing ordered field ion charges in units of e
	 * @param Tf vector containing ordered field ion temperatures in units of keV
	 * @param nf vector containing ordered field ion densities in units of 1/cm3
	 * @param Te the electron temperature in keV
 	 * @throws invalid_argument
	 */
	StopPow_Plasma(float mt, float Zt, std::vector<float> & mf, std::vector<float> & Zf, std::vector<float> & Tf, std::vector<float> & nf, float Te) throw(std::invalid_argument);

	/** Initialize the stopping power. Electrons should not be included - they will be added automatically!
	 * @param mt the test particle mass in AMU
	 * @param Zt the test particle in charge (units of e)
	 * @param field vector containing field ion info. Each row is one type of ion, then the array must contain:
	 * [mf,Zf,Tf,nf] in units of AMU, e, keV, and 1/cc
	 * @param Te the electron temperature in keV
 	 * @throws invalid_argument
	 */
	StopPow_Plasma(float mt, float Zt, std::vector< std::array<float,4> > & field, float Te) throw(std::invalid_argument);

	/** Destructor */
	~StopPow_Plasma();
	
	/** Modify the test particle used in the theory
	 * @param mt the test particle mass in AMU
	 * @param Zt the test particle in charge (units of e)
 	 * @throws invalid_argument
	 */
	void set_particle(float mt, float Zt) throw(std::invalid_argument);
	
	/** Modify the field ions used in the theory. Electrons should be included!
	 * @param mf vector containing ordered field ion masses in AMU
	 * @param Zf vector containing ordered field ion charges in units of e
	 * @param Tf vector containing ordered field ion temperatures in units of keV
	 * @param nf vector containing ordered field ion densities in units of 1/cm3
 	 * @throws invalid_argument
	 */
	void set_field(std::vector<float> & mf, std::vector<float> & Zf, std::vector<float> & Tf, std::vector<float> & nf) throw(std::invalid_argument);
	
	/** Modify the field ions used in the theory. Electrons should be included!
	 * @param field vector containing field ion info. Each row is one type of ion, then the array must contain:
	 * [mf,Zf,Tf,nf] in units of AMU, e, keV, and 1/cc
 	 * @throws invalid_argument
	 */
	void set_field(std::vector< std::array<float,4> > & field) throw(std::invalid_argument);
	
	/** Modify the field ions used in the theory. Electrons should not be included - they will be added automatically!
	 * @param mf vector containing ordered field ion masses in AMU
	 * @param Zf vector containing ordered field ion charges in units of e
	 * @param Tf vector containing ordered field ion temperatures in units of keV
	 * @param nf vector containing ordered field ion densities in units of 1/cm3
	 * @param Te the electron temperature in keV
 	 * @throws invalid_argument
	 */
	void set_field(std::vector<float> & mf, std::vector<float> & Zf, std::vector<float> & Tf, std::vector<float> & nf, float Te) throw(std::invalid_argument);
	
	/** Modify the field ions used in the theory. Electrons should not be included - they will be added automatically!
	 * @param field vector containing field ion info. Each row is one type of ion, then the array must contain:
	 * [mf,Zf,Tf,nf] in units of AMU, e, keV, and 1/cc
	 * @param Te the electron temperature in keV
 	 * @throws invalid_argument
	 */
	void set_field(std::vector< std::array<float,4> > & field, float Te) throw(std::invalid_argument);

protected:
	// data on the field particles:
	/** mass in atomic units */
	std::vector<float> mf; 
	/** charge in atomic units */
	std::vector<float> Zf; 
	/** Temperature in keV */
	std::vector<float> Tf; 
	/** particle density in 1/cc */
	std::vector<float> nf; 
	/** number of field particle species */
	int num; 
	/** mass density in g/cc */
	float rho; 

	// type of test particle:
	/** mass in atomic units */
	float mt; 
	/** charge in atomic units */
	float Zt; 
};

} // end namespace StopPow

#endif