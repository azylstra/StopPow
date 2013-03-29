/**
 * @brief Calculate Bethe-Bloch stopping power.
 * 
 * Implement a stopping-power calculator for arbitrary cold matter, using
 * the simple Bethe-Bloch theory.
 *
 * @class StopPow_BetheBloch
 * @author Alex Zylstra
 * @date 2013/03/29
 */

#ifndef STOPPOW_BETHEBLOCH_H
#define STOPPOW_BETHEBLOCH_H

#include "StopPow.h"
#include "StopPow_Constants.h"
#include <math.h>
#include <stdexcept>
#include <vector>

using namespace std;

class StopPow_BetheBloch : public StopPow
{
public:
	/** Initialize the Bethe-Bloch calculator.
	 * @param mt the test particle mass in AMU
	 * @param Zt the test particle in charge (units of e)
	 * @param mf vector containing ordered field particle masses in AMU
	 * @param Zt vector containing ordered field particle charges in units of e
	 * @param nf vector containing ordered field particle densities in units of 1/cm3
 	 * @throws invalid_argument
	 */
	StopPow_BetheBloch(float mt, float Zt, vector<float> mf , vector<float> Zf, vector<float> nf);

	/** Calculate the total stopping power
	 * @param E the test particle energy in MeV
	 * @return stopping power in units of MeV/um
 	 * @throws invalid_argument
	 */
	float dEdx_MeV_um(float E);

	/** Calculate the total stopping power
	 * @param E the test particle energy in MeV
	 * @return stopping power in units of MeV/(mg/cm2)
 	 * @throws invalid_argument
	 */
	float dEdx_MeV_mgcm2(float E);

private:
	/** Effecive ionization potential as a function of Z.
	 * @param Zf field particle charge in units of e
	 * @return Ibar in erg
	 */
	float Ibar(float Zf);

	// data on the field particles:
	vector<float> mf; /** mass in atomic units */
	vector<float> Zf; /** charge in atomic units */
	vector<float> nf; /** particle density in 1/cc */
	int num; /** number of field particle species */
	float rho; /** mass density in g/cc */
	// type of test particle:
	float mt; /** mass in atomic units */
	float Zt; /** charge in atomic units */

	// ionization data
	static const float IbarData[];
};

#endif