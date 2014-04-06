/*
 * Physical constants for stopping power calculators. 
 * 
 * @author Alex Zylstra
 * @date 2013/05/07
 * @copyright MIT / Alex Zylstra
 */

#ifndef STOPPOW_CONSTANTS_H
#define STOPPOW_CONSTANTS_H
 
namespace StopPow
{
/** Classical electron radius in cm, = e^2/mec^2 */
const double r0 = 2.82e-13;
/** Electron mass [keV/c^2] */
const double mec2 = 510.998928;
/** Proton mass [keV/c^2] */
const double mpc2 = 938272.046;
/** Speed of light [cm/s] */
const double c = 2.99792458e10;
/** Electron mass [g] */
const double me = 9.10938291e-28;
/** Proton mass [g] */
const double mp = 1.672621777e-24;
/** Boltzmann constant [erg/K] */
const double kB = 1.3806488e-16;
/** Avagadro's number */
const double Na = 6.02214129e23;
/** Planck constant [erg*s] */
const double h = 6.62606957e-27;
const double hbar = 1.054571726e-27;
/** Electron charge [statC] */
const double e = 4.80320451e-10;
const double esu = 4.80320451e-10;
/** Electron charge in Lorentz-Heaviside units */
const double e_LH = esu * sqrt(4.*M_PI);
/** Multiply by this to convert keV -> K */
const double keVtoK = 1.1604505e7;
/** Mass of one AMU in g */
const double amu = 1.66053892e-24;

// for some build systems, math.h does not define M_PI
#ifndef M_PI
/** Pi */
const double M_PI = 3.14159265359;
#endif
#ifndef M_E
/** E */
const double M_E = 2.718281828459;
#endif

} // end namespace StopPow

#endif