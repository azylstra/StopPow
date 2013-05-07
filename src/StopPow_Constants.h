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
const float r0 = 2.82e-13;
/** Electron mass [keV/c^2] */
const float mec2 = 511;
/** Proton mass [keV/c^2] */
const float mpc2 = 9.38e5;
/** Speed of light [cm/s] */
const float c = 2.997e10;
/** Electron mass [g] */
const float me = 9.109e-28;
/** Proton mass [g] */
const float mp = 1.6726e-24;
/** Boltzmann constant [erg/K] */
const float kB = 1.381e-16;
/** Avagadro's number */
const float Na = 6.022e23;
/** reduced Planck constant [erg*s] */
const float hbar = 1.054e-27;
/** Electron charge [statC] */
const float e = 4.8e-10;
/** Multiply by this to convert keV -> K */
const float keVtoK = 1.16e7;

// for some build systems, math.h does not define M_PI
#ifndef M_PI
/** Pi */
const float M_PI = 3.1415926;
#endif

} // end namespace StopPow

#endif