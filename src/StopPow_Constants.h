/** 
 * @brief Physical constants for stopping power calculators. 
 * 
 * @author Alex Zylstra
 * @date 2013/04/02
 * @copyright MIT / Alex Zylstra
 */


#ifndef STOPPOW_CONSTANTS_H
#define STOPPOW_CONSTANTS_H
 
namespace StopPow
{

const float r0 = 2.82e-13; //classical electron radius in cm
const float mec2 = 511; // keV/c2
const float mpc2 = 9.38e5; // keV/c2
const float c = 2.997e10;
const float me = 9.109e-28;
const float mp = 1.6726e-24;
const float kB = 1.381e-16;
const float Na = 6.022e23;
const float hbar = 1.054e-27;
const float e = 4.8e-10;
const float keVtoK = 1.16e7;
const float keVtoeV = 1.602e-9;

#ifndef M_PI
const float M_PI = 3.1415926;
#endif

} // end namespace StopPow

#endif