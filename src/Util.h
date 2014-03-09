/** 
 * @brief Misc utilities
 *  
 * @author Alex Zylstra
 * @date 2014/03/09
 * @copyright MIT / Alex Zylstra
 */

#ifndef UTIL_H
#define UTIL_H
 
#include <math.h>

namespace StopPow
{
	
template<class T>
inline bool approx(T a, T b, T tol)
{
	return (2*abs(a-b)/(a+b) < tol);
}
	
}

#endif