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