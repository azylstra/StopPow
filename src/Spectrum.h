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
 * @brief Spectrum utilities
 * 
 * This file defines several static methods withing
 * the StopPow namespace used as helper functions for
 * dealing with spectral analysis
 * 
 * @author Alex Zylstra
 * @date 2014/03/09
 * @copyright MIT / Alex Zylstra
 */

#ifndef SPECTRUM_H
#define SPECTRUM_H
 
#include <math.h>
#include <cmath>
#include <stdexcept>
#include <vector>
#include <array>
#include "StopPow.h"
#include "Util.h"

namespace StopPow
{
	/** Shift a spectrum using a stopping power model and a given thickness. Uses the given model. Result put in argument vectors
	* @param model the StopPow model to use
	* @param thickness the thickness to transmit the spectrum through. Note: uses `mode' that model is set to. Can be negative, in which case the spectrum is upshifted.
	* @param data_E the energy bin values in MeV
	* @param data_Y the yield values for each energy in Yield/MeV
	*/
	void shift(StopPow & model, double thickness, std::vector<double> & data_E, std::vector<double> & data_Y) throw(std::invalid_argument);

	/** Shift a spectrum using a stopping power model and a given thickness. Uses the given model. Result is put in argument vectors
	* @param model the StopPow model to use
	* @param thickness the thickness to transmit the spectrum through. Note: uses `mode' that model is set to. Can be negative, in which case the spectrum is upshifted.
	* @param data_E the energy bin values in MeV
	* @param data_Y the yield values for each energy in Yield/MeV
	* @param data_err the error bars on yield
	*/
	void shift(StopPow & model, double thickness, std::vector<double> & data_E, std::vector<double> & data_Y, std::vector<double> & data_err) throw(std::invalid_argument);
} // end of namespace StopPow

#endif