/** 
 * @brief Generate data series for plotting
 * 
 * This file defines several static methods withing
 * the StopPow namespace used for generating
 * data series appropriate for plotting. The methods
 * are consistently named as
 * get_(ABSCISSA)_vs_(ORDINATE)
 * 
 * @author Alex Zylstra
 * @date 2013/05/24
 * @copyright MIT / Alex Zylstra
 */

#ifndef PLOTGEN_H
#define PLOTGEN_H

#include <math.h>
#include <stdexcept>
#include <vector>
#include "StopPow.h"

namespace StopPow
{

static const float PLOT_DEFAULT_NUM_POINTS = 100;

/** Create a datset for plotting of stopping power versus energy.
  * Uses model's energy limits for bounds and default number of points.
  * @param model the StopPow model to use to calculate dE/dx
  * @param data a std::vector based container, the results will be stored
  * here. data will have two elements corresponding to the energy and dE/dx
  * values respectively. Each of those will have n elements. 
  * The units of data's abscissa will depend on the mode that model was
  * set to when this function is called.
  * @return true if the data was calculated successfully, false otherwise
  */
bool get_dEdx_vs_E(		StopPow & model ,
						std::vector< std::vector<float> > & data );

/** Create a datset for plotting of stopping power versus energy.
  * Uses model's energy limits for bounds.
  * @param model the StopPow model to use to calculate dE/dx
  * @param num_points the number of points to generate
  * @param data a std::vector based container, the results will be stored
  * here. data will have two elements corresponding to the energy and dE/dx
  * values respectively. Each of those will have n elements. 
  * The units of data's abscissa will depend on the mode that model was
  * set to when this function is called.
  * @return true if the data was calculated successfully, false otherwise
  */
bool get_dEdx_vs_E(		StopPow & model ,
						int num_points ,
						std::vector< std::vector<float> > & data );

/** Create a datset for plotting of stopping power versus energy
  * @param model the StopPow model to use to calculate dE/dx
  * @param Emin the minimum energy value to plot (MeV) [inclusive]
  * @param Emax the maximum energy value to plot (MeV) [inclusive]
  * @param num_points the number of points to generate
  * @param data a std::vector based container, the results will be stored
  * here. data will have two elements corresponding to the energy and dE/dx
  * values respectively. Each of those will have n elements. 
  * The units of data's abscissa will depend on the mode that model was
  * set to when this function is called.
  * @return true if the data was calculated successfully, false otherwise
  */
bool get_dEdx_vs_E( 	StopPow & model ,
						float Emin ,
						float Emax ,
						int num_points ,
						std::vector< std::vector<float> > & data );

/** Create a datset for plotting of range versus energy.
  * Uses model's energy limits for bound and default step size.
  * @param model the StopPow model to use to calculate range
  * @param data a std::vector based container, the results will be stored
  * here. data will have two elements corresponding to the ordinate and abscissa
  * values respectively. Each of those will have n elements. 
  * The units of data's abscissa will depend on the mode that model was
  * set to when this function is called.
  * @return true if the data was calculated successfully, false otherwise
  */
bool get_Range_vs_E(	StopPow & model ,
						std::vector< std::vector<float> > & data );

/** Create a datset for plotting of range versus energy.
  * Uses model's energy limits for bounds.
  * @param model the StopPow model to use to calculate range
  * @param num_points the number of points to generate
  * @param data a std::vector based container, the results will be stored
  * here. data will have two elements corresponding to the ordinate and abscissa
  * values respectively. Each of those will have n elements. 
  * The units of data's abscissa will depend on the mode that model was
  * set to when this function is called.
  * @return true if the data was calculated successfully, false otherwise
  */
bool get_Range_vs_E(	StopPow & model ,
						int num_points ,
						std::vector< std::vector<float> > & data );

/** Create a datset for plotting of range versus energy
  * @param model the StopPow model to use to calculate range
  * @param Emin the minimum energy value to plot (MeV) [inclusive]
  * @param Emax the maximum energy value to plot (MeV) [inclusive]
  * @param num_points the number of points to generate
  * @param data a std::vector based container, the results will be stored
  * here. data will have two elements corresponding to the ordinate and abscissa
  * values respectively. Each of those will have n elements. 
  * The units of data's abscissa will depend on the mode that model was
  * set to when this function is called.
  * @return true if the data was calculated successfully, false otherwise
  */
bool get_Range_vs_E(	StopPow & model ,
						float Emin ,
						float Emax ,
						int num_points ,
						std::vector< std::vector<float> > & data );

/** Create a dataset for Eout vs Ein for a given thickness
  * Uses default energy limits and step size.
  * @param model the StopPow model to use to calculate downshift
  * @param num_points the number of points to generate
  * @param Thickness the thickness of material (um or mg/cm2 depending on mode)
  * @param data a std::vector based container, the results will be stored
  * here. data will have two elements corresponding to the ordinate and abscissa
  * values respectively. Each of those will have n elements. 
  * @return true if the data was calculated successfully, false otherwise
  */
bool get_Eout_vs_Ein(	StopPow & model ,
						float Thickness ,
						std::vector< std::vector<float> > & data );

/** Create a dataset for Eout vs Ein for a given thickness
  * Uses default energy limits.
  * @param model the StopPow model to use to calculate downshift
  * @param num_points the number of points to generate
  * @param Thickness the thickness of material (um or mg/cm2 depending on mode)
  * @param data a std::vector based container, the results will be stored
  * here. data will have two elements corresponding to the ordinate and abscissa
  * values respectively. Each of those will have n elements. 
  * @return true if the data was calculated successfully, false otherwise
  */
bool get_Eout_vs_Ein(	StopPow & model ,
						int num_points ,
						float Thickness ,
						std::vector< std::vector<float> > & data );

/** Create a dataset for Eout vs Ein for a given thickness
  * @param model the StopPow model to use to calculate downshift
  * @param Emin the minimum energy value to plot (MeV) [inclusive]
  * @param Emax the maximum energy value to plot (MeV) [inclusive]
  * @param num_points the number of points to generate
  * @param Thickness the thickness of material (um or mg/cm2 depending on mode)
  * @param data a std::vector based container, the results will be stored
  * here. data will have two elements corresponding to the ordinate and abscissa
  * values respectively. Each of those will have n elements. 
  * @return true if the data was calculated successfully, false otherwise
  */
bool get_Eout_vs_Ein(	StopPow & model ,
						float Emin ,
						float Emax ,
						int num_points ,
						float Thickness ,
						std::vector< std::vector<float> > & data );

/** Create a dataset for Eout vs thickness for a given Ein
  * Uses default energy limits and step size.
  * @param model the StopPow model to use to calculate downshift
  * @param num_points the number of points to generate
  * @param Ein the incident particle energy (MeV)
  * @param data a std::vector based container, the results will be stored
  * here. data will have two elements corresponding to the ordinate and abscissa
  * values respectively. Each of those will have n elements. 
  * The units of data's ordinate will depend on the mode that model was
  * set to when this function is called.
  * @return true if the data was calculated successfully, false otherwise
  */
bool get_Eout_vs_Thickness(	StopPow & model ,
						float Ein ,
						std::vector< std::vector<float> > & data );

/** Create a dataset for Eout vs thickness for a given Ein
  * Uses default energy limits.
  * @param model the StopPow model to use to calculate downshift
  * @param num_points the number of points to generate
  * @param Ein the incident particle energy (MeV)
  * @param data a std::vector based container, the results will be stored
  * here. data will have two elements corresponding to the ordinate and abscissa
  * values respectively. Each of those will have n elements. 
  * The units of data's ordinate will depend on the mode that model was
  * set to when this function is called.
  * @return true if the data was calculated successfully, false otherwise
  */
bool get_Eout_vs_Thickness(	StopPow & model ,
						int num_points ,
						float Ein ,
						std::vector< std::vector<float> > & data );

/** Create a dataset for Eout vs thickness for a given Ein
  * @param model the StopPow model to use to calculate downshift
  * @param Tmin the minimum thickness to plot (um or mg/cm2 depending on mode) [inclusive]
  * @param Tmax the maximum thickness to plot (um or mg/cm2 depending on mode) [inclusive]
  * @param num_points the number of points to generate
  * @param Ein the incident particle energy (MeV)
  * @param data a std::vector based container, the results will be stored
  * here. data will have two elements corresponding to the ordinate and abscissa
  * values respectively. Each of those will have n elements. 
  * The units of data's ordinate will depend on the mode that model was
  * set to when this function is called.
  * @return true if the data was calculated successfully, false otherwise
  */
bool get_Eout_vs_Thickness(	StopPow & model ,
						float Tmin ,
						float Tmax ,
						int num_points ,
						float Ein ,
						std::vector< std::vector<float> > & data );

/** Create a dataset for Ein vs Eout for a given thickness.
  * Using default energy limits and number of points.
  * @param model the StopPow model to use to calculate downshift
  * @param data a std::vector based container, the results will be stored
  * here. data will have two elements corresponding to the ordinate and abscissa
  * values respectively. Each of those will have n elements. 
  * @return true if the data was calculated successfully, false otherwise
  */
bool get_Ein_vs_Eout(	StopPow & model ,
						float Thickness ,
						std::vector< std::vector<float> > & data );

/** Create a dataset for Ein vs Eout for a given thickness
  * Using default energy limits.
  * @param model the StopPow model to use to calculate downshift
  * @param num_points the number of points to generate
  * @param Thickness the thickness of material (um or mg/cm2 depending on mode)
  * @param data a std::vector based container, the results will be stored
  * here. data will have two elements corresponding to the ordinate and abscissa
  * values respectively. Each of those will have n elements. 
  * @return true if the data was calculated successfully, false otherwise
  */
bool get_Ein_vs_Eout(	StopPow & model ,
						int num_points ,
						float Thickness ,
						std::vector< std::vector<float> > & data );

/** Create a dataset for Ein vs Eout for a given thickness
  * @param model the StopPow model to use to calculate downshift
  * @param Emin the minimum energy value to plot (MeV) [inclusive]
  * @param Emax the maximum energy value to plot (MeV) [inclusive]
  * @param num_points the number of points to generate
  * @param Thickness the thickness of material (um or mg/cm2 depending on mode)
  * @param data a std::vector based container, the results will be stored
  * here. data will have two elements corresponding to the ordinate and abscissa
  * values respectively. Each of those will have n elements. 
  * @return true if the data was calculated successfully, false otherwise
  */
bool get_Ein_vs_Eout(	StopPow & model ,
						float Emin ,
						float Emax ,
						int num_points ,
						float Thickness ,
						std::vector< std::vector<float> > & data );

/** Create a dataset for Ein vs thickness for a given Eout
  * Uses default limits for thickness and default number of points.
  * @param model the StopPow model to use to calculate downshift
  * @param num_points the number of points to generate
  * @param Eout the outbound particle energy (MeV)
  * @param data a std::vector based container, the results will be stored
  * here. data will have two elements corresponding to the ordinate and abscissa
  * values respectively. Each of those will have n elements. 
  * The units of data's ordinate will depend on the mode that model was
  * set to when this function is called.
  * @return true if the data was calculated successfully, false otherwise
  */
bool get_Ein_vs_Thickness(	StopPow & model ,
							float Eout ,
							std::vector< std::vector<float> > & data );

/** Create a dataset for Ein vs thickness for a given Eout
  * Uses default limits for thickness.
  * @param model the StopPow model to use to calculate downshift
  * @param num_points the number of points to generate
  * @param Eout the outbound particle energy (MeV)
  * @param data a std::vector based container, the results will be stored
  * here. data will have two elements corresponding to the ordinate and abscissa
  * values respectively. Each of those will have n elements. 
  * The units of data's ordinate will depend on the mode that model was
  * set to when this function is called.
  * @return true if the data was calculated successfully, false otherwise
  */
bool get_Ein_vs_Thickness(	StopPow & model ,
							int num_points ,
							float Eout ,
							std::vector< std::vector<float> > & data );

/** Create a dataset for Ein vs thickness for a given Eout
  * @param model the StopPow model to use to calculate downshift
  * @param Tmin the minimum thickness to plot (um or mg/cm2 depending on mode) [inclusive]
  * @param Tmax the maximum thickness to plot (um or mg/cm2 depending on mode) [inclusive]
  * @param num_points the number of points to generate
  * @param Eout the outbound particle energy (MeV)
  * @param data a std::vector based container, the results will be stored
  * here. data will have two elements corresponding to the ordinate and abscissa
  * values respectively. Each of those will have n elements. 
  * The units of data's ordinate will depend on the mode that model was
  * set to when this function is called.
  * @return true if the data was calculated successfully, false otherwise
  */
bool get_Ein_vs_Thickness(	StopPow & model ,
							float Tmin ,
							float Tmax ,
							int num_points ,
							float Eout ,
							std::vector< std::vector<float> > & data );

/** Create a dataset for Thickness vs Eout for a given Ein.
  * Uses default energy limits and number of points.
  * @param model the StopPow model to use to calculate downshift
  * @param Ein the incident particle energy (MeV)
  * @param data a std::vector based container, the results will be stored
  * here. data will have two elements corresponding to the ordinate and abscissa
  * values respectively. Each of those will have n elements. 
  * The units of data's ordinate will depend on the mode that model was
  * set to when this function is called.
  * @return true if the data was calculated successfully, false otherwise
  */
bool get_Thickness_vs_Eout(	StopPow & model ,
							float Ein ,
							std::vector< std::vector<float> > & data );

/** Create a dataset for Thickness vs Eout for a given Ein.
  * Uses default energy limits.
  * @param model the StopPow model to use to calculate downshift
  * @param num_points the number of points to generate
  * @param Ein the incident particle energy (MeV)
  * @param data a std::vector based container, the results will be stored
  * here. data will have two elements corresponding to the ordinate and abscissa
  * values respectively. Each of those will have n elements. 
  * The units of data's ordinate will depend on the mode that model was
  * set to when this function is called.
  * @return true if the data was calculated successfully, false otherwise
  */
bool get_Thickness_vs_Eout(	StopPow & model ,
							int num_points ,
							float Ein ,
							std::vector< std::vector<float> > & data );

/** Create a dataset for Thickness vs Eout for a given Ein
  * @param model the StopPow model to use to calculate downshift
  * @param Emin the minimum energy to plot (MeV) [inclusive]
  * @param Emax the maximum energy to plot (MeV) [inclusive]
  * @param num_points the number of points to generate
  * @param Ein the incident particle energy (MeV)
  * @param data a std::vector based container, the results will be stored
  * here. data will have two elements corresponding to the ordinate and abscissa
  * values respectively. Each of those will have n elements. 
  * The units of data's ordinate will depend on the mode that model was
  * set to when this function is called.
  * @return true if the data was calculated successfully, false otherwise
  */
bool get_Thickness_vs_Eout(	StopPow & model ,
							float Emin ,
							float Emax ,
							int num_points ,
							float Ein ,
							std::vector< std::vector<float> > & data );

/** Create a dataset for Thickness vs Ein for a given Eout.
  * Uses default energy limits and number of points.
  * @param model the StopPow model to use to calculate downshift
  * @param Eout the outbound particle energy (MeV)
  * @param data a std::vector based container, the results will be stored
  * here. data will have two elements corresponding to the ordinate and abscissa
  * values respectively. Each of those will have n elements. 
  * The units of data's ordinate will depend on the mode that model was
  * set to when this function is called.
  * @return true if the data was calculated successfully, false otherwise
  */
bool get_Thickness_vs_Ein(	StopPow & model ,
							float Eout ,
							std::vector< std::vector<float> > & data );

/** Create a dataset for Thickness vs Ein for a given Eout.
  * Uses default energy limits.
  * @param model the StopPow model to use to calculate downshift
  * @param num_points the number of points to generate
  * @param Eout the outbound particle energy (MeV)
  * @param data a std::vector based container, the results will be stored
  * here. data will have two elements corresponding to the ordinate and abscissa
  * values respectively. Each of those will have n elements. 
  * The units of data's ordinate will depend on the mode that model was
  * set to when this function is called.
  * @return true if the data was calculated successfully, false otherwise
  */
bool get_Thickness_vs_Ein(	StopPow & model ,
							int num_points ,
							float Eout ,
							std::vector< std::vector<float> > & data );

/** Create a dataset for Thickness vs Ein for a given Eout
  * @param model the StopPow model to use to calculate downshift
  * @param Emin the minimum energy to plot (MeV) [inclusive]
  * @param Emax the maximum energy to plot (MeV) [inclusive]
  * @param num_points the number of points to generate
  * @param Eout the outbound particle energy (MeV)
  * @param data a std::vector based container, the results will be stored
  * here. data will have two elements corresponding to the ordinate and abscissa
  * values respectively. Each of those will have n elements. 
  * The units of data's ordinate will depend on the mode that model was
  * set to when this function is called.
  * @return true if the data was calculated successfully, false otherwise
  */
bool get_Thickness_vs_Ein(	StopPow & model ,
							float Emin ,
							float Emax ,
							int num_points ,
							float Eout ,
							std::vector< std::vector<float> > & data );
// Thickness vs Eout given Ein
// Thickness vs Ein given Eout

} // end of namespace StopPow

#endif