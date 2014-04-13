// StopPow.i - SWIG interface
%module StopPow
%{
	#include "../src/StopPow.h"
	#include "../src/StopPow_Plasma.h"
	#include "../src/StopPow_PartialIoniz.h"
	#include "../src/StopPow_SRIM.h"
	#include "../src/StopPow_LP.h"
	#include "../src/StopPow_BetheBloch.h"
	#include "../src/StopPow_AZ.h"
	#include "../src/StopPow_Mehlhorn.h"
	#include "../src/StopPow_Grabowski.h"
	#include "../src/StopPow_Zimmerman.h"
	#include "../src/StopPow_BPS.h"
	#include "../src/StopPow_Fit.h"
	#include "../src/AtomicData.h"
	#include "../src/PlotGen.h"
	#include "../src/Spectrum.h"
	#include "../src/Fit.h"
	#include "../src/Util.h"
%}

%include "cpointer.i"
%pointer_functions(int, intp);
%pointer_functions(float, floatp);
%pointer_functions(double, doublep);

%include "std_vector.i"
#include <vector>
// Instantiate templates
namespace std {
   %template(IntVector) vector<int>;
   %template(FloatVector) vector<float>;
   %template(FloatVector2D) vector< vector<float> >;
   %template(DoubleVector) vector<double>;
   %template(DoubleVector2D) vector< vector<double> >;
}

%include "std_string.i"
#include <string>

%include "exception.i"
#include <stdexcept>
#include <ios>
#include <iostream>
%exception {
    try {
        $action
    } catch (const std::exception &e) {
    	PyErr_SetString(PyExc_Exception, const_cast<char*>(e.what()));
    }
    catch(std::ios_base::failure &e) {
    	PyErr_SetString(PyExc_IOError, const_cast<char*>(e.what()));
    }
}

// Need to define the base class for SWIG:
namespace StopPow
{
class StopPow {
public:
	//StopPow();
	virtual double dEdx_MeV_um(double E) = 0;
	virtual double dEdx_MeV_mgcm2(double E) = 0;
	virtual double get_Emin() = 0;
	virtual double get_Emax() = 0;
	std::string get_type();
	std::string get_info();
	double dEdx(double E);
	double Eout(double E, double x);
	double Ein(double E, double x);
	double Thickness(double E1, double E2);
	double Range(double E);
	int get_mode();
	void set_mode(int new_mode);

	static const int MODE_LENGTH;
	static const int MODE_RHOR;
};
};

%include "../src/StopPow_Plasma.h"
%include "../src/StopPow_PartialIoniz.h"
%include "../src/StopPow_SRIM.h"
%include "../src/StopPow_LP.h"
%include "../src/StopPow_BetheBloch.h"
%include "../src/StopPow_AZ.h"
%include "../src/StopPow_Mehlhorn.h"
%include "../src/StopPow_Grabowski.h"
%include "../src/StopPow_Zimmerman.h"
%include "../src/StopPow_BPS.h"
%include "../src/StopPow_Fit.h"
%include "../src/AtomicData.h"
%include "../src/PlotGen.h"
%include "../src/Spectrum.h"
%include "../src/Fit.h"
%include "../src/Util.h"