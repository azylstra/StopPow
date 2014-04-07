// StopPow.i - SWIG interface
%module cStopPow
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
	#include "../src/PlotGen.h"
	#include "../src/AtomicData.h"
	#include "../src/Spectrum.h"
	#include "../src/Util.h"
%}

%include "cpointer.i"
%pointer_functions(int, intp);
%pointer_functions(float, floatp);

%include "std_vector.i"
#include <vector>
// Instantiate templates
namespace std {
   %template(IntVector) vector<int>;
   %template(FloatVector) vector<float>;
   %template(FloatVector2D) vector< vector<float> >;
}

%include "std_string.i"
#include <string>

%include "exception.i"
%typemap(throws, throws="java.io.IOException") std::ios_base::failure {
  jclass excep = jenv->FindClass("java/io/IOException");
  if (excep)
    jenv->ThrowNew(excep, $1.what());
  return $null;
}
%typemap(javabase) std::ios_base::failure "java.lang.Exception";
%typemap(throws, throws="java.lang.IllegalArgumentException") std::invalid_argument {
  jclass excep = jenv->FindClass("java/lang/IllegalArgumentException");
  if (excep)
    jenv->ThrowNew(excep, $1.what());
  
  return $null;
}
%typemap(javabase) std::invalid_argument "java.lang.Exception";


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
	double dEdx(double E) throw(std::invalid_argument);
	double Eout(double E, double x) throw(std::invalid_argument);
	double Ein(double E, double x) throw(std::invalid_argument);
	double Thickness(double E1, double E2) throw(std::invalid_argument);
	double Range(double E) throw(std::invalid_argument);
	double get_dx();
	void set_dx(double new_dx) throw(std::invalid_argument);
	int get_mode();
	void set_mode(int new_mode) throw(std::invalid_argument);

	static const double DEFAULT_DX;
	static const double DEFAULT_DRHOR;
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
%include "../src/PlotGen.h"
%include "../src/AtomicData.h"
%include "../src/Spectrum.h"
%include "../src/Util.h"