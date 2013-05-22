// StopPow.i - SWIG interface
%module cStopPow
%{
	#include "../src/StopPow.h"
	#include "../src/StopPow_SRIM.h"
	#include "../src/StopPow_LP.h"
	#include "../src/StopPow_BetheBloch.h"
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
	virtual float dEdx_MeV_um(float E) = 0;
	virtual float dEdx_MeV_mgcm2(float E) = 0;
	virtual float get_Emin() = 0;
	virtual float get_Emax() = 0;
	float dEdx(float E) throw(std::invalid_argument);
	float Eout(float E, float x) throw(std::invalid_argument);
	float Ein(float E, float x) throw(std::invalid_argument);
	float Thickness(float E1, float E2) throw(std::invalid_argument);
	float Range(float E) throw(std::invalid_argument);
	float get_dx();
	void set_dx(float new_dx) throw(std::invalid_argument);
	int get_mode();
	void set_mode(int new_mode) throw(std::invalid_argument);

	static const float DEFAULT_DX;
	static const float DEFAULT_DRHOR;
	static const int MODE_LENGTH;
	static const int MODE_RHOR;
};
};

%include "../src/StopPow_SRIM.h"
%include "../src/StopPow_LP.h"
%include "../src/StopPow_BetheBloch.h"
