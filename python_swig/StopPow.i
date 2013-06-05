// StopPow.i - SWIG interface
%module StopPow
%{
	#include "../src/StopPow.h"
	#include "../src/StopPow_SRIM.h"
	#include "../src/StopPow_LP.h"
	#include "../src/StopPow_BetheBloch.h"
	#include "../src/PlotGen.h"
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

//%nspace StopPow::StopPow;
//%nspace StopPow::StopPow_LP;

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
	std::string get_type();
	std::string get_info();
	float dEdx(float E);
	float Eout(float E, float x);
	float Ein(float E, float x);
	float Thickness(float E1, float E2);
	float Range(float E);
	float get_dx();
	void set_dx(float new_dx);
	int get_mode();
	void set_mode(int new_mode);

	static const float DEFAULT_DX;
	static const float DEFAULT_DRHOR;
	static const int MODE_LENGTH;
	static const int MODE_RHOR;
};
};

%include "../src/StopPow_SRIM.h"
%include "../src/StopPow_LP.h"
%include "../src/StopPow_BetheBloch.h"
%include "../src/PlotGen.h"