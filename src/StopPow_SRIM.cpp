
#include "StopPow_SRIM.h"

namespace StopPow
{

const char StopPow_SRIM::WHITESPACE = ' ';
const std::string StopPow_SRIM::header_sep = "--------------";
const std::string StopPow_SRIM::footer_sep = "--------------------";
const std::string StopPow_SRIM::KEY_DENSITY = "Target Density";

 // Constructor
StopPow_SRIM::StopPow_SRIM(std::string fname) throw(std::ios_base::failure)
{
	// default mode for SRIM:
	set_mode(MODE_LENGTH);

	// clear data:
	data.clear();
	
	// open file:
	std::string line;
	std::ifstream myfile ( fname.c_str() );
	// for splitting the file into 3 sections:
	bool header_complete = false;
	bool body_complete = false;
	std::stringstream header;
	std::stringstream body;
	std::stringstream footer;
	// make sure file is open
	if (myfile.is_open())
	{
		// read in every line sequentially:
		while ( getline(myfile,line) )
		{
			// check to see if we've reached the end of
			// the body section:
			if( line.find(footer_sep) != std::string::npos )
				body_complete = true;

			// add this line to the appropriate section:
			if( !header_complete )
				header << line << std::endl;
			else if( !body_complete )
				body << line << std::endl;
			else
				footer << line << std::endl;

			// check to see if we've reached the end of
			// the header:
			if( line.find(header_sep) != std::string::npos )
				header_complete = true;
		}
		myfile.close();	

		// parse three sections:
		parse_header(header);
		parse_body(body);
		parse_footer(footer);
	}

	// if we could not open the file:
	else
	{
		throw std::ios_base::failure("Could not read data from file.");
	}

	// Sort by particle energy, i.e. 1st element of data
	sort( data.begin() , data.end() , StopPow_SRIM::vector_compare );
}

//destructor
StopPow_SRIM::~StopPow_SRIM()
{
	data.clear();
}

// Calculate stopping power for an arbitrary energy (MeV). Returns MeV/um
float StopPow_SRIM::dEdx_MeV_um(float E) throw(std::invalid_argument)
{
	// check limits:
	if( E < data[0][0] || E > data[data.size()-1][0])
	{
		std::stringstream msg;
		msg << "Energy passed to StopPow_SRIM::dEdx is bad: " << E;
		throw std::invalid_argument(msg.str());
	}
	// check the upper bound to prevent interpolation errors
	if( E == data[data.size()-1][0] )
		return data[data.size()-1][1];

	// Find two data points which bracket the requested energy
	std::vector< std::vector<float> >::iterator val2 = lower_bound( data.begin() , data.end() , E, find_compare );
	std::vector< std::vector<float> >::iterator val1;
	if( val2 != data.begin() ) 
		val1 = val2-1;
	else
		val1 = data.begin();

	// linear interpolation:
	float slope = 0;
	if( val2 != val1)
		slope = ((*val2)[1] - (*val1)[1]) / ((*val2)[0]-(*val1)[0]);
	float ret = (*val1)[1] + slope*(E-(*val1)[0]);

	 // flip sign and convert to MeV/um
	return -1.0*scale_keV_um*1e-3*ret;
}

// Calculate stopping power for an arbitrary energy (MeV). Returns MeV/(mg/cm2)
float StopPow_SRIM::dEdx_MeV_mgcm2(float E) throw(std::invalid_argument)
{
	return dEdx_MeV_um(E) * (scale_Mev_mgcm2/(scale_keV_um*1e-3));
}


/**
 * Get the minimum energy that can be used for dE/dx calculations
 * @return Emin in MeV
 */
float StopPow_SRIM::get_Emin()
{
	return data[1][0];
}

/**
 * Get the maximum energy that can be used for dE/dx calculations
 * @return Emax in MeV
 */
float StopPow_SRIM::get_Emax()
{
	return data[data.size()-1][0];
}


// Compare two vectors of floats by first element
bool StopPow_SRIM::vector_compare(const std::vector<float>& v1, const std::vector<float>& v2)
{
	return v1[0] < v2[0];
}
// Function to find data in database based on the energy value
bool StopPow_SRIM::find_compare(const std::vector<float>& v, const float& E)
{
	return v[0] <= E;
}

/**
 * Parse utility for the SRIM file's header
 */
void StopPow_SRIM::parse_header(std::stringstream& header)
{
	// From the header, we only want to parse the density

	// read each line of header:
	std::string line;
	while( getline(header,line) )
	{
		if( line.find(KEY_DENSITY) != std::string::npos )
		{
			// split the line into two parts containing values:
			int i1 = line.find("=");
			int i2 = line.find("=",i1+1);
			std::string density1 = line.substr(i1+1,i2-i1-2);
			std::string density2 = line.substr(i2+1);

			// remove extra whitespace at beginning/end:
			i1 = density1.find_first_not_of(WHITESPACE);
			i2 = density1.find_last_not_of(WHITESPACE);
			density1 = density1.substr( i1 , i2-i1+1 );
			i1 = density2.find_first_not_of(WHITESPACE);
			i2 = density2.find_last_not_of(WHITESPACE);
			density2 = density2.substr( i1 , i2-i1+1 );

			// parse mass density value and units:
			float density1_val = atof( density1.substr(0,density1.find(WHITESPACE)).c_str() );
			std::string density1_units = density1.substr( density1.find(WHITESPACE)+1 );
			if(density1_units.find("g/cm3") != std::string::npos)
				density1_val = density1_val*1.0;
			else if(density1_units.find("kg/m3") != std::string::npos)
				density1_val = density1_val*1e-3;
			else
				throw std::ios_base::failure("Could not parse header from file.");

			// parse number density value and units:
			float density2_val = atof( density2.substr(0,density2.find(WHITESPACE)).c_str() );
			std::string density2_units = density2.substr( density2.find(WHITESPACE)+1 );
			if(density2_units.find("atoms/cm3") != std::string::npos)
				density2_val = density2_val*1.0;
			else if(density2_units.find("atoms/m3") != std::string::npos)
				density2_val = density2_val*1e-6;
			else
				throw std::ios_base::failure("Could not parse header from file.");

			// set class variables appropriately:
			rho = density1_val;
			ni = density2_val;
		}
	}
	return;
}

/**
 * Parse utility for the SRIM file's body
 */
void StopPow_SRIM::parse_body(std::stringstream& body)
{
	// loop over each line in the body:
	std::string line;
	while(getline(body,line))
	{
		// break up the line into elements
		// format is whitespace-separated values:
		std::stringstream linestream(line);
		std::string element;
		std::vector<std::string> line_elements;
		while( getline(linestream,element,WHITESPACE))
			if(element.size()>0)
				line_elements.push_back(element);

		// get the Energy from the first element:
		float Energy = atof( line_elements[0].c_str() );
		// parse units for energy:
		if( line_elements[1].find("keV") != std::string::npos)
			Energy = Energy*1e-3;
		else if( line_elements[1].find("MeV") != std::string::npos)
			Energy = Energy*1.0;
		else
			throw std::ios_base::failure("Could not parse data from file.");

		// stopping power is the second and third columns
		// (indices 2 and 3 to account for energy units)
		float dEdx = atof( line_elements[2].c_str() );
		dEdx += atof( line_elements[3].c_str() );

		// construct a new vector, and add it to the main vector<vector<float>>
		std::vector<float> new_data;
		new_data.push_back(Energy);
		new_data.push_back(dEdx);
		data.push_back(new_data);
	}

	return;
}

/**
 * Parse utility for the SRIM file's footer
 */
void StopPow_SRIM::parse_footer(std::stringstream& footer)
{
	// default initialization:
	scale_keV_um = scale_Mev_mgcm2 = 0;

	// loop over each line:
	std::string line;
	while(getline(footer,line))
	{
		// ignore separator lines, and the text header.
		// we search for things we don't want to find, and only
		// continue if they are not found:
		if( line.find("---") == std::string::npos 
			&& line.find("===") == std::string::npos
			&& line.find("Multiply") == std::string::npos 
			&& line.find("Ziegler") == std::string::npos )
		{
			// break up the line into elements
			// format is whitespace-separated values:
			std::stringstream linestream(line);
			std::string element;
			std::vector<std::string> line_elements;
			while( getline(linestream,element,WHITESPACE))
				if(element.size()>0)
					line_elements.push_back(element);

			// look for the line for keV / micron scale factor:
			if( line_elements[1].find("keV") != std::string::npos
				&& line_elements[3].find("micron") != std::string::npos )
				scale_keV_um = atof( line_elements[0].c_str() );
			// look for the line for MeV / (mg/cm2) scale factor:
			if( line_elements[1].find("MeV") != std::string::npos
				&& line_elements[3].find("mg/cm2") != std::string::npos )
				scale_Mev_mgcm2 = atof( line_elements[0].c_str() );
		}
	}

	// check if either scale factor is still zero
	// if it is, there was a problem and output cannot be trusted:
	if( scale_keV_um == 0 || scale_Mev_mgcm2 == 0)
		throw std::ios_base::failure("Could not read data from file.");
	return;
}

} // end namespace StopPow
