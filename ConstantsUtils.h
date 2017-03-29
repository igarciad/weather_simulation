#ifndef __CONSTANSTS_UTILS__ 
#define __CONSTANSTS_UTILS__

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>

///////////////////////////////////////
// WEATHER CONSTANTS

const float Kx = 500.; // diffusion coefficients
const float Ky = 500.; // diffusion coefficients
const float Kz = 100.;

const float Aw = 1.0f; // is there vertical advection
const float Ab = 1.0f; // is there buoyancy

const float Rd = 287.05f; // J kg-1 K-1
const float epsilon = 0.622f;

const float p_0 = 100000.0f;// PA
const float T0 = 273.15f;	// K
const float cpd = 1004.5f;	// PA

const float cvd = 717.5f;	// J kg-1 K-1
const float Llv = 2.501e6;	// J kg-1
const float Rv = 461.5f;	// J kg-1 K-1
const float cpv = 1850.0f;	// J kg-1 K-1
const float cvv = 1390.0f;	// J kg-1 K-1
const float cpL = 4186.;	// J kg-1 K-1

const float es_T0 = 610.7f; // PA
const float g = 9.81f;		// m s-2
const float phi = 0.785398163f;  // 45 degree in randians

const float ro_0 = 1.225f;  // kg m-3
const float cmax = 50.0f;   // anelastic speed of sound (300ms-1)

const float M_PI = 3.1415927410125732421875f;
const float DEG2RAD = M_PI / 180.0f;

///////////////////////////////////////////////////////
// INDEX TO VARIABLES

const int U = 0;		// u: wind component in the X direction
const int V = 1;		// v: wind component in the Y direction
const int W = 2;		// w: wind component in the Z direction
const int THETA = 3;	// Theta: Potential temperature
const int Pi = 4;		// Pi: Exener function
const int RO = 5;		// RO: Density
const int QV = 6;		// qv: Vapor mixing ratio
const int QC = 7;		// qc: Condensation mixing ratio
const int QR = 8;		// qr: Rain mixing ratio
const int VORT = 9;		// Vort: Vorticity

const int NUM_ELEM = 10;//NUMER OF ELEMENTS (last+1)


// INDEX TO GROUND VARIABLES

const int GR_TG = 0;		// Tg: Temperature of ground of first ds centimeters
const int GR_TA = 1;		// Ta: Temperayre air above z=0
const int GR_ALBEDO = 2;	// a: Albedo of the gridcell
const int GR_CGA = 3;		// Soildheat capacity per area of the gridcell

const int GR_TG_RESET = 4;	// Variables to reset after 24h
const int GR_TA_RESET = 5;
const int GR_TG_CORR = 6;
const int GR_TA_CORR = 7;

const int GR_BETA_INV = 8;		// Inverse of the Bowen ratio
const int GR_CLOUD_COVER = 9;	// Cloud coverage: Used in simulation and shadows

const int GR_NUM_ELEM = 10;

#include <stdio.h>
#ifdef __unix__ /* Just used to print the current path (useful to know where the system will search for the data files */
#include <unistd.h>
#define GetCurrentDir getcwd
#else
#include <direct.h>
#define GetCurrentDir _getcwd
#endif

class Utils {
	// We assume the file contains an initial line with the fields, followed by entries
	// The output is returned as a map of pairs field, value
public:

	/* Simple method to print where it is the current path. Used to see where to place the files. */
	static bool printCurrentDir() {
		char cCurrentPath[FILENAME_MAX];
		if (!GetCurrentDir(cCurrentPath, sizeof(cCurrentPath))) {
			return false;
		}
		printf("Current directory: %s\n", cCurrentPath);
		return true;
	}

	/* Read cvs file returning the values as an unsorted map of name and float values 
		@param fileName: Path of the cvs file to read
		@param cvsOutput: unordered map that contains pairs of field and value (the fields are extracted from the first line of text)
	*/
	static bool readCVSFile(std::string fileName, std::vector<std::unordered_map<std::string, float > >& cvsOutput) {
		std::ifstream  data(fileName);
		if (!data.good()) {
			printCurrentDir();
			printf("ERROR: reading file %s\n", fileName.c_str());
			return false;
		}
		std::string line;
		if (!std::getline(data, line)) {
			cvsOutput.clear();
			return false;
		}
		// extract header
		printf("readCVSFile: %s\n", line.c_str());
		std::vector<std::string> header;
		std::stringstream  lineStream(line);
		std::string cell;
		while (std::getline(lineStream, cell, ',')) {
			header.push_back(cell);
			printf("readCVSFile:   %s\n", cell.c_str());
		}
		// extract each line
		while (std::getline(data, line)) {
			// each line contains one entry
			std::unordered_map<std::string, float>  entry;

			std::stringstream  lineStream(line);
			std::string cell;
			int i = 0;
			while (std::getline(lineStream, cell, ',')) {
				// insert each pair of header to value
				if (i >= header.size()) {
					printf("ERROR: readCVSFile %s\n", line.c_str());
					cvsOutput.clear();
					return false;
				}
				try {
					entry.insert(std::pair<std::string, float>(header[i], std::stof(cell)));
				}
				catch (...) {};  // skip not float
				i++;
			}
			cvsOutput.push_back(entry);
		}
		return true;
	}
};

#endif  // __CONSTANSTS_UTILS__
