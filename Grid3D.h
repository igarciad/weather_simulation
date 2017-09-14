#ifndef __GRID3D__ 
#define __GRID3D__

#include "ConstantsUtils.h"
#include <algorithm>
#include <cmath>

/* Grid3D: Structure that contains all variables and applies the toroidal behavior */
struct Grid3D {
	std::vector<std::vector<float>> var; // vector of vectors that contains the variable values [VARIABLE][INDEX]
	std::vector<float> vecZero; // zero array to check wrong acces

	int gridX;  // number of gridcells in X direction
	int gridY;  // number of gridcells in Y direction
	int gridZ;  // number of gridcells in Z direction

	/* Operators to acces the data in a toroidal manner */
	float operator()(const int arr, const int x, const int y, const int z) const {
		int aI = x;  // Toroidal behavior: (x<0)-->(x=gridX-1-x) and (x>gridX-1)-->(x=x-gridX-1)
		if (aI < 0) aI = gridX + x;
		if (aI >= gridX)aI = (x%gridX);
		int aJ = y;  // Toroidal
		if (aJ < 0)aJ = gridY + y;
		if (aJ >= gridY)aJ = (y%gridY);
		int aZ = z;
		int gridK = gridZ;
		if (arr == W)
			gridK++;
		if ((z < 0 || z >= gridK) && (arr <= 2 || arr == 8)) {//for values (under z<0, or higher than the gridK) AND they are arr<=2 (u,v,w) or arr==8 (rL)
			return 0.0f;
		}
		aZ = std::max(0, aZ);//not under
		aZ = std::min(aZ, gridK - 2);//not over top
		return var[arr][aI + aJ*gridX + aZ*(gridX*gridY)];
	}

	/* Operators to acces the data in a toroidal manner */
	float &operator ()(const int arr, const int x, const int y, const int z) {
		int aI = x;  // Toroidal behavior: (x<0)-->(x=gridX-1-x) and (x>gridX-1)-->(x=x-gridX-1)
		if (aI < 0) aI = gridX + x;
		if (aI >= gridX)aI = (x%gridX);
		int aJ = y;  // Toroidal
		if (aJ < 0)aJ = gridY + y;
		if (aJ >= gridY)aJ = (y%gridY);
		int aZ = z;
		int gridK = gridZ;
		if (arr == W)
			gridK++;
		if ((z < 0 || z >= gridK) && (arr <= 2 || arr == 8)) {//for values (under z<0, or higher than the gridK) AND they are arr<=2 (u,v,w) or arr==8 (rL)
			if (vecZero[0] != 0) {//check to make sure we are not trying to change a value out of scope
				printf("Catastrofic change of vecZero: CAT arr %d i %d j %d k %d vecZero[0] %f\n", arr, x, y, z, vecZero[0]);
			}
			return vecZero[0];// return 0;
		}
		aZ = std::max(0, aZ);//not under
		aZ = std::min(aZ, gridK - 2);//not over top
		return var[arr][aI + aJ*gridX + aZ*(gridX*gridY)];
	}//

	/* Reset to zero a variable in the array 
		@param arr: Variable to be set to zero
		*/
	void resetVariable(const int arr) {
		std::fill(var[arr].begin(), var[arr].end(), 0.0f);
	}

	/* Initialize a Grid3D. It reads the file: "initValues3D.csv" that contains a sounding (see file as example)
		@param _gridI: Number of gridcells in the X direction (Note: The notation changed from X to be consisten with weather terminilogy)
		@param _gridJ: Number of gridcells in the Y direction
		@param _gridK: Number of gridcells in the Z direction
		@param gridSizeKAcc: It is the size (in meters) of each gridcell accumulated in Z direction (Note: this contains the accumulated values, not the heights)
	*/
	void init(const int _gridI, const int _gridJ, const int _gridK, std::vector<float> gridSizeKAcc) {
		printf("Init Grid3D %d %d %d\n", _gridI, _gridJ, _gridK);
		gridX = _gridI;
		gridY = _gridJ;
		gridZ = _gridK;

		const bool saveToFile = true;

		// clear files
		if (saveToFile == true) {
			std::ofstream oU, oU2;
			oU.open("OUT.csv", std::fstream::out | std::fstream::trunc);//trunc to clear file
			oU.close();
			oU2.open("initState.csv", std::fstream::out | std::fstream::trunc);//trunc to clear file
			oU2.close();
		}

		// 0. init size
		var.resize(NUM_ELEM);
		for (int v = 0; v < var.size(); v++) {
			var[v].resize(gridX*gridY*gridZ);
		}
		// 1. read raw data file
		std::ifstream input("initValues3D.csv");

		if (!input) {
			printf("Can't open file initValues3D.csv\n");
			exit(0);
		}
		std::vector<std::vector<float> > valuesT;
		float groundHeight = -1.0f;
		for (std::string line; getline(input, line);) {
			std::vector<float> values;
			std::stringstream ss(line);
			float i;
			while (ss >> i) {
				values.push_back(i);
				if (ss.peek() == ';')
					ss.ignore();
			}
			if (groundHeight == -1.0f) {
				groundHeight = values[1];//init height
				values[1] = 0;
			} else
				values[1] -= groundHeight;
			valuesT.push_back(values);
		}
		input.close();
		// 2. for each Z find the values
		std::vector<std::vector<float> > initValues;//init values for all Z
		initValues.resize(gridZ);//as many as heights
		std::vector<float> val;//init values for one Z
		val.resize(NUM_ELEM);
		std::vector<float> interpolateValues;
		int indexInArray1 = 0;
		int indexInArray2;
		for (int k = 0; k < gridZ; k++) {
			float currZ;
			if (k>0)
				currZ = (gridSizeKAcc[k - 1] + gridSizeKAcc[k]) / 2.0f;//midpoint
			else
				currZ = gridSizeKAcc[k] / 2.0f;//first just half
			// 2.1 find in raw data the closest elements
			while (indexInArray1 < valuesT.size() && valuesT[indexInArray1][1] < currZ) {
				indexInArray1++;
			}
			indexInArray1--;
			indexInArray2 = indexInArray1 + 1;
			if (indexInArray1 == -1) {
				interpolateValues = valuesT[indexInArray2];//first element
				indexInArray1 = 0;
			} else {

				if (indexInArray2 == valuesT.size()) {
					interpolateValues = valuesT[indexInArray1];//last element
				} else {
					// interpolate arrays
					interpolateValues.resize(valuesT[indexInArray1].size());
					float z1 = valuesT[indexInArray1][1];
					float z2 = valuesT[indexInArray2][1];

					float fact1 = 1.0f - ((currZ - z1) / (z2 - z1));
					float fact2 = 1.0f - ((z2 - currZ) / (z2 - z1));
					if (valuesT[indexInArray1][6] > 180.0f)valuesT[indexInArray1][6] -= 360.0;//wind make to be in the first cuadrant
					if (valuesT[indexInArray2][6] > 180.0f)valuesT[indexInArray2][6] -= 360.0;
					for (int inV = 0; inV < valuesT[indexInArray1].size(); inV++) {
						interpolateValues[inV] = fact1*valuesT[indexInArray1][inV] + fact2*valuesT[indexInArray2][inV];
					}
				}
			}
			// 2.2 set initial values
			val[U] = interpolateValues[7] * cos(interpolateValues[6] * 0.0174532925f)*0.514444444f;//wind speed*angle*knts/m/s
			val[V] = interpolateValues[7] * sin(interpolateValues[6] * 0.0174532925f)*0.514444444f;//wind speed*angle*knts/m/s
			val[W] = 0;
			val[THETA] = interpolateValues[10];//theta v
			float P = interpolateValues[0] * 100.0f;//preasure (in hPa)
			val[Pi] = pow((P / p_0), Rd / cpd);

			val[QV] = interpolateValues[5] * 0.001f;//rv kg/kg (0.001f: g->kg)
			val[QC] = 0;//rl g/kg

			float ro = P / ((Rd*(interpolateValues[2] + 273.15f))*(1 + val[QV] / epsilon));
			val[RO] = ro;//

			initValues[k] = val;//it could be done directly here
			// 2.2 set values in layer
			// init
			for (int v = 0; v < NUM_ELEM; v++) {
				std::fill(var[v].begin() + (k*gridX*gridY), var[v].begin() + ((k + 1)*gridX*gridY), val[v]);
			}
			//save to file
			if (saveToFile == true) {
				std::ofstream oU;
				oU.open("initState.csv", std::fstream::out | std::fstream::app);
				printf("    CurrZ %f\t", currZ);
				oU << currZ;
				for (int v = 0; v < NUM_ELEM; v++) {
					printf("%f\t", val[v]);//c_str()
					oU << "," << val[v];
				}
				oU << "\n";
				printf("\n");
			}//saveToFile
		}
	}//

	/* Basic function to read the file. Remove extra spaces (In the actual code we use regular expresions */
	static inline bool removeSpaces(char lhs, char rhs) {
		return (lhs == rhs) && (lhs == ' ');
	}

	/* Basic function to read the file. Remove extra spaces (In the actual code we use regular expresions */
	static inline std::string trim(std::string& str) {
		size_t first = str.find_first_not_of(' ');
		size_t last = str.find_last_not_of(' ');
		return str.substr(first, (last - first + 1));
	}

	/* Basic function to read the file. Remove extra spaces (In the actual code we use regular expresions */
	static inline std::string cleanString(std::string& in) {
		std::string::iterator new_end = std::unique(in.begin(), in.end(), removeSpaces); // remove multiple whitespaces
		in.erase(new_end, in.end());
		std::transform(in.begin(), in.end(), in.begin(), ::tolower); // to lower case
		return trim(in); // trim
	}

	/* Initialize a Grid3D. It reads the file: "sounding_XXX.csv" that contains a sounding (see file as example)
		@param soundingNumber: Number of the sounding to read, it defines the number as "sounding_<soundingNumber>.csv"
		@param _gridX: Number of gridcells in the X direction
		@param _gridY: Number of gridcells in the Y direction
		@param _gridZ: Number of gridcells in the Z direction
		@param gridSizeKAcc: It is the size (in meters) of each gridcell accumulated in Z direction (Note: this contains the accumulated values, not the heights)
	*/
	bool initFromSounding(const int soundingNumber, const int _gridX, const int _gridY, const int _gridZ, std::vector<float> gridSizeKAcc) {
		printf("Init Grid3D %d %d %d From Sounding %d\n", _gridX, _gridY, _gridZ, soundingNumber);
		gridX = _gridX;
		gridY = _gridY;
		gridZ = _gridZ;

		// 0. init size
		var.resize(NUM_ELEM);
		for (int v = 0; v < var.size(); v++) {
			var[v].resize(gridX*gridY*gridZ);
		}

		const bool saveComputedSoundingToFile = true;
		std::string OUTSOUNDING3 = "computed_sounding.csv";//file to see the computed sounding

		//////////////////////////////////////
		// 1. read raw data file
		std::string fileName = "sounding_" + std::to_string(soundingNumber) + ".txt";
		std::ifstream input(fileName);

		if (!input) {
			printf("Can't open file %s\n", fileName.c_str());
			return false;
		}
		// READ Header
		std::vector<std::string> fields;
		char splitV = ' ';//find if it is ; or , used as CVS
		std::string line;

		std::getline(input, line);// skip ----
		{
			std::getline(input, line);
			line = cleanString(line);
			printf("Fields Line: %s\n", line.c_str());
			std::stringstream  lineStream(line);
			std::string cell;
			while (std::getline(lineStream, cell, splitV)) {
				fields.push_back(cell);
				// printf("field: <%s>\n", cell.c_str());
			}
		}
		getline(input, line);// skip  units
		getline(input, line);// skip ----

		std::vector<std::unordered_map<std::string, float>> valuesT;
		float groundHeight = -1.0f;
		for (std::string line; getline(input, line);) {
			line = cleanString(line);
			// printf("Line: <%s>\n", line.c_str());
			std::unordered_map<std::string, float> values;
			std::stringstream  lineStream(line);
			std::string cell;
			int i = 0;
			while (std::getline(lineStream, cell, splitV)) {
				if (i >= fields.size()) {
					printf("SKIPPED: line %s\n", line.c_str());
					continue;
				}
				try {
					values.insert(std::pair<std::string, float>(fields[i], std::stof(cell)));
					// printf("Insert %s value %f\n", fields[i], std::stof(cell));
					i++;
				}
				catch (...) {
					printf("SKIPPED: reading line %s\n", line.c_str());
					continue;
				};//skip not float
			}

			if (groundHeight == -1.0f) {//remove z0
				groundHeight = values["hght"];
				values["hght"] = 0;
			} else {
				values["hght"] -= groundHeight;
			}
			valuesT.push_back(values);
		}
		input.close();
		///////////////////////////////////////////////
		// 2. for each Z find the values
		std::vector<std::vector<float> > initValues;//init values for all Z
		initValues.resize(gridZ);//as many as heights
		std::vector<float> val;//init values for one Z
		val.resize(NUM_ELEM);
		std::unordered_map<std::string, float> interpolateValues;
		int indexInArray1 = 0;
		int indexInArray2;
		for (int k = 0; k < gridZ; k++) {
			float currZ;
			if (k>0)
				currZ = (gridSizeKAcc[k - 1] + gridSizeKAcc[k]) / 2.0f;//midpoint
			else
				currZ = gridSizeKAcc[k] / 2.0f;//first just half
			// 2.1 find in raw data the closest elements
			while (indexInArray1 < valuesT.size() && valuesT[indexInArray1]["hght"] < currZ) {
				indexInArray1++;
			}
			indexInArray1--;
			indexInArray2 = indexInArray1 + 1;
			if (indexInArray1 == -1) {
				interpolateValues = valuesT[indexInArray2];//first element
				indexInArray1 = 0;
			} else {

				if (indexInArray2 == valuesT.size()) {
					interpolateValues = valuesT[indexInArray1];//last element
				} else {
					// interpolate arrays
					float z1 = valuesT[indexInArray1]["hght"];
					float z2 = valuesT[indexInArray2]["hght"];

					float fact1 = 1.0f - ((currZ - z1) / (z2 - z1));
					float fact2 = 1.0f - ((z2 - currZ) / (z2 - z1));
					if (valuesT[indexInArray1]["drct"]>180.0f)valuesT[indexInArray1]["drct"] -= 360.0;//wind make to be in the first cuadrant
					if (valuesT[indexInArray2]["drct"] > 180.0f)valuesT[indexInArray2]["drct"] -= 360.0;
					for (int fN = 0; fN < fields.size(); fN++) {
						interpolateValues[fields[fN]] = fact1*valuesT[indexInArray1][fields[fN]] + fact2*valuesT[indexInArray2][fields[fN]];
					}
				}
			}
			// 2.2 set initial values

			val[U] = interpolateValues["sknt"] * cos(interpolateValues["drct"] * 0.0174532925f)*0.514444444f*0.1f;//wind speed*angle*knts/m/s
			val[V] = interpolateValues["sknt"] * sin(interpolateValues["drct"] * 0.0174532925f)*0.514444444f*0.1f;//wind speed*angle*knts/m/s
			if (k == 0 || k == 1) {
				val[U] = val[V] = 0.0f;
			}
			val[W] = 0;// interpolateValues["w"];
			val[THETA] = interpolateValues["thta"];//theta v
			float P = interpolateValues["pres"] * 100.0f;//preasure (in hPa)
			val[Pi] = std::pow((P / p_0), Rd / cpd);

			val[QV] = interpolateValues["mixr"] * 0.001f;// (interpolateValues["mixr"] / (1.0f + interpolateValues["mixr"])) * 0.001f;//qv=(rv/1+rv) kg/kg (0.001f: g->kg)

			val[RO] = P / ((Rd*(interpolateValues["temp"] + 273.15f))*(1 + val[QV] / epsilon));

			initValues[k] = val;//it could be done directly here
			// 2.2 set values in layer
			// init
			for (int v = 0; v < NUM_ELEM; v++) {
				std::fill(var[v].begin() + (k*gridX*gridY), var[v].begin() + ((k + 1)*gridX*gridY), val[v]);
			}
			//save to file
			if (saveComputedSoundingToFile == true) {
				std::ofstream oU;
				oU.open(OUTSOUNDING3, std::fstream::out | std::fstream::app);
				printf("    CurrZ %f\t", currZ);
				oU << currZ;
				for (int v = 0; v < NUM_ELEM; v++) {
					printf(" %f\t", val[v]);//c_str()
					oU << "," << val[v];
				}
				oU << "\n";
				printf("\n");
				oU.close();
			}
		}
		// init no-initlialized variables to zero
		resetVariable(W);
		resetVariable(QC);
		resetVariable(QR);
		resetVariable(VORT);
		printf("-- -- LOAD sounding %d\n", soundingNumber);
		return true;
	}//

	/* Basic funtion to initialize all the variable needed for the simulation. It will initialize the variables using the "sounding_0.cvs"
		@param grid0Var: Output Grid3D where the values are initialized (the perturbation variables will be set to zero)
		@param gridInitVar: Output Grid3D where the values are initialized
		@param _gridX: Number of gridcells in the X direction
		@param _gridY: Number of gridcells in the Y direction
		@param _gridZ: Number of gridcells in the Z direction
		@param gridSizeKAcc: It is the size (in meters) of each gridcell accumulated in Z direction (Note: this contains the accumulated values, not the heights)
	*/
	static bool initGrid3D(Grid3D& grid0Var, Grid3D& gridInitVar, const int _gridX, const int _gridY, const int _gridZ, const std::vector<float>& _gridSizeKAcc) {
		std::vector<float> vecZero(1, 0.0f);
		grid0Var.vecZero = vecZero;

		grid0Var.gridX = _gridX;
		grid0Var.gridY = _gridY;
		grid0Var.gridZ = _gridZ;

		// Load sounding
		grid0Var.initFromSounding(0, _gridX, _gridY, _gridZ, _gridSizeKAcc);
		gridInitVar = grid0Var;
		// Set variables to zero when they will use their prime form
		grid0Var.resetVariable(QV);		//it is QV prime
		grid0Var.resetVariable(Pi);		//it is Pi prime
		grid0Var.resetVariable(THETA);	//it is theta prime
		return true;
	}

};

#endif  // __GRID3D__
