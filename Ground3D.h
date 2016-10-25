#ifndef __GROUND3D__ 
#define __GROUND3D__

#include "ConstantsUtils.h"

/* Ground3D: Structure that contains the variables of each bottomost layer of the ground used in radiation model */
struct Ground3D {
public:

	std::vector<std::vector<float>> var;  // vector of vectors that contains the variable values [GROUND_VARIABLE][INDEX]
	int gridX;  // number of gridcells in X direction
	int gridY;  // number of gridcells in Y direction

	Ground3D() {}

	/* Basic constructor. */
	Ground3D(const int _gridX, const int _gridY, const std::vector<std::vector<float>>& _var) {
		gridX = _gridX;
		gridY = _gridY;
		var = _var;
	}

	/* Operators to acces the data in a toroidal manner */
	float operator()(const int arr, const int x, const int y) const {
		int aY = y;  // toroidal behavior
		if (y < 0)aY = gridY + y;
		if (y >= gridY)aY = (y % gridY);

		int aX = x;  // toroidal behavior
		if (x < 0)aX = gridX + x;
		if (x >= gridX)aX = (x % gridX);
		return var[arr][aX + aY*gridY];
	}

	/* Operators to acces the data in a toroidal manner */
	float &operator ()(const int arr, const int x, const int y) {
		int aY = y;  // toroidal behavior
		if (y < 0)aY = gridY + y;
		if (y >= gridY)aY = (y % gridY);

		int aX = x;  // toroidal behavior
		if (x < 0)aX = gridX + x;
		if (x >= gridX)aX = (x % gridX);

		return var[arr][aX + aY*gridY];//return that value
	}//

	/* Initialize the ground with random land use variable values. In the main
		project we load a png file that contains the ditribution of each land
		use type for each gridcell and interpolate the values. Here, we randomly
		asign land use to each gridcell and asign the values of the constants from
		a cvs file.
		*/
	int randBetween(const int start, const int end) {
		return (rand() % (end - start)) + start;
	}

	bool initRandomlyGroundVariables(const int _gridX, const int _gridY) {
		printf("initRandomlyGroundVariables %d %d\n", _gridX, _gridY);
		gridX = _gridX;
		gridY = _gridY;
		// resize to contain all variables
		var.resize(GR_NUM_ELEM);
		for (int i = 0; i < GR_NUM_ELEM; i++) {
			var[i].resize(_gridX*_gridY);
		}
		printf("readCVSFile\n");
		std::vector<std::unordered_map<std::string, float>> cvsOutput;
		if (!Utils::readCVSFile("landuse_constants.csv", cvsOutput)) {
			printf("ERROR: landuse_constants\n");
			return false;
		}
		printf("Select random landuse and set constants\n");
		for (int i = 0; i < gridX*gridY; i++) {
			int randLandUse = randBetween(0, (int)cvsOutput.size());
			var[GR_ALBEDO][i] = cvsOutput[randLandUse]["ALBEDO"];
			var[GR_CGA][i] = cvsOutput[randLandUse]["CGA"];
			var[GR_BETA_INV][i] = 1.0f / cvsOutput[randLandUse]["BETA"];
			printf("grid[%d]-> alb %f cga %f beta_inv %f\n", i, var[GR_ALBEDO][i], var[GR_CGA][i], var[GR_BETA_INV][i]);
		}
		return true;
	}

};

#endif  // __GROUND3D__

