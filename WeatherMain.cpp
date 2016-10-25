// WeatherSimple.cpp : Defines the entry point for the console application.
//
#include "Weather3D.h"
#include <stdio.h>
#include <algorithm> // for std::min/max

int main() {
	printf("main start\n");
	Utils::printCurrentDir();

	float dT = 1.0f;
	float gridXSize = 1000.0f; // size of the grid
	float gridYSize = 1000.0f;
	float side = 10000.f; // world will be side x side meters

	printf("////////////////////////////\nINIT\n////////////////////////////\n");

	// GridX GridY
	int gridX = int(side / gridXSize);
	int gridY = int(side / gridYSize);

	// GridZ
	std::vector<float> gridSizeK;
	std::vector<float> gridSizeKAcc;
	int gridZ = 56;
	float zLEVEL_sr = 1.025f;
	float zLEVEL_dz = 200.0f;
	gridSizeK.resize(gridZ + 1);
	gridSizeKAcc.resize(gridZ + 1);
	gridSizeK[0] = zLEVEL_dz / 2.0f;
	gridSizeKAcc[0] = zLEVEL_dz / 2.0f;
	for (int k = 1; k < gridZ + 1; k++) {
		gridSizeK[k] = std::min<float>(zLEVEL_dz*pow(zLEVEL_sr, k), 1000.0f);
		gridSizeKAcc[k] = gridSizeKAcc[k - 1] + std::min<float>(zLEVEL_dz*pow(zLEVEL_sr, k), 1000.0);
	}

	// Radiation parameters (book values)
	float initDayInYearUTC = 100.0f; // out of 355.5f
	float initTimeUTC_hours = 12.0f; // 0-24h
	float timeZone = -6.0f;
	float latitudeRad = 37.0f*DEG2RAD; // in radians 
	float longitudeRad = -122.0f*DEG2RAD; // in radians
	float rainProbability = 1.0f;
	RadStruct radStruct(initDayInYearUTC, initTimeUTC_hours, timeZone, longitudeRad, latitudeRad, rainProbability);


	// Load Ground
	Ground3D ground;
	if (!ground.initRandomlyGroundVariables(gridX, gridY)) {
		printf("ERROR: initRandomlyGroundVariables\n");
		return -1;
	}

	// Load sounding as initial state
	Grid3D grid0Var;
	Grid3D gridInitVar;
	Grid3D::initGrid3D(grid0Var, gridInitVar, gridX, gridY, gridZ, gridSizeKAcc);  // reads

	// Init simulation
	initSimulation(dT, gridX, gridY, gridZ, gridXSize, gridYSize, gridSizeK, grid0Var, gridInitVar, ground, radStruct);

	printf("////////////////////////////\nSIMULATE\n////////////////////////////\n");

	float simulationTime = 0.0f;
	const int numIterPerStep = 100;
	const int numIter = 10;
	Grid3D gridCurrState; Ground3D groundCurrState;
	for (int i = 0; i < numIter; i++) {
		printf("////////////////////////////\nSIMULATE STEP %d\n", i*numIterPerStep);
		simulateStep(numIterPerStep, simulationTime);
		// Save to hdd (see function to define what to save)
		saveToHDDCurrentState();
	}

	return 0;
}//

