#ifndef __WEATHER3D__ 
#define __WEATHER3D__

#include <vector>
#include <string>
#include "ConstantsUtils.h"
#include "Ground3D.h"
#include "Grid3D.h"
#include "RadStruct.h"


/* Initialize the simulation. This method is used to set all the necessary initial variables to start the simulation.
	@param _dT: It is the delta time, time step, usually defined to 1.0f seconds
	@param _gridX: It is the number of gridcells in the X direction
	@param _gridY: It is the number of gridcells in the Y direction
	@param _gridZ: It is the number of gridcells in the Z direction
	@param _gridSizeX: It is the size (in meters) of each gridcell in X direction
	@param _gridSizeY: It is the size (in meters) of each gridcell in Y direction
	@param _gridSizeZ: It is the size (in meters) of each gridcell in Z direction (Note: It has more accuracy close to the surface)
	@param _grid0Var: It is the initial grid values (usually load from a sounding) with some variables computed as perturbations (so set to zero)
	@param _gridInitVar: It is the initial grid values (usually load from a sounding)
	@param _ground: It is the initial ground variables and values (usually load from file with land use)
	@param _radStruct: Structure that contains info for the radiation model (e.g., latitude and longitude)
*/
void initSimulation(const float _dT, const int _gridX, const int _gridY, const int _gridZ, const float _gridSizeI, const float _gridSizeJ, std::vector<float>& _gridSizeK,
	Grid3D& _grid0Var, Grid3D& _gridInitVar,
	Ground3D& _ground, RadStruct& _radStruct);

/* Run the simulation. Define the number of steps to run, and the simulation time (usually starts at 0.0f).
   This methods calls the four simulation steps and updates all simulation variables. (Note: initSimulation method shoyld be called fist)
	@param numSteps: Number of steps to run in the simulation
	@param simulationTime: float that contains the current sumulationTime
*/
void simulateStep(const int numSteps, float& simulationTime);

/* Copy to the structs the current state of the simulation 
	@param _currGrid: Struct where to copy the current Grid3D
	@param _currGround: Struct where to copythe current Ground3D
*/
void copyCurrentState(Grid3D& _currGrid, Ground3D& _currGround);

/* Saves to hard drive the current state of the selected variables at different heights to be ploted */
void saveToHDDCurrentState();

#endif  // __WEATHER3D__