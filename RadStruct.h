#ifndef __RADSTRUCT__ 
#define __RADSTRUCT__

/* RadStruct: Structure that contains information for the radiation model.*/
struct RadStruct {
	float initDayInYearUTC;  // Initial day in the year (0-366.0f)
	float initTimeUTC_hours; // Time in UTC (0-24.0f)
	float timeZone;			 // Time zone (-12.0f-12.0f)
	float longitudeRad;		 // Longitude of the simulation
	float latitudeRad;		 // Laltitude of the simulation
	float rainProbability;	 // Probability of rain

	/* Basic contructor to inialize variables. */
	RadStruct() {
		initDayInYearUTC = 100.0f; 
		initTimeUTC_hours = 12.0f; 
		timeZone = -6.0f; 
		longitudeRad = 1.564f; 
		latitudeRad = 0.784f;  
		rainProbability=1.0f;//book values
	}
	RadStruct(float _initDayInYearUTC, float _initTimeUTC_hours, float _timeZone, float _longitudeRad, float _latitudeRad, float _rainProbability) {
		initDayInYearUTC = _initDayInYearUTC; 
		initTimeUTC_hours = _initTimeUTC_hours; 
		timeZone = _timeZone; 
		longitudeRad = _longitudeRad; 
		latitudeRad = _latitudeRad;
		rainProbability = _rainProbability;
	}
};

#endif  // __RADSTRUCT__