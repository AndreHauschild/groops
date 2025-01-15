/***********************************************/
/**
* @file gnssTroposphere2GriddedDataTimeSeries.cpp
*
* @brief Evaluate a troposphere class by computing slant delays and mapping function values.
*
* @author Andreas Kvas
* @date 2024-06-21
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Evaluate a troposphere class by computing slant delays and mapping function values.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileGriddedDataTimeSeries.h"
#include "classes/troposphere/troposphere.h"
#include "classes/timeSeries/timeSeries.h"

/***** CLASS ***********************************/

/** @brief Evaluate a troposphere class by computing slant delays and mapping function values.
* Evaluate a troposphere class by computing slant delays and mapping function values.
* @ingroup programsGroup */
class GnssTroposphere2GriddedDataTimeSeries
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GnssTroposphere2GriddedDataTimeSeries, SINGLEPROCESS, "Evaluate a troposphere class by computing slant delays and mapping function values.", Gnss)

/***********************************************/

void GnssTroposphere2GriddedDataTimeSeries::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName fileNameOut;
    TropospherePtr troposphere;

    Angle elevationMin, elevationMax, elevationSampling;
    Angle azimuthMin, azimuthMax, azimuthSampling;
    TimeSeriesPtr timeSeries;

    Angle stationLongitude, stationLatitude;
    Double stationHeight;
    Double a, f;

    readConfig(config, "outputfileGriddedDataTimeSeries", fileNameOut,  Config::MUSTSET,  "", "output file name for the gridded data time series");
    readConfig(config, "troposphere",                     troposphere,  Config::MUSTSET,  "", "troposphere class to be evaluated");
    readConfig(config, "timeSeries",                      timeSeries,   Config::MUSTSET,  "", "time series at which the troposphere class is evaluated");

    readConfig(config, "stationLongitude",                stationLongitude,   Config::MUSTSET,  "", "station longitude in degrees");
    readConfig(config, "stationLatitude",                 stationLatitude,   Config::MUSTSET,  "", "station latitude in degrees");
    readConfig(config, "stationHeight",                   stationHeight,   Config::MUSTSET,  "", "station height in meters");


    readConfig(config, "elevationMin",                    elevationMin, Config::DEFAULT,  "0.5", "minimum elevation (degrees) at which the troposphere class is evaluated");
    readConfig(config, "elevationMax",                    elevationMax, Config::DEFAULT,  "89.5", "maximum elevation (degrees) at which the troposphere class is evaluated");
    readConfig(config, "elevationSampling",               elevationSampling, Config::DEFAULT,  "1", "elevation sampling (degrees)");
    readConfig(config, "azimuthMin",                      azimuthMin, Config::DEFAULT,  "-179.5", "minimum azimuth (degrees) at which the troposphere class is evaluated");
    readConfig(config, "azimuthMax",                      azimuthMax, Config::DEFAULT,  "179.5", "maximum azimuth (degrees) at which the troposphere class is evaluated");
    readConfig(config, "azimuthSampling",                 azimuthSampling, Config::DEFAULT,  "1", "elevation sampling (degrees)");

    readConfig(config, "R",                               a,          Config::DEFAULT, STRING_DEFAULT_GRS80_a, "reference radius for ellipsoidal coordinates");
    readConfig(config, "inverseFlattening",               f,          Config::DEFAULT, STRING_DEFAULT_GRS80_f, "reference flattening for ellipsoidal coordinates");
    if(isCreateSchema(config)) return;


    // create station position and init troposphere
    Ellipsoid trfEll(a, f);
    std::vector<Vector3d> stationPositions(1, trfEll(stationLongitude, stationLatitude, stationHeight));
    troposphere->init(stationPositions);

    // set up evaluation grid (time, elevation, azimuth)
    std::vector<Time> times = timeSeries->times();

    std::vector<Angle> elevation(1, elevationMin);
    while(elevation.back() <= elevationMax)
      elevation.push_back(elevation.back() + elevationSampling);
    std::reverse(elevation.begin(), elevation.end());

    std::vector<Angle> azimuth(1, azimuthMin);
    while(azimuth.back() <= azimuthMax)
      azimuth.push_back(azimuth.back() + azimuthSampling);

    // evaluate troposphere and create gridded data time series
    GriddedData grid;
    grid.ellipsoid = trfEll;
    std::vector<Matrix> data;
    for(const Time& time : times)
    {
      UInt stationId = 0;
      Matrix newData(elevation.size() * azimuth.size(), 5, NAN_EXPR);
      UInt pointIdx = 0;

      for(const Angle& el : elevation)
      {
        for(const Angle& az : azimuth)
        {
          newData(pointIdx, 0) = troposphere->slantDelay(time, stationId, az, el);
          newData(pointIdx, 1) = troposphere->mappingFunctionWet(time, stationId, az, el);
          newData(pointIdx, 2) = troposphere->mappingFunctionHydrostatic(time, stationId, az, el);
          troposphere->mappingFunctionGradient(time, stationId, az, el, newData(pointIdx, 3), newData(pointIdx, 4));

          if(grid.points.size() < newData.rows())
          {
            Vector3d p = grid.ellipsoid(az, el, 0);
            grid.points.push_back(p);
          }

          pointIdx++;
        }
      }
      data.push_back(newData);
    }

    // write data to file
    writeFileGriddedDataTimeSeries(fileNameOut, 1, times, grid, data);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
