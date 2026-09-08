/***********************************************/
/**
* @file simulateStarCameraMetopSG.cpp
*
* @brief Simulate star camera data for Metop-SG.
*
* @author Andre Hauschild
* @date 2026-02-10
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program simulates \file{star camera}{instrument} measurements at each satellite's position.
The satellite's orientation follows a local orbit frame with the z-axis pointing downwards normal 
to the Earth's ellipsoid, the y-axis perpendicular to the Earth-fixed velocity and the x-axis 
oriented in the direction of the Earth-fixed velocity.
The resulting rotation matrices rotate from satellite frame to inertial frame.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"
#include "classes/earthRotation/earthRotation.h"

/***** CLASS ***********************************/

/** @brief Simulate star camera data for Metop-SG.
* @ingroup programsGroup */
class SimulateStarCameraMetopSG
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(SimulateStarCameraMetopSG, PARALLEL, "Simulate star camera data for Metop-SG.", Simulation, Instrument)

/***********************************************/

void SimulateStarCameraMetopSG::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName         orbitName, starCameraName;
    EarthRotationPtr earthRotation;
    Bool     isNadirPointing;

    readConfig(config, "outputfileStarCamera", starCameraName,  Config::MUSTSET, "",  "rotation from satellite to inertial frame (x: along Earth-fixed velocity, y: cross, z: normal to Earth ellipsoid)");
    readConfig(config, "inputfileOrbit",       orbitName,       Config::MUSTSET, "",  "position and velocity defines the orientation of the satellite at each epoch");
    readConfig(config, "earthRotation",        earthRotation,   Config::MUSTSET,  "", "transformation from CRF to TRF");

    if(isCreateSchema(config)) return;

    logStatus<<"read orbit and generate star camera data"<<Log::endl;
    InstrumentFile  orbitFile(orbitName);
    std::vector<Arc> arcList(orbitFile.arcCount());

    Parallel::forEach(arcList, [&](UInt arcNo)
    {
      OrbitArc orbit = orbitFile.readArc(arcNo);
      StarCameraArc arc;
      for(UInt i=0; i<orbit.size(); i++)
      {

        // CRF position and velocity
        Vector3d position = orbit.at(i).position;
        Vector3d velocity = orbit.at(i).velocity; // velocity vector
        if(velocity.r()==0)
        {
          if(i<orbit.size()-1)
            velocity = orbit.at(i+1).position - orbit.at(i).position;
          else
            velocity = orbit.at(i).position - orbit.at(i-1).position;
        }

        // CRF to TRF rotation
        const Rotary3d rot = earthRotation->rotaryMatrix(orbit.at(i).time);
        const Vector3d omega = rot.rotate(earthRotation->rotaryAxis(orbit.at(i).time));

        // TRF position and velocity
        position = rot.rotate(orbit.at(i).position);
        velocity = rot.rotate(orbit.at(i).velocity) - crossProduct(omega, position);

        position = localNorthEastUp(position, Ellipsoid()).transform(Vector3d(0,0,1)); // up component

        Vector3d y = normalize(crossProduct(velocity, position)); // cross
        Vector3d x = normalize(crossProduct(position, y));        // along

        StarCameraEpoch epoch;
        epoch.time   = orbit.at(i).time;
        epoch.rotary = inverse(rot)*Rotary3d(x, y); // rotation crf <- trf * trf <- sat
        arc.push_back(epoch);
      }
      return arc;
    }, comm);

    if(Parallel::isMaster(comm))
    {
      logStatus<<"write star camera data to file <"<<starCameraName<<">"<<Log::endl;
      InstrumentFile::write(starCameraName, arcList);
      Arc::printStatistics(arcList);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
