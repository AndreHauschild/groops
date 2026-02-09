/***********************************************/
/**
* @file simulateStarCameraMetopSG.cpp
*
* @brief Rotation from satellite (x: along, y: cross, z: nadir) to inertial frame.
*
* @author Torsten Mayer-Guerr
* @date 2005-01-21
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program simulates \file{star camera}{instrument} measurements at each satellite's position.
The satellite's orientation follows a local orbit frame with the x-axis in along track (along velocity),
y-axis is cross track (normal to position and velocity vector) and z-axis pointing nadir (negative position vector).
As for non circular orbit the position and velocity are not exact normal, the default is the x-axis to be exact
along velocity and the z-axis forms a right hand system (not exact nadir) or with \config{nadirPointing} the z-axis
is exact nadir and x-axis approximates along.
The resulting rotation matrices rotate from satellite frame to inertial frame.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"
#include "classes/earthRotation/earthRotation.h"

/***** CLASS ***********************************/

/** @brief Rotation from satellite (x: along, y: cross, z: nadir) to inertial frame.
* @ingroup programsGroup */
class SimulateStarCameraMetopSG
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(SimulateStarCameraMetopSG, PARALLEL, "Rotation from satellite (x: along, y: cross, z: nadir) to inertial frame.", Simulation, Instrument)

/***********************************************/

void SimulateStarCameraMetopSG::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName         orbitName, starCameraName;
    EarthRotationPtr earthRotation;
    Bool     isNadirPointing;

    readConfig(config, "outputfileStarCamera", starCameraName,  Config::MUSTSET, "",  "rotation from satellite to inertial frame (x: along, y: cross, z: nadir)");
    readConfig(config, "inputfileOrbit",       orbitName,       Config::MUSTSET, "",  "position and velocity defines the orientation of the satellite at each epoch");
    readConfig(config, "earthRotation",        earthRotation,   Config::MUSTSET,  "", "transformation from CRF to TRF");
    readConfig(config, "nadirPointing",        isNadirPointing, Config::DEFAULT, "1", "false: exact along and nearly nadir, true: nearly along and exact nadir");

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

        const Vector3d position = rot.rotate(orbit.at(i).position);
        velocity = rot.rotate(orbit.at(i).velocity) - crossProduct(omega, position);

        Vector3d y = normalize(crossProduct(velocity, position)); // cross
        Vector3d x = normalize((isNadirPointing) ? crossProduct(position, y) : velocity);  // along

        StarCameraEpoch epoch;
        epoch.time   = orbit.at(i).time;
        epoch.rotary = inverse(rot)*Rotary3d(x, y); // rotation crf <- trf + trf <- sat
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
