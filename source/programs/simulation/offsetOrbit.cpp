/***********************************************/
/**
* @file offsetOrbit.cpp
*
* @brief Add along-track, cross-track and radial offset to orbit positions.
*
* @author Andre Hauschild
* @date 2024-111-06
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program adds along-track, cross-track and radial offsets to \file{satellite}{instrument}'s positions.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief Add along-track, cross-track and radial offset to orbit positions.
* @ingroup programsGroup */
class OffsetOrbit
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(OffsetOrbit, PARALLEL, "add offsets to orbit positions", Simulation, Orbit, Instrument)

/***********************************************/

void OffsetOrbit::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName          inName, outName;
    Double            offsetAlong, offsetCross, offsetRadial;
    Vector3d          offsetPosition = Vector3d();

    readConfig(config, "outputfileOrbit",       outName,                Config::MUSTSET,  "", "");
    readConfig(config, "inputfileOrbit",        inName,                 Config::MUSTSET,  "", "");
    readConfig(config, "offsetPositionAlong",   offsetAlong,            Config::DEFAULT,  "0.0", "along [m]");
    readConfig(config, "offsetPositionCross",   offsetCross,            Config::DEFAULT,  "0.0", "cross [m]");
    readConfig(config, "offsetPositionRadial",  offsetRadial,           Config::DEFAULT,  "0.0", "radial [m]");
    if(isCreateSchema(config)) return;

    offsetPosition = Vector3d(offsetAlong,offsetCross,offsetRadial);

    // Read satellite
    // -----------------
    logStatus<<"add offsets to orbit data <"<<inName<<">"<<Log::endl;
    InstrumentFile  orbitFile(inName);
    std::vector<Arc> arcList(orbitFile.arcCount());

    Parallel::forEach(arcList, [&](UInt arcNo)
    {
      OrbitArc orbit  = orbitFile.readArc(arcNo);

      for(UInt i=0; i<orbit.size(); i++)
      {
        // rotate into satellite system
        // ----------------------------
        Rotary3d rot;
        if(orbit.size()>1)
        {
          Vector3d x;
          if(i==0)
            x = orbit.at(i+1).position - orbit.at(i).position;
          else
            x = orbit.at(i).position - orbit.at(i-1).position;
          Vector3d z = normalize(orbit.at(i).position);
          Vector3d y = normalize(crossProduct(z, x));
          x = crossProduct(y, z);
          rot = Rotary3d(x,y);
        }

        orbit.at(i).position += rot.rotate(offsetPosition);

      }
      return orbit;
    }, comm);

    // Save
    // ----
    if(Parallel::isMaster(comm))
    {
      logStatus<<"write orbit data to file <"<<outName<<">"<<Log::endl;
      InstrumentFile::write(outName, arcList);
      Arc::printStatistics(arcList);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
