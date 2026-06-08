/***********************************************/
/**
* @file gnssProcessingStepWriteResidualsIsl.h
*
* @brief GNSS processing step: WriteResidualsIsl.
*
* @author Andre Hauschild
* @date 2026-03-26
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPROCESSINGSTEPWRITERESIDUALSISL__
#define __GROOPS_GNSSPROCESSINGSTEPWRITERESIDUALSISL__

// Latex documentation
#ifdef DOCSTRING_GnssProcessingStep
static const char *docstringGnssProcessingStepWriteResidualsIsl = R"(
\subsection{WriteResidualsIsl}\label{gnssProcessingStepType:writeResidualsIsl}
Writes the ISL \file{observation residuals}{instrument} for all
\configClass{selectTransmitters}{platformSelectorType}.
One file is written for each transmitter. The file name is interpreted as
a template with the variable \verb|{prn}| being replaced by the satellite PRN.
)";
#endif

/***********************************************/

#include "config/config.h"
#include "classes/platformSelector/platformSelector.h"
#include "gnss/gnssProcessingStep/gnssProcessingStep.h"

/***** CLASS ***********************************/

/** @brief GNSS processing step: WriteResidualsIsl.
* @ingroup gnssProcessingStepGroup
* @see GnssProcessingStep */
class GnssProcessingStepWriteResidualsIsl : public GnssProcessingStepBase
{
  PlatformSelectorPtr selectTransmitters;
  FileName            fileNameResiduals;

public:
  GnssProcessingStepWriteResidualsIsl(Config &config);
  void process(GnssProcessingStep::State &state) override;
};

/***********************************************/

inline GnssProcessingStepWriteResidualsIsl::GnssProcessingStepWriteResidualsIsl(Config &config)
{
  try
  {
    readConfig(config, "selectTransmitter",   selectTransmitters, Config::OPTIONAL, "", "subset of used transmitters");
    readConfig(config, "outputfileResiduals", fileNameResiduals,  Config::OPTIONAL, "output/residuals_{loopTime:%D}.{prn}.dat", "variable {prn} available");
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void GnssProcessingStepWriteResidualsIsl::process(GnssProcessingStep::State &state)
{
  try
  {
    auto selectedTransmitters = state.gnss->selectTransmitters(selectTransmitters);
    VariableList fileNameVariableList;
    fileNameVariableList.setVariable("prn", "***");
    logStatus<<"write ISL residuals to file <"<<fileNameResiduals(fileNameVariableList)<<">"<<Log::endl;
    for(auto recv : state.gnss->transmitters)
      if(selectedTransmitters.at(recv->idTrans()))
      {
        GnssReceiverArc arc;
        for(UInt idEpoch : state.normalEquationInfo.idEpochs)
          if(recv->useable(idEpoch))
          {
            // get types
            std::vector<GnssType> types;
            for(UInt idTrans=0; idTrans<recv->idTransmitterSize(idEpoch); idTrans++)
              if(recv->observationIsl(idTrans, idEpoch) && state.gnss->transmitters.at(idTrans)->useable(idEpoch))
              {
                const GnssType typeIsl = GnssType("C1C") + state.gnss->transmitters.at(idTrans)->PRN();
                if(!typeIsl.isInList(types))
                  types.push_back(typeIsl & ~(GnssType::PRN+GnssType::FREQ_NO));
              }
            std::sort(types.begin(), types.end());
            if(!types.size())
              continue;

            GnssReceiverEpoch epoch;
            GnssType system = GnssType::SYSTEM;
            for(UInt idType=0; idType<types.size(); idType++)
            {
              if(types.at(idType) != system)
              {
                system = types.at(idType) & GnssType::SYSTEM;
                epoch.obsType.push_back( GnssType::AZIMUT    + GnssType::L1 + system );
                epoch.obsType.push_back( GnssType::ELEVATION + GnssType::L1 + system );
                epoch.obsType.push_back( GnssType::AZIMUT    + GnssType::L2 + system );
                epoch.obsType.push_back( GnssType::ELEVATION + GnssType::L2 + system );
              }
              // residuals, redundancy, sigma
              epoch.obsType.insert(epoch.obsType.end(), {types.at(idType), types.at(idType), types.at(idType)});
            }

            for(UInt idTrans=0; idTrans<recv->idTransmitterSize(idEpoch); idTrans++)
              if(recv->observationIsl(idTrans, idEpoch) && state.gnss->transmitters.at(idTrans)->useable(idEpoch))
              {
                const GnssObservationIsl &obs = *recv->observationIsl(idTrans, idEpoch);
                const GnssObservationEquationIsl eqn(obs, *recv, *state.gnss->transmitters.at(idTrans),
                                                     state.gnss->funcReduceModelsIsl, idEpoch, FALSE);
                const GnssType prn = state.gnss->transmitters.at(idTrans)->PRN();
                UInt idType = std::distance(epoch.obsType.begin(), std::find(epoch.obsType.begin(), epoch.obsType.end(), prn));
                if(idType >= epoch.obsType.size())
                  continue;

                epoch.time = eqn.timeRecv;
                epoch.satellite.push_back(prn);
                epoch.observation.insert(epoch.observation.end(), {eqn.azimutRecvAnt, eqn.elevationRecvAnt, eqn.azimutTrans, eqn.elevationTrans});
                idType += 4;

                for(; (idType<epoch.obsType.size()) && (epoch.obsType.at(idType) == prn); idType+=3)
                {
                  epoch.observation.insert(epoch.observation.end(), {0., 0., 1.});
                  const GnssType typeIsl = GnssType("C1C") + state.gnss->transmitters.at(idTrans)->PRN();
                  if(typeIsl == epoch.obsType.at(idType))
                  {
                    epoch.observation.at(epoch.observation.size()-3) = obs.residual;
                    epoch.observation.at(epoch.observation.size()-2) = obs.redundancy;
                    epoch.observation.at(epoch.observation.size()-1) = obs.sigma/obs.sigma0;
                    break;
                  }
                }
              } // for(idTrans)

            if(epoch.satellite.size())
              arc.push_back(epoch);
          } // for(idEpoch)

        if(!arc.size())
          continue;

        VariableList fileNameVariableList;
        fileNameVariableList.setVariable("prn", recv->name());
        InstrumentFile::write(fileNameResiduals(fileNameVariableList), arc);
      } // for(trans)
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
