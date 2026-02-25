/***********************************************/
/**
* @file gnssTransmitter.cpp
*
* @brief GNSS transmitter.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2013-06-28
*
*/
/***********************************************/

#include "base/import.h"
#include "base/string.h"
#include "files/fileInstrument.h"
#include "inputOutput/logging.h"
#include "gnss/gnssObservationIsl.h"
#include "gnss/gnssTransmitter.h"
#include "classes/noiseGenerator/noiseGenerator.h"

/***********************************************/

GnssTransmitter::~GnssTransmitter()
{
  for(auto &obsEpoch : observations_)
    for(auto &obs : obsEpoch)
      delete obs;
}

/***********************************************/

void GnssTransmitter::disable(UInt idEpoch, const std::string &reason)
{
  try
  {
    GnssTransceiver::disable(idEpoch, reason);
    if(idEpoch < observations_.size())
      observations_.at(idEpoch).clear();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssTransmitter::disable(const std::string &reason)
{
  try
  {
    GnssTransceiver::disable(reason);
    if(!reason.empty() && isMyRank())
      disableReason = reason;
    isMyRank_ = FALSE;
    observations_.clear();
    observations_.shrink_to_fit();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

GnssObservationIsl *GnssTransmitter::observationIsl(UInt idTrans, UInt idEpoch) const
{
  if((idEpoch < observations_.size()) && (idTrans < observations_.at(idEpoch).size()))
    return observations_[idEpoch][idTrans];
  return nullptr;
}

/***********************************************/

void GnssTransmitter::deleteObservationIsl(UInt idTrans, UInt idEpoch)
{
  try
  {
    if(!observationIsl(idTrans, idEpoch))
      return;
    observations_[idEpoch][idTrans] = nullptr;
    if(std::all_of(observations_[idEpoch].begin(), observations_[idEpoch].end(), [](auto obs) {return obs == nullptr;}))
      disable(idEpoch, "no valid epochs left");
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

// C1C  range
// X1A  terminalSend
// X1B  terminalRecv
// S1C  sigma, accuracy
void GnssTransmitter::readObservationsIsl(const FileName &fileName, const std::vector<GnssTransmitterPtr> &transmitters, const std::vector<Time> &times, const Time &timeMargin)
{
  try
  {
    GnssReceiverArc arc = InstrumentFile::read(fileName);

    UInt idEpoch = 0;
    for(UInt arcEpoch=0; arcEpoch<arc.size(); arcEpoch++)
    {
      // search time slot
      while((idEpoch < times.size()) && (times.at(idEpoch)+timeMargin < arc.at(arcEpoch).time))
        idEpoch++;
      if(idEpoch >= times.size())
        break;
      if((arc.at(arcEpoch).time+timeMargin < times.at(idEpoch)) || !useable(idEpoch))
        continue;

      // create observation class for each satellite
      UInt idObs  = 0;
      for(UInt k=0; k<arc.at(arcEpoch).satellite.size(); k++)
      {
        // find list of observation types for this satellite
        GnssType satType = arc.at(arcEpoch).satellite.at(k);
        UInt idType = 0;
        while(arc.at(arcEpoch).obsType.at(idType) != satType)
          idType++;

        // search transmitter index for satellite number (PRN)
        const UInt idTrans = std::distance(transmitters.begin(), std::find_if(transmitters.begin(), transmitters.end(),
                                                                              [&](auto t) {return t->PRN() == satType;}));

        GnssObservationIsl *obs = new GnssObservationIsl();
        obs->time = arc.at(arcEpoch).time;
        for(; (idType<arc.at(arcEpoch).obsType.size()) && (arc.at(arcEpoch).obsType.at(idType)==satType); idType++, idObs++)
          if((idTrans < transmitters.size()) && arc.at(arcEpoch).observation.at(idObs)  && !std::isnan(arc.at(arcEpoch).observation.at(idObs)))
          {
            GnssType type = arc.at(arcEpoch).obsType.at(idType);
            if(type == GnssType::RANGE)
              obs->observation = arc.at(arcEpoch).observation.at(idObs);
            else if(type == GnssType::SNR)
              obs->sigma = obs->sigma0 = arc.at(arcEpoch).observation.at(idObs);
            else if(type == GnssType::CHANNEL+GnssType::B)
              obs->terminalRecv = static_cast<UInt>(arc.at(arcEpoch).observation.at(idObs));
            else if(type == GnssType::CHANNEL+GnssType::A)
              obs->terminalSend = static_cast<UInt>(arc.at(arcEpoch).observation.at(idObs));
          }

        if(observations_.size() <= idEpoch)
          observations_.resize(idEpoch+1);
        if(observations_.at(idEpoch).size() <= idTrans)
          observations_.at(idEpoch).resize(idTrans+1, nullptr);
        if(observations_[idEpoch][idTrans])
          logWarning<<name()<<" -> "<<transmitters.at(idTrans)->name()<<" at "<<times.at(idEpoch).dateTimeStr()<<": observation already exists"<<Log::endl;
        std::swap(observations_[idEpoch][idTrans], obs);
        delete obs;
      } // for(satellite)

      idEpoch++;
    } // for(arcEpoch)
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssTransmitter::simulateObservationsIsl(NoiseGeneratorPtr noiseObs,
                                              const std::vector<GnssTransmitterPtr> &transmitters,
                                              const std::vector<Time> &times,
                                              const GnssReceiverArc   &scheduleIsl,
                                              const std::function<void(GnssObservationEquationIsl &eqn)> &reduceModels)
{
  try
  {

    // GNSS signal types for transmit and receive terminal IDs
    // -------------------------------------------------------
    const GnssType typeTx = GnssType::CHANNEL + GnssType::E1  + GnssType::A;
    const GnssType typeRx = GnssType::CHANNEL + GnssType::E1  + GnssType::B;

    // Simulate zero observations
    // --------------------------
    const Vector eps = noiseObs->noise(times.size()); // obs noise
    for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
    {

      GnssReceiverEpoch epochSchedule = scheduleIsl.at(idEpoch);

      // create observation class for each satellite
      for(UInt idTrans=0; idTrans<transmitters.size(); idTrans++)
      {

        // Filter for transmitter PRN from the ISL schedule
        // ------------------------------------------------
        if(!transmitters.at(idTrans)->PRN().isInList(scheduleIsl.at(idEpoch).satellite))
          continue;

        // Skip receiver PRN
        // -----------------
        if (this->name()==transmitters.at(idTrans)->name())
        {
          logError<<"Identical transmitter and receiver PRN in ISL schedule for "
                  <<transmitters.at(idTrans)->PRN().str()<<" at "
                  <<times.at(idEpoch).dateTimeStr()<<Log::endl;
          continue;
        }

        UInt idChanTx = std::distance(epochSchedule.obsType.begin(), std::find(epochSchedule.obsType.begin(), epochSchedule.obsType.end(),
                                      typeTx + transmitters.at(idTrans)->PRN()));
        UInt idChanRx = std::distance(epochSchedule.obsType.begin(), std::find(epochSchedule.obsType.begin(), epochSchedule.obsType.end(),
                                      typeRx + transmitters.at(idTrans)->PRN()));

        GnssObservationIsl *obs = new GnssObservationIsl();
        obs->time = times.at(idEpoch);
        obs->terminalRecv = scheduleIsl.at(idEpoch).observation.at(idChanRx);
        obs->terminalSend = scheduleIsl.at(idEpoch).observation.at(idChanTx);

        GnssObservationEquationIsl *eqn = new GnssObservationEquationIsl(*obs, *this, *transmitters.at(idTrans), reduceModels, idEpoch, FALSE);

        obs->observation = -eqn->l(0) + eps.at(idEpoch);
        obs->sigma0 = sqrt(noiseObs->covarianceFunction(1)(0,1));

        if(observations_.size() <= idEpoch)
          observations_.resize(idEpoch+1);
        if(observations_.at(idEpoch).size() <= idTrans)
          observations_.at(idEpoch).resize(idTrans+1, nullptr);
        if(observations_[idEpoch][idTrans])
          logWarning<<name()<<" -> "<<transmitters.at(idTrans)->name()<<" at "<<times.at(idEpoch).dateTimeStr()<<": observation already exists"<<Log::endl;
        std::swap(observations_[idEpoch][idTrans], obs);

        delete obs;
        delete eqn;
      } // for(satellite)

      if((observations_.size() <= idEpoch) || (observations_[idEpoch].size() == 0))
        disable(idEpoch, "no observations simulated (elevationCutOff, use/ignoreTypes, defined receiver/transmitter types, missing antenna patterns)");
    } // for(arcEpoch)

  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}
