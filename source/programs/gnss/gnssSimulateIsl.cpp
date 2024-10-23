/***********************************************/
/**
* @file gnssSimulateIsl.cpp
*
* @brief GNSS inter-satellite link simulation
*
* @author Andre Hauschild
* @date 2024-10-01
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program simulates inter-satellite link (ISL) observations from receivers to GNSS satellites.
These simulated observations can then be used in \program{GnssProcessing}, for example to conduct closed-loop simulations.

One or more GNSS constellations must be defined via \configClass{transmitter}{gnssTransmitterGeneratorType}.

The \configClass{parametrization}{gnssParametrizationType} are used to simulate a priori models (e.g. ISL biases).
Parameter settings and output files are ignored.
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "files/fileInstrument.h"
#include "classes/earthRotation/earthRotation.h"
#include "classes/noiseGenerator/noiseGenerator.h"
#include "classes/timeSeries/timeSeries.h"
#include "gnss/gnss.h"
#include "gnss/gnssObservation.h"
#include "gnss/gnssTransmitterGenerator/gnssTransmitterGenerator.h"
#include "gnss/gnssReceiverGenerator/gnssReceiverGenerator.h"

/***** CLASS ***********************************/

/** @brief GNSS inter-satellite link simulation.
* @ingroup programsGroup */
class GnssSimulateIsl
{
  void readFile(const FileName &fileName, const std::vector<Time> &times, const Time &timeMargin, GnssReceiverArc &scheduleIsl) const;
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GnssSimulateIsl, PARALLEL, "GNSS ISL simulation", Gnss, Simulation)

/***********************************************/

void GnssSimulateIsl::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName                    fileNameIsl, fileNameClock, fileNameSchedule;
    TimeSeriesPtr               timeSeries;
    Double                      marginSeconds;
    GnssTransmitterGeneratorPtr transmitterGenerator;
    GnssReceiverGeneratorPtr    receiverGenerator;
    GnssParametrizationPtr      gnssParametrization;
    EarthRotationPtr            earthRotation;
    NoiseGeneratorPtr           noiseClock, noiseObs;
    std::vector<GnssType>       obsTypes;

    readConfig(config, "outputfileGnssIsl",       fileNameIsl,          Config::MUSTSET,  "gnssIsl_{loopTime:%D}.{station}.dat.gz", "variable {station} available, simulated observations");
    readConfig(config, "outputfileClock",         fileNameClock,        Config::OPTIONAL, "clock_{loopTime:%D}.{station}.dat", "variable {station} available, simulated receiver clock errors");
    readConfig(config, "timeSeries",              timeSeries,           Config::MUSTSET,  "",    "defines observation epochs");
    readConfig(config, "timeMargin",              marginSeconds,        Config::DEFAULT,  "0.1", "[seconds] margin to consider two times identical");
    readConfig(config, "transmitter",             transmitterGenerator, Config::MUSTSET,  "",    "constellation of GNSS satellites");
    readConfig(config, "receiver",                receiverGenerator,    Config::MUSTSET,  "",    "ground station network or LEO satellite");
    readConfig(config, "islSchedule",             fileNameSchedule,     Config::MUSTSET,  "islSchedule_{loopTime:%D}.{station}.txt", "variable {station} available, schedule with ISL connections");
    readConfig(config, "earthRotation",           earthRotation,        Config::MUSTSET,  "",    "apriori earth rotation");
    readConfig(config, "parametrization",         gnssParametrization,  Config::DEFAULT,  R"(["signalBiasesIsl"])", "models and parameters");
    readConfig(config, "noiseObservation",        noiseObs,             Config::DEFAULT,  "",    "[-] noise is multiplied with type accuracy pattern of receiver");
    readConfig(config, "noiseClockReceiver",      noiseClock,           Config::DEFAULT,  "",    "[m] noise added to the simulated receiver clock");
    if(isCreateSchema(config)) return;

    // ============================

    // Set fixed list of observation types for all simulated GNSSs.
    //
    // NOTE: only C1X and X1X will alter be returned as observations, the rest
    //       is required internally for validation checks and is filtered prior
    //       to output.
    //
    // TODO: check if the X-type is the best choice due to the bias handling of
    //       combined observation types!
    // -------------------------------------------------------------------------

    obsTypes.push_back(GnssType::RANGE   + GnssType::E1  + GnssType::X);
    obsTypes.push_back(GnssType::RANGE   + GnssType::E5a + GnssType::X);
    obsTypes.push_back(GnssType::PHASE   + GnssType::E1  + GnssType::X);
    obsTypes.push_back(GnssType::PHASE   + GnssType::E5a + GnssType::X);
    obsTypes.push_back(GnssType::CHANNEL + GnssType::E1  + GnssType::X);

    // init the GNSS system
    // --------------------
    logInfo<<"Init GNSS"<<Log::endl;
    std::sort(obsTypes.begin(), obsTypes.end());
    std::vector<Time> times = timeSeries->times();
    Gnss gnss;
    gnss.init(times, seconds2time(marginSeconds), transmitterGenerator, receiverGenerator, earthRotation, gnssParametrization, comm);
    receiverGenerator->simulation(obsTypes, noiseClock, noiseObs, &gnss, comm);
    gnss.synchronizeTransceivers(comm);
    logInfo<<"  transmitter: "<<std::count_if(gnss.transmitters.begin(), gnss.transmitters.end(), [](auto t) {return t->useable();})<<Log::endl;
    logInfo<<"  receiver:    "<<std::count_if(gnss.receivers.begin(),    gnss.receivers.end(),    [](auto r) {return r->useable();})<<Log::endl;
    if(!std::any_of(gnss.transmitters.begin(), gnss.transmitters.end(), [](auto trans){return trans->useable();}))
    {
      logWarningOnce<<times.front().dateTimeStr()<<" - "<<times.back().dateTimeStr()<<": no useable transmitters"<<Log::endl;
      return;
    }
    if(!std::any_of(gnss.receivers.begin(), gnss.receivers.end(), [](auto recv){return recv->useable();}))
    {
      logWarningOnce<<times.front().dateTimeStr()<<" - "<<times.back().dateTimeStr()<<": no useable receivers"<<Log::endl;
      return;
    }

    // count observation types
    // -----------------------
    logInfo<<"types and number of observations:"<<Log::endl;
    std::vector<GnssType> types = gnss.types(~(GnssType::PRN + GnssType::FREQ_NO));
    Vector countTypes(types.size());
    for(auto recv : gnss.receivers)
      if(recv->isMyRank())
        for(UInt idEpoch=0; idEpoch<recv->idEpochSize(); idEpoch++)
          for(UInt idTrans=0; idTrans<recv->idTransmitterSize(idEpoch); idTrans++)
          {
            auto obs = recv->observation(idTrans, idEpoch);
            if(obs)
              for(UInt idType=0; idType<obs->size(); idType++)
              {
                const UInt idx = GnssType::index(types, obs->at(idType).type);
                if(idx != NULLINDEX)
                  countTypes(idx)++;
              }
          }
    Parallel::reduceSum(countTypes, 0, comm);

    for(UInt idType=0; idType<types.size(); idType++)
      logInfo<<"  "<<types.at(idType).str()<<":"<<countTypes(idType)%"%10i"s<<Log::endl;
    logInfo<<"        + ========="<<Log::endl;
    logInfo<<"  total:"<<sum(countTypes)%"%11i"s<<Log::endl;

    UInt countTracks = 0;
    for(auto recv : gnss.receivers)
      if(recv->isMyRank())
        countTracks += recv->tracks.size();
    Parallel::reduceSum(countTracks, 0, comm);
    logInfo<<"  number of tracks: "<<countTracks<<Log::endl;

    // ============================

    // Write observations
    // ------------------
    if(!fileNameIsl.empty())
    {
      VariableList fileNameVariableList;
      fileNameVariableList.setVariable("station", "****");
      logStatus<<"write receiver observations to files <"<<fileNameIsl(fileNameVariableList)<<">"<<Log::endl;
      for(auto recv : gnss.receivers)
        if(recv->isMyRank())
        {

          // Loading ISL schedule file
          // -------------------------
          fileNameVariableList.setVariable("station", recv->markerNumber());
          logInfo<<"Load ISL schedule "<<fileNameSchedule(fileNameVariableList)<<Log::endl;

          GnssReceiverArc scheduleIsl;
          readFile(fileNameSchedule(fileNameVariableList),gnss.times,seconds2time(marginSeconds),scheduleIsl);

          GnssReceiverArc arc;
          for(UInt idEpoch=0; idEpoch<gnss.times.size(); idEpoch++)
            if(recv->useable(idEpoch))
            {
              GnssReceiverEpoch epoch;
              epoch.time       = gnss.times.at(idEpoch);
              epoch.clockError = recv->clockError(idEpoch);

              // get types
              for(UInt idTrans=0; idTrans<recv->idTransmitterSize(idEpoch); idTrans++)
                if(recv->observation(idTrans, idEpoch) && gnss.transmitters.at(idTrans)->useable(idEpoch))
                  for(UInt idType=0; idType<recv->observation(idTrans, idEpoch)->size(); idType++)
                  {
                    // Filter all observation types except for C1X and X1X (works for all GNSSs!)
                    GnssType obsType = recv->observation(idTrans, idEpoch)->at(idType).type;
                    if(obsType != GnssType(GnssType::RANGE   + GnssType::E1 + GnssType::X) &&
                       obsType != GnssType(GnssType::CHANNEL + GnssType::E1 + GnssType::X))
                      continue;
                    if(!recv->observation(idTrans, idEpoch)->at(idType).type.isInList(epoch.obsType))
                      epoch.obsType.push_back(recv->observation(idTrans, idEpoch)->at(idType).type & ~(GnssType::PRN+GnssType::FREQ_NO));
                  }
              std::sort(epoch.obsType.begin(), epoch.obsType.end());
              if(!epoch.obsType.size())
                continue;

              for(UInt idTrans=0; idTrans<recv->idTransmitterSize(idEpoch); idTrans++)
              {

                // Filter for transmitter PRN from the ISL schedule
                // ------------------------------------------------
                if (!gnss.transmitters.at(idTrans)->PRN().isInList(scheduleIsl.at(idEpoch).satellite)) {
                  continue;
                }

                if(recv->observation(idTrans, idEpoch) && gnss.transmitters.at(idTrans)->useable(idEpoch))
                {
                  const GnssObservation &obs = *recv->observation(idTrans, idEpoch);
                  const GnssType prn = obs.at(0).type & (GnssType::SYSTEM + GnssType::PRN + GnssType::FREQ_NO);
                  UInt idType = std::distance(epoch.obsType.begin(), std::find(epoch.obsType.begin(), epoch.obsType.end(), prn));

                  epoch.satellite.push_back(prn);
                  for(; (idType<epoch.obsType.size()) && (epoch.obsType.at(idType) == prn); idType++)
                  {
                    epoch.observation.push_back(NAN_EXPR);
                    for(UInt i=0; i<obs.size(); i++)
                      if(obs.at(i).type == epoch.obsType.at(idType))
                      {
                        epoch.observation.back() = obs.at(i).observation;
                        break;
                      }
                  }
                } // for(idTrans)
              }

              if(epoch.satellite.size())
                arc.push_back(epoch);
            } // for(idEpoch)

          fileNameVariableList.setVariable("station", recv->markerNumber());
          InstrumentFile::write(fileNameIsl(fileNameVariableList), arc);
        } // for(recv)
    } // if(fileNameReceiver)

    // ============================

    // write clock errors (TODO: remove this?)
    // ------------------
    if(!fileNameClock.empty())
    {
      VariableList fileNameVariableList;
      fileNameVariableList.setVariable("station", "****");
      logStatus<<"write receiver clocks to files <"<<fileNameClock(fileNameVariableList)<<">"<<Log::endl;
      for(auto recv : gnss.receivers)
        if(recv->isMyRank())
        {
          MiscValueArc arc;
          for(UInt idEpoch=0; idEpoch<gnss.times.size(); idEpoch++)
            if(recv->useable(idEpoch))
            {
              MiscValueEpoch epoch;
              epoch.time  = gnss.times.at(idEpoch);
              epoch.value = recv->clockError(idEpoch);
              arc.push_back(epoch);
            }
          fileNameVariableList.setVariable("station", recv->markerNumber());
          InstrumentFile::write(fileNameClock(fileNameVariableList), arc);
        } // for(recv)
    } // if(fileNameClock)
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssSimulateIsl::readFile(const FileName &fileName, const std::vector<Time> &times, const Time &timeMargin, GnssReceiverArc &scheduleIsl) const
{
  try
  {
    InFile file(fileName);

    GnssReceiverArc arc;

    // read data
    // ---------
    std::string line;
    while(std::getline(file, line))
    {
      UInt   mjdInt = String::toInt(line.substr(0, 5));
      Double mjdMod = String::toDouble("0"+line.substr(5, 21));
      Time   time   = Time(mjdInt,mjdMod);

      GnssType prn = GnssType("***"+line.substr(44, 3));
      Double   cha = String::toInt(line.substr(48, 1));

      if(arc.size() && (time<arc.back().time))
        throw(Exception("epochs not in increasing order"));

      if((arc.size()==0) || (time>arc.back().time))
      {
        GnssReceiverEpoch epoch;
        epoch.time = time;
        arc.push_back(epoch);
      }
      arc.back().satellite.push_back(prn);
      arc.back().obsType.push_back(GnssType::CHANNEL + GnssType::E1 + GnssType::X);
      arc.back().observation.push_back(cha);
    }

    // Generate lists consistent with the epochs in input time series

    UInt idEpoch = 0;
    for(UInt arcEpoch=0; arcEpoch<arc.size(); arcEpoch++)
    {
      while((idEpoch < times.size()) && (times.at(idEpoch)+timeMargin < arc.at(arcEpoch).time))
      {
        if(times.at(idEpoch)+timeMargin < arc.at(0).time)
        {
          GnssReceiverEpoch epoch;
          epoch.time = times.at(idEpoch++);
          scheduleIsl.push_back(epoch);
          logWarning<<"missing ISL scheduler data at "<<epoch.time.dateTimeStr()<<" (front)"<<Log::endl;
        }
        else
        {
          GnssReceiverEpoch epoch;
          epoch.time        = times.at(idEpoch++);
          epoch.satellite   = arc.at(arcEpoch).satellite;
          epoch.observation = arc.at(arcEpoch).observation;
          epoch.obsType     = arc.at(arcEpoch).obsType;
          scheduleIsl.push_back(epoch);
        }
      }
      if(idEpoch >= times.size())
        break;
      if((arc.at(arcEpoch).time+timeMargin < times.at(idEpoch)))
        continue;
      GnssReceiverEpoch epoch;
      epoch.time        = times.at(idEpoch++);
      epoch.satellite   = arc.at(arcEpoch).satellite;
      epoch.observation = arc.at(arcEpoch).observation;
      epoch.obsType     = arc.at(arcEpoch).obsType;
      scheduleIsl.push_back(epoch);
    }
    for(; idEpoch<times.size(); idEpoch++)
    {
      GnssReceiverEpoch epoch;
      epoch.time = times.at(idEpoch);
      scheduleIsl.push_back(epoch);
      logWarning<<"missing ISL scheduler data at "<<epoch.time.dateTimeStr()<<" (end)"<<Log::endl;
    }

    GnssReceiverArc::printStatistics(scheduleIsl);

  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
