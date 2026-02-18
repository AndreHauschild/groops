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
This program simulates inter-satellite link (ISL) observations between terminals on GNSS satellites.
These simulated observations can then be used in \program{GnssProcessing}, for example to conduct closed-loop simulations.

One or more GNSS constellations must be defined via \configClass{transmitter}{gnssTransmitterGeneratorType}.

The \configClass{parametrization}{gnssParametrizationType} are used to simulate a priori models (e.g. ISL biases).
Parameter settings and output files are ignored.
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
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
    FileName                    fileNameIsl, fileNameSchedule;
    TimeSeriesPtr               timeSeries;
    Double                      marginSeconds;
    GnssTransmitterGeneratorPtr transmitterGenerator;
    GnssReceiverGeneratorPtr    receiverGenerator;
    GnssParametrizationPtr      gnssParametrization;
    EarthRotationPtr            earthRotation;
    NoiseGeneratorPtr           noiseObs;

    readConfig(config, "outputfileGnssIsl",       fileNameIsl,          Config::MUSTSET,  "gnssIsl_{loopTime:%D}.{prn}.dat.gz", "variable {prn} available, simulated observations");
    readConfig(config, "timeSeries",              timeSeries,           Config::MUSTSET,  "",    "defines observation epochs");
    readConfig(config, "timeMargin",              marginSeconds,        Config::DEFAULT,  "0.1", "[seconds] margin to consider two times identical");
    readConfig(config, "transmitter",             transmitterGenerator, Config::MUSTSET,  "",    "constellation of GNSS satellites");
    readConfig(config, "islSchedule",             fileNameSchedule,     Config::MUSTSET,  "islSchedule_{loopTime:%D}.{prn}.txt", "variable {prn} available, schedule with ISL connections");
    readConfig(config, "earthRotation",           earthRotation,        Config::MUSTSET,  "",    "apriori earth rotation");
    readConfig(config, "parametrization",         gnssParametrization,  Config::DEFAULT,  R"(["signalBiasesIsl"])", "models and parameters");
    readConfig(config, "noiseObservation",        noiseObs,             Config::DEFAULT,  "",    "[m] ISL observation noise");
    if(isCreateSchema(config)) return;

    // ============================

    // init the GNSS system
    // --------------------
    logInfo<<"Init GNSS"<<Log::endl;
    std::vector<Time> times = timeSeries->times();
    Gnss gnss;
    gnss.init({}, times, seconds2time(marginSeconds), transmitterGenerator, receiverGenerator, earthRotation, gnssParametrization, comm);

    // distribute transmitters to nodes
    // --------------------------------
    for(UInt idTrans=0; idTrans<gnss.transmitters.size(); idTrans++)
      if(idTrans%Parallel::size(comm) == Parallel::myRank(comm)) // distribute to nodes
        gnss.transmitters.at(idTrans)->isMyRank_ = TRUE;

    VariableList fileNameVariableList;
    fileNameVariableList.setVariable("prn", "***");
    logInfo<<"Load ISL schedule files <"<<fileNameSchedule(fileNameVariableList)<<">"<<Log::endl;

    // inter satellite links
    // ---------------------
    for(auto recv : gnss.transmitters)
      if (recv->isMyRank())
      {
        // Loading ISL schedule file
        // -------------------------
        fileNameVariableList.setVariable("prn", recv->name());
        GnssReceiverArc scheduleIsl;
        try
        {
          readFile(fileNameSchedule(fileNameVariableList),gnss.times,seconds2time(marginSeconds),scheduleIsl);
        }
        catch(std::exception &/*e*/)
        {
          logWarning<<"Unable to read ISL schedule <"<<fileNameSchedule(fileNameVariableList)<<">, skipping satellite."<<Log::endl;
          continue;
        }
        if(!scheduleIsl.size())
        {
          logWarning<<"no scheduled ISLs found"<<Log::endl;
          continue;
        }

        // Simulate ISL observations
        // -------------------------
        recv->simulateObservationsIsl(noiseObs,gnss.transmitters,times,scheduleIsl,gnss.funcReduceModelsIsl);
      }

    Parallel::barrier(comm);
    gnss.synchronizeTransceiversIsl(comm);
    Parallel::barrier(comm);
    logInfo<<"  transmitter: "<<std::count_if(gnss.transmitters.begin(), gnss.transmitters.end(), [](auto t) {return t->useable();})<<Log::endl;
    if(!std::any_of(gnss.transmitters.begin(), gnss.transmitters.end(), [](auto trans){return trans->useable();}))
    {
      logWarningOnce<<times.front().dateTimeStr()<<" - "<<times.back().dateTimeStr()<<": no useable transmitters"<<Log::endl;
      return;
    }

    // Count ISL observations
    // ----------------------
    logInfo<<"number of ISL observations:"<<Log::endl;
    Vector countTypes(1);
    for(auto recv : gnss.transmitters)
      if(recv->isMyRank())
        for(UInt idEpoch=0; idEpoch<recv->idEpochSize(); idEpoch++)
          for(UInt idTrans=0; idTrans<recv->idTransmitterSize(idEpoch); idTrans++)
          {
            auto obs = recv->observationIsl(idTrans, idEpoch);
            if(obs)
            {
              countTypes(0)++;
            }
          }
    Parallel::reduceSum(countTypes, 0, comm);
    logInfo<<"  total:"<<countTypes(0)%"%11i"s<<Log::endl;

    // Write observations
    // ------------------
    if(!fileNameIsl.empty())
    {
      VariableList fileNameVariableList;
      fileNameVariableList.setVariable("prn", "***");
      logStatus<<"write ISL observations to files <"<<fileNameIsl(fileNameVariableList)<<">"<<Log::endl;

      for(auto recv : gnss.transmitters)
        if(recv->useable() && recv->isMyRank())
        {

          GnssReceiverArc arc;
          for(UInt idEpoch=0; idEpoch<gnss.times.size(); idEpoch++)
            if(recv->useable(idEpoch))
            {
              GnssReceiverEpoch epoch;
              epoch.time       = gnss.times.at(idEpoch);
              epoch.clockError = recv->clockError(idEpoch);

              // get types
              for(UInt idTrans=0; idTrans<recv->idTransmitterSize(idEpoch); idTrans++)
                if(recv->observationIsl(idTrans, idEpoch) && gnss.transmitters.at(idTrans)->useable(idEpoch)) {
                  epoch.obsType.push_back(GnssType::RANGE + GnssType::E1 + GnssType::C);
                  epoch.obsType.push_back(GnssType::CHANNEL + GnssType::E1 + GnssType::A);
                  epoch.obsType.push_back(GnssType::CHANNEL + GnssType::E1 + GnssType::B);
                  epoch.obsType.push_back(GnssType::SNR + GnssType::E1 + GnssType::C);
                }
              std::sort(epoch.obsType.begin(), epoch.obsType.end());
              if(!epoch.obsType.size())
                continue;

              for(UInt idTrans=0; idTrans<recv->idTransmitterSize(idEpoch); idTrans++)
              {

                if(recv->observationIsl(idTrans, idEpoch) && gnss.transmitters.at(idTrans)->useable(idEpoch))
                {
                  const GnssObservationIsl &obs = *recv->observationIsl(idTrans, idEpoch);
                  const GnssType prn = gnss.transmitters.at(idTrans)->PRN();
                  epoch.satellite.push_back(prn);

                  UInt idType = std::distance(epoch.obsType.begin(), std::find(epoch.obsType.begin(), epoch.obsType.end(), prn));
                  for(; (idType<epoch.obsType.size()) && (epoch.obsType.at(idType) == prn); idType++)
                  {
                    epoch.observation.push_back(NAN_EXPR);
                    if(epoch.obsType.at(idType) == GnssType(GnssType::SNR + GnssType::E1 + GnssType::C))
                      epoch.observation.back() = obs.sigma0;
                    else if(epoch.obsType.at(idType) == GnssType(GnssType::CHANNEL + GnssType::E1 + GnssType::A))
                      epoch.observation.back() = obs.terminalSend;
                    else if(epoch.obsType.at(idType) == GnssType(GnssType::CHANNEL + GnssType::E1 + GnssType::B))
                      epoch.observation.back() = obs.terminalRecv;
                    else
                      epoch.observation.back() = obs.observation;
                  }
                } // for(idTrans)
              }

              if(epoch.satellite.size())
                arc.push_back(epoch);
            } // for(idEpoch)

          fileNameVariableList.setVariable("prn", recv->name());
          if(arc.size())
            InstrumentFile::write(fileNameIsl(fileNameVariableList), arc);
        } // for(recv)
    } // if(fileNameIsl)

    // ============================

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
void GnssSimulateIsl::readFile(const FileName &fileName, const std::vector<Time> &times, const Time &timeMargin, GnssReceiverArc &scheduleIsl) const
{
  try
  {
    InFile file(fileName);

    scheduleIsl = GnssReceiverArc();

    // read data
    // ---------
    GnssReceiverArc arc;
    std::string     line, prn_, chaTx_,chaRx_;
    while(std::getline(file, line))
    {
      // skip comment lines
      // ------------------
      if (line[0]=='#')
          continue;

      UInt   mjdInt = String::toInt(line.substr(0, 5));
      Double mjdMod = String::toDouble("0"+line.substr(5, 21));
      Time   time   = Time(mjdInt,mjdMod);

      std::stringstream ss(line.substr(27));
      ss>>prn_>>chaTx_>>chaRx_;

      if(arc.size() && (time<arc.back().time))
        throw(Exception("epochs not in increasing order"));

      if((arc.size()==0) || (time>arc.back().time))
      {
        GnssReceiverEpoch epoch;
        epoch.time = time;
        arc.push_back(epoch);
      }
      arc.back().satellite.push_back(GnssType("***"+prn_));
      arc.back().obsType.push_back(GnssType("X1A"+prn_));
      arc.back().observation.push_back(String::toDouble(chaTx_));
      arc.back().obsType.push_back(GnssType("X1B"+prn_));
      arc.back().observation.push_back(String::toDouble(chaRx_));

    }
    if(!arc.size())
      return;

    // generate lists consistent with the epochs in input time series
    // --------------------------------------------------------------
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
          // fill-in epochs missing in scheduler file
          // ----------------------------------------
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
