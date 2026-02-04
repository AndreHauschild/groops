/***********************************************/
/**
* @file gnss.cpp
*
* @brief global navigation satellite system.
*
* @author Torsten Mayer-Guerr
* @date 2010-08-03
*
*/
/***********************************************/

#define DEBUG_SYNC_ISL 1
#define DEBUG          0

#include <vector>
#include <queue>
#include <unordered_set>

#include "base/import.h"
#include "base/planets.h"
#include "config/config.h"
#include "inputOutput/logging.h"
#include "classes/earthRotation/earthRotation.h"
#include "classes/platformSelector/platformSelector.h"
#include "gnss.h"
#include "gnss/gnssObservation.h"
#include "gnss/gnssDesignMatrix.h"
#include "gnss/gnssTransmitter.h"
#include "gnss/gnssReceiver.h"
#include "gnss/gnssParametrization/gnssParametrization.h"
#include "gnss/gnssTransmitterGenerator/gnssTransmitterGenerator.h"
#include "gnss/gnssReceiverGenerator/gnssReceiverGenerator.h"

/***********************************************/

// Breadth-First Search (BFS)

std::vector<UInt> bfs(UInt start, const std::vector<std::vector<UInt>>& graph) {

  std::unordered_set<UInt> visited;
  std::queue<UInt> q;
  std::vector<UInt> result;

  q.push(start);
  visited.insert(start);

  while (!q.empty())
  {
    int node = q.front();
    q.pop();
    result.push_back(node);

    for (int neighbor : graph[node])
      if (!visited.count(neighbor))
      {
        visited.insert(neighbor);
        q.push(neighbor);
      }
  }
  return result;
}

/***********************************************/

bool isInList(std::vector<UInt> v, UInt i) {
  return std::find(v.begin(), v.end(), i)!=v.end();
};

/***********************************************/

void Gnss::init(std::vector<GnssType> simulationTypes, const std::vector<Time> &times, const Time &timeMargin,
                GnssTransmitterGeneratorPtr transmitterGenerator, GnssReceiverGeneratorPtr receiverGenerator,
                EarthRotationPtr earthRotation, GnssParametrizationPtr parametrization, Parallel::CommunicatorPtr comm)
{
  try
  {
    this->times = times;

    // init earth rotation
    // -------------------
    eop = Matrix(times.size(), 8); // Matrix eop columns: xp, yp, sp, deltaUT, LOD, X, Y, S
    for(UInt i=0; i<times.size(); i++)
      earthRotation->earthOrientationParameter(times.at(i), eop(i,0), eop(i,1), eop(i,2), eop(i,3), eop(i,4), eop(i,5), eop(i,6), eop(i,7));
    // UT1-UTC => UT1-GPS (avoid leap seconds jumps for interpolation)
    for(UInt i=0; i<times.size(); i++)
      eop(i,3) -= (times.at(i)-timeGPS2UTC(times.at(i))).seconds();

    funcRotationCrf2Trf = std::bind(&Gnss::rotationCrf2Trf, this, std::placeholders::_1);

    // init transmitters
    // -----------------
    transmitters = transmitterGenerator->transmitters(times, timeMargin, comm);
    for(UInt idTrans=0; idTrans<transmitters.size(); idTrans++)
      transmitters.at(idTrans)->id_ = idTrans;

    // init receivers
    // --------------
    if(receiverGenerator)
    {
      receivers = receiverGenerator->receivers(simulationTypes, times, timeMargin, transmitters, earthRotation, comm);
      for(UInt idRecv=0; idRecv<receivers.size(); idRecv++)
        receivers.at(idRecv)->id_ = idRecv;
    }
    synchronizeTransceivers(comm);

    // init parametrization
    // --------------------
    this->parametrization = parametrization;
    if(parametrization)
    {
      parametrization->init(this, comm);
      funcReduceModels    = std::bind(&GnssParametrization::observationCorrections,    parametrization, std::placeholders::_1);
      funcReduceModelsIsl = std::bind(&GnssParametrization::observationCorrectionsIsl, parametrization, std::placeholders::_1);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Rotary3d Gnss::rotationCrf2Trf(const Time &time) const
{
  try
  {
    const UInt idx = std::min((times.size()-1), static_cast<UInt>(std::distance(times.begin(),
                     std::upper_bound(times.begin(), times.end(), time, [](const Time &t, const Time &s) {return (t-s).seconds() < 0.5;}))));
    const Double xp      = eop(idx, 0);
    const Double yp      = eop(idx, 1);
    const Double sp      = eop(idx, 2);
    const Double deltaUT = eop(idx, 3) + (time-timeGPS2UTC(time)).seconds();
    const Double X       = eop(idx, 5);
    const Double Y       = eop(idx, 6);
    const Double S       = eop(idx, 7);

    const Double ERA = Planets::ERA(timeGPS2UTC(time) + seconds2time(deltaUT));
    const Double r2  = X*X + Y*Y;
    const Double E   = (r2!=0.) ? std::atan2(Y, X) : 0.;
    const Double D   = std::atan(std::sqrt(r2/(1-r2)));

    return  rotaryX(Angle(-yp)) * rotaryY(Angle(-xp)) *
            rotaryZ(Angle(sp+ERA-S-E)) *
            rotaryY(Angle(D)) * rotaryZ(Angle(E));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Gnss::synchronizeTransceivers(Parallel::CommunicatorPtr comm)
{
  try
  {
    // distribute process id of receivers
    // ----------------------------------
    Vector recvProcess(receivers.size());
    for(UInt idRecv=0; idRecv<receivers.size(); idRecv++)
      if(receivers.at(idRecv)->isMyRank())
        recvProcess(idRecv) = Parallel::myRank(comm)+1; // process number for each receiver
    Parallel::reduceSum(recvProcess, 0, comm); // get process number ( all others have zero )
    Parallel::broadCast(recvProcess, 0, comm); // distribute

    // synchronize transceivers
    // ------------------------
    for(UInt idRecv=0; idRecv<receivers.size(); idRecv++)
      if(recvProcess(idRecv))
        Parallel::broadCast(static_cast<GnssTransceiver&>(*receivers.at(idRecv)), static_cast<UInt>(recvProcess(idRecv)-1), comm);
      else if(receivers.at(idRecv)->useable())
        receivers.at(idRecv)->disable("");

    // collect observation types
    // -------------------------
    typesRecvTrans.clear();
    typesRecvTrans.resize(receivers.size(), std::vector<std::vector<GnssType>>(transmitters.size()));
    for(auto recv : receivers)
    {
      if(recv->isMyRank())
        for(UInt idTrans=0; idTrans<transmitters.size(); idTrans++)
        {
          for(UInt idEpoch=0; idEpoch<recv->idEpochSize(); idEpoch++)
          {
            auto obs = recv->observation(idTrans, idEpoch);
            if(obs)
              for(UInt idType=0; idType<obs->size(); idType++)
                if(!obs->at(idType).type.isInList(typesRecvTrans.at(recv->idRecv()).at(idTrans)))
                  typesRecvTrans.at(recv->idRecv()).at(idTrans).push_back(obs->at(idType).type);
          }
          std::sort(typesRecvTrans.at(recv->idRecv()).at(idTrans).begin(), typesRecvTrans.at(recv->idRecv()).at(idTrans).end());
        }
      if(recv->useable())
        Parallel::broadCast(typesRecvTrans.at(recv->idRecv()), static_cast<UInt>(recvProcess(recv->idRecv())-1), comm); // synchronize types of this process to all others
    }

    // adjust signal biases to available observation types
    // NOTE: if no observations are found, the biases are not adjusted and
    //       a-priori values are retained!!
    // ---------------------------------------------------
    for(auto trans : transmitters)
    {
      std::vector<GnssType> types;
      for(auto &typesTrans : typesRecvTrans)
        for(GnssType type : typesTrans.at(trans->idTrans()))
          if((type == GnssType::PHASE) || (type == GnssType::RANGE))
            if(!type.isInList(types))
              types.push_back(type + trans->PRN());
      types = GnssType::replaceCompositeSignals(types);

      if(types.size())
      {
        trans->signalBias.biases = trans->signalBias.compute(types); // apriori signal bias
        trans->signalBias.types  = types;
      }
    }

    for(auto recv : receivers)
    {
      std::vector<GnssType> types;
      for(auto &typesTrans : typesRecvTrans.at(recv->idRecv()))
        for(GnssType type : typesTrans)
          if((type == GnssType::PHASE) || (type == GnssType::RANGE))
            if(!type.isInList(types))
              types.push_back(type & ~GnssType::PRN);
      std::sort(types.begin(), types.end());

      if(types.size())
      {
        recv->signalBias.biases = recv->signalBias.compute(types); // apriori signal bias
        recv->signalBias.types  = types;
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Gnss::synchronizeTransceiversIsl(Parallel::CommunicatorPtr comm)
{
  try
  {

#if DEBUG_SYNC_ISL > 0
    logWarning<<"synchronizeTransceiversIsl() start"
              <<Log::endl;
#endif

    // distribute process id of transmitters
    // -------------------------------------
    Vector recvProcess(transmitters.size());
    for(UInt idTrans=0; idTrans<transmitters.size(); idTrans++)
      if(transmitters.at(idTrans)->isMyRank())
        recvProcess(idTrans) = Parallel::myRank(comm)+1; // process number for each transmitter
    Parallel::reduceSum(recvProcess, 0, comm); // get process number ( all others have zero )
    Parallel::broadCast(recvProcess, 0, comm); // distribute

    // collect ISL observations
    // ------------------------
    /* start obsolete */
    typesRecvTransIsl.clear();
    typesRecvTransIsl.resize(transmitters.size(), std::vector<std::vector<GnssType>>(transmitters.size()));
    /* end obsolete */
    islTerminalRecv.clear();
    islTerminalRecv.resize(transmitters.size(), std::vector<std::vector<UInt>>(transmitters.size()));
    islTerminalTrans.clear();
    islTerminalTrans.resize(transmitters.size(), std::vector<std::vector<UInt>>(transmitters.size()));
    for(auto recvSatellite : transmitters)
    {
      if(recvSatellite->isMyRank())
      {
        for(UInt idTrans=0; idTrans<transmitters.size(); idTrans++)
        {
          for(UInt idEpoch=0; idEpoch<recvSatellite->idEpochSize(); idEpoch++)
          {
            auto obs = recvSatellite->observationIsl(idTrans, idEpoch);
            if(obs)
            {
              if(!isInList(islTerminalRecv.at(recvSatellite->idTrans()).at(idTrans),obs->terminalRecv))
                islTerminalRecv.at(recvSatellite->idTrans()).at(idTrans).push_back(obs->terminalRecv);
              if(!isInList(islTerminalTrans.at(recvSatellite->idTrans()).at(idTrans),obs->terminalSend))
                islTerminalTrans.at(recvSatellite->idTrans()).at(idTrans).push_back(obs->terminalSend);
              /* start obsolete */
              const GnssType typeIsl = GnssType("C1C") + transmitters.at(idTrans)->PRN();
              if(!typeIsl.isInList(typesRecvTransIsl.at(recvSatellite->idTrans()).at(idTrans)))
                typesRecvTransIsl.at(recvSatellite->idTrans()).at(idTrans).push_back(typeIsl);
              /* end obsolete */
            }
          }
        }
      }
      if(recvSatellite->useable())
      {
        Parallel::broadCast(islTerminalRecv.at(recvSatellite->idTrans()), static_cast<UInt>(recvProcess(recvSatellite->idTrans())-1), comm); // synchronize data of this process to all others
        Parallel::broadCast(islTerminalTrans.at(recvSatellite->idTrans()), static_cast<UInt>(recvProcess(recvSatellite->idTrans())-1), comm); // synchronize data of this process to all others        
        Parallel::broadCast(typesRecvTransIsl.at(recvSatellite->idTrans()), static_cast<UInt>(recvProcess(recvSatellite->idTrans())-1), comm); // synchronize types of this process to all others
      }
    }

#if DEBUG_SYNC_ISL > 0
    logWarning<<"synchronizeTransceiversIsl() mid"
              <<Log::endl;
#endif

    // adjust ISL biases to available terminals
    // NOTE: if no observations are found, a-priori biases are removed!
    // ----------------------------------------------------------------
    for(auto trans : transmitters)
    {
      std::vector<UInt> terminals;
      for(auto &termTrans : islTerminalTrans)
        for(UInt terminal : termTrans.at(trans->idTrans()))
          if(!isInList(terminals,terminal))
            terminals.push_back(terminal);
      std::sort(terminals.begin(), terminals.end());

#if DEBUG_SYNC_ISL > 0
      if(Parallel::isMaster(comm))
        logWarning<<"synchronizeTransceiversIsl() send ISL terminal "<<trans->name()<<" "
                  <<terminals.size()%"# terminals %i"s
                  <<Log::endl;
#endif
      // NOTE: a-priori ISL biases are NOT retained!
      trans->signalBiasIslTx.biases    = trans->signalBiasIslTx.compute(terminals); // apriori ISL bias
      trans->signalBiasIslTx.terminals = terminals;
      /*
      if(terminals.size())
      {
        trans->signalBiasIslTx.biases    = trans->signalBiasIslTx.compute(terminals); // apriori ISL bias
        trans->signalBiasIslTx.terminals = terminals;
      }
      */
    }

    for(auto recv : transmitters)
    {
      std::vector<UInt> terminals;
      for(auto &termRecv : islTerminalRecv.at(recv->idTrans()))
        for(UInt idTerm=0; idTerm<termRecv.size(); idTerm++)
          if(!isInList(terminals, termRecv.at(idTerm)))
            terminals.push_back(termRecv.at(idTerm));
      std::sort(terminals.begin(), terminals.end());

#if DEBUG_SYNC_ISL > 0
      if(Parallel::isMaster(comm))
        logWarning<<"synchronizeTransceiversIsl() recv ISL terminal "<<recv->name()<<" "
                  <<terminals.size()%"# terminals %i"s
                  <<Log::endl;
#endif

      // NOTE: a-priori ISL biases are NOT retained!
      recv->signalBiasIslRx.biases    = recv->signalBiasIslRx.compute(terminals); // apriori ISL bias
      recv->signalBiasIslRx.terminals = terminals;
      /*
      if(terminals.size())
      {
        recv->signalBiasIslRx.biases    = recv->signalBiasIslRx.compute(terminals); // apriori ISL bias
        recv->signalBiasIslRx.terminals = terminals;
      }
      */
    }

#if DEBUG_SYNC_ISL > 0
    if(Parallel::isMaster(comm))
    {
      for(auto transmitter : transmitters)
        for(UInt i=0; i<transmitter->signalBiasIslTx.terminals.size(); i++)
          logWarning<<"synchronizeTransceiversIsl() send ISL terminal bias "<<transmitter->name()<<" "
                    <<transmitter->signalBiasIslTx.terminals.at(i)%"%i"s<< " : "
                    <<transmitter->signalBiasIslTx.biases.at(i)%" %6.2f"s
                    <<Log::endl;

      for(auto transmitter : transmitters)
        for(UInt i=0; i<transmitter->signalBiasIslRx.terminals.size(); i++)
          logWarning<<"synchronizeTransceiversIsl() recv ISL terminal bias "<<transmitter->name()<<" "
                    <<transmitter->signalBiasIslRx.terminals.at(i)%"%i"s<< " : "
                    <<transmitter->signalBiasIslRx.biases.at(i)%" %6.2f"s
                    <<Log::endl;
    }
#endif

#if DEBUG_SYNC_ISL > 0
    logWarning<<"synchronizeTransceiversIsl() end"
              <<Log::endl;
#endif

  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Gnss::initParameter(GnssNormalEquationInfo &normalEquationInfo)
{
  try
  {
    logStatus<<"setup parameters"<<Log::endl;
    if(!parametrization)
      throw(Exception("no parametrization given"));

    // distribute process id of receivers and transmitters
    // ---------------------------------------------------
    Vector recvProcess(receivers.size());
    for(UInt idRecv=0; idRecv<receivers.size(); idRecv++)
      if(receivers.at(idRecv)->isMyRank())
        recvProcess(idRecv) = Parallel::myRank(normalEquationInfo.comm)+1;
    Parallel::reduceSum(recvProcess, 0, normalEquationInfo.comm);
    Parallel::broadCast(recvProcess, 0, normalEquationInfo.comm);

    Vector transProcess(transmitters.size());
    for(UInt idTrans=0; idTrans<transmitters.size(); idTrans++)
      if(transmitters.at(idTrans)->isMyRank())
        transProcess(idTrans) = Parallel::myRank(normalEquationInfo.comm)+1;
    Parallel::reduceSum(transProcess, 0, normalEquationInfo.comm);
    Parallel::broadCast(transProcess, 0, normalEquationInfo.comm);

    // Build connection matrix for receivers and transmitters
    // ------------------------------------------------------
    UInt nRecv  = receivers.size();
    UInt nTrans = transmitters.size();
    UInt nTotal = nRecv+nTrans;

    for(UInt idEpoch : normalEquationInfo.idEpochs)
    {
#if DEBUG > 10
      logStatus<<"setup links "<<times.at(idEpoch).dateTimeStr()<<Log::endl;
#endif
      links.clear();
      links.resize(nTotal);

      // GNSS observations transmitter -> receiver
      // -----------------------------------------
      for(const auto &recv : receivers)
      {
        if(normalEquationInfo.estimateReceiver.at(recv->idRecv()) && recv->useable(idEpoch) && recv->isMyRank())
          for(const auto &trans : transmitters)
            if(trans->useable(idEpoch) && recv->observation(trans->idTrans(),idEpoch))
            {
              links.at(recv->idRecv()).push_back(nRecv+trans->idTrans());
#if DEBUG > 10
              logWarning<<"Link  "<<times.at(idEpoch).dateTimeStr()<<" "<<recv->name()<<"<-"<<trans->name()<<Log::endl;
#endif
            }
        if(recv->useable(idEpoch))
          Parallel::broadCast(links.at(recv->idRecv()), static_cast<UInt>(recvProcess(recv->idRecv())-1), normalEquationInfo.comm);
      }
      // ISL observations transmitter -> transmitter
      // -------------------------------------------
      for(const auto &recv : transmitters)
      {
        if(recv->useable(idEpoch) && recv->isMyRank())
          for(const auto &trans : transmitters)
            if(trans->idTrans()!=recv->idTrans() && trans->useable(idEpoch) && recv->observationIsl(trans->idTrans(),idEpoch))
            {
              links.at(nRecv+recv->idTrans()).push_back(nRecv+trans->idTrans());
#if DEBUG > 10
              logWarning<<"Link  "<<times.at(idEpoch).dateTimeStr()<<" "<<recv->name()<<"<-"<<trans->name()<<Log::endl;
#endif
            }
        if(recv->useable(idEpoch))
          Parallel::broadCast(links.at(nRecv+recv->idTrans()), static_cast<UInt>(transProcess(recv->idTrans())-1), normalEquationInfo.comm);
      }
      Parallel::barrier(normalEquationInfo.comm);

      std::vector<UInt> Q;
      if(Parallel::isMaster(normalEquationInfo.comm))
      {

        // Select reference
        // ----------------
        int reference = -1;
        if(receivers.size())
        {
          for(const auto &recv : receivers)
            if(normalEquationInfo.estimateReceiver.at(recv->idRecv()) && recv->useable(idEpoch))
            {
              reference = recv->idRecv();
#if DEBUG > 10
              logWarningOnce<<"Pivot "<<times.at(idEpoch).dateTimeStr()<<" "<<recv->name()<<Log::endl;
#endif
              break;
            }
        }
        else if(transmitters.size())
        {
          for(const auto &trans : transmitters)
            if(trans->useable(idEpoch))
            {
              reference = nRecv+trans->idTrans();
#if DEBUG > 10
              logWarningOnce<<"Pivot "<<times.at(idEpoch).dateTimeStr()<<" "<<trans->name()<<Log::endl;
#endif
              break;
            }
        }
        else
          logWarningOnce<<"no receiver/satellite found as reference at "<<times.at(idEpoch).dateTimeStr()<<Log::endl;

        if (reference!=-1)
        {
          std::vector<std::vector<UInt>> graph(nTotal);
          if(reference != -1)
            for (UInt i=0; i<links.size(); i++)
              for (UInt j : links.at(i))
              {
                graph[i].push_back(j);
                graph[j].push_back(i);
              }

          Q = bfs(reference, graph);
        }

      } // if(Parallel::isMaster(normalEquationInfo.comm))

      Parallel::broadCast(Q, 0, normalEquationInfo.comm);

      for(const auto &recv : receivers)
        if(std::find(Q.begin(),Q.end(),recv->idRecv())==Q.end())
        {
          recv->disable(idEpoch, "insufficient observations");
#if DEBUG > 0
          logStatus<<"Disable "<<times.at(idEpoch).dateTimeStr()<<" "<<recv->name()<<Log::endl;
#endif
        }
#if DEBUG > 1
        else
          logStatus<<"Enable  "<<times.at(idEpoch).dateTimeStr()<<" "<<recv->name()<<Log::endl;
#endif

      for(const auto &trans : transmitters)
        if(std::find(Q.begin(),Q.end(),nRecv+trans->idTrans())==Q.end())
        {
          trans->disable(idEpoch, "insufficient observations");
#if DEBUG > 0
          logStatus<<"Disable "<<times.at(idEpoch).dateTimeStr()<<" "<<trans->name()<<Log::endl;
#endif
        }
#if DEBUG > 1
        else
          logStatus<<"Enable  "<<times.at(idEpoch).dateTimeStr()<<" "<<trans->name()<<Log::endl;
#endif

    } // for(UInt idEpoch : normalEquationInfo.idEpochs)

    // synchronize transceivers
    synchronizeTransceivers(normalEquationInfo.comm);
    synchronizeTransceiversIsl(normalEquationInfo.comm);

    // disable un-useable transmitters/receivers/epochs
    // ------------------------------------------------
    // check number of required observations
    std::vector<UInt> transCount(transmitters.size(), 0), transCountEpoch(transmitters.size(), 0);
    std::vector<UInt> recvCount(receivers.size(), 0),     recvCountEpoch(receivers.size(), 0);
    parametrization->requirements(normalEquationInfo, transCount, transCountEpoch, recvCount, recvCountEpoch);
    UInt disabledEpochsTrans = 0;
    UInt disabledEpochsRecv = 0;
    for(;;)
    {
      UInt mustSync = 0;

      // disable transmitters
      for(auto trans : transmitters)
        if(trans->useable() && (transCount.at(trans->idTrans()) || transCountEpoch.at(trans->idTrans())))
        {
          Vector countEpoch(times.size());
          for(const auto &recv : receivers)
            if(normalEquationInfo.estimateReceiver.at(recv->idRecv()) && recv->isMyRank())
              for(UInt idEpoch : normalEquationInfo.idEpochs)
                if(trans->useable(idEpoch) && recv->useable(idEpoch) && recv->observation(trans->idTrans(), idEpoch))
                  countEpoch(idEpoch)++;

          // Count ISL observations to other transmitters
          // TODO: this does not cover all necessary cases!
          for(const auto &recv : transmitters)
            if(trans->idTrans()!=recv->idTrans() && recv->isMyRank())
              for(UInt idEpoch : normalEquationInfo.idEpochs)
                if(trans->useable(idEpoch) && recv->useable(idEpoch) &&
                    (recv->observationIsl(trans->idTrans(), idEpoch) ||
                     trans->observationIsl(recv->idTrans(), idEpoch)))
                  countEpoch(idEpoch)++;
          Parallel::reduceSum(countEpoch, 0, normalEquationInfo.comm);
          Parallel::broadCast(countEpoch, 0, normalEquationInfo.comm);

          for(UInt idEpoch : normalEquationInfo.idEpochs)
            if(trans->useable(idEpoch) && (countEpoch(idEpoch) < transCountEpoch.at(trans->idTrans())))
            {
              disabledEpochsTrans++;
              trans->disable(idEpoch, "failed parametrization requirements");
              mustSync = TRUE;
              for(const auto &recv : receivers)
                if(recv->isMyRank() && recv->observation(trans->idTrans(), idEpoch))
                  recv->deleteObservation(trans->idTrans(), idEpoch);
              for(const auto &recv : transmitters)
                if(recv->isMyRank() && recv->observationIsl(trans->idTrans(), idEpoch))
                  recv->deleteObservationIsl(trans->idTrans(), idEpoch);
            }

          UInt epochCount = 0;
          for(UInt idEpoch : normalEquationInfo.idEpochs)
            if(countEpoch(idEpoch) > transCountEpoch.at(trans->idTrans()))
              epochCount++;
          if(epochCount < transCount.at(trans->idTrans()))
          {
            logWarningOnce<<trans->name()<<" disabled: not enough estimable epochs ("<<epochCount<<")"<<Log::endl;
            trans->disable("not enough estimable epochs ("+epochCount%"%i)"s);
            mustSync = TRUE;
            for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
            {
              for(const auto &recv : receivers)
                if(recv->isMyRank() && recv->observation(trans->idTrans(), idEpoch))
                  recv->deleteObservation(trans->idTrans(), idEpoch);
              for(const auto &recv : transmitters)
                if(recv->isMyRank() && recv->observationIsl(trans->idTrans(), idEpoch))
                  recv->deleteObservationIsl(trans->idTrans(), idEpoch);
            }
          }
        }

      // disable receivers
      for(auto recv : receivers)
        if(normalEquationInfo.estimateReceiver.at(recv->idRecv()) && recv->isMyRank() && (recvCount.at(recv->idRecv()) || recvCountEpoch.at(recv->idRecv())))
        {
          UInt epochCount = 0;
          for(UInt idEpoch : normalEquationInfo.idEpochs)
            if(recv->useable(idEpoch))
            {
              UInt count = 0;
              for(auto trans : transmitters)
                if(trans->useable(idEpoch) && recv->observation(trans->idTrans(), idEpoch))
                  count++;
              if(count < recvCountEpoch.at(recv->idRecv()))
              {
                disabledEpochsRecv++;
                recv->disable(idEpoch, "failed parametrization requirements");
                mustSync = TRUE;
              }
              else if(count > recvCountEpoch.at(recv->idRecv()))
                epochCount++;
            }

          if(epochCount < recvCount.at(recv->idRecv()))
          {
            logWarning<<recv->name()<<" disabled: not enough estimable epochs ("<<epochCount<<")"<<Log::endl;
            recv->disable("not enough estimable epochs ("+epochCount%"%i)"s);
            mustSync = TRUE;
          }
        }

      // is something disabled?
      Parallel::reduceSum(mustSync, 0, normalEquationInfo.comm);
      Parallel::broadCast(mustSync, 0, normalEquationInfo.comm);
      if(!mustSync)
        break;

      // synchronize transceivers
      synchronizeTransceivers(normalEquationInfo.comm);
      synchronizeTransceiversIsl(normalEquationInfo.comm);
      Parallel::reduceSum(disabledEpochsRecv, 0, normalEquationInfo.comm);
      if(!Parallel::isMaster(normalEquationInfo.comm))
        disabledEpochsRecv = 0;
    } // for(;;)

    if(disabledEpochsTrans) logWarningOnce<<disabledEpochsTrans<<" disabled transmitter epochs"<<Log::endl;
    if(disabledEpochsRecv)  logWarningOnce<<disabledEpochsRecv <<" disabled receiver epochs"<<Log::endl;

    for(auto recv : receivers)
      if(!recv->useable())
        normalEquationInfo.estimateReceiver.at(recv->idRecv()) = FALSE;

    // init parameters
    // ---------------
    normalEquationInfo.initNewParameterNames();
    parametrization->initParameter(normalEquationInfo);
    normalEquationInfo.calculateIndex(recvProcess);
    logInfo<<"+ ======="<<Log::endl;
    logInfo<<normalEquationInfo.parameterCount()%"%9i parameters in "s<<normalEquationInfo.blockCount()<<" normal equation matrix blocks"<<Log::endl;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector Gnss::aprioriParameter(const GnssNormalEquationInfo &normalEquationInfo) const
{
  try
  {
    Vector x0 = parametrization->aprioriParameter(normalEquationInfo);
    Parallel::reduceSum(x0, 0, normalEquationInfo.comm);
    return x0;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool Gnss::basicObservationEquations(const GnssNormalEquationInfo &/*normalEquationInfo*/, UInt idRecv, UInt idTrans, UInt idEpoch, GnssObservationEquation &eqn) const
{
  try
  {
    std::vector<GnssType> type;
    if(!receivers.at(idRecv)->observation(idTrans, idEpoch) ||
       !receivers.at(idRecv)->observation(idTrans, idEpoch)->observationList(GnssObservation::RANGE | GnssObservation::PHASE, type))
      return FALSE;
    eqn.compute(*receivers.at(idRecv)->observation(idTrans, idEpoch), *receivers.at(idRecv), *transmitters.at(idTrans),
                funcRotationCrf2Trf, funcReduceModels, idEpoch, TRUE, type);
    return TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool Gnss::basicObservationEquationsIsl(const GnssNormalEquationInfo &/*normalEquationInfo*/, UInt idRecv, UInt idTrans, UInt idEpoch, GnssObservationEquationIsl &eqn) const
{
  try
  {
    if(!transmitters.at(idRecv)->observationIsl(idTrans, idEpoch))
      return FALSE;
    eqn.compute(*transmitters.at(idRecv)->observationIsl(idTrans, idEpoch), *transmitters.at(idRecv), *transmitters.at(idTrans),
                funcReduceModelsIsl, idEpoch, TRUE);
    return TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Gnss::designMatrix(const GnssNormalEquationInfo &normalEquationInfo, const GnssObservationEquation &eqn, GnssDesignMatrix &A) const
{
  try
  {
    if(eqn.l.rows())
      parametrization->designMatrix(normalEquationInfo, eqn, A);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Gnss::designMatrixIsl(const GnssNormalEquationInfo &normalEquationInfo, const GnssObservationEquationIsl &eqn, GnssDesignMatrix &A) const
{
  try
  {
    if(eqn.l.rows())
      parametrization->designMatrixIsl(normalEquationInfo, eqn, A);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Gnss::constraintsEpoch(const GnssNormalEquationInfo &normalEquationInfo, UInt idEpoch, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const
{
  try
  {
    parametrization->constraintsEpoch(normalEquationInfo, idEpoch, normals, n, lPl, obsCount);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Gnss::constraints(const GnssNormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const
{
  try
  {
    parametrization->constraints(normalEquationInfo, normals, n, lPl, obsCount);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double Gnss::ambiguityResolve(const GnssNormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount,
                              const std::vector<Byte> &selectedTransmitters, const std::vector<Byte> &selectedReceivers,
                              const std::function<Vector(const_MatrixSliceRef xFloat, MatrixSliceRef W, const_MatrixSliceRef d, Vector &xInt, Double &sigma)> &searchInteger)
{
  try
  {
    return parametrization->ambiguityResolve(normalEquationInfo, normals, n, lPl, obsCount,
                                             selectedTransmitters, selectedReceivers, searchInteger);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double Gnss::updateParameter(const GnssNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef Wz)
{
  try
  {
    return parametrization->updateParameter(normalEquationInfo, x, Wz);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Gnss::updateCovariance(const GnssNormalEquationInfo &normalEquationInfo, const MatrixDistributed &covariance)
{
  try
  {
    parametrization->updateCovariance(normalEquationInfo, covariance);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Gnss::writeResults(const GnssNormalEquationInfo &normalEquationInfo, const std::string &suffix)
{
  try
  {
    parametrization->writeResults(normalEquationInfo, suffix);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

std::vector<GnssType> Gnss::types(const GnssType mask) const
{
  try
  {
    std::vector<GnssType> types;
    for(UInt idRecv=0; idRecv<receivers.size(); idRecv++)
      for(UInt idTrans=0; idTrans<transmitters.size(); idTrans++)
        for(GnssType type : typesRecvTrans.at(idRecv).at(idTrans))
          if(!type.isInList(types))
             types.push_back(type & mask);
    std::sort(types.begin(), types.end());
    return types;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::vector<GnssType> Gnss::typesIsl(const GnssType mask) const
{
  try
  {
    std::vector<GnssType> types;
    for(UInt idRecv=0; idRecv<transmitters.size(); idRecv++)
      for(UInt idTrans=0; idTrans<transmitters.size(); idTrans++)
        for(GnssType type : typesRecvTransIsl.at(idRecv).at(idTrans))
          if(!type.isInList(types))
             types.push_back(type & mask);
    std::sort(types.begin(), types.end());
    return types;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

std::vector<Byte> Gnss::selectTransmitters(PlatformSelectorPtr selector)
{
  try
  {
    std::vector<const Platform*> platforms(transmitters.size(), nullptr);
    for(UInt idTrans=0; idTrans<transmitters.size(); idTrans++)
      if(transmitters.at(idTrans)->useable())
        platforms.at(idTrans) = &transmitters.at(idTrans)->platform;
    return selector->select(times.front(), times.back(), platforms);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::vector<Byte> Gnss::selectReceivers(PlatformSelectorPtr selector)
{
  try
  {
    std::vector<const Platform*> platforms(receivers.size(), nullptr);
    for(UInt idRecv=0; idRecv<receivers.size(); idRecv++)
      if(receivers.at(idRecv)->useable())
        platforms.at(idRecv) = &receivers.at(idRecv)->platform;
    return selector->select(times.front(), times.back(), platforms);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

Bool Gnss::InfoParameterChange::update(Double change)
{
  try
  {
    count++;
    rms += change*change;
    if(std::fabs(change) > std::fabs(maxChange))
    {
      maxChange = change;
      return TRUE;
    }
    return FALSE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Gnss::InfoParameterChange::synchronizeAndPrint(Parallel::CommunicatorPtr comm, Double convertToMeter, Double &maxChangeTotal)
{
  try
  {
    Vector change(Parallel::size(comm));
    change(Parallel::myRank(comm)) = maxChange;
    Parallel::reduceSum(change, 0, comm);
    Parallel::broadCast(change, 0, comm);
    UInt idProcess = 0;
    maxChange = 0;
    for(UInt i=0; i<change.rows(); i++)
      if(std::fabs(change(i)) > std::fabs(maxChange))
      {
        maxChange = change(i);
        idProcess = i;
      }
    maxChangeTotal = std::max(maxChangeTotal, std::fabs(convertToMeter*maxChange));

    Parallel::reduceSum(count, 0, comm);
    Parallel::reduceSum(rms,   0, comm);

    if((idProcess != 0) && (idProcess == Parallel::myRank(comm)))
      Parallel::send(info, 0, comm);
    else if(Parallel::isMaster(comm))
    {
      if(idProcess != 0)
        Parallel::receive(info, idProcess, comm);
      rms = std::sqrt(rms/count);
      if(!info.empty())
      {
        std::string space(5-std::min(unit.size(), std::size_t(4)), ' ');
        logInfo<<"  rms ="<<rms%"%7.1f "s<<unit<<","<<space<<"max ="<<maxChange%"%8.1f "s<<unit<<","<<space<<info<<Log::endl;
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
