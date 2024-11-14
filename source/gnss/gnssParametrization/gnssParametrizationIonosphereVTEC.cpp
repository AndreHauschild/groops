/***********************************************/
/**
* @file gnssParametrizationIonosphereVTEC.cpp
*
* @brief IonosphereVTEC.
* @see GnssParametrization
*
* @author Torsten Mayer-Guerr
* @date 2021-01-23
*
*/
/***********************************************/

#include "base/import.h"
#include "base/planets.h"
#include "config/config.h"
#include "files/fileInstrument.h"
#include "classes/magnetosphere/magnetosphere.h"
#include "gnss/gnssParametrization/gnssParametrizationIonosphereVTEC.h"

/***********************************************/

GnssParametrizationIonosphereVTEC::GnssParametrizationIonosphereVTEC(Config &config)
{
  try
  {
    readConfig(config, "name",                   name,                    Config::OPTIONAL, "parameter.VTEC", "");
    readConfig(config, "selectReceivers",        selectReceivers,         Config::MUSTSET,   R"(["all"])", "");
    readConfig(config, "outputfileVTEC",         fileNameVTEC,            Config::OPTIONAL, "output/vtec_{loopTime:%D}.{station}.dat", "variable {station} available");
    readConfig(config, "mapR",                   mapR,                    Config::DEFAULT,  "6371e3",     "constant of MSLM mapping function");
    readConfig(config, "mapH",                   mapH,                    Config::DEFAULT,  "506.7e3",    "constant of MSLM mapping function");
    readConfig(config, "mapAlpha",               mapAlpha,                Config::DEFAULT,  "0.9782",     "constant of MSLM mapping function");
    readConfig(config, "vtecGradientEstimation", parametrizationGradient, Config::DEFAULT,  "",  "[degree] parametrization of north and east gradients");
    readConfig(config, "magnetosphere",          magnetosphere,           Config::MUSTSET,  "",  "");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationIonosphereVTEC::init(Gnss *gnss, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    this->gnss = gnss;
    selectedReceivers = gnss->selectReceivers(selectReceivers);
    index.resize(gnss->receivers.size());
    VTEC.resize(gnss->receivers.size());

    indexGradient.resize(gnss->receivers.size());
    xGradient.resize(gnss->receivers.size());
    gradientX.resize(gnss->receivers.size());
    gradientY.resize(gnss->receivers.size());

    for(auto recv : gnss->receivers)
    {
      const UInt idRecv = recv->idRecv();
      if(recv->useable() && recv->isMyRank())
      {
        VTEC.at(idRecv).resize(gnss->times.size(), 0);
        xGradient.at(idRecv) = Vector(2*parametrizationGradient->parameterCount());
        gradientX.at(idRecv).resize(gnss->times.size(), 0);
        gradientY.at(idRecv).resize(gnss->times.size(), 0);
      }
    }

  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationIonosphereVTEC::initParameter(GnssNormalEquationInfo &normalEquationInfo)
{
  try
  {
    index.clear();
    index.resize(gnss->receivers.size());
    indexGradient.clear();
    indexGradient.resize(gnss->receivers.size(),GnssParameterIndex());
    if(!isEnabled(normalEquationInfo, name))
      return;

    UInt countPara = 0;
    for(auto recv : gnss->receivers)
    {
      const UInt idRecv = recv->idRecv();
      if(recv->useable() && normalEquationInfo.estimateReceiver.at(idRecv) && selectedReceivers.at(idRecv))
      {
        index.at(idRecv).resize(gnss->times.size());
        for(UInt idEpoch : normalEquationInfo.idEpochs)
          if(recv->useable(idEpoch))
          {
            index.at(idRecv).at(idEpoch) = normalEquationInfo.parameterNamesEpochReceiver(idEpoch, idRecv, {ParameterName(recv->name(), "VTEC", "", gnss->times.at(idEpoch))});
            countPara++;
          }
      }
    }
    if(countPara)
      logInfo<<countPara%"%9i VTEC epoch parameters"s<<Log::endl;

    // VTEC gradient
    UInt countParaGradient = 0;
    if(parametrizationGradient->parameterCount())
      for(auto recv : gnss->receivers)
      {
        const UInt idRecv = recv->idRecv();
        if(recv->useable() && normalEquationInfo.estimateReceiver.at(idRecv) && selectedReceivers.at(idRecv))
        {
          std::vector<ParameterName> parameterNames;
          std::vector<ParameterName> name({{gnss->receivers.at(idRecv)->name(), "vtecGradient.x"}, {gnss->receivers.at(idRecv)->name(), "vtecGradient.y"}});
          parametrizationGradient->parameterName(name, parameterNames);
          indexGradient.at(idRecv) = normalEquationInfo.parameterNamesReceiver(idRecv, parameterNames);
          countParaGradient += parameterNames.size();
        }
      }
    if(countParaGradient)
      logInfo<<countParaGradient%"%9i VTEC gradient parameters"s<<Log::endl;

  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationIonosphereVTEC::aprioriParameter(const GnssNormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const
{
  try
  {
    for(auto recv : gnss->receivers)
      if(recv->isMyRank())
      {
        const UInt idRecv = recv->idRecv();
        for(UInt idEpoch : normalEquationInfo.idEpochs)
          if(index.at(idRecv).size() && index.at(idRecv).at(idEpoch))
            x0(normalEquationInfo.index(index.at(idRecv).at(idEpoch)), 0) = VTEC.at(idRecv).at(idEpoch);
        if(indexGradient.at(idRecv))
          copy(xGradient.at(idRecv), x0.row(normalEquationInfo.index(indexGradient.at(idRecv)), xGradient.at(idRecv).rows()));
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

// Mapping function
Double GnssParametrizationIonosphereVTEC::mapping(Angle elevation) const
{
  return 1./std::cos(std::asin(mapR/(mapR+mapH) * std::sin(mapAlpha*(PI/2-elevation))));
}

/***********************************************/

// intersection point in ionosphere height
Vector3d GnssParametrizationIonosphereVTEC::intersection(const Double radiusIono, const Vector3d &posRecv, const Vector3d &posTrans) const
{
  try
  {
    const Double   rRecv = std::min(posRecv.r(), radiusIono); // LEO satellites flying possibly higher
    const Vector3d k     = normalize(posRecv-posTrans);       // direction from transmitter
    const Double   rk    = inner(posRecv, k);
    return posRecv - (std::sqrt(rk*rk+radiusIono*radiusIono-rRecv*rRecv)+rk) * k;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

// Mapping function for gradient
void GnssParametrizationIonosphereVTEC::mappingGradient(const GnssObservationEquation &eqn, Double &dx, Double &dy) const
{
  try
  {
    /*
    const Double M_h = 1.0/(sin(eqn.elevationRecvLocal)*tan(eqn.elevationRecvLocal) + 3*(mapH+mapR)/(2*(mapH+mapR)+eqn.posRecv.r()));
    dx = M_h*sin(eqn.azimutRecvLocal);  // East
    dy = M_h*cos(eqn.azimutRecvLocal);  // North
     */
    dx = sin(eqn.azimutRecvLocal);  // East
    dy = cos(eqn.azimutRecvLocal);  // North

  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationIonosphereVTEC::designMatrix(const GnssNormalEquationInfo &/*normalEquationInfo*/, const GnssObservationEquation &eqn, GnssDesignMatrix &A) const
{
  try
  {
    if(!index.at(eqn.receiver->idRecv()).size() || !index.at(eqn.receiver->idRecv()).at(eqn.idEpoch))
      return;

    // VTEC at station per epoch
    const Double map_vtec = mapping(eqn.elevationRecvLocal);
    axpy(map_vtec, eqn.A.column(GnssObservationEquation::idxSTEC), A.column(index.at(eqn.receiver->idRecv()).at(eqn.idEpoch)));

    // temporal parametrization
    auto designMatrixTemporal = [&](ParametrizationTemporalPtr parametrization, const_MatrixSliceRef B, const GnssParameterIndex &index)
    {
      std::vector<UInt>   idx;
      std::vector<Double> factor;
      parametrization->factors(std::max(eqn.timeRecv, gnss->times.at(0)), idx, factor);
      MatrixSlice Design(A.column(index));
      for(UInt i=0; i<factor.size(); i++)
        axpy(factor.at(i), B, Design.column(B.columns()*idx.at(i), B.columns()));
    };

    // VTEC gradient at station
    if(indexGradient.at(eqn.receiver->idRecv()))
    {
      Double dx, dy;
      mappingGradient(eqn, dx, dy);
      Matrix B(eqn.A.rows(), 2);
      axpy(map_vtec*dx, eqn.A.column(GnssObservationEquation::idxSTEC), B.column(0));
      axpy(map_vtec*dy, eqn.A.column(GnssObservationEquation::idxSTEC), B.column(1));
      designMatrixTemporal(parametrizationGradient, B, indexGradient.at(eqn.receiver->idRecv()));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double GnssParametrizationIonosphereVTEC::updateParameter(const GnssNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef /*Wz*/)
{
  try
  {
    if(!isEnabled(normalEquationInfo, name))
      return 0;

    Double minVTEC   = 1e+99;
    Double maxVTEC   = 1e-99;
    Double meanVTEC  = 0;
    Double stdVTEC   = 0;
    UInt   countVTEC = 0;

    Gnss::InfoParameterChange info("tec");
    for(auto recv : gnss->receivers)
      if(recv->isMyRank())
      {
        const UInt idRecv = recv->idRecv();
        for(UInt idEpoch : normalEquationInfo.idEpochs)
          if(index.at(idRecv).size() && index.at(idRecv).at(idEpoch))
          {
            // update VTEC
            const Double dVTEC = x(normalEquationInfo.index(index.at(idRecv).at(idEpoch)), 0);
            VTEC.at(idRecv).at(idEpoch) += dVTEC;
            if(info.update(dVTEC))
              info.info = "VTEC ("+recv->name()+", "+gnss->times.at(idEpoch).dateTimeStr()+")";

            minVTEC   = std::min(VTEC.at(idRecv).at(idEpoch), minVTEC);
            maxVTEC   = std::max(VTEC.at(idRecv).at(idEpoch), maxVTEC);
            meanVTEC += VTEC.at(idRecv).at(idEpoch);
            stdVTEC  += VTEC.at(idRecv).at(idEpoch)*VTEC.at(idRecv).at(idEpoch);
            countVTEC++;

            // update STEC
            for(UInt idTrans=0; idTrans<recv->idTransmitterSize(idEpoch); idTrans++)
              if(recv->observation(idTrans, idEpoch))
              {
                GnssObservationEquation eqn(*recv->observation(idTrans, idEpoch), *recv, *gnss->transmitters.at(idTrans), gnss->funcRotationCrf2Trf,
                                            nullptr/*reduceModels*/, idEpoch, FALSE/*decorrelate*/, {}/*types*/);
                recv->observation(idTrans, idEpoch)->STEC += mapping(eqn.elevationRecvLocal) * dVTEC;
              }
          }
      }

    // VTEC statistics
    Parallel::reduceMin(minVTEC,   0, normalEquationInfo.comm);
    Parallel::reduceMax(maxVTEC,   0, normalEquationInfo.comm);
    Parallel::reduceSum(meanVTEC,  0, normalEquationInfo.comm);
    Parallel::reduceSum(stdVTEC,   0, normalEquationInfo.comm);
    Parallel::reduceSum(countVTEC, 0, normalEquationInfo.comm);
    stdVTEC   = std::sqrt((stdVTEC-meanVTEC*meanVTEC/countVTEC)/(countVTEC-1));
    meanVTEC /= countVTEC;
    std::string infoStr = " (total: "+meanVTEC%"%.2f +- "s+stdVTEC%"%.2f ["s+minVTEC%"%.2f -- "s+maxVTEC%"%.2f])"s;
    Parallel::broadCast(infoStr, 0, normalEquationInfo.comm);
    info.info += infoStr;

    // Gradient update

    Double maxChange = 0;
    info.synchronizeAndPrint(normalEquationInfo.comm, 0, maxChange);

    // update VTEC gradient
    Gnss::InfoParameterChange infoGradient("tec");
    for(auto recv : gnss->receivers)
      if(recv->isMyRank())
      {
        const UInt idRecv = recv->idRecv();
        if(indexGradient.at(idRecv))
        {
          xGradient.at(idRecv) += x.row(normalEquationInfo.index(indexGradient.at(idRecv)), 2*parametrizationGradient->parameterCount());
          std::vector<UInt>   index;
          std::vector<Double> factor;
          for(UInt idEpoch=0; idEpoch<gnss->times.size(); idEpoch++)
          {
            parametrizationGradient->factors(std::max(recv->timeCorrected(idEpoch), gnss->times.at(0)), index, factor);
            Double gx=0, gy=0;
            for(UInt k=0; k<factor.size(); k++)
            {
              gx += factor.at(k) * xGradient.at(idRecv)(2*index.at(k)+0);
              gy += factor.at(k) * xGradient.at(idRecv)(2*index.at(k)+1);
            }
            const Double dgx = gx - gradientX.at(idRecv).at(idEpoch);
            const Double dgy = gy - gradientY.at(idRecv).at(idEpoch);
            gradientX.at(idRecv).at(idEpoch) = gx;
            gradientY.at(idRecv).at(idEpoch) = gy;
            const Double dxdy = std::sqrt(dgx*dgx+dgy*dgy);
            if(infoGradient.update(dxdy))
              infoGradient.info = "VTEC gradient ("+recv->name()+", "+gnss->times.at(idEpoch).dateTimeStr()+")";

            // update STEC
            for(UInt idTrans=0; idTrans<recv->idTransmitterSize(idEpoch); idTrans++)
              if(recv->observation(idTrans, idEpoch))
              {
                GnssObservationEquation eqn(*recv->observation(idTrans, idEpoch), *recv, *gnss->transmitters.at(idTrans), gnss->funcRotationCrf2Trf,
                    nullptr/*reduceModels*/, idEpoch, FALSE/*decorrelate*/, {}/*types*/);
                Double dx,dy;
                mappingGradient(eqn, dx, dy);
                recv->observation(idTrans, idEpoch)->STEC += mapping(eqn.elevationRecvLocal) * (dx*dgx + dy*dgy);
              }
          }
        }
      }
    infoGradient.synchronizeAndPrint(normalEquationInfo.comm, 0, maxChange);
    return maxChange;

  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationIonosphereVTEC::writeResults(const GnssNormalEquationInfo &normalEquationInfo, const std::string &suffix) const
{
  try
  {
    if(!isEnabled(normalEquationInfo, name))
      return;

    if(!fileNameVTEC.empty() && VTEC.size())
    {
      VariableList fileNameVariableList;
      fileNameVariableList.setVariable("station", "****");
      logStatus<<"write VTEC to files <"<<fileNameVTEC(fileNameVariableList).appendBaseName(suffix)<<">"<<Log::endl;

      for(auto &recv : gnss->receivers)
        if(normalEquationInfo.estimateReceiver.at(recv->idRecv()) && recv->isMyRank())
        {
          MiscValuesArc arc;
          for(UInt idEpoch : normalEquationInfo.idEpochs)
          {
            if(index.at(recv->idRecv()).size() && index.at(recv->idRecv()).at(idEpoch))
            {
              Vector data(3);
              data.at(0) = VTEC.at(recv->idRecv()).at(idEpoch);
              data.at(1) = gradientX.at(recv->idRecv()).at(idEpoch);
              data.at(2) = gradientY.at(recv->idRecv()).at(idEpoch);
              MiscValuesEpoch epoch(data.size());
              epoch.time  = gnss->times.at(idEpoch);
              epoch.setData(data);
              arc.push_back(epoch);
            }
          }
          fileNameVariableList.setVariable("station", recv->name());
          if(arc.size())
            InstrumentFile::write(fileNameVTEC(fileNameVariableList).appendBaseName(suffix), arc);
        }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
