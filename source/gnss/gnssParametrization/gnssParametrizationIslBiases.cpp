/***********************************************/
/**
* @file gnssParametrizationIslBiases.cpp
*
* @brief Inter Satellite Link biases.
* @see GnssParametrization
*
* @author Torsten Mayer-Guerr
* @author Andre Hauschild
* @date 2024-01-30
*
*/
/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "classes/platformSelector/platformSelector.h"
#include "gnss/gnss.h"
#include "gnss/gnssParametrization/gnssParametrizationIslBiases.h"

/***********************************************/

GnssParametrizationIslBiases::GnssParametrizationIslBiases(Config &config)
{
  try
  {
    readConfig(config, "name",                          name,                       Config::OPTIONAL, "parameter.islBiases", "used for parameter selection");
    readConfig(config, "selectSendTerminal",            selectSendTerminal,         Config::DEFAULT,  R"(["all"])", "");
    readConfig(config, "selectRecvTerminal",            selectRecvTerminal,         Config::DEFAULT,  R"(["all"])", "");
    readConfig(config, "nameConstraint",                nameConstraint,             Config::OPTIONAL, "constraint.islBiases", "used for parameter selection");
    readConfig(config, "selectSendTerminalZeroMean",    selectSendTerminalZeroMean, Config::DEFAULT,  R"(["all"])", "");
    readConfig(config, "selectRecvTerminalZeroMean",    selectRecvTerminalZeroMean, Config::DEFAULT,  "", "");
    readConfig(config, "sigmaZeroMeanConstraint",       sigmaZeroMean,              Config::DEFAULT,  "0.0001", "(0 = unconstrained) sigma [m] for constraint over all selected biases");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

GnssParametrizationIslBiases::~GnssParametrizationIslBiases()
{
  for(auto para : paraSendTerminal)
    delete para;
  for(auto para : paraRecvTerminal)
    delete para;
}

/***********************************************/

void GnssParametrizationIslBiases::init(Gnss *gnss, Parallel::CommunicatorPtr comm)
{
  try
  {
    this->gnss = gnss;

    logInfo<<"GnssParametrizationIslBiases::init() init a-priori values"<<Log::endl;

    auto selectedSendTerminal = gnss->selectTransmitters(selectSendTerminal);
    paraSendTerminal.resize(gnss->transmitters.size(), nullptr);
    for(UInt idTrans=0; idTrans<gnss->transmitters.size(); idTrans++)
      if(selectedSendTerminal.at(idTrans) && gnss->transmitters.at(idTrans)->useable())
      {
        auto para = new Parameter();
        paraSendTerminal.at(idTrans) = para;
        para->trans = gnss->transmitters.at(idTrans);
        logInfo<<"Init send ISL terminal bias parameter "<<para->trans->name()<<Log::endl;
      }

    auto selectedRecvTerminal = gnss->selectTransmitters(selectRecvTerminal);
    paraRecvTerminal.resize(gnss->transmitters.size(), nullptr);
    for(UInt idRecv=0; idRecv<gnss->transmitters.size(); idRecv++)
      if(selectedRecvTerminal.at(idRecv) && gnss->transmitters.at(idRecv)->useable())
      {
        auto para = new Parameter();
        paraRecvTerminal.at(idRecv) = para;
        para->trans = gnss->transmitters.at(idRecv);
        logInfo<<"Init recv ISL terminal bias parameter "<<para->trans->name()<<Log::endl;
      }

  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationIslBiases::initParameter(GnssNormalEquationInfo &normalEquationInfo)
{
  try
  {
    for(auto para : paraSendTerminal)
      if(para)
        para->index = GnssParameterIndex();
    for(auto para : paraRecvTerminal)
      if(para)
        para->index = GnssParameterIndex();

    applyConstraint = FALSE;
    if(!isEnabled(normalEquationInfo, name) || normalEquationInfo.isEachReceiverSeparately)
      return;

    // transmitter ISL terminal parameters
    // -----------------------------------
    UInt countparaSendTerminal = 0;
    for(auto para : paraSendTerminal)
      if(para && para->trans->useable() && para->trans->signalBiasIslTx.biases.size())
      {
        logInfo<<"initParameter() send ISL terminal bias parameter "<<para->trans->name()<<para->trans->signalBiasIslTx.biases.size()%" (%02i)"s<<Log::endl;
        // determine parameter names
        std::vector<ParameterName> parameterNames;
        parameterNames.push_back(ParameterName(para->trans->name(), "sendIslTerminalBias"));
        para->index = normalEquationInfo.parameterNamesTransmitter(para->trans->idTrans(), parameterNames);
        countparaSendTerminal += parameterNames.size();
      }
    if(countparaSendTerminal)
      logInfo<<countparaSendTerminal%"%9i send ISL terminal bias parameters"s<<Log::endl;

    // receiving ISL terminal parameters
    // ---------------------------------
    UInt countparaRecvTerminal = 0;
    for(auto para : paraRecvTerminal)
      if(para && para->trans->useable() && para->trans->signalBiasIslRx.biases.size())
      {
        logInfo<<"initParameter() recv ISL terminal bias parameter "<<para->trans->name()<<Log::endl;
        // determine parameter names
        std::vector<ParameterName> parameterNames;
        parameterNames.push_back(ParameterName(para->trans->name(), "recvIslTerminalBias"));
        para->index = normalEquationInfo.parameterNamesTransmitter(para->trans->idTrans(), parameterNames);
        countparaRecvTerminal += parameterNames.size();
      }
    if(countparaRecvTerminal)
      logInfo<<countparaRecvTerminal%"%9i recv ISL terminal bias parameters"s<<Log::endl;

    // Zero-mean constraint

    selectedSendTerminalZeroMean = gnss->selectTransmitters(selectSendTerminalZeroMean);
    selectedRecvTerminalZeroMean = gnss->selectTransmitters(selectRecvTerminalZeroMean);

    // Initialize a-priori values

    // NOTE: the a-priori values for the zero-mean constrained are stored in
    //       this function and used in the constraints() function. It cannot be
    //       done in the init() function since the a-priori values may not yet
    //       initialized at that point depending on the parametrization order
    //       of gnssParametrizationislBiases and gnssParametrizationSignalBiasesIsl.

    UInt countZeroMean = 0;
    x0SendTerminal.resize(gnss->transmitters.size());
    x0RecvTerminal.resize(gnss->transmitters.size());
    for(auto trans : gnss->transmitters)
    {
      const UInt idTrans = trans->idTrans();
      if(trans->useable() && selectedSendTerminalZeroMean.at(idTrans))
      {
        x0SendTerminal.at(idTrans) = trans->signalBiasesIslTx();
        countZeroMean++;
        logInfo<<"initParameter() store send ISL terminal bias a-priori "<<trans->name()<<x0SendTerminal.at(idTrans)%" %6.2f m"s<<Log::endl;
      }
      if(trans->useable() && selectedRecvTerminalZeroMean.at(idTrans))
      {
        x0RecvTerminal.at(idTrans) = trans->signalBiasesIslRx();
        countZeroMean++;
        logInfo<<"initParameter() store recv ISL terminal bias a-priori "<<trans->name()<<x0RecvTerminal.at(idTrans)%" %6.2f m"s<<Log::endl;
      }
    }

    applyConstraint = isEnabled(normalEquationInfo, nameConstraint) && sigmaZeroMean
                      && (countparaSendTerminal+countparaRecvTerminal)
                      && countZeroMean;

  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationIslBiases::aprioriParameter(const GnssNormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const
{
  try
  {
    if(Parallel::isMaster(normalEquationInfo.comm))
    {
      for(auto para : paraSendTerminal)
        if(para && para->index)
        {
          logInfo<<"aprioriParameter() send ISL terminal bias  "<<para->trans->name() << para->trans->signalBiasesIslTx()%" %6.2f"s <<Log::endl;
          x0(normalEquationInfo.index(para->index),0) = para->trans->signalBiasesIslTx();
        }
      for(auto para : paraRecvTerminal)
        if(para && para->index)
        {
          logInfo<<"aprioriParameter() recv ISL terminal bias  "<<para->trans->name() << para->trans->signalBiasesIslRx()%" %6.2f"s <<Log::endl;
          x0(normalEquationInfo.index(para->index),0) = para->trans->signalBiasesIslRx();
        }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationIslBiases::designMatrixIsl(const GnssNormalEquationInfo &/*normalEquationInfo*/, const GnssObservationEquationIsl &eqn, GnssDesignMatrix &A) const
{
  try
  {
#if DEBUG >0
    logInfo << "GnssParametrizationIslBiases::designMatrixIsl() "
            << "idRecv "  << eqn.receiver->idTrans()
            <<" idTrans " << eqn.transmitter->idTrans()
            <<Log::endl;
#endif
    auto paraSendTerminal = this->paraSendTerminal.at(eqn.transmitter->idTrans());
    if(paraSendTerminal && paraSendTerminal->index)
      axpy(1., eqn.A.column(GnssObservationEquationIsl::idxRange), A.column(paraSendTerminal->index));

    auto paraRecvTerminal = this->paraRecvTerminal.at(eqn.receiver->idTrans());
    if(paraRecvTerminal && paraRecvTerminal->index)
      axpy(1., eqn.A.column(GnssObservationEquationIsl::idxRange), A.column(paraRecvTerminal->index));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationIslBiases::constraints(const GnssNormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const
{
  try
  {
    if(!applyConstraint)
      return;

    logStatus<<"apply zero mean equations for ISL terminal bias parameters"<<Log::endl;

    // zero-mean constraint of ISL terminal biases
    // -------------------------------------------

    UInt countTrans = 0, countRecv = 0;
    for(UInt idTrans=0; idTrans<paraSendTerminal.size(); idTrans++)
      if(selectedSendTerminalZeroMean.at(idTrans) && paraSendTerminal.at(idTrans) && paraSendTerminal.at(idTrans)->index)
        countTrans++;
    for(UInt idTrans=0; idTrans<paraRecvTerminal.size(); idTrans++)
      if(selectedRecvTerminalZeroMean.at(idTrans) && paraRecvTerminal.at(idTrans) && paraRecvTerminal.at(idTrans)->index)
        countRecv++;
    UInt count = countTrans + countRecv;
    if(!count)
      return;

    // collect apriori bias values
    // ---------------------------

    if(Parallel::isMaster(normalEquationInfo.comm))
    {
      GnssDesignMatrix A(normalEquationInfo, Vector(1));
      for(auto para : paraSendTerminal)
        if(para && para->index && selectedSendTerminalZeroMean.at(para->trans->idTrans()))
        {
          A.l(0) += (x0SendTerminal.at(para->trans->idTrans()) - para->trans->signalBiasesIslTx())/count/sigmaZeroMean; // remove apriori value -> regularization towards 0
          A.column(para->index)(0,0) = 1./count/sigmaZeroMean;
          logInfo<<"constraint() send ISL terminal bias "<<para->trans->name()
                 <<x0SendTerminal.at(para->trans->idTrans())%" a-priori %6.2f m"s
                 <<para->trans->signalBiasesIslTx()%" estimate %6.2f m"s
                 <<Log::endl;
        }
      for(auto para : paraRecvTerminal)
        if(para && para->index && selectedRecvTerminalZeroMean.at(para->trans->idTrans()))
        {
          A.l(0) += (x0RecvTerminal.at(para->trans->idTrans()) - para->trans->signalBiasesIslRx())/count/sigmaZeroMean; // remove apriori value -> regularization towards 0
          A.column(para->index)(0,0) = 1./count/sigmaZeroMean;
          logInfo<<"constraint() recv ISL terminal bias "<<para->trans->name()
                 <<x0RecvTerminal.at(para->trans->idTrans())%" a-priori %6.2f m"s
                 <<para->trans->signalBiasesIslRx()%" estimate %6.2f m"s
                 <<Log::endl;
        }
      A.accumulateNormals(normals, n, lPl, obsCount);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double GnssParametrizationIslBiases::updateParameter(const GnssNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef /*Wz*/)
{
  try
  {
    Double maxChange = 0;
    Gnss::InfoParameterChange infoTrans("mm");
    for(auto para : paraSendTerminal)
      if(para && para->index)
      {
        const Vector dBias = x.row(normalEquationInfo.index(para->index), 1);
        for(UInt idType=0; idType<dBias.size(); idType++)
          para->trans->signalBiasIslTx.biases.at(idType) += dBias(idType);
        for(UInt idType=0; idType<dBias.size(); idType++)
          if(infoTrans.update(1e3*dBias(idType)))
            infoTrans.info = "send ISL terminal bias ("+para->trans->signalBiasIslTx.types.at(idType).str()+")";
      }
    infoTrans.synchronizeAndPrint(normalEquationInfo.comm, 1e-3, maxChange);

    Gnss::InfoParameterChange infoRecv("mm");
    for(auto para : paraRecvTerminal)
      if(para && para->index)
      {
        const Vector dBias = x.row(normalEquationInfo.index(para->index), 1);
        for(UInt idType=0; idType<dBias.size(); idType++)
          para->trans->signalBiasIslRx.biases.at(idType) += dBias(idType);
        for(UInt idType=0; idType<dBias.size(); idType++)
          if(infoRecv.update(1e3*dBias(idType)))
            infoTrans.info = "recv ISL terminal bias ("+para->trans->signalBiasIslRx.types.at(idType).str()+")";
      }
    infoRecv.synchronizeAndPrint(normalEquationInfo.comm, 1e-3, maxChange);

    return maxChange;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
