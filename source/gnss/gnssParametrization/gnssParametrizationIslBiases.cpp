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

#define DEBUG 1

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
    readConfig(config, "name",                           name,                           Config::OPTIONAL, "parameter.islTerminalBiases", "used for parameter selection");
    readConfig(config, "selectTransmitIslTerminal",      selectTransmitTerminal,         Config::DEFAULT,  R"(["all"])", "");
    readConfig(config, "selectReceiveIslTerminal",       selectReceiveTerminal,          Config::DEFAULT,  R"(["all"])", "");
    readConfig(config, "nameConstraint",                 nameConstraint,                 Config::OPTIONAL, "constraint.islBiases", "used for parameter selection");
    readConfig(config, "selectTransmitTerminalZeroMean", selectTransmitTerminalZeroMean, Config::DEFAULT,  R"(["all"])", "");
    readConfig(config, "selectReceiveTerminalZeroMean",  selectReceiveTerminalZeroMean,  Config::DEFAULT,  "", "");
    readConfig(config, "sigmaZeroMeanConstraint",        sigmaZeroMean,                  Config::DEFAULT,  "0.0001", "(0 = unconstrained) sigma [m] for constraint over all selected biases");
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
  for(auto para : paraTransmitTerminal)
    delete para;
  for(auto para : paraReceiveTerminal)
    delete para;
}

/***********************************************/

// NOTE: ISL terminal biases can only be estimated if ISL observations are
//       available. This is checked in Gnss::synchronizeTransceiversIsl() and
//       if no observations are available, the ISL bias lists are set to zero
//       length. This happens after init() and before initParameter() is called.
//       Therefore the bias parameters for all transmit and receive terminals
//       are created here, but the parameter indices are created in the call of
//       initParameter() depending on the availability of ISL observations.

void GnssParametrizationIslBiases::init(Gnss *gnss, Parallel::CommunicatorPtr comm)
{
  try
  {
    this->gnss = gnss;

    auto selectedTransmitTerminal = gnss->selectTransmitters(selectTransmitTerminal);
    paraTransmitTerminal.resize(gnss->transmitters.size(), nullptr);
    for(UInt idTrans=0; idTrans<gnss->transmitters.size(); idTrans++)
      if(selectedTransmitTerminal.at(idTrans) && gnss->transmitters.at(idTrans)->useable())
      {
        auto para = new Parameter();
        paraTransmitTerminal.at(idTrans) = para;
        para->trans = gnss->transmitters.at(idTrans);
      }

    auto selectedReceiveTerminal = gnss->selectTransmitters(selectReceiveTerminal);
    paraReceiveTerminal.resize(gnss->transmitters.size(), nullptr);
    for(UInt idRecv=0; idRecv<gnss->transmitters.size(); idRecv++)
      if(selectedReceiveTerminal.at(idRecv) && gnss->transmitters.at(idRecv)->useable())
      {
        auto para = new Parameter();
        paraReceiveTerminal.at(idRecv) = para;
        para->trans = gnss->transmitters.at(idRecv);
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
    for(auto para : paraTransmitTerminal)
      if(para)
        para->index = GnssParameterIndex();
    for(auto para : paraReceiveTerminal)
      if(para)
        para->index = GnssParameterIndex();

    applyConstraint = FALSE;
    if(!isEnabled(normalEquationInfo, name) || normalEquationInfo.isEachReceiverSeparately)
      return;

    // transmitting ISL terminal parameters
    // ------------------------------------
    UInt countparaTransmitTerminal = 0;
    for(auto para : paraTransmitTerminal)
      if(para && para->trans->useable() && para->trans->signalBiasIslTx.terminals.size())
      {

        if(!para->trans->signalBiasIslTx.biases.size())
          continue;

        // determine parameter names
        std::vector<ParameterName> parameterNames;
        for(UInt i=0; i<para->trans->signalBiasIslTx.biases.size(); i++)
        {
          parameterNames.push_back(ParameterName(para->trans->name(), "islBiasTx"+i%"%02i"s));
#if DEBUG > 0
          logInfo<<"initParameter() transmit ISL terminal bias parameter "<<para->trans->name()<<":"<<"islBiasTx"+i%"%02i"s<<Log::endl;
#endif
        }
        para->index = normalEquationInfo.parameterNamesTransmitter(para->trans->idTrans(), parameterNames);
        countparaTransmitTerminal += parameterNames.size();
      }
    if(countparaTransmitTerminal)
      logInfo<<countparaTransmitTerminal%"%9i ISL transmit terminal bias parameters"s<<Log::endl;

    // receiving ISL terminal parameters
    // ---------------------------------
    UInt countparaReceiveTerminal = 0;
    for(auto para : paraReceiveTerminal)
      if(para && para->trans->useable() && para->trans->signalBiasIslRx.terminals.size())
      {

        if(!para->trans->signalBiasIslRx.terminals.size())
          continue;

        // determine parameter names
        std::vector<ParameterName> parameterNames;
        for(UInt i=0; i<para->trans->signalBiasIslRx.terminals.size(); i++)
        {
          parameterNames.push_back(ParameterName(para->trans->name(), "islBiasRx"+i%"%02i"s));
#if DEBUG > 0
          logInfo<<"initParameter() receive ISL terminal bias parameter "<<para->trans->name()<<":"<<"islBiasRx"+i%"%02i"s<<Log::endl;
#endif
        }
        para->index = normalEquationInfo.parameterNamesTransmitter(para->trans->idTrans(), parameterNames);
        countparaReceiveTerminal += parameterNames.size();
      }
    if(countparaReceiveTerminal)
      logInfo<<countparaReceiveTerminal%"%9i ISL receive terminal bias parameters"s<<Log::endl;

    // Zero-mean constraint

    selectedTransmitTerminalZeroMean = gnss->selectTransmitters(selectTransmitTerminalZeroMean);
    selectedReceiveTerminalZeroMean = gnss->selectTransmitters(selectReceiveTerminalZeroMean);

    // Initialize a-priori values

    // NOTE: the a-priori values for the zero-mean constraint are stored in
    //       this function and used in the constraints() function. It cannot be
    //       done in the init() function, since the a-priori values may not yet
    //       initialized at that point depending on the parametrization order of
    //       gnssParametrizationislBiases and gnssParametrizationSignalBiasesIsl

    UInt countZeroMean = 0;
    x0TransmitTerminal.resize(gnss->transmitters.size());
    for(auto para : paraTransmitTerminal)
      if(para && para->index && selectedTransmitTerminalZeroMean.at(para->trans->idTrans()))
      {
        x0TransmitTerminal.at(para->trans->idTrans()) = Vector(para->trans->signalBiasIslTx.biases);
        countZeroMean++;
#if DEBUG >0
        for(UInt i=0; i<para->trans->signalBiasIslTx.biases.size(); i++)
          logInfo<<"initParameter() store initial transmit ISL terminal bias parameter "
                 <<para->trans->name()<<para->trans->signalBiasIslTx.biases.at(i)%" %5.2f m"s<<Log::endl;
#endif
      }

    x0ReceiveTerminal.resize(gnss->transmitters.size());
    for(auto para : paraReceiveTerminal)
      if(para && para->index && selectedReceiveTerminalZeroMean.at(para->trans->idTrans()))
      {
        x0ReceiveTerminal.at(para->trans->idTrans()) = Vector(para->trans->signalBiasIslRx.biases);
        countZeroMean++;
#if DEBUG >0
        for(UInt i=0; i<para->trans->signalBiasIslRx.biases.size(); i++)
          logInfo<<"initParameter() store initial receive  ISL terminal bias parameter "
                 <<para->trans->name()<<para->trans->signalBiasIslRx.biases.at(i)%" %5.2f m"s<<Log::endl;
#endif
      }
    applyConstraint = isEnabled(normalEquationInfo, nameConstraint) && sigmaZeroMean
                      && (countparaTransmitTerminal+countparaReceiveTerminal)
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
      for(auto para : paraTransmitTerminal)
        if(para && para->index)
        {
          copy(Vector(para->trans->signalBiasIslTx.biases), x0.row(normalEquationInfo.index(para->index), para->trans->signalBiasIslRx.biases.size()));
#if DEBUG > 0
          for(UInt idTerm=0; idTerm<para->trans->signalBiasIslRx.biases.size(); idTerm++)
            logInfo<<"aprioriParameter() transmit ISL terminal bias  "<<para->trans->name()<<idTerm%"TX%2i: "s<<para->trans->signalBiasIslTx.biases.at(idTerm)%" %6.2f"s <<Log::endl;
#endif
        }
      for(auto para : paraReceiveTerminal)
        if(para && para->index)
        {
          copy(Vector(para->trans->signalBiasIslRx.biases), x0.row(normalEquationInfo.index(para->index), para->trans->signalBiasIslRx.biases.size()));
#if DEBUG > 0
          for(UInt idTerm=0; idTerm<para->trans->signalBiasIslRx.biases.size(); idTerm++)
            logInfo<<"aprioriParameter() receive  ISL terminal bias  "<<para->trans->name()<<idTerm%"RX%2i: "s<< para->trans->signalBiasIslRx.biases.at(idTerm)%" %6.2f"s <<Log::endl;
#endif
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

#if DEBUG > 0
    logInfo <<"GnssParametrizationIslBiases::designMatrixIsl() "
            <<" sat rx "<<eqn.receiver->name()<<eqn.terminalRecv%" %2i)"s
            <<" idx "<<eqn.receiver->signalBiasIslRx.index(eqn.terminalRecv)%"%2i"s
            <<" sat tx "<<eqn.transmitter->name()<<eqn.terminalSend%" (%2i)"s
            <<" idx "<<eqn.transmitter->signalBiasIslTx.index(eqn.terminalSend)%"%2i"s
            <<Log::endl;
#endif

    // TODO: check the indexing here!

    auto paraTransmitTerminal = this->paraTransmitTerminal.at(eqn.transmitter->idTrans());
    if(paraTransmitTerminal && paraTransmitTerminal->index)
    {
      UInt idx = eqn.transmitter->signalBiasIslTx.index(eqn.terminalSend);
      axpy(1., eqn.A.column(GnssObservationEquationIsl::idxRange), A.column(paraTransmitTerminal->index));
    }

    auto paraReceiveTerminal = this->paraReceiveTerminal.at(eqn.receiver->idTrans());
    if(paraReceiveTerminal && paraReceiveTerminal->index)
    {
      UInt idx = eqn.receiver->signalBiasIslRx.index(eqn.terminalRecv);
      axpy(1., eqn.A.column(GnssObservationEquationIsl::idxRange), A.column(paraReceiveTerminal->index));
    }

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

    // zero-mean constraint of ISL terminal biases
    // -------------------------------------------

    UInt countTrans = 0, countRecv = 0;
    for(UInt idTrans=0; idTrans<paraTransmitTerminal.size(); idTrans++)
      if(selectedTransmitTerminalZeroMean.at(idTrans) && paraTransmitTerminal.at(idTrans) && paraTransmitTerminal.at(idTrans)->index)
        countTrans++;
    for(UInt idTrans=0; idTrans<paraReceiveTerminal.size(); idTrans++)
      if(selectedReceiveTerminalZeroMean.at(idTrans) && paraReceiveTerminal.at(idTrans) && paraReceiveTerminal.at(idTrans)->index)
        countRecv++;
    UInt count = countTrans + countRecv;
    if(!count)
      return;

    // collect apriori bias values
    // ---------------------------

    if(Parallel::isMaster(normalEquationInfo.comm))
    {
      logWarning<<"ISL bias constraints temporarily disabled!"<<Log::endl;
      /*
      Vector           l(1);
      GnssDesignMatrix A(normalEquationInfo, 1);
      const std::vector<UInt> terminals = { 0 }; // TODO: use variable terminal ID!
      for(auto para : paraTransmitTerminal)
        if(para && para->index && selectedTransmitTerminalZeroMean.at(para->trans->idTrans()))
        {
          l(0) += (x0TransmitTerminal.at(para->trans->idTrans()) - para->trans->signalBiasesIslTx(terminals).at(0))/count/sigmaZeroMean; // remove apriori value -> regularization towards 0
          A.column(para->index)(0,0) = 1./count/sigmaZeroMean;
#if DEBUG > 0
          logInfo<<"constraint() transmit ISL terminal bias "<<para->trans->name()
                 <<x0TransmitTerminal.at(para->trans->idTrans())%" a-priori %6.2f m"s
                 <<para->trans->signalBiasesIslTx(terminals).at(0)%" estimate %6.2f m"s
                 <<Log::endl;
#endif
        }
      for(auto para : paraReceiveTerminal)
        if(para && para->index && selectedReceiveTerminalZeroMean.at(para->trans->idTrans()))
        {
          l(0) += (x0ReceiveTerminal.at(para->trans->idTrans()) - para->trans->signalBiasesIslRx(terminals).at(0))/count/sigmaZeroMean; // remove apriori value -> regularization towards 0
          A.column(para->index)(0,0) = 1./count/sigmaZeroMean;
#if DEBUG > 0
          logInfo<<"constraint() receive  ISL terminal bias "<<para->trans->name()
                 <<x0ReceiveTerminal.at(para->trans->idTrans())%" a-priori %6.2f m"s
                 <<para->trans->signalBiasesIslRx(terminals).at(0)%" estimate %6.2f m"s
                 <<Log::endl;
#endif
        }
      GnssDesignMatrix::accumulateNormals(A, l, normals, n, lPl, obsCount);
      */
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
    for(auto para : paraTransmitTerminal)
      if(para && para->index)
      {
        const Vector dBias = x.row(normalEquationInfo.index(para->index), para->trans->signalBiasIslTx.biases.size());
        for(UInt idType=0; idType<dBias.size(); idType++)
          para->trans->signalBiasIslTx.biases.at(idType) += dBias(idType);
        for(UInt idType=0; idType<dBias.size(); idType++)
          if(infoTrans.update(1e3*dBias(idType)))
            infoTrans.info = "transmit ISL terminal bias ("+para->trans->name()+para->trans->signalBiasIslTx.terminals.at(idType)%"TX %i"s+")";
      }
    infoTrans.synchronizeAndPrint(normalEquationInfo.comm, 1e-3, maxChange);

    Gnss::InfoParameterChange infoRecv("mm");
    for(auto para : paraReceiveTerminal)
      if(para && para->index)
      {
        const Vector dBias = x.row(normalEquationInfo.index(para->index), para->trans->signalBiasIslRx.biases.size());
        for(UInt idType=0; idType<dBias.size(); idType++)
          para->trans->signalBiasIslRx.biases.at(idType) += dBias(idType);
        for(UInt idType=0; idType<dBias.size(); idType++)
          if(infoRecv.update(1e3*dBias(idType)))
            infoTrans.info = "receive ISL terminal bias ("+para->trans->name()+para->trans->signalBiasIslRx.terminals.at(idType)%"RX %i"s+")";
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
