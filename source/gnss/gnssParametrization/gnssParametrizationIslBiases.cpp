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
    readConfig(config, "name",                           name,                           Config::OPTIONAL, "parameter.islBiases", "used for parameter selection");
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
  for(UInt idTrans=0; idTrans<gnss->transmitters.size(); idTrans++)
  {
    for(auto para : paraTransmitTerminal.at(idTrans))
      delete para;
    for(auto para : paraReceiveTerminal.at(idTrans))
      delete para;
  }
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
    paraTransmitTerminal.resize(gnss->transmitters.size());
    for(UInt idTrans=0; idTrans<gnss->transmitters.size(); idTrans++)
      if(selectedTransmitTerminal.at(idTrans) && gnss->transmitters.at(idTrans)->useable())
      {
        UInt nTerm=gnss->transmitters.at(idTrans)->islBiasSend.terminals.size();
        for(UInt idTerm=0; idTerm<nTerm; idTerm++)
        {
          auto para = new Parameter();
          paraTransmitTerminal.at(idTrans).push_back(para);
          para->trans = gnss->transmitters.at(idTrans);
        }
      }

    auto selectedReceiveTerminal = gnss->selectTransmitters(selectReceiveTerminal);
    paraReceiveTerminal.resize(gnss->transmitters.size());
    for(UInt idRecv=0; idRecv<gnss->transmitters.size(); idRecv++)
      if(selectedReceiveTerminal.at(idRecv) && gnss->transmitters.at(idRecv)->useable())
      {
        UInt nTerm=gnss->transmitters.at(idRecv)->islBiasRecv.terminals.size();
        for(UInt idTerm=0; idTerm<nTerm; idTerm++)
        {
          auto para = new Parameter();
          paraReceiveTerminal.at(idRecv).push_back(para);
          para->trans = gnss->transmitters.at(idRecv);
        }
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
    for(UInt idTrans=0; idTrans<gnss->transmitters.size(); idTrans++)
    {
      for(auto para : paraTransmitTerminal.at(idTrans))
        if(para)
          para->index = GnssParameterIndex();
      for(auto para : paraReceiveTerminal.at(idTrans))
        if(para)
          para->index = GnssParameterIndex();
    }

    applyConstraint = FALSE;
    if(!isEnabled(normalEquationInfo, name) || normalEquationInfo.isEachReceiverSeparately)
      return;

    // transmit ISL terminal parameters
    // --------------------------------
    UInt countparaTransmitTerminal = 0;
    for(UInt idTrans=0; idTrans<gnss->transmitters.size(); idTrans++)
      for(UInt idTerm=0; idTerm<gnss->transmitters.at(idTrans)->islBiasSend.terminals.size(); idTerm++)
        if(paraTransmitTerminal.at(idTrans).size() && idTerm<paraTransmitTerminal.at(idTrans).size())
        {
          auto para = paraTransmitTerminal.at(idTrans).at(idTerm);
          std::vector<ParameterName> parameterNames;
          parameterNames.push_back(ParameterName(para->trans->name(), "islBiasSend"+idTerm%"(%i)"s));
          para->index = normalEquationInfo.parameterNamesTransmitter(para->trans->idTrans(), parameterNames);
          countparaTransmitTerminal += parameterNames.size();
#if DEBUG > 0
          logInfo<<"initParameter() transmit ISL terminal bias parameter "<<para->trans->name()<<":"<<"islBiasSend"+idTerm%"%i"s<<Log::endl;
#endif
        }
    if(countparaTransmitTerminal)
      logInfo<<countparaTransmitTerminal%"%9i ISL transmit terminal bias parameters"s<<Log::endl;

    // receive ISL terminal parameters
    // -------------------------------
    UInt countparaReceiveTerminal = 0;
    for(UInt idTrans=0; idTrans<gnss->transmitters.size(); idTrans++)
      for(UInt idTerm=0; idTerm<gnss->transmitters.at(idTrans)->islBiasRecv.terminals.size(); idTerm++)
        if(paraReceiveTerminal.at(idTrans).size() && idTerm<paraReceiveTerminal.at(idTrans).size())
        {
          auto para = paraReceiveTerminal.at(idTrans).at(idTerm);
          std::vector<ParameterName> parameterNames;
          parameterNames.push_back(ParameterName(para->trans->name(), "islBiasRecv"+idTerm%"(%i)"s));
          para->index = normalEquationInfo.parameterNamesTransmitter(para->trans->idTrans(), parameterNames);
          countparaReceiveTerminal += parameterNames.size();
#if DEBUG > 0
          logInfo<<"initParameter() receive ISL terminal bias parameter "<<para->trans->name()<<":"<<"islBiasRecv"+idTerm%"%i"s<<Log::endl;
#endif
        }
    if(countparaReceiveTerminal)
      logInfo<<countparaReceiveTerminal%"%9i ISL receive  terminal bias parameters"s<<Log::endl;

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
    for(UInt idTrans=0; idTrans<gnss->transmitters.size(); idTrans++)
    {
      x0TransmitTerminal.at(idTrans).resize(gnss->transmitters.at(idTrans)->islBiasSend.terminals.size());
      for(UInt idTerm=0; idTerm<x0TransmitTerminal.at(idTrans).size(); idTerm++)
        if(paraTransmitTerminal.at(idTrans).size() && idTerm<paraTransmitTerminal.at(idTrans).size() && \
           selectedTransmitTerminalZeroMean.at(idTrans))
        {
          auto para = paraTransmitTerminal.at(idTrans).at(idTerm);
          x0TransmitTerminal.at(idTrans).at(idTerm) = para->trans->islBiasSend.biases.at(idTerm);
          countZeroMean++;
#if DEBUG >0
          for(UInt i=0; i<para->trans->islBiasSend.biases.size(); i++)
            logInfo<<"initParameter() store initial transmit ISL terminal bias parameter "
                   <<para->trans->name()<<para->trans->islBiasSend.biases.at(i)%" %5.2f m"s<<Log::endl;
#endif
        }
      }

    x0ReceiveTerminal.resize(gnss->transmitters.size());
    for(UInt idTrans=0; idTrans<gnss->transmitters.size(); idTrans++)
    {
      x0ReceiveTerminal.at(idTrans).resize(gnss->transmitters.at(idTrans)->islBiasRecv.terminals.size());
      for(UInt idTerm=0; idTerm<x0ReceiveTerminal.at(idTrans).size(); idTerm++)
        if(paraReceiveTerminal.at(idTrans).size() && idTerm<paraReceiveTerminal.at(idTrans).size() && \
           selectedReceiveTerminalZeroMean.at(idTrans))
        {
          auto para = paraReceiveTerminal.at(idTrans).at(idTerm);
          x0ReceiveTerminal.at(idTrans).at(idTerm) = para->trans->islBiasRecv.biases.at(idTerm);
          countZeroMean++;
#if DEBUG >0
          for(UInt i=0; i<para->trans->islBiasRecv.biases.size(); i++)
            logInfo<<"initParameter() store initial receive  ISL terminal bias parameter "
                   <<para->trans->name()<<para->trans->islBiasRecv.biases.at(i)%" %5.2f m"s<<Log::endl;
#endif
        }
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
      for(UInt idTrans=0; idTrans<gnss->transmitters.size(); idTrans++)
        for(auto para : paraTransmitTerminal.at(idTrans))
        {
          if(para && para->index)
          {
            copy(Vector(para->trans->islBiasSend.biases), x0.row(normalEquationInfo.index(para->index), para->trans->islBiasRecv.biases.size()));
#if DEBUG > 0
            for(UInt idTerm=0; idTerm<para->trans->islBiasRecv.biases.size(); idTerm++)
              logInfo<<"aprioriParameter() transmit ISL terminal bias  "<<para->trans->name()<<idTerm%"TX%2i: "s<<para->trans->islBiasSend.biases.at(idTerm)%" %6.2f"s <<Log::endl;
#endif
          }
          for(auto para : paraReceiveTerminal.at(idTrans))
            if(para && para->index)
            {
              copy(Vector(para->trans->islBiasRecv.biases), x0.row(normalEquationInfo.index(para->index), para->trans->islBiasRecv.biases.size()));
#if DEBUG > 0
              for(UInt idTerm=0; idTerm<para->trans->islBiasRecv.biases.size(); idTerm++)
                logInfo<<"aprioriParameter() receive  ISL terminal bias  "<<para->trans->name()<<idTerm%"RX%2i: "s<< para->trans->islBiasRecv.biases.at(idTerm)%" %6.2f"s <<Log::endl;
#endif
            }
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
#if DEBUG > 1
    logInfo << "GnssParametrizationIslBiases::designMatrixIsl() "
            << eqn.receiver->name() << eqn.receiver->idTrans()%" (%2i)"s
            << " <- "
            << eqn.transmitter->name() << eqn.transmitter->idTrans()%" (%2i)"s
            << Log::endl;
#endif
    UInt idTerm;

    // transmitter terminal bias
    idTerm = eqn.transmitter->islBiasSend.index(eqn.terminalSend);
    auto paraTransmitTerminal = this->paraTransmitTerminal.at(eqn.transmitter->idTrans()).at(idTerm);
    if(paraTransmitTerminal && paraTransmitTerminal->index)
      copy(eqn.A.column(GnssObservationEquationIsl::idxRange,1), A.column(paraTransmitTerminal->index));

    // receiver terminal bias
    idTerm = eqn.transmitter->islBiasRecv.index(eqn.terminalRecv);
    auto paraReceiveTerminal = this->paraReceiveTerminal.at(eqn.receiver->idTrans()).at(idTerm);
    if(paraReceiveTerminal && paraReceiveTerminal->index)
      copy(eqn.A.column(GnssObservationEquationIsl::idxRange,1), A.column(paraReceiveTerminal->index));
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

    /*
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
    }
     */
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
    for(UInt idTrans=0; idTrans<paraTransmitTerminal.size(); idTrans++)
      for(auto para : paraTransmitTerminal.at(idTrans))
        if(para && para->index)
        {
          const Vector dBias = x.row(normalEquationInfo.index(para->index), para->trans->islBiasSend.biases.size());
          for(UInt idType=0; idType<dBias.size(); idType++)
            para->trans->islBiasSend.biases.at(idType) += dBias(idType);
          for(UInt idType=0; idType<dBias.size(); idType++)
            if(infoTrans.update(1e3*dBias(idType)))
              infoTrans.info = "transmit ISL terminal bias ("+para->trans->name()+para->trans->islBiasSend.terminals.at(idType)%"Send(%i)"s+")";
        }
    infoTrans.synchronizeAndPrint(normalEquationInfo.comm, 1e-3, maxChange);

    Gnss::InfoParameterChange infoRecv("mm");
    for(UInt idTrans=0; idTrans<paraTransmitTerminal.size(); idTrans++)
      for(auto para : paraReceiveTerminal.at(idTrans))
        if(para && para->index)
        {
          const Vector dBias = x.row(normalEquationInfo.index(para->index), para->trans->islBiasRecv.biases.size());
          for(UInt idType=0; idType<dBias.size(); idType++)
            para->trans->islBiasRecv.biases.at(idType) += dBias(idType);
          for(UInt idType=0; idType<dBias.size(); idType++)
            if(infoRecv.update(1e3*dBias(idType)))
              infoTrans.info = "receive  ISL terminal bias ("+para->trans->name()+para->trans->islBiasRecv.terminals.at(idType)%"Recv(%i)"s+")";
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
