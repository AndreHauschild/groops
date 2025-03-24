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

#define DEBUG 0

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
// TODO: check if the parameters can also be created in initParameter().

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
#if DEBUG > 0
        logInfo<<"init() transmit ISL terminal bias parameter "<<para->trans->name()<<Log::endl;
#endif
      }

    auto selectedReceiveTerminal = gnss->selectTransmitters(selectReceiveTerminal);
    paraReceiveTerminal.resize(gnss->transmitters.size(), nullptr);
    for(UInt idRecv=0; idRecv<gnss->transmitters.size(); idRecv++)
      if(selectedReceiveTerminal.at(idRecv) && gnss->transmitters.at(idRecv)->useable())
      {
        auto para = new Parameter();
        paraReceiveTerminal.at(idRecv) = para;
        para->trans = gnss->transmitters.at(idRecv);
#if DEBUG > 0
        logInfo<<"init() receive ISL terminal bias parameter "<<para->trans->name()<<Log::endl;
#endif
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

    // transmitter ISL terminal parameters
    // -----------------------------------
    UInt countparaTransmitTerminal = 0;
    for(auto para : paraTransmitTerminal)
      if(para && para->trans->useable() && para->trans->signalBiasIslTx.biases.size())
      {
#if DEBUG >0
        logInfo<<"initParameter() transmit ISL terminal bias parameter "<<para->trans->name()<<Log::endl;
#endif
        // determine parameter names
        std::vector<ParameterName> parameterNames;
        parameterNames.push_back(ParameterName(para->trans->name(), "transmitIslTerminalBias"));
        para->index = normalEquationInfo.parameterNamesTransmitter(para->trans->idTrans(), parameterNames);
        countparaTransmitTerminal += parameterNames.size();
      }
    if(countparaTransmitTerminal)
      logInfo<<countparaTransmitTerminal%"%9i ISL transmit terminal bias parameters"s<<Log::endl;

    // receiving ISL terminal parameters
    // ---------------------------------
    UInt countparaReceiveTerminal = 0;
    for(auto para : paraReceiveTerminal)
      if(para && para->trans->useable() && para->trans->signalBiasIslRx.biases.size())
      {
#if DEBUG >0
        logInfo<<"initParameter() receive ISL terminal bias parameter "<<para->trans->name()<<Log::endl;
#endif
        // determine parameter names
        std::vector<ParameterName> parameterNames;
        parameterNames.push_back(ParameterName(para->trans->name(), "receiveIslTerminalBias"));
        para->index = normalEquationInfo.parameterNamesTransmitter(para->trans->idTrans(), parameterNames);
        countparaReceiveTerminal += parameterNames.size();
      }
    if(countparaReceiveTerminal)
      logInfo<<countparaReceiveTerminal%"%9i ISL receive terminal bias parameters"s<<Log::endl;

    // Zero-mean constraint

    selectedTransmitTerminalZeroMean = gnss->selectTransmitters(selectTransmitTerminalZeroMean);
    selectedReceiveTerminalZeroMean = gnss->selectTransmitters(selectReceiveTerminalZeroMean);

    // Initialize a-priori values

    // NOTE: the a-priori values for the zero-mean constrained are stored in
    //       this function and used in the constraints() function. It cannot be
    //       done in the init() function, since the a-priori values may not yet
    //       initialized at that point depending on the parametrization order of
    //       gnssParametrizationislBiases and gnssParametrizationSignalBiasesIsl

    UInt countZeroMean = 0;
    x0TransmitTerminal.resize(gnss->transmitters.size());
    for(auto para : paraTransmitTerminal)
      if(para && para->index && selectedTransmitTerminalZeroMean.at(para->trans->idTrans()))
      {
        x0TransmitTerminal.at(para->trans->idTrans()) = para->trans->signalBiasesIslTx();
        countZeroMean++;
#if DEBUG >0
        logInfo<<"initParameter() store initial send ISL terminal bias parameter "
               <<para->trans->name()<<para->trans->signalBiasesIslTx()%" %5.2fm"s<<Log::endl;
#endif
      }

    x0ReceiveTerminal.resize(gnss->transmitters.size());
    for(auto para : paraReceiveTerminal)
      if(para && para->index && selectedReceiveTerminalZeroMean.at(para->trans->idTrans()))
      {
        x0ReceiveTerminal.at(para->trans->idTrans()) = para->trans->signalBiasesIslRx();
        countZeroMean++;
#if DEBUG >0
        logInfo<<"initParameter() store initial recv ISL terminal bias parameter "
               <<para->trans->name()<<para->trans->signalBiasesIslRx()%" %5.2fm"s<<Log::endl;
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
          x0(normalEquationInfo.index(para->index),0) = para->trans->signalBiasesIslTx();
#if DEBUG > 0
          logInfo<<"aprioriParameter() transmit ISL terminal bias  "<<para->trans->name() << para->trans->signalBiasesIslTx()%" %6.2f"s <<Log::endl;
#endif
        }
      for(auto para : paraReceiveTerminal)
        if(para && para->index)
        {
          x0(normalEquationInfo.index(para->index),0) = para->trans->signalBiasesIslRx();
#if DEBUG > 0
          logInfo<<"aprioriParameter() receiver ISL terminal bias  "<<para->trans->name() << para->trans->signalBiasesIslRx()%" %6.2f"s <<Log::endl;
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
#if DEBUG >0
    logInfo << "GnssParametrizationIslBiases::designMatrixIsl() "
            << "idRecv "  << eqn.receiver->idTrans()
            <<" idTrans " << eqn.transmitter->idTrans()
            <<Log::endl;
#endif
    auto paraTransmitTerminal = this->paraTransmitTerminal.at(eqn.transmitter->idTrans());
    if(paraTransmitTerminal && paraTransmitTerminal->index)
      axpy(1., eqn.A.column(GnssObservationEquationIsl::idxRange), A.column(paraTransmitTerminal->index));

    auto paraReceiveTerminal = this->paraReceiveTerminal.at(eqn.receiver->idTrans());
    if(paraReceiveTerminal && paraReceiveTerminal->index)
      axpy(1., eqn.A.column(GnssObservationEquationIsl::idxRange), A.column(paraReceiveTerminal->index));
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
      GnssDesignMatrix A(normalEquationInfo, Vector(1));
      for(auto para : paraTransmitTerminal)
        if(para && para->index && selectedTransmitTerminalZeroMean.at(para->trans->idTrans()))
        {
          A.l(0) += (x0TransmitTerminal.at(para->trans->idTrans()) - para->trans->signalBiasesIslTx())/count/sigmaZeroMean; // remove apriori value -> regularization towards 0
          A.column(para->index)(0,0) = 1./count/sigmaZeroMean;
#if DEBUG > 0
          logInfo<<"constraint() transmit ISL terminal bias "<<para->trans->name()
                 <<x0TransmitTerminal.at(para->trans->idTrans())%" a-priori %6.2f m"s
                 <<para->trans->signalBiasesIslTx()%" estimate %6.2f m"s
                 <<Log::endl;
#endif
        }
      for(auto para : paraReceiveTerminal)
        if(para && para->index && selectedReceiveTerminalZeroMean.at(para->trans->idTrans()))
        {
          A.l(0) += (x0ReceiveTerminal.at(para->trans->idTrans()) - para->trans->signalBiasesIslRx())/count/sigmaZeroMean; // remove apriori value -> regularization towards 0
          A.column(para->index)(0,0) = 1./count/sigmaZeroMean;
#if DEBUG > 0
          logInfo<<"constraint() receive ISL terminal bias "<<para->trans->name()
                 <<x0ReceiveTerminal.at(para->trans->idTrans())%" a-priori %6.2f m"s
                 <<para->trans->signalBiasesIslRx()%" estimate %6.2f m"s
                 <<Log::endl;
#endif
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
    for(auto para : paraTransmitTerminal)
      if(para && para->index)
      {
        const Vector dBias = x.row(normalEquationInfo.index(para->index), 1);
        for(UInt idType=0; idType<dBias.size(); idType++)
          para->trans->signalBiasIslTx.biases.at(idType) += dBias(idType);
        for(UInt idType=0; idType<dBias.size(); idType++)
          if(infoTrans.update(1e3*dBias(idType)))
            infoTrans.info = "transmit ISL terminal bias ("+para->trans->signalBiasIslTx.types.at(idType).str()+")";
      }
    infoTrans.synchronizeAndPrint(normalEquationInfo.comm, 1e-3, maxChange);

    Gnss::InfoParameterChange infoRecv("mm");
    for(auto para : paraReceiveTerminal)
      if(para && para->index)
      {
        const Vector dBias = x.row(normalEquationInfo.index(para->index), 1);
        for(UInt idType=0; idType<dBias.size(); idType++)
          para->trans->signalBiasIslRx.biases.at(idType) += dBias(idType);
        for(UInt idType=0; idType<dBias.size(); idType++)
          if(infoRecv.update(1e3*dBias(idType)))
            infoTrans.info = "receive ISL terminal bias ("+para->trans->signalBiasIslRx.types.at(idType).str()+")";
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
