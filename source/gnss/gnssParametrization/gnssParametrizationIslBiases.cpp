/***********************************************/
/**
* @file gnssParametrizationIslBiases.cpp
*
* @brief Inter Satellite Link biases.
* @see GnssParametrization
*
* @author Torsten Mayer-Guerr
* @date 2021-01-23
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
    readConfig(config, "name",                            name,                Config::OPTIONAL, "parameter.islBiases", "used for parameter selection");
    readConfig(config, "selectTransmitters",              selectTransmitters,  Config::DEFAULT,  R"(["all"])", "");
    readConfig(config, "selectReceivers",                 selectReceivers,     Config::DEFAULT,  R"(["all"])", "");
    readConfig(config, "nameConstraint",                  nameConstraint,      Config::OPTIONAL, "constraint.islBiases", "used for parameter selection");
    readConfig(config, "sigmaZeroMeanConstraint",         sigmaZeroMean,       Config::DEFAULT,  "0.0001", "(0 = unconstrained) sigma [m] for null space constraint");
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
  for(auto para : paraTrans)
    delete para;
  for(auto para : paraRecv)
    delete para;
}

/***********************************************/

void GnssParametrizationIslBiases::init(Gnss *gnss, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    this->gnss = gnss;

    auto selectedTransmitters = gnss->selectTransmitters(selectTransmitters);
    paraTrans.resize(gnss->transmitters.size(), nullptr);
    for(UInt idTrans=0; idTrans<gnss->transmitters.size(); idTrans++)
      if(selectedTransmitters.at(idTrans) && gnss->transmitters.at(idTrans)->useable())
      {
        auto para = new Parameter();
        paraTrans.at(idTrans) = para;
        para->trans = gnss->transmitters.at(idTrans);
        logInfo<<"Init transmitter ISL bias parameter "<<para->trans->name()<<Log::endl;
      }

    auto selectedReceivers = gnss->selectTransmitters(selectReceivers);
    paraRecv.resize(gnss->transmitters.size(), nullptr);
    for(UInt idRecv=0; idRecv<gnss->transmitters.size(); idRecv++)
      if(selectedReceivers.at(idRecv) && gnss->transmitters.at(idRecv)->useable())
      {
        auto para = new Parameter();
        paraRecv.at(idRecv) = para;
        para->trans = gnss->transmitters.at(idRecv);
        logInfo<<"Init receiver ISL bias parameter "<<para->trans->name()<<Log::endl;
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
    for(auto para : paraTrans)
      if(para)
        para->index = GnssParameterIndex();
    for(auto para : paraRecv)
      if(para)
        para->index = GnssParameterIndex();

    applyConstraint = FALSE;
    if(!isEnabled(normalEquationInfo, name))
      return;

    // transmitter parameters
    // ----------------------
    UInt countParaTrans = 0;
    if(!normalEquationInfo.isEachReceiverSeparately) // TODO: what does this mean?
      for(auto para : paraTrans)
        if(para && para->trans->useable())
        {
          logInfo<<"initParameter() transmitter ISL bias parameter "<<para->trans->name()<<Log::endl;
          // determine parameter names
          std::vector<ParameterName> parameterNames;
          parameterNames.push_back(ParameterName(para->trans->name(), "IslTxBias"));
          para->index = normalEquationInfo.parameterNamesTransmitter(para->trans->idTrans(), parameterNames);
          countParaTrans += parameterNames.size();
        }
    if(countParaTrans)
      logInfo<<countParaTrans%"%9i transmitter ISL bias parameters"s<<Log::endl;

    // receiver parameters
    // -------------------
    UInt countParaRecv = 0;
    for(auto para : paraRecv)
      if(para && para->trans->useable())
      {
        logInfo<<"initParameter() receiver ISL bias parameter "<<para->trans->name()<<Log::endl;
        // determine parameter names
        std::vector<ParameterName> parameterNames;
        parameterNames.push_back(ParameterName(para->trans->name(), "IslRxBias"));
        para->index = normalEquationInfo.parameterNamesTransmitter(para->trans->idTrans(), parameterNames);
        countParaRecv += parameterNames.size();
      }
    if(countParaRecv)
      logInfo<<countParaRecv%"%9i receiver ISL bias parameters"s<<Log::endl;

    applyConstraint = isEnabled(normalEquationInfo, nameConstraint) && sigmaZeroMean
                    && (countParaTrans+countParaRecv) &&  !normalEquationInfo.isEachReceiverSeparately;

    // TODO: constraints?

  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
// TODO: make sure correct bias value is selected!
void GnssParametrizationIslBiases::aprioriParameter(const GnssNormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const
{
  try
  {
    if(Parallel::isMaster(normalEquationInfo.comm))
      for(auto para : paraTrans)
        if(para && para->index)
        {
          logInfo<<"aprioriParameter() transmitter ISL bias  "<<para->trans->name() << " " << para->trans->signalBiasIslTx.biases.at(0) <<Log::endl;
          x0(normalEquationInfo.index(para->index),0) = para->trans->signalBiasIslTx.biases.at(0);
        }
    for(auto para : paraRecv)
      if(para && para->index)
      {
        logInfo<<"aprioriParameter() receiver ISL bias  "<<para->trans->name() << " " << para->trans->signalBiasIslRx.biases.at(0) <<Log::endl;
        x0(normalEquationInfo.index(para->index),0) = para->trans->signalBiasIslRx.biases.at(0);
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
    auto paraTrans = this->paraTrans.at(eqn.transmitter->idTrans());
    if(paraTrans && paraTrans->index)
      copy(eqn.A.column(GnssObservationEquationIsl::idxRange,1), A.column(paraTrans->index));

    auto paraRecv = this->paraRecv.at(eqn.receiver->idTrans());
    if(paraRecv && paraRecv->index)
    {
      copy(eqn.A.column(GnssObservationEquationIsl::idxRange,1), A.column(paraRecv->index));
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
    if(Parallel::isMaster(normalEquationInfo.comm) && applyConstraint)
    {
      logStatus<<"apply "<<zeroMeanDesign.rows()<<" zero mean equations for ISL bias parameters"<<Log::endl;
      if(zeroMeanDesign.size())
      {
        // NOTE: axpy(c,A,B) -> B += c* A
        GnssDesignMatrix A(normalEquationInfo, Vector(zeroMeanDesign.rows()));
        for(auto para : paraTrans)
          if(para && para->index && (idxBiasTrans.at(para->trans->idTrans()) != NULLINDEX))
            axpy(1./sigmaZeroMean, zeroMeanDesign.column(idxBiasTrans.at(para->trans->idTrans()), 1), A.column(para->index));
        for(auto para : paraRecv)
          if(para && para->index && (idxBiasRecv.at(para->trans->idTrans()) != NULLINDEX))
            axpy(1./sigmaZeroMean, zeroMeanDesign.column(idxBiasRecv.at(para->trans->idTrans()), 1), A.column(para->index));
        A.accumulateNormals(normals, n, lPl, obsCount);
      }
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
    for(auto para : paraTrans)
      if(para && para->index)
      {
        const Vector dBias = x.row(normalEquationInfo.index(para->index), 1);
        for(UInt idType=0; idType<dBias.size(); idType++)
          para->trans->signalBiasIslTx.biases.at(idType) += dBias(idType);
        for(UInt idType=0; idType<dBias.size(); idType++)
          if(infoTrans.update(1e3*dBias(idType)))
            infoTrans.info = "ISL bias transmitter ("+para->trans->signalBiasIslTx.types.at(idType).str()+")";
      }
    infoTrans.synchronizeAndPrint(normalEquationInfo.comm, 1e-3, maxChange);

    Gnss::InfoParameterChange infoRecv("mm");
    for(auto para : paraRecv)
      if(para && para->index)
      {
        const Vector dBias = x.row(normalEquationInfo.index(para->index), 1);
        for(UInt idType=0; idType<dBias.size(); idType++)
          para->trans->signalBiasIslRx.biases.at(idType) += dBias(idType);
        for(UInt idType=0; idType<dBias.size(); idType++)
          if(infoRecv.update(1e3*dBias(idType)))
            infoTrans.info = "ISL bias receiver ("+para->trans->signalBiasIslRx.types.at(idType).str()+")";
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
