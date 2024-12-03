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
  for(auto para : parameter)
    delete para;
}

/***********************************************/

void GnssParametrizationIslBiases::init(Gnss *gnss, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    this->gnss = gnss;

    auto selectedTransmitters = gnss->selectTransmitters(selectTransmitters);
    parameter.resize(gnss->transmitters.size(), nullptr);
    for(UInt idTrans=0; idTrans<gnss->transmitters.size(); idTrans++)
      if(selectedTransmitters.at(idTrans) && gnss->transmitters.at(idTrans)->useable())
      {
        auto para = new Parameter();
        parameter.at(idTrans) = para;
        para->trans = gnss->transmitters.at(idTrans);
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
    for(auto para : parameter)
      if(para)
        para->index = GnssParameterIndex();
    applyConstraint = FALSE;
    if(!isEnabled(normalEquationInfo, name))
      return;

    auto nullSpace = [](Matrix &N, Bool remove)
    {
      Vector eigen = eigenValueDecomposition(N, TRUE);
      eigen *= (eigen(eigen.rows()-1) > 1e-4) ? 1./eigen(eigen.rows()-1) : 0.;
      UInt countZeros = 0;
      while((countZeros < eigen.rows()) && (eigen(countZeros) < 1e-8))
        countZeros++;
      return (remove) ? N.column(countZeros, N.columns()-countZeros) : N.column(0, countZeros);
    };

    // transmitter parameters
    // ----------------------
    UInt countParaTrans = 0;
    if(!normalEquationInfo.isEachReceiverSeparately)
      for(auto para : parameter)
        if(para && para->trans->useable())
        {
          std::vector<GnssType> biasTypes;
          for(GnssType type : para->trans->signalBias.types)
            if((type == GnssType::RANGE) && !type.isInList(biasTypes))
              biasTypes.push_back(type);
          if(biasTypes.size() <= 2)
            continue;

          // transformation matrix
          Matrix T(para->trans->signalBias.types.size(), biasTypes.size());
          for(UInt i=0; i<para->trans->signalBias.types.size(); i++)
            for(UInt k=0; k<biasTypes.size(); k++)
              if(para->trans->signalBias.types.at(i) == biasTypes.at(k))
                T(i, k) = 1.;

          // determine nullspace
          Matrix N(1+T.columns(), Matrix::SYMMETRIC);
          for(UInt idRecv=0; idRecv<gnss->receivers.size(); idRecv++)
          {
            std::vector<GnssType> types;
            for(GnssType type : gnss->typesRecvTrans.at(idRecv).at(para->trans->idTrans()))
              if((type == GnssType::RANGE) && !type.isInList(types))
                types.push_back(type);
            if(!types.size())
              continue;

            // Composed signals (e.g. C2DG)
            std::vector<GnssType> typesTrans;
            Matrix Compose;
            gnss->receivers.at(idRecv)->signalComposition(NULLINDEX, types, typesTrans, Compose);

            Matrix A(types.size(), 1+T.columns());
            for(UInt i=0; i<types.size(); i++)
              if(types.at(i) == GnssType::RANGE)
                A(i, 0) = -1.;  // clock
            UInt idx;
            for(UInt i=0; i<typesTrans.size(); i++)
              if(typesTrans.at(i).isInList(para->trans->signalBias.types, idx))
                matMult(1., Compose.column(i), T.row(idx), A.column(1, T.columns()));
            Vector STEC(types.size());
            for(UInt i=0; i<types.size(); i++)
              STEC(i) = types.at(i).ionosphericFactor();
            eliminationParameter(STEC, {A});
            rankKUpdate(1., A, N);
          }

          // eliminate clock
          rankKUpdate(-1./N(0,0), N.slice(0, 1, 1, N.rows()-1), N.slice(1, 1, N.rows()-1, N.rows()-1));
          N = N.slice(1, 1, N.rows()-1, N.rows()-1); // without clock

          para->Bias = T * nullSpace(N, TRUE);
          if(!para->Bias.size())
            continue;

          // determine parameter names
          std::vector<ParameterName> parameterNames;
          for(UInt i=0; i<para->Bias.columns(); i++)
          {
            std::string typeStr;
            for(UInt idType=0; idType<para->Bias.rows(); idType++)
              if(std::fabs(para->Bias(idType, i)) > 1e-4)
                typeStr += ((para->Bias(idType, i) > 0) ? "+" : "") + para->Bias(idType, i)%"%.2f"s + para->trans->signalBias.types.at(idType).str();
            parameterNames.push_back(ParameterName(para->trans->name(), "codeBias"+(i+1)%"%02i("s+typeStr+")"));
          }
          para->index = normalEquationInfo.parameterNamesTransmitter(para->trans->idTrans(), parameterNames);
          countParaTrans += parameterNames.size();
        }
    if(countParaTrans)
      logInfo<<countParaTrans%"%9i transmitter code bias parameters"s<<Log::endl;


    applyConstraint = isEnabled(normalEquationInfo, nameConstraint) && sigmaZeroMean
                    && (countParaTrans) && !normalEquationInfo.isEachReceiverSeparately;

    // TODO: fix this!!

    // calculate constraint equations (zero mean of signal biases)
    // -----------------------------------------------------------
    if(applyConstraint)
    {
      // parameter indices
      UInt parameterCount = 0;
      // indices of transmitter biases
      idxBias.clear();
      idxBias.resize(gnss->transmitters.size(), NULLINDEX);
      for(auto para : parameter)
        if(para && para->index)
        {
          idxBias.at(para->trans->idTrans()) = parameterCount;
          parameterCount += para->Bias.columns();
        }

      if(parameterCount)
      {
        UInt idxClocks = parameterCount;
        // indices of transmitter clocks
        std::vector<UInt> idxClockTrans(gnss->transmitters.size(), NULLINDEX);
        for(auto trans : gnss->transmitters)
          if(trans->useable())
            idxClockTrans.at(trans->idTrans()) = parameterCount++;
        // indices of receiver clocks
        std::vector<UInt> idxClockRecv(gnss->receivers.size(), NULLINDEX);
        // TODO: fix this!!
        //typesRecvTrans.clear();
        //typesRecvTrans.resize(gnss->receivers.size());
        //for(auto recv : gnss->receivers)
        //  if(normalEquationInfo.estimateReceiver.at(recv->idRecv()))
        //  {
        //    const UInt idRecvOld = std::distance(typesRecvTrans.begin(), std::find(typesRecvTrans.begin(), typesRecvTrans.end(), gnss->typesRecvTrans.at(recv->idRecv())));
        //    if(idRecvOld >= typesRecvTrans.size())
        //    {
        //      typesRecvTrans.at(recv->idRecv()) = gnss->typesRecvTrans.at(recv->idRecv());
        //      idxClockRecv.at(recv->idRecv())   = parameterCount++;
        //    }
        //    else
        //      idxClockRecv.at(recv->idRecv()) = idxClockRecv.at(idRecvOld);
        //  }

        // normals of pseudo observations
        Matrix N(parameterCount, Matrix::SYMMETRIC);
        for(auto recv : gnss->receivers)
          if(recv->isMyRank() && normalEquationInfo.estimateReceiver.at(recv->idRecv()))
            for(auto trans : gnss->transmitters)
            {
              // observation types
              std::vector<GnssType> types;
              for(GnssType type : gnss->typesRecvTrans.at(recv->idRecv()).at(trans->idTrans()))
                if((type == GnssType::RANGE) && !type.isInList(types))
                  types.push_back(type);
              if(!types.size())
                continue;

              // transmitted types
              std::vector<GnssType> typesTrans;
              Matrix T;
              recv->signalComposition(NULLINDEX/*idEpoch*/, types, typesTrans, T);

              // design matrix
              Matrix A(types.size(), parameterCount);
              if(idxClockRecv.at(recv->idRecv()) != NULLINDEX)    // clock recv
                for(UInt i=0; i<types.size(); i++)
                  A(i, idxClockRecv.at(recv->idRecv())) = 1.;
              if(idxClockTrans.at(trans->idTrans()) != NULLINDEX) // clock trans
                for(UInt i=0; i<types.size(); i++)
                  A(i, idxClockTrans.at(trans->idTrans())) = -1.;
              if(idxBias.at(trans->idTrans()) != NULLINDEX)  // bias trans
                for(UInt k=0; k<typesTrans.size(); k++)
                  matMult(1., T.column(k), parameter.at(trans->idTrans())->Bias.row(GnssType::index(trans->signalBias.types, typesTrans.at(k))),
                          A.column(idxBias.at(trans->idTrans()), parameter.at(trans->idTrans())->Bias.columns()));
              // eliminate STEC
              Vector STEC(types.size());
              for(UInt i=0; i<types.size(); i++)
                STEC(i) = types.at(i).ionosphericFactor();
              eliminationParameter(STEC, {A});
              rankKUpdate(1., A, N);
            }
        Parallel::reduceSum(N, 0, normalEquationInfo.comm);

        if(Parallel::isMaster(normalEquationInfo.comm))
        {
          // add zero mean of clocks and eliminate clocks
          Matrix N11 = N.slice(idxClocks, idxClocks, parameterCount-idxClocks, parameterCount-idxClocks);
          for(UInt i=0; i<N11.rows(); i++)
            if(N11(i,i) == 0)
              N11(i,i) = 1.;
          rankKUpdate(1., Matrix(1, N11.rows(), 1.), N11); // add zero mean
          cholesky(N11);
          triangularSolve(1., N11.trans(), N.slice(0, idxClocks, idxClocks, N11.rows()).trans());
          rankKUpdate(-1., N.slice(0, idxClocks, idxClocks, N11.rows()).trans(), N.slice(0, 0, idxClocks, idxClocks));
          N = N.slice(0, 0, idxClocks, idxClocks);

          zeroMeanDesign = nullSpace(N, FALSE).trans(); // null space defines the constraint equations
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

void GnssParametrizationIslBiases::aprioriParameter(const GnssNormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const
{
  try
  {
    if(Parallel::isMaster(normalEquationInfo.comm))
      for(auto para : parameter)
        if(para && para->index)
          copy(leastSquares(Matrix(para->Bias), Vector(para->trans->signalBias.biases)), x0.row(normalEquationInfo.index(para->index), para->Bias.columns()));
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
    auto parameter = this->parameter.at(eqn.transmitter->idTrans());
    if(parameter && parameter->index)
    {
      UInt idx;
      MatrixSlice Design(A.column(parameter->index));
      // TODO: fix this!!
      //for(UInt idType=0; idType<eqn.typesTransmitted.size(); idType++)
      //  if((eqn.typesTransmitted.at(idType) == GnssType::RANGE) && eqn.typesTransmitted.at(idType).isInList(eqn.transmitter->signalBias.types, idx))
      //   matMult(1., eqn.A.column(GnssObservationEquationIsl::idxUnit + eqn.types.size() + idType), parameter->Bias.row(idx), Design);
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
      logStatus<<"apply "<<zeroMeanDesign.rows()<<" zero mean equations for code bias parameters"<<Log::endl;
      if(zeroMeanDesign.size())
      {
        GnssDesignMatrix A(normalEquationInfo, Vector(zeroMeanDesign.rows()));
        for(auto para : parameter)
          if(para && para->index && (idxBias.at(para->trans->idTrans()) != NULLINDEX))
            axpy(1./sigmaZeroMean, zeroMeanDesign.column(idxBias.at(para->trans->idTrans()), para->Bias.columns()), A.column(para->index));
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
    for(auto para : parameter)
      if(para && para->index)
      {
        const Vector dBias = para->Bias * x.row(normalEquationInfo.index(para->index), para->Bias.columns());
        for(UInt idType=0; idType<dBias.size(); idType++)
          para->trans->signalBias.biases.at(idType) += dBias(idType);
        for(UInt idType=0; idType<dBias.size(); idType++)
          if(infoTrans.update(1e3*dBias(idType)))
            infoTrans.info = "ISL bias transmitter ("+para->trans->signalBias.types.at(idType).str()+")";
      }
    infoTrans.synchronizeAndPrint(normalEquationInfo.comm, 1e-3, maxChange);

    return maxChange;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
