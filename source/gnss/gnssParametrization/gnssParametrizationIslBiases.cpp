/***********************************************/
/**
* @file gnssParametrizationIslBiases.cpp
*
* @brief Inter-satellite link biases.
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
    readConfig(config, "name",                           name,                           Config::OPTIONAL, "parameter.islBiases", "used for parameter selection");
    readConfig(config, "selectTransmitIslTerminal",      selectTransmitTerminal,         Config::DEFAULT,  R"(["all"])", "");
    readConfig(config, "selectReceiveIslTerminal",       selectReceiveTerminal,          Config::DEFAULT,  R"(["all"])", "");
    readConfig(config, "outputfileTransmitIslTerminal",  fileNameOutTransmitter,         Config::OPTIONAL, "", "variable {prn} available");
    readConfig(config, "outputfileReceiveIslTerminal",   fileNameOutReceiver,            Config::OPTIONAL, "", "variable {prn} available");
    readConfig(config, "inputfileTransmitIslTerminal",   fileNameInTransmitter,          Config::OPTIONAL, "", "variable {prn} available");
    readConfig(config, "inputfileReceiveIslTerminal",    fileNameInReceiver,             Config::OPTIONAL, "", "variable {prn} available");
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

    // Load a-priori values
    //---------------------

    if(!fileNameInTransmitter.empty())
    {
      VariableList fileNameVariableList;
      auto selectedTransmitters = gnss->selectTransmitters(selectTransmitTerminal);
      for(UInt idTrans=0; idTrans<gnss->transmitters.size(); idTrans++)
        if(selectedTransmitters.at(idTrans) && gnss->transmitters.at(idTrans)->useable())
        {
          fileNameVariableList.setVariable("prn", gnss->transmitters.at(idTrans)->name());
          try
          {
            readFileIslBias(fileNameInTransmitter(fileNameVariableList), gnss->transmitters.at(idTrans)->islBiasSend);
          }
          catch(std::exception &/*e*/)
          {
            logWarningOnce<<"Unable to read transmit ISL terminal bias file <"<<fileNameInTransmitter(fileNameVariableList)<<">, disabling transmitter."<<Log::endl;
            gnss->transmitters.at(idTrans)->disable("Unable to read transmit ISL terminal bias file <"+fileNameInTransmitter(fileNameVariableList).str()+">");
          }
        }
    }

    if(!fileNameInReceiver.empty())
    {
      VariableList fileNameVariableList;
      auto selectedTransmitters = gnss->selectTransmitters(selectReceiveTerminal);
      for(UInt idTrans=0; idTrans<gnss->transmitters.size(); idTrans++)
        if(selectedTransmitters.at(idTrans) && gnss->transmitters.at(idTrans)->useable())
        {
          fileNameVariableList.setVariable("prn", gnss->transmitters.at(idTrans)->name());
          try
          {
            readFileIslBias(fileNameInReceiver(fileNameVariableList), gnss->transmitters.at(idTrans)->islBiasRecv);
          }
          catch(std::exception &/*e*/)
          {
            logWarningOnce<<"Unable to read receive ISL terminal bias file <"<<fileNameInReceiver(fileNameVariableList)<<">, disabling transmitter."<<Log::endl;
            gnss->transmitters.at(idTrans)->disable("Unable to read receive ISL terminal bias file <"+fileNameInReceiver(fileNameVariableList).str()+">");
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

void GnssParametrizationIslBiases::initParameter(GnssNormalEquationInfo &normalEquationInfo)
{
  try
  {

    // Setup parameters
    //-----------------

    auto selectedTransmitTerminal = gnss->selectTransmitters(selectTransmitTerminal);
    paraTransmitTerminal.resize(gnss->transmitters.size());
    for(UInt idTrans=0; idTrans<gnss->transmitters.size(); idTrans++)
      if(selectedTransmitTerminal.at(idTrans) && gnss->transmitters.at(idTrans)->useable())
      {
        UInt nTerm=gnss->transmitters.at(idTrans)->islBiasSend.terminals.size();
#if DEBUG > 0
          logInfo<<"init() send ISL terminal bias parameters "<<gnss->transmitters.at(idTrans)->name()<<nTerm%" %i"s<<Log::endl;
#endif
        for(UInt idTerm=0; idTerm<nTerm; idTerm++)
        {
          auto para = new Parameter();
          paraTransmitTerminal.at(idTrans).push_back(para);
          para->trans = gnss->transmitters.at(idTrans);
        /*para->terminal = gnss->transmitters.at(idTrans)->islBiasSend.terminals.at(idTerm);*/
          para->index = GnssParameterIndex();
#if DEBUG > 0
          logInfo<<"init() send ISL terminal bias parameter "<<para->trans->name()<<Log::endl;
#endif
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
        /*para->terminal = gnss->transmitters.at(idRecv)->islBiasRecv.terminals.at(idTerm);*/
          para->index = GnssParameterIndex();
#if DEBUG > 0
          logInfo<<"init() recv ISL terminal bias parameter "<<para->trans->name()<<Log::endl;
#endif
        }
      }

    applyConstraint = FALSE;
    if(!isEnabled(normalEquationInfo, name))
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
          parameterNames.push_back(ParameterName(para->trans->name(), idTerm%"islBiasSend.%i"s));
          para->index = normalEquationInfo.parameterNamesTransmitter(para->trans->idTrans(), parameterNames);
          countparaTransmitTerminal += parameterNames.size();
#if DEBUG > 0
          logInfo<<"initParameter() send ISL terminal bias parameter "<<parameterNames.at(0).str()<<Log::endl;
#endif
        }
    if(countparaTransmitTerminal)
      logInfo<<countparaTransmitTerminal%"%9i ISL send terminal bias parameters"s<<Log::endl;

    // receive ISL terminal parameters
    // -------------------------------
    UInt countparaReceiveTerminal = 0;
    for(UInt idTrans=0; idTrans<gnss->transmitters.size(); idTrans++)
      for(UInt idTerm=0; idTerm<gnss->transmitters.at(idTrans)->islBiasRecv.terminals.size(); idTerm++)
        if(paraReceiveTerminal.at(idTrans).size() && idTerm<paraReceiveTerminal.at(idTrans).size())
        {
          auto para = paraReceiveTerminal.at(idTrans).at(idTerm);
          std::vector<ParameterName> parameterNames;
          parameterNames.push_back(ParameterName(para->trans->name(), idTerm%"islBiasRecv.%i"s));
          para->index = normalEquationInfo.parameterNamesTransmitter(para->trans->idTrans(), parameterNames);
          countparaReceiveTerminal += parameterNames.size();
#if DEBUG > 0
          logInfo<<"initParameter() recv ISL terminal bias parameter "<<parameterNames.at(0).str()<<Log::endl;
#endif
        }
    if(countparaReceiveTerminal)
      logInfo<<countparaReceiveTerminal%"%9i ISL recv terminal bias parameters"s<<Log::endl;

    // Zero-mean constraint

    selectedTransmitTerminalZeroMean = gnss->selectTransmitters(selectTransmitTerminalZeroMean);
    selectedReceiveTerminalZeroMean = gnss->selectTransmitters(selectReceiveTerminalZeroMean);

    // Initialize a-priori values

    // NOTE: the a-priori values for the zero-mean constraint are stored in
    //       this function and used in the constraints() function.

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
            logInfo<<"initParameter() store initial send ISL terminal bias parameter "
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
            logInfo<<"initParameter() store initial recv ISL terminal bias parameter "
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
      for(UInt idTrans=0; idTrans<gnss->transmitters.size(); idTrans++)
      {
        for(UInt idTerm=0; idTerm<gnss->transmitters.at(idTrans)->islBiasSend.terminals.size(); idTerm++)
          if(paraTransmitTerminal.at(idTrans).size() && idTerm<paraTransmitTerminal.at(idTrans).size())
          {
            auto para = paraTransmitTerminal.at(idTrans).at(idTerm);
            if(para && para->index)
            {
              x0(normalEquationInfo.index(para->index),0) = para->trans->islBiasSend.biases.at(idTerm);
#if DEBUG > 0
              logInfo<<"aprioriParameter() send ISL terminal bias "
                     <<normalEquationInfo.parameterNames().at(normalEquationInfo.index(para->index)).str()
                     <<para->trans->islBiasSend.biases.at(idTerm)%" %6.2f"s <<Log::endl;
#endif
              }
          }
        for(UInt idTerm=0; idTerm<gnss->transmitters.at(idTrans)->islBiasRecv.terminals.size(); idTerm++)
          if(paraReceiveTerminal.at(idTrans).size() && idTerm<paraReceiveTerminal.at(idTrans).size())
          {
            auto para = paraReceiveTerminal.at(idTrans).at(idTerm);
            if(para && para->index)
            {
              x0(normalEquationInfo.index(para->index),0) = para->trans->islBiasRecv.biases.at(idTerm);
#if DEBUG > 0
              logInfo<<"aprioriParameter() recv ISL terminal bias "
                     <<normalEquationInfo.parameterNames().at(normalEquationInfo.index(para->index)).str()
                     <<para->trans->islBiasRecv.biases.at(idTerm)%" %6.2f"s <<Log::endl;
#endif
              }
          }
      } // end for(idTrans...
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
    logInfo<<"GnssParametrizationIslBiases::designMatrixIsl() "
           <<eqn.transmitter->name()<< eqn.terminalSend%"(%2i)"s<<" -> "
           <<eqn.receiver->name()<< eqn.terminalRecv%"(%2i)"s
           <<Log::endl;
#endif
    UInt idTerm;

    // transmitter terminal bias
    // -------------------------
    if(eqn.transmitter->islBiasSend.isInList(eqn.terminalSend, idTerm))
    {
#if DEBUG > 1
      logInfo<<"GnssParametrizationIslBiases::designMatrixIsl() send "
             <<eqn.transmitter->name()<<eqn.terminalSend%"(%2i)"s<<idTerm%" (idx %i)"s
             <<Log::endl;
#endif
      auto para = paraTransmitTerminal.at(eqn.transmitter->idTrans()).at(idTerm);
      if(para && para->index)
        copy(eqn.A.column(GnssObservationEquationIsl::idxRange,1), A.column(para->index));
    }

    // receiver terminal bias
    // -------------------------
    if(eqn.receiver->islBiasRecv.isInList(eqn.terminalRecv, idTerm))
    {
#if DEBUG > 1
      logInfo<<"GnssParametrizationIslBiases::designMatrixIsl() recv "
             <<eqn.receiver->name()<<eqn.terminalRecv%"(%2i)"s<<idTerm%" (idx %i)"s
             <<Log::endl;
#endif
      auto para =paraReceiveTerminal.at(eqn.receiver->idTrans()).at(idTerm);
      if(para && para->index)
        copy(eqn.A.column(GnssObservationEquationIsl::idxRange,1), A.column(para->index));
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
    for(UInt idTrans=0; idTrans<gnss->transmitters.size(); idTrans++)
    {
      for(auto para : paraTransmitTerminal.at(idTrans))
        if(para && para->index)
          countTrans++;
      for(auto para : paraReceiveTerminal.at(idTrans))
        if(para && para->index)
          countRecv++;
    }
    UInt count = countTrans + countRecv;
    if(!count)
      return;

    logWarning<<"ISL bias zer-mean constraints temporarily disabled!"<<Log::endl;

    /*
    // collect apriori bias values
    // ---------------------------

    if(Parallel::isMaster(normalEquationInfo.comm))
    {
      Vector           l(1);
      GnssDesignMatrix A(normalEquationInfo, 1);
      const std::vector<UInt> terminals = { 0 }; // TODO: use variable terminal ID!
      for(auto para : paraTransmitTerminal)
        if(para && para->index && selectedTransmitTerminalZeroMean.at(para->trans->idTrans()))
        {
          l(0) += (x0TransmitTerminal.at(para->trans->idTrans()) - para->trans->sendIslBias(terminals).at(0))/count/sigmaZeroMean; // remove apriori value -> regularization towards 0
          A.column(para->index)(0,0) = 1./count/sigmaZeroMean;
#if DEBUG > 0
          logInfo<<"constraint() send ISL terminal bias "<<para->trans->name()
                 <<x0TransmitTerminal.at(para->trans->idTrans())%" a-priori %6.2f m"s
                 <<para->trans->sendIslBias(terminals).at(0)%" estimate %6.2f m"s
                 <<Log::endl;
#endif
        }
      for(auto para : paraReceiveTerminal)
        if(para && para->index && selectedReceiveTerminalZeroMean.at(para->trans->idTrans()))
        {
          l(0) += (x0ReceiveTerminal.at(para->trans->idTrans()) - para->trans->recvIslBias(terminals).at(0))/count/sigmaZeroMean; // remove apriori value -> regularization towards 0
          A.column(para->index)(0,0) = 1./count/sigmaZeroMean;
#if DEBUG > 0
          logInfo<<"constraint() recv ISL terminal bias "<<para->trans->name()
                 <<x0ReceiveTerminal.at(para->trans->idTrans())%" a-priori %6.2f m"s
                 <<para->trans->recvIslBias(terminals).at(0)%" estimate %6.2f m"s
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
    Gnss::InfoParameterChange infoSend("mm");
    Gnss::InfoParameterChange infoRecv("mm");

    for(UInt idTrans=0; idTrans<gnss->transmitters.size(); idTrans++)
    {
      for(UInt idTerm=0; idTerm<gnss->transmitters.at(idTrans)->islBiasSend.terminals.size(); idTerm++)
        if(paraTransmitTerminal.at(idTrans).size() && idTerm<paraTransmitTerminal.at(idTrans).size())
        {
          auto para = paraTransmitTerminal.at(idTrans).at(idTerm);
          if(para && para->index)
          {
            double dBias = x(normalEquationInfo.index(para->index),0);
            para->trans->islBiasSend.biases.at(idTerm) += dBias;
            if(infoSend.update(1e3*dBias))
              infoSend.info = "send ISL terminal bias ("+normalEquationInfo.parameterNames().at(normalEquationInfo.index(para->index)).str()+")";
          }
        }

      for(UInt idTerm=0; idTerm<gnss->transmitters.at(idTrans)->islBiasRecv.terminals.size(); idTerm++)
        if(paraReceiveTerminal.at(idTrans).size() && idTerm<paraReceiveTerminal.at(idTrans).size())
        {
          auto para = paraReceiveTerminal.at(idTrans).at(idTerm);
          if(para && para->index)
          {
            double dBias = x(normalEquationInfo.index(para->index),0);
            para->trans->islBiasRecv.biases.at(idTerm) += dBias;
            if(infoRecv.update(1e3*dBias))
              infoRecv.info = "recv ISL terminal bias ("+normalEquationInfo.parameterNames().at(normalEquationInfo.index(para->index)).str()+")";
          }
        }
    }
    infoSend.synchronizeAndPrint(normalEquationInfo.comm, 1e-3, maxChange);
    infoRecv.synchronizeAndPrint(normalEquationInfo.comm, 1e-3, maxChange);

    return maxChange;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationIslBiases::writeResults(const GnssNormalEquationInfo &normalEquationInfo, const std::string &suffix) const
{
  try
  {
    if(!isEnabled(normalEquationInfo, name))
      return;

    if(!fileNameOutTransmitter.empty() && Parallel::isMaster(normalEquationInfo.comm))
    {
      VariableList fileNameVariableList;
      fileNameVariableList.setVariable("prn", "***");
      logStatus<<"write transmit ISL terminal biases to files <"<<fileNameOutTransmitter(fileNameVariableList).appendBaseName(suffix)<<">"<<Log::endl;
      auto selectedTransmitters = gnss->selectTransmitters(selectTransmitTerminal);
      for(auto trans : gnss->transmitters)
        if(trans->useable() && selectedTransmitters.at(trans->idTrans()))
        {
          IslBias islBias = trans->islBiasSend;
          fileNameVariableList.setVariable("prn", trans->name());
          writeFileIslBias(fileNameOutTransmitter(fileNameVariableList).appendBaseName(suffix), islBias);
        }
    }

    if(!fileNameOutReceiver.empty() && Parallel::isMaster(normalEquationInfo.comm))
    {
      VariableList fileNameVariableList;
      fileNameVariableList.setVariable("prn", "***");
      logStatus<<"write receive ISL terminal biases to files <"<<fileNameOutReceiver(fileNameVariableList).appendBaseName(suffix)<<">"<<Log::endl;
      auto selectedTransmitters = gnss->selectTransmitters(selectReceiveTerminal);
      for(auto trans : gnss->transmitters)
        if(trans->useable() && selectedTransmitters.at(trans->idTrans()))
        {
          IslBias islBias = trans->islBiasRecv;
          fileNameVariableList.setVariable("prn", trans->name());
          writeFileIslBias(fileNameOutReceiver(fileNameVariableList).appendBaseName(suffix), islBias);
        }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
