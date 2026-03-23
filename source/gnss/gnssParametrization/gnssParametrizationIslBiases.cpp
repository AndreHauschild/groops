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
    readConfig(config, "selectTransmitterSendTerminal",  selectTransmitterSendTerminal,  Config::DEFAULT,  R"(["all"])", "");
    readConfig(config, "selectTransmitterRecvTerminal",  selectTransmitterRecvTerminal,  Config::DEFAULT,  R"(["all"])", "");
    readConfig(config, "outputfileSendTerminal",         fileNameOutSendTerminal,        Config::OPTIONAL, "", "variable {prn} available");
    readConfig(config, "outputfileRecvTerminal",         fileNameOutRecvTerminal,        Config::OPTIONAL, "", "variable {prn} available");
    readConfig(config, "inputfileSendTerminal",          fileNameInSendTerminal,         Config::OPTIONAL, "", "variable {prn} available");
    readConfig(config, "inputfileRecvTerminal",          fileNameInRecvTerminal,         Config::OPTIONAL, "", "variable {prn} available");
    readConfig(config, "nameConstraint",                 nameConstraint,                 Config::OPTIONAL, "constraint.islBiases", "used for parameter selection");
    readConfig(config, "selectSendTerminalZeroMean",     selectSendTerminalZeroMean,     Config::DEFAULT,  R"(["all"])", "");
    readConfig(config, "selectRecvTerminalZeroMean",     selectRecvTerminalZeroMean,     Config::DEFAULT,  "", "");
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
//       Therefore the a priori bias values are loaded here,, but the parameter
//       indices are created in the call of initParameter() depending on the
//       availability of ISL observations.

void GnssParametrizationIslBiases::init(Gnss *gnss, Parallel::CommunicatorPtr comm)
{
  try
  {
    this->gnss = gnss;

    // Load a priori values
    //---------------------

    if(!fileNameInSendTerminal.empty())
    {
      VariableList fileNameVariableList;
      auto selectedTransmitters = gnss->selectTransmitters(selectTransmitterSendTerminal);
      for(UInt idTrans=0; idTrans<gnss->transmitters.size(); idTrans++)
        if(selectedTransmitters.at(idTrans) && gnss->transmitters.at(idTrans)->useable())
        {
          fileNameVariableList.setVariable("prn", gnss->transmitters.at(idTrans)->name());
          try
          {
            readFileIslBias(fileNameInSendTerminal(fileNameVariableList), gnss->transmitters.at(idTrans)->islBiasSend);
          }
          catch(std::exception &/*e*/)
          {
            logWarningOnce<<"Unable to read transmit ISL terminal bias file <"<<fileNameInSendTerminal(fileNameVariableList)<<">, disabling transmitter."<<Log::endl;
            gnss->transmitters.at(idTrans)->disable("Unable to read transmit ISL terminal bias file <"+fileNameInSendTerminal(fileNameVariableList).str()+">");
          }
        }
    }

    if(!fileNameInRecvTerminal.empty())
    {
      VariableList fileNameVariableList;
      auto selectedTransmitters = gnss->selectTransmitters(selectTransmitterRecvTerminal);
      for(UInt idTrans=0; idTrans<gnss->transmitters.size(); idTrans++)
        if(selectedTransmitters.at(idTrans) && gnss->transmitters.at(idTrans)->useable())
        {
          fileNameVariableList.setVariable("prn", gnss->transmitters.at(idTrans)->name());
          try
          {
            readFileIslBias(fileNameInRecvTerminal(fileNameVariableList), gnss->transmitters.at(idTrans)->islBiasRecv);
          }
          catch(std::exception &/*e*/)
          {
            logWarningOnce<<"Unable to read receive ISL terminal bias file <"<<fileNameInRecvTerminal(fileNameVariableList)<<">, disabling transmitter."<<Log::endl;
            gnss->transmitters.at(idTrans)->disable("Unable to read receive ISL terminal bias file <"+fileNameInRecvTerminal(fileNameVariableList).str()+">");
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

    auto selectedTransmitterSendTerminal = gnss->selectTransmitters(selectTransmitterSendTerminal);
    paraSendTerminal.clear();
    paraSendTerminal.resize(gnss->transmitters.size());
    for(UInt idTrans=0; idTrans<gnss->transmitters.size(); idTrans++)
      if(selectedTransmitterSendTerminal.at(idTrans) && gnss->transmitters.at(idTrans)->useable())
        for(UInt idTerm=0; idTerm<gnss->transmitters.at(idTrans)->islBiasSend.terminals.size(); idTerm++)
        {
          auto para = new Parameter();
          paraSendTerminal.at(idTrans).push_back(para);
          para->trans = gnss->transmitters.at(idTrans);
          para->terminal = gnss->transmitters.at(idTrans)->islBiasSend.terminals.at(idTerm);
          para->index = GnssParameterIndex();
#if DEBUG > 0
          logInfo<<"initParameter() send ISL terminal bias parameter "
                 <<para->trans->name()<<para->terminal%"(%i)"s<<Log::endl;
#endif
        }

    auto selectedTransmitterRecvTerminal = gnss->selectTransmitters(selectTransmitterRecvTerminal);
    paraRecvTerminal.clear();
    paraRecvTerminal.resize(gnss->transmitters.size());
    for(UInt idTrans=0; idTrans<gnss->transmitters.size(); idTrans++)
      if(selectedTransmitterRecvTerminal.at(idTrans) && gnss->transmitters.at(idTrans)->useable())
        for(UInt idTerm=0; idTerm<gnss->transmitters.at(idTrans)->islBiasRecv.terminals.size(); idTerm++)
        {
          auto para = new Parameter();
          paraRecvTerminal.at(idTrans).push_back(para);
          para->trans = gnss->transmitters.at(idTrans);
          para->terminal = gnss->transmitters.at(idTrans)->islBiasRecv.terminals.at(idTerm);
          para->index = GnssParameterIndex();
#if DEBUG > 0
          logInfo<<"initParameter() recv ISL terminal bias parameter "
                 <<para->trans->name()<<para->terminal%"(%i)"s<<Log::endl;
#endif
        }

    applyConstraint = FALSE;
    if(!isEnabled(normalEquationInfo, name))
      return;

    // transmit ISL terminal parameters
    // --------------------------------
    UInt countparaSendTerminal = 0;
    for(UInt idTrans=0; idTrans<gnss->transmitters.size(); idTrans++)
      for(auto para : paraSendTerminal.at(idTrans))
      {
        std::vector<ParameterName> parameterNames;
        parameterNames.push_back(ParameterName(para->trans->name(), para->terminal%"islBiasSend.%i"s));
        para->index = normalEquationInfo.parameterNamesTransmitter(para->trans->idTrans(), parameterNames);
        countparaSendTerminal += parameterNames.size();
#if DEBUG > 0
        logInfo<<"initParameter() send ISL terminal bias parameter "<<parameterNames.at(0).str()<<Log::endl;
#endif
      }
    if(countparaSendTerminal)
      logInfo<<countparaSendTerminal%"%9i ISL send terminal bias parameters"s<<Log::endl;

    // receive ISL terminal parameters
    // -------------------------------
    UInt countparaRecvTerminal = 0;
    for(UInt idTrans=0; idTrans<gnss->transmitters.size(); idTrans++)
      for(auto para : paraRecvTerminal.at(idTrans))
      {
        std::vector<ParameterName> parameterNames;
        parameterNames.push_back(ParameterName(para->trans->name(), para->terminal%"islBiasRecv.%i"s));
        para->index = normalEquationInfo.parameterNamesTransmitter(para->trans->idTrans(), parameterNames);
        countparaRecvTerminal += parameterNames.size();
#if DEBUG > 0
        logInfo<<"initParameter() recv ISL terminal bias parameter "<<parameterNames.at(0).str()<<Log::endl;
#endif
      }
    if(countparaRecvTerminal)
      logInfo<<countparaRecvTerminal%"%9i ISL recv terminal bias parameters"s<<Log::endl;

    // Zero-mean constraint

    selectedSendTerminalZeroMean = gnss->selectTransmitters(selectSendTerminalZeroMean);
    selectedRecvTerminalZeroMean = gnss->selectTransmitters(selectRecvTerminalZeroMean);

    // Initialize a priori values

    // NOTE: the a priori values for the zero-mean constraint are stored in
    //       this function and used in the constraints() function.

    UInt countZeroMean = 0;
    x0SendTerminal.clear();
    x0RecvTerminal.clear();
    x0SendTerminal.resize(gnss->transmitters.size());
    x0RecvTerminal.resize(gnss->transmitters.size());
    for(UInt idTrans=0; idTrans<gnss->transmitters.size(); idTrans++)
    {
      for(auto para : paraSendTerminal.at(idTrans))
        if(para && para->index && selectedSendTerminalZeroMean.at(idTrans))
        {
          x0SendTerminal.at(idTrans).push_back(para->trans->sendIslBias(para->terminal));
          countZeroMean++;
#if DEBUG >0
          logInfo<<"initParameter() store initial send ISL terminal bias parameter "
                 <<para->trans->name()<<para->terminal%"(%i)"s
                 <<para->trans->sendIslBias(para->terminal)%": %5.2f m"s<<Log::endl;
#endif
        }
      for(auto para : paraRecvTerminal.at(idTrans))
        if(para && para->index && selectedRecvTerminalZeroMean.at(idTrans))
        {
          x0RecvTerminal.at(idTrans).push_back(para->trans->recvIslBias(para->terminal));
          countZeroMean++;
#if DEBUG >0
          logInfo<<"initParameter() store initial recv ISL terminal bias parameter "
                 <<para->trans->name()<<para->terminal%"(%i)"s
                 <<para->trans->recvIslBias(para->terminal)%": %5.2f m"s<<Log::endl;
#endif
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
      for(UInt idTrans=0; idTrans<gnss->transmitters.size(); idTrans++)
      {
        for(auto para : paraSendTerminal.at(idTrans))
          if(para && para->index)
          {
            x0(normalEquationInfo.index(para->index),0) = para->trans->sendIslBias(para->terminal);
#if DEBUG > 1
            logInfo<<"aprioriParameter() send ISL terminal bias "
                <<normalEquationInfo.parameterNames().at(normalEquationInfo.index(para->index)).str()
                <<para->trans->sendIslBias(para->terminal)%" %6.2f m"s <<Log::endl;
#endif
          }
        for(auto para : paraRecvTerminal.at(idTrans))
          if(para && para->index)
          {
            x0(normalEquationInfo.index(para->index),0) = para->trans->recvIslBias(para->terminal);
#if DEBUG > 1
            logInfo<<"aprioriParameter() recv ISL terminal bias "
                <<normalEquationInfo.parameterNames().at(normalEquationInfo.index(para->index)).str()
                <<para->trans->islBiasRecv.biases.at(idTerm)%" %6.2f m"s <<Log::endl;
#endif
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
#if DEBUG > 1
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
#if DEBUG > 2
      logInfo<<"GnssParametrizationIslBiases::designMatrixIsl() send "
             <<eqn.transmitter->name()<<eqn.terminalSend%"(%2i)"s<<idTerm%" (idx %i)"s
             <<Log::endl;
#endif
      auto para = paraSendTerminal.at(eqn.transmitter->idTrans()).at(idTerm);
      if(para && para->index)
        copy(eqn.A.column(GnssObservationEquationIsl::idxRange,1), A.column(para->index));
    }

    // receiver terminal bias
    // ----------------------
    if(eqn.receiver->islBiasRecv.isInList(eqn.terminalRecv, idTerm))
    {
#if DEBUG > 2
      logInfo<<"GnssParametrizationIslBiases::designMatrixIsl() recv "
             <<eqn.receiver->name()<<eqn.terminalRecv%"(%2i)"s<<idTerm%" (idx %i)"s
             <<Log::endl;
#endif
      auto para =paraRecvTerminal.at(eqn.receiver->idTrans()).at(idTerm);
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

    UInt countSend = 0, countRecv = 0;
    for(UInt idTrans=0; idTrans<gnss->transmitters.size(); idTrans++)
    {
      for(auto para : paraSendTerminal.at(idTrans))
        if(para && para->index)
          countSend++;
      for(auto para : paraRecvTerminal.at(idTrans))
        if(para && para->index)
          countRecv++;
    }
    UInt count = countSend + countRecv;
    if(!count)
      return;

    // collect a priori bias values
    // ----------------------------

    if(Parallel::isMaster(normalEquationInfo.comm))
    {
      Vector           l(1);
      GnssDesignMatrix A(normalEquationInfo, 1);
      for(UInt idTrans=0; idTrans<gnss->transmitters.size(); idTrans++)
      {
        for(auto para : paraSendTerminal.at(idTrans))
          if(para && para->index && selectedSendTerminalZeroMean.at(idTrans))
          {
            UInt idTerm = para->trans->islBiasSend.index(para->terminal);
            l(0) += (x0SendTerminal.at(idTrans).at(idTerm) - para->trans->sendIslBias(para->terminal))/count/sigmaZeroMean; // remove a priori value -> regularization towards 0
            A.column(para->index)(0,0) = 1./count/sigmaZeroMean;
#if DEBUG > 0
            logInfo<<"constraints() send ISL terminal bias "<<para->trans->name()<<para->terminal%"(%i)"s
                   <<x0SendTerminal.at(idTrans).at(idTerm)%" a priori %6.2f m"s
                   <<para->trans->sendIslBias(para->terminal)%" estimate %6.2f m"s
                   <<Log::endl;
#endif
          }
        for(auto para : paraRecvTerminal.at(idTrans))
          if(para && para->index && selectedRecvTerminalZeroMean.at(idTrans))
          {
            UInt idTerm = para->trans->islBiasRecv.index(para->terminal);
            l(0) += (x0RecvTerminal.at(idTrans).at(idTerm) - para->trans->recvIslBias(para->terminal))/count/sigmaZeroMean; // remove a priori value -> regularization towards 0
            A.column(para->index)(0,0) = 1./count/sigmaZeroMean;
#if DEBUG > 0
            logInfo<<"constraints() recv ISL terminal bias "<<para->trans->name()<<para->terminal%"(%i)"s
                   <<x0RecvTerminal.at(idTrans).at(idTerm)%" a priori %6.2f m"s
                   <<para->trans->recvIslBias(para->terminal)%" estimate %6.2f m"s
                   <<Log::endl;
#endif
          }
      }
      GnssDesignMatrix::accumulateNormals(A, l, normals, n, lPl, obsCount);
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
    Gnss::InfoParameterChange infoSend("mm");
    Gnss::InfoParameterChange infoRecv("mm");

    for(UInt idTrans=0; idTrans<gnss->transmitters.size(); idTrans++)
    {
      for(auto para : paraSendTerminal.at(idTrans))
        if(para && para->index)
        {
          UInt idTerm = para->trans->islBiasSend.index(para->terminal);
          double dBias = x(normalEquationInfo.index(para->index),0);
          para->trans->islBiasSend.biases.at(idTerm) += dBias;
          if(infoSend.update(1e3*dBias))
            infoSend.info = "send ISL terminal bias ("+normalEquationInfo.parameterNames().at(normalEquationInfo.index(para->index)).str()+")";
        }
      for(auto para : paraRecvTerminal.at(idTrans))
        if(para && para->index)
        {
          UInt idTerm = para->trans->islBiasRecv.index(para->terminal);
          double dBias = x(normalEquationInfo.index(para->index),0);
          para->trans->islBiasRecv.biases.at(idTerm) += dBias;
          if(infoRecv.update(1e3*dBias))
            infoRecv.info = "recv ISL terminal bias ("+normalEquationInfo.parameterNames().at(normalEquationInfo.index(para->index)).str()+")";
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

    if(!fileNameOutSendTerminal.empty() && Parallel::isMaster(normalEquationInfo.comm))
    {
      VariableList fileNameVariableList;
      fileNameVariableList.setVariable("prn", "***");
      logStatus<<"write transmit ISL terminal biases to files <"<<fileNameOutSendTerminal(fileNameVariableList).appendBaseName(suffix)<<">"<<Log::endl;
      auto selectedTransmitters = gnss->selectTransmitters(selectTransmitterSendTerminal);
      for(auto trans : gnss->transmitters)
        if(trans->useable() && selectedTransmitters.at(trans->idTrans()))
        {
          fileNameVariableList.setVariable("prn", trans->name());
          writeFileIslBias(fileNameOutSendTerminal(fileNameVariableList).appendBaseName(suffix), trans->islBiasSend);
        }
    }

    if(!fileNameOutRecvTerminal.empty() && Parallel::isMaster(normalEquationInfo.comm))
    {
      VariableList fileNameVariableList;
      fileNameVariableList.setVariable("prn", "***");
      logStatus<<"write receive ISL terminal biases to files <"<<fileNameOutRecvTerminal(fileNameVariableList).appendBaseName(suffix)<<">"<<Log::endl;
      auto selectedTransmitters = gnss->selectTransmitters(selectTransmitterRecvTerminal);
      for(auto trans : gnss->transmitters)
        if(trans->useable() && selectedTransmitters.at(trans->idTrans()))
        {
          fileNameVariableList.setVariable("prn", trans->name());
          writeFileIslBias(fileNameOutRecvTerminal(fileNameVariableList).appendBaseName(suffix), trans->islBiasRecv);
        }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
