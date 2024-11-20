/***********************************************/
/**
* @file gnssParametrizationSignalBiasesIsl.cpp
*
* @brief Inter-satellite link signal biases.
* @see GnssParametrization
*
* @author Andre Hauschild
* @date 2024-10-21
*
*/
/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "classes/platformSelector/platformSelector.h"
#include "gnss/gnss.h"
#include "gnss/gnssParametrization/gnssParametrizationSignalBiasesIsl.h"

/***********************************************/

GnssParametrizationSignalBiasesIsl::GnssParametrizationSignalBiasesIsl(Config &config)
{
  try
  {
    readConfig(config, "name",                            name,                   Config::OPTIONAL, "parameter.signalBiasesIsl", "used for parameter selection");
    readConfig(config, "selectTransmitters",              selectTransmitters,     Config::DEFAULT,  R"(["all"])", "");
    readConfig(config, "selectReceivers",                 selectReceivers,        Config::DEFAULT,  R"(["all"])", "");
    readConfig(config, "outputfileSignalBiasTransmitter", fileNameOutTransmitter, Config::OPTIONAL, "", "variable {prn} available");
    readConfig(config, "outputfileSignalBiasReceiver",    fileNameOutReceiver,    Config::OPTIONAL, "", "variable {prn} available");
    readConfig(config, "inputfileSignalBiasTransmitter",  fileNameInTransmitter,  Config::OPTIONAL, "", "variable {prn} available");
    readConfig(config, "inputfileSignalBiasReceiver",     fileNameInReceiver,     Config::OPTIONAL, "", "variable {prn} available");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationSignalBiasesIsl::init(Gnss *gnss, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    this->gnss = gnss;

    if(!fileNameInTransmitter.empty())
    {
      VariableList fileNameVariableList;
      auto selectedTransmitters = gnss->selectTransmitters(selectTransmitters);
      for(UInt idTrans=0; idTrans<gnss->transmitters.size(); idTrans++)
        if(selectedTransmitters.at(idTrans) && gnss->transmitters.at(idTrans)->useable())
        {
          fileNameVariableList.setVariable("prn", gnss->transmitters.at(idTrans)->name());
          try
          {
            readFileGnssSignalBias(fileNameInTransmitter(fileNameVariableList), gnss->transmitters.at(idTrans)->signalBiasIslTx);
          }
          catch(std::exception &/*e*/)
          {
            logWarningOnce<<"Unable to read ISL transmitter signal bias file <"<<fileNameInTransmitter(fileNameVariableList)<<">, disabling transmitter."<<Log::endl;
            gnss->transmitters.at(idTrans)->disable("Unable to read signal bias file <"+fileNameInTransmitter(fileNameVariableList).str()+">");
          }
        }
    }

    if(!fileNameInReceiver.empty())
    {
      VariableList fileNameVariableList;
      auto selectedTransmitters = gnss->selectTransmitters(selectTransmitters);
      for(UInt idTrans=0; idTrans<gnss->transmitters.size(); idTrans++)
        if(selectedTransmitters.at(idTrans) && gnss->transmitters.at(idTrans)->useable())
        {
          fileNameVariableList.setVariable("prn", gnss->transmitters.at(idTrans)->name());
          try
          {
            readFileGnssSignalBias(fileNameInReceiver(fileNameVariableList), gnss->transmitters.at(idTrans)->signalBiasIslRx);
          }
          catch(std::exception &/*e*/)
          {
            logWarningOnce<<"Unable to read ISL receiver signal bias file <"<<fileNameInReceiver(fileNameVariableList)<<">, disabling transmitter."<<Log::endl;
            gnss->transmitters.at(idTrans)->disable("Unable to read signal bias file <"+fileNameInReceiver(fileNameVariableList).str()+">");
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

void GnssParametrizationSignalBiasesIsl::writeResults(const GnssNormalEquationInfo &normalEquationInfo, const std::string &suffix) const
{
  try
  {
    if(!isEnabled(normalEquationInfo, name))
      return;

    if(!fileNameOutTransmitter.empty() && !normalEquationInfo.isEachReceiverSeparately && Parallel::isMaster(normalEquationInfo.comm))
    {
      VariableList fileNameVariableList;
      fileNameVariableList.setVariable("prn", "***");
      logStatus<<"write ISL transmitter signal biases to files <"<<fileNameOutTransmitter(fileNameVariableList).appendBaseName(suffix)<<">"<<Log::endl;
      auto selectedTransmitters = gnss->selectTransmitters(selectTransmitters);
      for(auto trans : gnss->transmitters)
        if(trans->useable() && selectedTransmitters.at(trans->idTrans()))
        {
          GnssSignalBias signalBias = trans->signalBiasIslTx;
          for(UInt idType=0; idType<signalBias.types.size(); idType++)
            if(signalBias.types.at(idType) == GnssType::PHASE)
              signalBias.biases.at(idType) = std::remainder(signalBias.biases.at(idType), signalBias.types.at(idType).wavelength());
          fileNameVariableList.setVariable("prn", trans->name());
          writeFileGnssSignalBias(fileNameOutTransmitter(fileNameVariableList).appendBaseName(suffix), signalBias);
        }
    }

    if(!fileNameOutReceiver.empty())
    {
      VariableList fileNameVariableList;
      fileNameVariableList.setVariable("prn", "***");
      logStatus<<"write ISL receiver signal biases to files <"<<fileNameOutReceiver(fileNameVariableList).appendBaseName(suffix)<<">"<<Log::endl;
      auto selectedTransmitters = gnss->selectTransmitters(selectTransmitters);
      for(auto trans : gnss->transmitters)
        if(trans->useable() && selectedTransmitters.at(trans->idTrans()))
        {
          GnssSignalBias signalBias = trans->signalBiasIslRx;
          for(UInt idType=0; idType<signalBias.types.size(); idType++)
            if(signalBias.types.at(idType) == GnssType::PHASE)
              signalBias.biases.at(idType) = std::remainder(signalBias.biases.at(idType), signalBias.types.at(idType).wavelength());
          fileNameVariableList.setVariable("prn", trans->name());
          writeFileGnssSignalBias(fileNameOutReceiver(fileNameVariableList).appendBaseName(suffix), signalBias);
        }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
