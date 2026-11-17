/***********************************************/
/**
* @file gnssMatrix2SignalBias.cpp
*
* @brief Converts values for GNSS signal biases from matrix file to GnssSignalBias format.
*
* @author Andre Hauschild
* @date 2026-11-16
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Converts values for GNSS signal biases from matrix format to \file{GnssSignalBias format}{gnssSignalBias}.

Combined \verb|X| signals are transformed to actually transmitted signals.

See also \program{GnssSignalBias2Matrix}.
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "files/fileMatrix.h"
#include "files/fileStringTable.h"
#include "files/fileGnssSignalBias.h"
#include "gnss/gnssReceiver.h"

/***** CLASS ***********************************/

/** @brief Converts GNSS signal biases from matrix format to GnssSignalBias format.
* @ingroup programsGroup */
class GnssMatrix2SignalBias
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GnssMatrix2SignalBias, SINGLEPROCESS, "Converts GNSS signal biases from  matrix format to GnssSignalBias format.", Conversion, Gnss)

/***********************************************/

void GnssMatrix2SignalBias::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNameSignalBias, fileNameMatrix, fileNameTypes;
    std::vector<GnssType> types;

    readConfig(config, "outputfileSignalBias", fileNameSignalBias, Config::MUSTSET,  "", "");
    readConfig(config, "inputfileMatrix",      fileNameMatrix,     Config::MUSTSET,  "", "");
    readConfig(config, "inputfileTypes",       fileNameTypes,      Config::MUSTSET,  "", "ASCII list of types");
    readConfig(config, "types",                types,              Config::OPTIONAL, "", "if not set, all types in the file are used");
    if(isCreateSchema(config)) return;

    // ==================================================

    logStatus<<"Load signal bias types list <"<<fileNameTypes<<">"<<Log::endl;

    std::vector<std::string> typesIn;
    readFileStringList(fileNameTypes, typesIn);

    logStatus<<"Load signal bias matrix <"<<fileNameMatrix<<">"<<Log::endl;

    Matrix biasesIn;
    readFileMatrix(fileNameMatrix, biasesIn);

    if(biasesIn.columns() != 1)
      throw(Exception("Signal bias matrix must have exactly one column, found "+biasesIn.columns()%"%i"s));
    if(typesIn.size() != biasesIn.rows())
      throw(Exception("Inconsistent number of types ("+typesIn.size()%"%i"s+") and biases ("+biasesIn.rows()%"%i"s+") in input files!"));

    // select types (all if no selection given)
    GnssSignalBias signalBias;
    for(UInt i=0; i<typesIn.size(); i++)
    {
      const GnssType type = GnssType(typesIn.at(i));
      if(types.size() && !type.isInList(types))
        continue;

      signalBias.types.push_back(type);
      signalBias.biases.push_back(biasesIn(i, 0));
      logStatus<<"add signal bias type <"<<type.str()<<">"<<Log::endl;
    }

    if(!signalBias.types.size())
      throw(Exception("No signal bias types selected."));

    // ==================================================

    logStatus<<"write signal bias files <"<<fileNameSignalBias<<">"<<Log::endl;

    // Transform combined X signals to actually transmitted signals.
    Matrix T;
    GnssReceiver::signalCompositionDefault(signalBias.types, signalBias.types, T);
    T = pseudoInverse(T);
    signalBias.biases = Vector(T * Vector(signalBias.biases));

    writeFileGnssSignalBias(fileNameSignalBias, signalBias);

  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
