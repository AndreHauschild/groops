/***********************************************/
/**
* @file gnssMatrix2IslBias.cpp
*
* @brief Converts values for ISL terminal biases from matrix file to IslBias format.
*
* @author Andre Hauschild
* @date 2026-11-17
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Converts values for ISL terminal biases from matrix format to \file{IslBias format}{islBias}.

See also \program{GnssIslBias2Matrix}.
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "files/fileMatrix.h"
#include "files/fileStringTable.h"
#include "files/fileIslBias.h"

/***** CLASS ***********************************/

/** @brief Converts ISLs terminal biases from matrix format to IslBias format.
* @ingroup programsGroup */
class GnssMatrix2IslBias
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GnssMatrix2IslBias, SINGLEPROCESS, "Converts ISL terminal biases from  matrix format to IslBias format.", Conversion, Gnss)

/***********************************************/

void GnssMatrix2IslBias::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNameIslBias, fileNameMatrix, fileNameTerminals;
    std::vector<UInt> terminals;

    readConfig(config, "outputfileIslBias",  fileNameIslBias,   Config::MUSTSET,  "", "");
    readConfig(config, "inputfileMatrix",    fileNameMatrix,    Config::MUSTSET,  "", "");
    readConfig(config, "inputfileTerminals", fileNameTerminals, Config::MUSTSET,  "", "ASCII list of terminal IDs");
    readConfig(config, "terminals",          terminals,         Config::OPTIONAL, "", "if not set, all terminals in the file are used");
    if(isCreateSchema(config)) return;

    // ==================================================

    logStatus<<"Load ISL terminal list <"<<fileNameTerminals<<">"<<Log::endl;

    Matrix terminalsIn;
    readFileMatrix(fileNameTerminals, terminalsIn);

    logStatus<<"Load signal bias matrix <"<<fileNameMatrix<<">"<<Log::endl;

    Matrix biasesIn;
    readFileMatrix(fileNameMatrix, biasesIn);

    if(biasesIn.columns() != 1)
      throw(Exception("ISL bias matrix must have exactly one column, found "+biasesIn.columns()%"%i"s));
    if(terminalsIn.size() != biasesIn.rows())
      throw(Exception("Inconsistent number of terminals ("+terminalsIn.size()%"%i"s+") and biases ("+biasesIn.rows()%"%i"s+") in input files!"));

    // select terminals (all if no selection given)
    IslBias islBias;
    for(UInt i=0; i<terminalsIn.size(); i++)
    {
      UInt terminal = terminalsIn(i,0);
      if(terminals.size() /*&& !terminal.isInList(terminals)*/)
        continue;

      islBias.terminals.push_back(terminal);
      islBias.biases.push_back(biasesIn(i, 0));
      logStatus<<"add ISL terminal <"<<terminal%"%i"s<<">"<<Log::endl;
    }

    if(!islBias.terminals.size())
      throw(Exception("No ISL bias terminals selected."));

    // ==================================================

    logStatus<<"write ISL bias file <"<<fileNameIslBias<<">"<<Log::endl;

    writeFileIslBias(fileNameIslBias, islBias);

  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
