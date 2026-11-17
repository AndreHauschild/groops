/***********************************************/
/**
* @file gnssIslBias2Matrix.cpp
*
* @brief Computes ISL terminal biases for a given type list.
*
* @author Andre Hauschild
* @date 2026-11-17
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Computes ISL terminal biases for a given list of terminal IDs.
If the terminal list is empty, all terminals contained in \configFile{inputfileIslBias}{islBias} are used.
The resulting \configFile{outputfileMatrix}{matrix} contains a vector with an entry for each terminal.

See also \program{GnssMatrix2IslBias}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileIslBias.h"
#include "files/fileMatrix.h"
#include "files/fileStringTable.h"

/***** CLASS ***********************************/

/** @brief Computes ISL terminal biases for a given type list.
* @ingroup programsGroup */
class GnssIslBias2Matrix
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GnssIslBias2Matrix, SINGLEPROCESS, "Computes ISL terminal biases for a given terminal list", Gnss)


/***********************************************/

void GnssIslBias2Matrix::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNameIslBias, fileNameMatrix, fileNameTerminals;
    std::vector<UInt> terminals;

    readConfig(config, "outputfileMatrix",    fileNameMatrix,    Config::MUSTSET,  "", "");
    readConfig(config, "outputfileTerminals", fileNameTerminals, Config::OPTIONAL, "", "ASCII list of terminal IDs");
    readConfig(config, "inputfileIslBias",    fileNameIslBias,   Config::MUSTSET,  "", "");
    readConfig(config, "terminals",           terminals,         Config::OPTIONAL, "", "if not set, all terminals in the file are used");
    if(isCreateSchema(config)) return;

    logStatus<<"read ISL terminal bias file <"<<fileNameIslBias<<">"<<Log::endl;

    IslBias islBias;
    readFileIslBias(fileNameIslBias, islBias);

    if(terminals.size() == 0)
      terminals = islBias.terminals;
    Vector bias = islBias.compute(terminals);

    if(!fileNameTerminals.empty())
    {
      logStatus<<"write terminals to file <"<<fileNameTerminals<<">"<<Log::endl;
      std::vector<std::string> terminalNames;
      for(const auto &terminal : terminals)
        terminalNames.push_back(terminal%"%i"s);
      writeFileStringList(fileNameTerminals, terminalNames);
    }

    logStatus<<"write matrix to file <"<<fileNameMatrix<<">"<<Log::endl;
    writeFileMatrix(fileNameMatrix, bias);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
