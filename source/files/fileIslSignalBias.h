/***********************************************/
/**
* @file fileIslSignalBias.h
*
* @brief ISL biases.
*
* @author Andre Hauschild
* @date 2026-02-02
*
*/
/***********************************************/

#ifndef __GROOPS_ISLSIGNALBIAS__
#define __GROOPS_ISLSIGNALBIAS__

// Latex documentation
#ifdef DOCSTRING_FILEFORMAT_IslSignalBias
static const char *docstringIslSignalBias = R"(
Signal biases of GNSS transmitters or receivers for different \configClass{gnssType}{gnssType}.

\begin{verbatim}
groops IslSignalBias version=20200123
          1 # number of terminals/biases
# terminal   bias [m]
# ====================================
  0          -1.752461109688110974e-01
 \end{verbatim}

NOTE: currently only one bias per satellite is supported. It must have terminal ID 0. 

See also \program{GnssProcessing}, \program{GnssSimulateIsl}.
)";
#endif

/***********************************************/

#include "base/gnssType.h"
#include "inputOutput/fileName.h"
#include "inputOutput/fileArchive.h"

/** @addtogroup filesGroup */
/// @{

/***** CONSTANTS ********************************/

const char *const FILE_ISLSIGNALBIAS_TYPE    = "islSignalBias";
constexpr UInt    FILE_ISLSIGNALBIAS_VERSION = std::max(UInt(20200123), FILE_BASE_VERSION);

/***** TYPES ***********************************/

class IslSignalBias;
typedef std::shared_ptr<IslSignalBias> IslSignalBiasPtr;

/***** CLASS ***********************************/

/** @brief Code/Phase biases. */
class IslSignalBias
{
  public:
  std::vector<UInt>   terminals;
  std::vector<Double> biases;

  Vector compute(std::vector<UInt> terminals) const;
  Bool   isInList(UInt terminal, UInt &idx) const;
  UInt   index(UInt terminal) const;
};

/***** FUNCTIONS *******************************/

template<> void save(OutArchive &ar, const IslSignalBias &x);
template<> void load(InArchive  &ar, IslSignalBias &x);

/** @brief Write into a IslSignalBias file. */
void writeFileIslSignalBias(const FileName &fileName, const IslSignalBias &x);

/** @brief Read from a IslSignalBias file. */
void readFileIslSignalBias(const FileName &fileName, IslSignalBias &x);

/// @}

/***********************************************/

#endif /* __GROOPS___ */
