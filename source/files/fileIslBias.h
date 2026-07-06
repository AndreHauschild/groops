/***********************************************/
/**
* @file fileIslBias.h
*
* @brief ISL biases.
*
* @author Andre Hauschild
* @date 2026-02-02
*
*/
/***********************************************/

#ifndef __GROOPS_ISLBIAS__
#define __GROOPS_ISLBIAS__

// Latex documentation
#ifdef DOCSTRING_FILEFORMAT_IslBias
static const char *docstringIslBias = R"(
Biases of ISL terminals on GNSS satellites.

\begin{verbatim}
groops IslBias version=20200123
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

#include "inputOutput/fileName.h"
#include "inputOutput/fileArchive.h"

/** @addtogroup filesGroup */
/// @{

/***** CONSTANTS ********************************/

const char *const FILE_ISLBIAS_TYPE    = "islBias";
constexpr UInt    FILE_ISLBIAS_VERSION = std::max(UInt(20200123), FILE_BASE_VERSION);

/***** TYPES ***********************************/

class IslBias;
typedef std::shared_ptr<IslBias> IslBiasPtr;

/***** CLASS ***********************************/

/** @brief ISL biases. */
class IslBias
{
  public:
  std::vector<UInt>   terminals;
  std::vector<Double> biases;

  Vector compute(const std::vector<UInt> &terminals) const;
  Bool   isInList(UInt terminal, UInt &idx) const;
  UInt   index(UInt terminal) const;
};

/***** FUNCTIONS *******************************/

template<> void save(OutArchive &ar, const IslBias &x);
template<> void load(InArchive  &ar, IslBias &x);

/** @brief Write into a IslBias file. */
void writeFileIslBias(const FileName &fileName, const IslBias &x);

/** @brief Read from a IslBias file. */
void readFileIslBias(const FileName &fileName, IslBias &x);

/// @}

/***********************************************/

#endif /* __GROOPS___ */
