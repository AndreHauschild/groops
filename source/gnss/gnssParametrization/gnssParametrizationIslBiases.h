/***********************************************/
/**
* @file gnssParametrizationIslBiases.h
*
* @brief Inter Satellite Link biases.
* @see GnssParametrization
*
* @author Torsten Mayer-Guerr
* @author Andre Hauschild
* @date 2024-01-30
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPARAMETRIZATIONISLBIASES__
#define __GROOPS_GNSSPARAMETRIZATIONISLBIASES__

// Latex documentation
#ifdef DOCSTRING_GnssParametrization
static const char *docstringGnssParametrizationIslBiases = R"(
\subsection{IslBiases}\label{gnssParametrizationType:islBiases}
Each inter-satellite-link observation contains a bias for the ISL terminals of the transmitting and receiving satellite
\begin{equation}
  []_{r,i}^{s,k}(t) = \dots + \text{bias}[]^{s,k} + \text{bias}[]_{r,i} + \dots
\end{equation}

The \configFile{inputfileIslBiasTransmitter/Receiver}{islBias} are read
for each satellite. The biases of all transmitting and receiving terminals must be stored in separate files.
The corresponding terminal must be identified by an integer number.
The file name is interpreted as a template with the variable \verb|{prn}| being replaced by the satellite PRN.
(Infos regarding the variable \verb|{prn}| can be found in
\configClass{gnssTransmitterGeneratorType}{gnssTransmitterGeneratorType}).
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "classes/platformSelector/platformSelector.h"
#include "gnss/gnss.h"
#include "gnss/gnssParametrization/gnssParametrization.h"

/***** CLASS ***********************************/

/** @brief Code biases.
* @ingroup gnssParametrizationGroup
* @see GnssParametrization */
class GnssParametrizationIslBiases : public GnssParametrizationBase
{
  class Parameter
  {
  public:
    GnssTransmitterPtr trans;
    GnssParameterIndex index;
  };

  Gnss                     *gnss;
  std::string               name, nameConstraint;
  PlatformSelectorPtr       selectTransmitTerminal, selectReceiveTerminal;
  PlatformSelectorPtr       selectTransmitTerminalZeroMean, selectReceiveTerminalZeroMean;
  std::vector<Byte>         selectedTransmitTerminalZeroMean, selectedReceiveTerminalZeroMean;
  Bool                      applyConstraint;
  Double                    sigmaZeroMean;
  std::vector<std::vector<Parameter*>>  paraTransmitTerminal, paraReceiveTerminal;
  std::vector<std::vector<Double>>      x0TransmitTerminal, x0ReceiveTerminal;      // a-priori values for each trans/recv terminal

  FileName                  fileNameOutTransmitter, fileNameOutReceiver;
  FileName                  fileNameInTransmitter, fileNameInReceiver;

public:
  GnssParametrizationIslBiases(Config &config);

  void   init(Gnss *gnss, Parallel::CommunicatorPtr comm) override;
  void   initParameter(GnssNormalEquationInfo &normalEquationInfo) override;
  void   aprioriParameter(const GnssNormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const override;
  void   designMatrixIsl(const GnssNormalEquationInfo &normalEquationInfo, const GnssObservationEquationIsl &eqn, GnssDesignMatrix &A) const override;
  void   constraints(const GnssNormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const override;
  Double updateParameter(const GnssNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef Wz) override;
  void   writeResults(const GnssNormalEquationInfo &normalEquationInfo, const std::string &suffix) const override;
};

/***********************************************/

#endif
