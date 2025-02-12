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
TODO!
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
  PlatformSelectorPtr       selectSendTerminal, selectRecvTerminal;
  PlatformSelectorPtr       selectSendTerminalZeroMean, selectRecvTerminalZeroMean;
  std::vector<Byte>         selectedSendTerminalZeroMean, selectedRecvTerminalZeroMean;
  Bool                      applyConstraint;
  Double                    sigmaZeroMean;
  std::vector<Parameter*>   paraSendTerminal, paraRecvTerminal;
  std::vector<Double>       x0SendTerminal, x0RecvTerminal; // a-priori values for each send/recv terminal

public:
  GnssParametrizationIslBiases(Config &config);
 ~GnssParametrizationIslBiases();

  void   init(Gnss *gnss, Parallel::CommunicatorPtr comm) override;
  void   initParameter(GnssNormalEquationInfo &normalEquationInfo) override;
  void   aprioriParameter(const GnssNormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const override;
  void   designMatrixIsl(const GnssNormalEquationInfo &normalEquationInfo, const GnssObservationEquationIsl &eqn, GnssDesignMatrix &A) const override;
  void   constraints(const GnssNormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const override;
  Double updateParameter(const GnssNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef Wz) override;
};

/***********************************************/

#endif
