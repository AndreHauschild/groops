/***********************************************/
/**
* @file gnssParametrizationIslBiases.h
*
* @brief Inter Satellite Link biases.
* @see GnssParametrization
*
* @author Torsten Mayer-Guerr
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
// TODO: properly handle transmitter and receiver biases for ISL!!
class GnssParametrizationIslBiases : public GnssParametrizationBase
{
  class Parameter
  {
  public:
    GnssTransmitterPtr trans;
    GnssParameterIndex index;
    Matrix             Bias;
  };

  Gnss                      *gnss;
  std::string                name, nameConstraint;
  PlatformSelectorPtr        selectTransmitters;
  Bool                       applyConstraint;
  Double                     sigmaZeroMean;
  std::vector<Parameter*>    parameter;
  std::vector<UInt>          idxBias;        // indices in zeroMean matrix
  Matrix                     zeroMeanDesign; // zero mean observation equations

public:
  GnssParametrizationIslBiases(Config &config);
 ~GnssParametrizationIslBiases();

  void   init(Gnss *gnss, Parallel::CommunicatorPtr comm) override;
//void   observationCorrections(GnssObservationEquationIsl &eqn) const override;
  void   initParameter(GnssNormalEquationInfo &normalEquationInfo) override;
  void   aprioriParameter(const GnssNormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const override;
  void   designMatrixIsl(const GnssNormalEquationInfo &normalEquationInfo, const GnssObservationEquationIsl &eqn, GnssDesignMatrix &A) const override;
  void   constraints(const GnssNormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const override;
  Double updateParameter(const GnssNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef Wz) override;
};

/***********************************************/

#endif
