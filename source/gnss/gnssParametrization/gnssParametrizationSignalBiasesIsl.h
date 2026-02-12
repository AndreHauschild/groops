/***********************************************/
/**
* @file gnssParametrizationSignalBiasesIsl.h
*
* @brief Inter-satellite link signal biases.
* @see GnssParametrization
*
* @author Andre Hauschild
* @date 2024-10-21
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPARAMETRIZATIONSIGNALBIASESISL__
#define __GROOPS_GNSSPARAMETRIZATIONSIGNALBIASESISL__

// Latex documentation
#ifdef DOCSTRING_GnssParametrization
static const char *docstringGnssParametrizationSignalBiasesIsl = R"(
\subsection{SignalBiasesIsl}\label{gnssParametrizationType:signalBiasesIsl}
Each inter-satellite-link observation contains a bias for ISL terminal of the transmitting and receiving satellite
\begin{equation}
  [\tau\nu a]_r^s(t) = \dots + \text{bias}[\tau\nu a]^s + \text{bias}[\tau\nu a]_r + \dots
\end{equation}
This class provides the apriori model $\M f(\M x_0)$ of eq. \eqref{gnssParametrizationType:model} only.

The \configFile{inputfileSignalBiasTransmitter/Receiver}{gnssSignalBias} are read
for each satellite. The transmitter and receiver biases are stored in separate files.
The file name is interpreted as a template with the variable \verb|{prn}| being replaced by the satellite PRN.
(Infos regarding the variable \verb|{prn}| can be found in
\configClass{gnssTransmitterGeneratorType}{gnssTransmitterGeneratorType}).
The bias files must contain one bias with the \verb|{C1C}| observation code.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "classes/platformSelector/platformSelector.h"
#include "gnss/gnss.h"
#include "gnss/gnssParametrization/gnssParametrization.h"

/***** CLASS ***********************************/

/** @brief Inter-satellite link signal biases.
* @ingroup gnssParametrizationGroup
* @see GnssParametrization */
class GnssParametrizationSignalBiasesIsl : public GnssParametrizationBase
{
  Gnss                *gnss;
  std::string          name;
  PlatformSelectorPtr  selectTransmitters, selectReceivers;
  FileName             fileNameOutTransmitter, fileNameOutReceiver;
  FileName             fileNameInTransmitter, fileNameInReceiver;

public:
  GnssParametrizationSignalBiasesIsl(Config &config);

  void init(Gnss *gnss, Parallel::CommunicatorPtr comm) override;
  void writeResults(const GnssNormalEquationInfo &normalEquationInfo, const std::string &suffix) const override;
};

/***********************************************/

#endif
