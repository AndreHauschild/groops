/***********************************************/
/**
* @file gnssProcessingStepPrintResidualStatisticsIsl.h
*
* @brief GNSS processing step: PrintResidualStatisticsIsl.
*
* @author Andre Hauschild
* @date 2026-03-26
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPROCESSINGSTEPPRINTRESIDUALSTATISTICSISL__
#define __GROOPS_GNSSPROCESSINGSTEPPRINTRESIDUALSTATISTICSISL__

// Latex documentation
#ifdef DOCSTRING_GnssProcessingStep
static const char *docstringGnssProcessingStepPrintResidualStatisticsIsl = R"(
\subsection{PrintResidualStatisticsIsl}\label{gnssProcessingStepType:PrintResidualStatisticsIsl}
Print residual statistics for ISL observations.
\begin{verbatim}
  E02 : ISL   : factor =  0.00, sigma0 = 1.00, count =  5758, outliers =     0 (0.00 \%)
  E03 : ISL   : factor =  0.00, sigma0 = 1.00, count =  5758, outliers =     0 (0.00 \%)
  E04 : ISL   : factor =  0.00, sigma0 = 1.00, count =  5758, outliers =     0 (0.00 \%)
  E05 : ISL   : factor =  0.00, sigma0 = 0.99, count =  5758, outliers =     0 (0.00 \%)
  E06 : ISL   : factor =  0.00, sigma0 = 0.97, count =  5758, outliers =     0 (0.00 \%)
  E07 : ISL   : factor =  0.00, sigma0 = 0.99, count =  5758, outliers =     0 (0.00 \%)
  E08 : ISL   : factor =  0.00, sigma0 = 0.99, count =  5758, outliers =     0 (0.00 \%)
  E09 : ISL   : factor =  0.00, sigma0 = 0.99, count =  5758, outliers =     0 (0.00 \%)
  E10 : ISL   : factor =  0.00, sigma0 = 1.02, count =  5758, outliers =     0 (0.00 \%)
  E11 : ISL   : factor =  0.00, sigma0 = 1.02, count =  5758, outliers =     0 (0.00 \%)
  E12 : ISL   : factor =  0.00, sigma0 = 1.01, count =  5758, outliers =     0 (0.00 \%)
  E13 : ISL   : factor =  0.00, sigma0 = 0.97, count =  5758, outliers =     0 (0.00 \%)
  E15 : ISL   : factor =  0.00, sigma0 = 0.99, count =  5758, outliers =     0 (0.00 \%)
  E16 : ISL   : factor =  0.00, sigma0 = 0.98, count =  5758, outliers =     0 (0.00 \%)
  E19 : ISL   : factor =  0.00, sigma0 = 1.01, count =  5758, outliers =     0 (0.00 \%)
  E21 : ISL   : factor =  0.00, sigma0 = 1.00, count =  5758, outliers =     0 (0.00 \%)
  E23 : ISL   : factor =  0.00, sigma0 = 1.00, count =  5758, outliers =     0 (0.00 \%)
  E24 : ISL   : factor =  0.00, sigma0 = 1.01, count =  5758, outliers =     0 (0.00 \%)
  E25 : ISL   : factor =  0.00, sigma0 = 0.99, count =  5758, outliers =     0 (0.00 \%)
  E26 : ISL   : factor =  0.00, sigma0 = 1.00, count =  5758, outliers =     0 (0.00 \%)
  E27 : ISL   : factor =  0.00, sigma0 = 0.99, count =  5758, outliers =     0 (0.00 \%)
  E29 : ISL   : factor =  0.00, sigma0 = 1.00, count =  5758, outliers =     0 (0.00 \%)
  E30 : ISL   : factor =  0.00, sigma0 = 1.00, count =  5758, outliers =     0 (0.00 \%)
  E31 : ISL   : factor =  0.00, sigma0 = 0.99, count =  5758, outliers =     0 (0.00 \%)
  E33 : ISL   : factor =  0.00, sigma0 = 0.99, count =  5758, outliers =     0 (0.00 \%)
  E34 : ISL   : factor =  0.00, sigma0 = 0.98, count =  5758, outliers =     0 (0.00 \%)
  E36 : ISL   : factor =  0.00, sigma0 = 1.00, count =  5758, outliers =     0 (0.00 \%)
  ...
\end{verbatim}
)";
#endif

/***********************************************/

#include "config/config.h"
#include "gnss/gnssProcessingStep/gnssProcessingStep.h"

/***** CLASS ***********************************/

/** @brief GNSS processing step: PrintResidualStatisticsIsl.
* @ingroup gnssProcessingStepGroup
* @see GnssProcessingStep */
class GnssProcessingStepPrintResidualStatisticsIsl : public GnssProcessingStepBase
{
public:
  GnssProcessingStepPrintResidualStatisticsIsl(Config &/*config*/) {}
  void process(GnssProcessingStep::State &state) override;
};

/***********************************************/

inline void GnssProcessingStepPrintResidualStatisticsIsl::process(GnssProcessingStep::State &state)
{
  try
  {
    logStatus<<"=== print ISL residual statistics  =========================="<<Log::endl;
    for(UInt idTrans=0; idTrans<state.gnss->transmitters.size(); idTrans++)
    {
      Double   ePe=0, redundancy=0;
      UInt     obsCount=0, outlierCount=0;
      state.residualsStatisticsIsl(idTrans, ePe, redundancy, obsCount, outlierCount);
      Double factor = state.transmitters.at(idTrans).sigmaFactor;
	    Parallel::reduceSum(factor, 0, state.normalEquationInfo.comm);
      if(Parallel::isMaster(state.normalEquationInfo.comm))
        if(obsCount>0)
          logInfo<<"  "<<state.gnss->transmitters.at(idTrans)->name()<<" : ISL   "
                 <<": factor = "    <<factor%"%5.2f"s
                 <<", sigma0 = "    <<Vce::standardDeviation(ePe, redundancy, 2.5/*huber*/, 1.5/*huberPower*/)%"%4.2f"s
                 <<", count = "     <<obsCount%"%5i"s
                 <<", outliers = "  <<outlierCount%"%5i"s<<" ("<<(100.*outlierCount/obsCount)%"%4.2f"s<<" %)"<<Log::endl;
    } // for(idTrans)

  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
