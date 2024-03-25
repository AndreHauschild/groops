/***********************************************/
/**
* @file gnssObservationIsl.h
*
* @brief Intersatellite links.
*
* @author Torsten Mayer-Guerr
* @date 2024-01-30
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSOBSERVATIONISL__
#define __GROOPS_GNSSOBSERVATIONISL__

#include "base/gnssType.h"

/** @addtogroup gnssGroup */
/// @{

/***** TYPES ***********************************/

class GnssTransmitter;

/***** CLASS ***********************************/

/** @brief Observations.
* Between one receiver and one transmitter at one epoch. */
class GnssObservationIsl
{
public:
  Time     time;
  Double   observation; ///< original observations
  Double   residual;    ///< estimated postfit residuals
  Double   redundancy;  ///< partial redundancies of the least squares adjustment
  Double   sigma0;      ///< expected (apriori) accuracies
  Double   sigma;       ///< modified accuracies (downweighted outliers)
  UInt     terminalRecv, terminalSend;

  GnssObservationIsl() : observation(0), residual(0), redundancy(0), sigma0(1), sigma(1), terminalRecv(0), terminalSend(0) {}

  void setDecorrelatedResiduals(Double residual, Double redundancy);
};

/***** CLASS ***********************************/

/** @brief Reduced observations (obs - computed) and design matrix.
* Between one receiver and one transmitter at one epoch. */
class GnssObservationEquationIsl
{
public:
  enum {idxPosRecv    = 0, // x,y,z (CRF)
        idxClockRecv  = 3,
        idxPosTrans   = 4, // x,y,z (CRF)
        idxClockTrans = 7,
        idxRange      = 8};

  UInt  idEpoch;
  const GnssTransmitter *receiver;
  const GnssTransmitter *transmitter;
  UInt  terminalRecv, terminalSend;

  // weighted observations (with 1/sigma)
  Vector l;            ///< weighted reduced observations
  Vector sigma;
  Vector sigma0;

  // design matrix
  Matrix A;      ///< columns: dl/dx, dl/dy, dl/dz, dl/dClock, unit matrix, transformation matrix (typesTransmitter->types)

  // approximate values (Taylor point)
  Time     timeRecv, timeTrans;
  Vector3d posRecv,  posTrans;
  Vector3d velocityRecv, velocityTrans;
  Angle    azimutRecvAnt,   elevationRecvAnt;
  Angle    azimutTrans,     elevationTrans;

  GnssObservationEquationIsl() : idEpoch(NULLINDEX), receiver(nullptr), transmitter(nullptr) {}

  GnssObservationEquationIsl(const GnssObservationIsl &observation, const GnssTransmitter &receiver, const GnssTransmitter &transmitter,
                             const std::function<void(GnssObservationEquationIsl &eqn)> &reduceModels,
                             UInt idEpoch, Bool decorrelate)
          {compute(observation, receiver, transmitter, reduceModels, idEpoch, decorrelate);}

  void compute(const GnssObservationIsl &observation, const GnssTransmitter &receiver, const GnssTransmitter &transmitter,
               const std::function<void(GnssObservationEquationIsl &eqn)> &reduceModels,
               UInt idEpoch, Bool decorrelate);
};

/***********************************************/

/// @}

#endif /* __GROOPS___ */
