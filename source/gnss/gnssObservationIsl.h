/***********************************************/
/**
* @file gnssObservationIsl.h
*
* @brief Inter-satellite links.
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

/** @brief ISL observations.
* Between one receiver and one transmitter terminal at one epoch. */
class GnssObservationIsl
{
public:
  Double   observation; ///< original observations
  Double   residual;    ///< estimated postfit residuals
  Double   redundancy;  ///< partial redundancies of the least squares adjustment
  Double   sigma0;      ///< expected (a-priori) accuracies
  Double   sigma;       ///< modified accuracies (downweighted outliers)
  UInt     terminalRecv, terminalSend;

  GnssObservationIsl() : observation(0), residual(0), redundancy(0), sigma0(0), sigma(0), terminalRecv(NULLINDEX), terminalSend(NULLINDEX) {}

  Bool init(const GnssTransmitter &receiver, const GnssTransmitter &transmitter, const Time &time, UInt idEpoch, UInt terminalRecv_, UInt terminalSend_);

  void setHomogenizedResiduals(Double residual, Double redundancy);
};

/***** CLASS ***********************************/

/** @brief Reduced ISL observations (obs - computed) and design matrix.
* Between one receiver and one transmitter terminal at one epoch. */
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
  Matrix A;      ///< columns: dl/dx, dl/dy, dl/dz, dl/dClock, unit matrix

  // approximate values (Taylor point)
  Time     timeRecv, timeTrans;
  Vector3d posRecv,  posTrans;
  Vector3d velocityRecv, velocityTrans;
  Angle    azimutRecvAnt,   elevationRecvAnt;
  Angle    azimutTrans,     elevationTrans;

  GnssObservationEquationIsl() : idEpoch(NULLINDEX), receiver(nullptr), transmitter(nullptr), terminalRecv(NULLINDEX), terminalSend(NULLINDEX) {}

  GnssObservationEquationIsl(const GnssObservationIsl &observation, const GnssTransmitter &receiver, const GnssTransmitter &transmitter,
                             const std::function<void(GnssObservationEquationIsl &eqn)> &reduceModels,
                             Time time, UInt idEpoch, Bool homogenize)
          {compute(observation, receiver, transmitter, reduceModels, time, idEpoch, homogenize);}

  void compute(const GnssObservationIsl &observation, const GnssTransmitter &receiver, const GnssTransmitter &transmitter,
               const std::function<void(GnssObservationEquationIsl &eqn)> &reduceModels,
               Time time, UInt idEpoch, Bool homogenize);
};

/***********************************************/

/// @}

#endif /* __GROOPS___ */
