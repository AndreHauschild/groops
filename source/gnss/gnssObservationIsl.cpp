/***********************************************/
/**
 * @file gnssObservationIsl.cpp
 *
 * @brief Intersatellite links.
 *
 * @author Torsten Mayer-Guerr
 * @date 2024-01-30
 *
 */
/***********************************************/

#include "base/import.h"
#include "gnss/gnssTransmitter.h"
#include "gnss/gnssReceiver.h"
#include "gnss/gnssObservationIsl.h"

/***********************************************/

static void positionVelocityTime(const GnssTransmitter &receiver, const GnssTransmitter &transmitter, const Time &time, UInt idEpoch,
                                 Time &timeRecv,  Vector3d &posRecv,  Vector3d &velRecv,  Angle &azimutRecv,  Angle &elevationRecv,
                                 Time &timeTrans, Vector3d &posTrans, Vector3d &velTrans, Angle &azimutTrans, Angle &elevationTrans,
                                 Vector3d &k, Vector3d &kRecv, Vector3d &kTrans)
{
  try
  {
    // receiver in celestial reference frame
    timeRecv = time - seconds2time(receiver.clockError(idEpoch));
    posRecv  = receiver.position(idEpoch, timeRecv);
    velRecv  = receiver.velocity(timeRecv);

    // transmitter position and time
    posTrans = transmitter.position(idEpoch, timeRecv);
    Vector3d posOld;
    for(UInt i=0; (i<10) && ((posTrans-posOld).r() > 0.0001); i++) // iteration
    {
      timeTrans = timeRecv - seconds2time((posTrans-posRecv).r()/LIGHT_VELOCITY);
      posOld    = posTrans;
      posTrans  = transmitter.position(idEpoch, timeTrans);
    }
    velTrans = transmitter.velocity(timeTrans);

    // line of sight from transmitter to receiver
    k              = normalize(posRecv - posTrans);
    kRecv          = receiver.celestial2antennaFrame(idEpoch, timeRecv).transform(-k); // line of sight in receiver antenna system
    azimutRecv     = kRecv.lambda();
    elevationRecv  = kRecv.phi();
    kTrans         = transmitter.celestial2antennaFrame(idEpoch, timeTrans).transform(k);
    azimutTrans    = kTrans.lambda();
    elevationTrans = kTrans.phi();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssObservationIsl::setDecorrelatedResiduals(Double residual, Double redundancy)
{
  this->residual   = residual * sigma;
  this->redundancy = redundancy;
}

/***********************************************/
/***********************************************/

void GnssObservationEquationIsl::compute(const GnssObservationIsl &observation, const GnssTransmitter &receiver_, const GnssTransmitter &transmitter_,
                                         const std::function<void(GnssObservationEquationIsl &eqn)> &reduceModels,
                                         UInt idEpoch_, Bool decorrelate)
{
  try
  {
    const UInt obsCount = 1;
    std::vector<GnssType> types;

    idEpoch      = idEpoch_;
    receiver     = &receiver_;
    transmitter  = &transmitter_;
    terminalRecv = observation.terminalRecv;
    terminalSend = observation.terminalSend;
    types.resize(obsCount);
    l            = Vector(obsCount);
    sigma        = Vector(obsCount);
    sigma0       = Vector(obsCount);

    types.at(0) = GnssType("C1C***"); // TODO: avoid hard-coding this here!
    l(0)        = observation.observation;
    sigma(0)    = observation.sigma;
    sigma0(0)   = observation.sigma0;

    // position, time of transmitter & receiver
    // ----------------------------------------
    Vector3d k, kRecvAnt, kTrans;
    positionVelocityTime(receiver_, transmitter_, observation.time, idEpoch_,
                         timeRecv,  posRecv,  velocityRecv,  azimutRecvAnt, elevationRecvAnt,
                         timeTrans, posTrans, velocityTrans, azimutTrans,   elevationTrans,
                         k, kRecvAnt, kTrans);
    const Double rDotTrans = inner(k, velocityTrans)/LIGHT_VELOCITY;
    const Double rDotRecv  = inner(k, velocityRecv) /LIGHT_VELOCITY;

    // Corrected range
    // ---------------
    const Double r1  = posTrans.r();
    const Double r2  = posRecv.r();
    Double       r12 = (posRecv - posTrans).r();
    r12 += 2*DEFAULT_GM/pow(LIGHT_VELOCITY,2)*log((r1+r2+r12)/(r1+r2-r12)); // curved space-time
    r12 += 2*inner(posTrans, velocityTrans)/LIGHT_VELOCITY;                 // relativistic clock correction
    r12 -= 2*inner(posRecv, velocityRecv)/LIGHT_VELOCITY;                   // relativistic clock correction
    r12 -= LIGHT_VELOCITY * transmitter->clockError(idEpoch);
    r12 += LIGHT_VELOCITY * receiver->clockError(idEpoch);

    // approximate range
    // -----------------
    for(UInt i=0; i<obsCount; i++)
      l(i) -= r12;

    // design matrix
    // -------------
    A = Matrix(obsCount, 9);
    for(UInt i=0; i<obsCount; i++)
    {
      A(i, idxPosRecv+0)  = k.x() * (1+rDotTrans);   // receiver coord x
      A(i, idxPosRecv+1)  = k.y() * (1+rDotTrans);   // receiver coord y
      A(i, idxPosRecv+2)  = k.z() * (1+rDotTrans);   // receiver coord z
      A(i, idxClockRecv)  = 1.0+rDotTrans-rDotRecv;  // receiver clock correction
      A(i, idxPosTrans+0) = k.x() * (-1+rDotTrans);  // transmitter coord x
      A(i, idxPosTrans+1) = k.y() * (-1+rDotTrans);  // transmitter coord y
      A(i, idxPosTrans+2) = k.z() * (-1+rDotTrans);  // transmitter coord z
      A(i, idxClockTrans) = -1.0;                    // transmitter clock correction
      A(i, idxRange)      =  1.0;                    // range correction (???, ...)
    }  // for(i=0..obsCount)

    // antenna correction and other corrections
    // ----------------------------------------
    l -= receiver->antennaVariations(timeRecv, azimutRecvAnt,  elevationRecvAnt,  types);
    l -= receiver->signalBiasesIslRx(types);
    l -= transmitter->antennaVariations(timeTrans, azimutTrans, elevationTrans, types);
    l -= transmitter->signalBiasesIslTx(types);
    if(reduceModels)
      reduceModels(*this);

    // Decorrelate
    // -----------
    if(decorrelate)
      for(UInt i=0; i<obsCount; i++)
      {
        if(l.size()) l.row(i) *= 1/sigma(i);
        if(A.size()) A.row(i) *= 1/sigma(i);
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
