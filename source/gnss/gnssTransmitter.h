/***********************************************/
/**
* @file gnssTransmitter.h
*
* @brief GNSS transmitter.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2013-06-28
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSTRANSMITTER__
#define __GROOPS_GNSSTRANSMITTER__

#include "base/polynomial.h"
#include "base/gnssType.h"
#include "files/fileInstrument.h"
#include "gnss/gnssObservationIsl.h"
#include "gnss/gnssTransceiver.h"
#include "classes/noiseGenerator/noiseGenerator.h"

/** @addtogroup gnssGroup */
/// @{

/***** TYPES ***********************************/

class GnssTransmitter;
typedef std::shared_ptr<GnssTransmitter> GnssTransmitterPtr;

/***** CLASS ***********************************/

/** @brief Abstract class for GNSS transmitter.
* eg. GPS satellites. */
class GnssTransmitter : public GnssTransceiver
{
  GnssType                 type; // system + PRN
  Polynomial               polynomial;
  std::vector<Double>      clk;
  std::vector<Vector3d>    offset;      // between CoM and ARF in SRF
  std::vector<Vector3d>    offsetIsl;   // between CoM and ISL RF in SRF
  std::vector<Transform3d> crf2srf, srf2arf;
  std::vector<Transform3d> srf2arfIsl;

public:
  Bool              isMyRank_;
  std::vector<Time> timesPosVel;
  Matrix            pos, vel; // CoM in CRF (epoch times (x,y,z))
  std::string       disableReason;

  GnssTransmitter(GnssType prn, const Platform &platform,
                  GnssAntennaDefinition::NoPatternFoundAction noPatternFoundAction,
                  const Vector &useableEpochs, const std::vector<Double> &clock, const std::vector<Vector3d> &offset,
                  const std::vector<Transform3d> &crf2srf, const std::vector<Transform3d> &srf2arf,
                  const std::vector<Vector3d> &offsetIsl, const std::vector<Transform3d> &srf2arfIsl,
                  const std::vector<Time> &timesPosVel, const_MatrixSliceRef position, const_MatrixSliceRef velocity, UInt interpolationDegree)
  : GnssTransceiver(platform, noPatternFoundAction, useableEpochs),
    type(prn), polynomial(timesPosVel, interpolationDegree, TRUE/*throwException*/, FALSE/*leastSquares*/, -(interpolationDegree+1.1), -1.1, 1e-7),
    clk(clock), offset(offset), crf2srf(crf2srf), srf2arf(srf2arf), offsetIsl(offsetIsl), srf2arfIsl(srf2arfIsl), timesPosVel(timesPosVel), pos(position), vel(velocity),
    isMyRank_(FALSE) {}

  /// Destructor.
  virtual ~GnssTransmitter();

  /** @brief Identify number in the GNSS system. */
  UInt idTrans() const {return id_;}

  /** @brief Disable given epoch (or all epochs). */
  void disable(UInt idEpoch, const std::string &reason) override;

  /** @brief Disable transmitter completely. */
  void disable(const std::string &reason) override;

  /** @brief Returns true if transmitter is assigned to current node. */
  Bool isMyRank() const {return isMyRank_;}

  /** @brief PRN number of satellite.
  *  = prn + GnssType::SYSTEM. */
  GnssType PRN() const {return type;}

  /** @brief Clock error.
  * error = clock time - system time [s] */
  Double clockError(UInt idEpoch) const {return clk.at(idEpoch);}

  /** @brief set clock error.
  * error = observed clock time - system time [s] */
  void updateClockError(UInt idEpoch, Double deltaClock) {clk.at(idEpoch) += deltaClock;}

  /** @brief center of mass in celestial reference frame (CRF). */
  Vector3d positionCoM(const Time &time) const;

  /** @brief antenna reference point in celestial reference frame (CRF). */
  Vector3d position(UInt idEpoch, const Time &time) const {return positionCoM(time) + crf2srf.at(idEpoch).inverseTransform(offset.at(idEpoch));}

  /** @brief terminal reference point in celestial reference frame (CRF). */
  Vector3d positionIsl(UInt idEpoch, const Time &time) const {return positionCoM(time) + crf2srf.at(idEpoch).inverseTransform(offsetIsl.at(idEpoch));}

  /** @brief velocity in CRF [m/s]. */
  Vector3d velocity(const Time &time) const;

  /** @brief Rotation from celestial reference frame (CRF) to left-handed antenna system. */
  Transform3d celestial2antennaFrame(UInt idEpoch, const Time &/*time*/) const {return srf2arf.at(idEpoch) * crf2srf.at(idEpoch);}

  /** @brief Rotation from celestial reference frame (CRF) to left-handed ISL terminal system. */
  Transform3d celestial2islTerminalFrame(UInt idEpoch, const Time &/*time*/) const {return srf2arfIsl.at(idEpoch) * crf2srf.at(idEpoch);}

  // Inter satellite links (ISL)
  // ---------------------------
private:
  std::vector<std::vector<GnssObservationIsl*>> observations_; // observations at receiver (for each epoch, for each transmitter)

public:
  /** @brief ISL observation between receiver and transmitter at one epoch. */
  GnssObservationIsl *observationIsl(UInt idTrans, UInt idEpoch) const;

  /** @brief Delete ISL observation. */
  void deleteObservationIsl(UInt idTrans, UInt idEpoch);

  /** @brief Max. observed epoch id+1. */
  UInt idEpochSize() const {return observations_.size();}

  /** @brief Max. observed transmitter id+1 at @a idEpoch. */
  UInt idTransmitterSize(UInt idEpoch) const {return (idEpoch < observations_.size()) ? observations_.at(idEpoch).size() : 0;}

  void readObservationsIsl(const FileName &fileName, const std::vector<GnssTransmitterPtr> &transmitters, const std::vector<Time> &times, const Time &timeMargin);
  void simulateObservationsIsl(NoiseGeneratorPtr noiseObs, const std::vector<GnssTransmitterPtr> &transmitters,
                               const std::vector<Time> &times, const GnssReceiverArc   &scheduleIsl,
                               const std::function<void(GnssObservationEquationIsl &eqn)> &reduceModels);
};

/***********************************************/

inline Vector3d GnssTransmitter::positionCoM(const Time &time) const
{
  try
  {
    return Vector3d(polynomial.interpolate({time}, pos));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Vector3d GnssTransmitter::velocity(const Time &time) const
{
  try
  {
    return Vector3d(polynomial.interpolate({time}, vel));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

/// @}

#endif /* __GROOPS___ */
