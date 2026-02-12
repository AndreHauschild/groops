/***********************************************/
/**
* @file gnssTransceiver.h
*
* @brief GNSS receiver or transmitter.
*
* @author Torsten Mayer-Guerr
* @date 2021-01-30
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSTRANSCEIVER__
#define __GROOPS_GNSSTRANSCEIVER__

#include "base/gnssType.h"
#include "files/fileGnssSignalBias.h"
#include "files/fileIslSignalBias.h"
#include "files/filePlatform.h"

/** @addtogroup gnssGroup */
/// @{

/***** TYPES ***********************************/

class GnssTransceiver;
typedef std::shared_ptr<GnssTransceiver> GnssTransceiverPtr;

/***** CLASS ***********************************/

/** @brief Abstract class for GNSS receiver or transmitter.
* eg. GPS satellites. */
class GnssTransceiver
{
   Vector      useableEpochs;
   UInt        countUseableEpochs;
   GnssAntennaDefinition::NoPatternFoundAction noPatternFoundAction;

public:
   UInt           id_; // set by Gnss::init()
   Platform       platform;
   GnssSignalBias signalBias;
   IslSignalBias  signalBiasIslTx;
   IslSignalBias  signalBiasIslRx;

public:
  /// Constructor.
  GnssTransceiver(const Platform &platform, GnssAntennaDefinition::NoPatternFoundAction noPatternFoundAction, const Vector &useableEpochs);

  /// Destructor.
  virtual ~GnssTransceiver() {}

  /** @brief name. */
  std::string name() const {return platform.name;}

  /** @brief marker number. */
  std::string markerNumber() const {return platform.markerNumber;}

  /** @brief Is the platform usable at given epoch (or all epochs). */
  Bool useable(UInt idEpoch=NULLINDEX) const {return countUseableEpochs && ((idEpoch == NULLINDEX) || useableEpochs(idEpoch));}

  /** @brief Disable given epoch. */
  virtual void disable(UInt idEpoch, const std::string &reason);

  /** @brief Disable transceiver. */
  virtual void disable(const std::string &reason);

  /** @brief Allowed signal types. Empty if no  definition was provided. */
  std::vector<GnssType> definedTypes(const Time &time) const;

  /** @brief Signal bias corrections.
  * observed range = range + bias. */
  Vector signalBiases(const std::vector<GnssType> &type) const;

  /** @brief Direction dependent corrections.
  * observed range = range (ARPs of transmitter and receiver) + antennaVariations. */
  Vector antennaVariations(const Time &time, Angle azimut, Angle elevation, const std::vector<GnssType> &type) const;

  /** @brief Direction (and other parameters) dependent standard deviation.
  * @a azimuth and @a elevation must be given in the antenna frame (left-handed). */
  Vector accuracy(const Time &time, Angle azimut, Angle elevation, const std::vector<GnssType> &type) const;

  /** @brief ISL terminal bias corrections.
  * observed range = range + bias. */
  Vector signalBiasesIslTx(std::vector<UInt> terminals) const;
  Vector signalBiasesIslRx(std::vector<UInt> terminals) const;

  /** @brief Direction dependent corrections for ISL observations.
  * observed range = range (ARPs of ISL transmit and receive terminal) + islTerminalVariations. */
  Double islTerminalVariations(const Time &time, Angle azimut, Angle elevation) const;


  void save(OutArchive &oa) const;
  void load(InArchive  &ia);
};

/***********************************************/

inline GnssTransceiver::GnssTransceiver(const Platform &platform, GnssAntennaDefinition::NoPatternFoundAction noPatternFoundAction, const Vector &useableEpochs)
  : useableEpochs(useableEpochs), countUseableEpochs(sum(useableEpochs)), noPatternFoundAction(noPatternFoundAction), platform(platform) {}

/***********************************************/

inline void GnssTransceiver::disable(UInt idEpoch, const std::string &reason)
{
  try
  {
    if(useableEpochs(idEpoch))
      countUseableEpochs--;
    useableEpochs(idEpoch) = FALSE;
    if(countUseableEpochs == 0)
      disable(reason);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void GnssTransceiver::disable(const std::string &/*reason*/)
{
  try
  {
    countUseableEpochs = 0;
    useableEpochs.setNull();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline std::vector<GnssType> GnssTransceiver::definedTypes(const Time &time) const
{
  try
  {
    auto receiver = platform.findEquipment<PlatformGnssReceiver>(time);
    if(!receiver || !receiver->receiverDef)
      return std::vector<GnssType>();
    return receiver->receiverDef->types;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Vector GnssTransceiver::signalBiases(const std::vector<GnssType> &types) const
{
  try
  {
    return signalBias.compute(types);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Vector GnssTransceiver::antennaVariations(const Time &time, Angle azimut, Angle elevation, const std::vector<GnssType> &types) const
{
  try
  {
    auto antenna = platform.findEquipment<PlatformGnssAntenna>(time);
    if(!antenna)
      throw(Exception(platform.markerName+"."+platform.markerNumber+": no antenna definition found at "+time.dateTimeStr()));
    if(!antenna->antennaDef)
      throw(Exception("no antenna definition for "+antenna->str()));
    return antenna->antennaDef->antennaVariations(azimut, elevation, types, noPatternFoundAction);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Vector GnssTransceiver::accuracy(const Time &time, Angle azimut, Angle elevation, const std::vector<GnssType> &types) const
{
  try
  {
    auto antenna = platform.findEquipment<PlatformGnssAntenna>(time);
    if(!antenna)
      throw(Exception(platform.markerName+"."+platform.markerNumber+": no antenna accuracy found at "+time.dateTimeStr()));
    if(!antenna->accuracyDef)
      throw(Exception("no accuracy definition for "+antenna->str()));
    return antenna->accuracyDef->antennaVariations(azimut, elevation, types, noPatternFoundAction);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Vector GnssTransceiver::signalBiasesIslTx(std::vector<UInt> terminals) const
{
  try
  {
    return signalBiasIslTx.compute(terminals);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Vector GnssTransceiver::signalBiasesIslRx(std::vector<UInt> terminals) const
{
  try
  {
    return signalBiasIslRx.compute(terminals);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Double GnssTransceiver::islTerminalVariations(const Time &time, Angle azimut, Angle elevation) const
{
  try
  {
    // ISL observation code for antenna corrections
    // TODO: avoid using observation types here!
    // --------------------------------------------
    const std::vector<GnssType> types = { GnssType("C1C***") };
    auto terminal = platform.findEquipment<PlatformIslTerminal>(time);
    if(!terminal)
      throw(Exception(platform.markerName+"."+platform.markerNumber+": no ISL terminal definition found at "+time.dateTimeStr()));
    if(!terminal->antennaDef)
      throw(Exception("no ISL terminal definition for "+terminal->str()));
    return terminal->antennaDef->antennaVariations(azimut, elevation, types, noPatternFoundAction).at(0);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void GnssTransceiver::save(OutArchive &oa) const
{
  oa<<nameValue("useableEpochs",      useableEpochs);
  oa<<nameValue("countUseableEpochs", countUseableEpochs);
  oa<<nameValue("signalBias",         signalBias);
  oa<<nameValue("signalBiasIsl",      signalBiasIslRx);
  oa<<nameValue("signalBiasIsl",      signalBiasIslTx);
}

/***********************************************/

inline void GnssTransceiver::load(InArchive  &ia)
{
  ia>>nameValue("useableEpochs",      useableEpochs);
  ia>>nameValue("countUseableEpochs", countUseableEpochs);
  ia>>nameValue("signalBias",         signalBias);
  ia>>nameValue("signalBiasIsl",      signalBiasIslRx);
  ia>>nameValue("signalBiasIsl",      signalBiasIslTx);
}

/***********************************************/

/// @}

#endif /* __GROOPS___ */
