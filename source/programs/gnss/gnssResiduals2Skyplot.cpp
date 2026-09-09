/***********************************************/
/**
* @file gnssResiduals2Skyplot.cpp
*
* @brief Convert residuals into griddedData format for plotting.
*
* @author Torsten Mayer-Guerr
* @date 2013-07-12
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Write GNSS residuals together with azimuth and elevation to be plotted with \program{PlotMap}.
Azimuth and elevation are written as ellipsoidal longitude and latitude in a \file{griddedData file}{griddedData}.
The choosen ellipsoid parameters \config{R} and \config{inverseFlattening} are arbitrary but should be the same
as in \program{PlotMap}. If with \configClass{typeTransmitter}{gnssType} (e.g. '\verb|***G18|')
a single transmitter is selected the azimuth and elevation are computed from the transmitter point of view.

For each GNSS \configClass{type}{gnssType} an extra data column is created.

A \file{GNSS residual file}{instrument} includes additional information
besides the residuals, which can also be selected with \configClass{type}{gnssType}
\begin{itemize}
\item \verb|A1*|, \verb|E1*|: azimuth and elevation at receiver
\item \verb|A2*|, \verb|E2*|: azimuth and elevation at transmitter
\item \verb|I**|: Estimated slant total electron content (STEC)
\end{itemize}

Furthermore these files may include for each residual \configClass{type}{gnssType}
information about the redundancy and the accuracy $\sigma$ from the least squares adjustment.
The 3 values (residuals, redundancy, $\sigma$) are coded with the same type.
To get access to all values the corresponding type must be repeated in \configClass{type}{gnssType}.

\fig{!hb}{0.5}{gnssResiduals2Skyplot}{fig:gnssResiduals2Skyplot}{GPS C2W residuals of GRAZ station at 2012-01-01}
)";

/***********************************************/

#include "programs/program.h"
#include "inputOutput/system.h"
#include "files/fileInstrument.h"
#include "files/fileGriddedData.h"
#include "misc/miscGriddedData.h"

/***** CLASS ***********************************/

/** @brief Convert residuals into griddedData format for plotting.
* @ingroup programsGroup */
class GnssResiduals2Skyplot
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GnssResiduals2Skyplot, SINGLEPROCESS, "Convert residuals into griddedData format for plotting", Gnss, Grid)
GROOPS_RENAMED_PROGRAM(GnssResiduals2GriddedData, GnssResiduals2Skyplot, date2time(2019, 9, 9))

/***********************************************/

void GnssResiduals2Skyplot::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName              fileNameGriddedData;
    std::vector<FileName> fileNameResiduals;
    std::vector<GnssType> types;
    GnssType              typeTransmitter;
    Double                a, f;

    readConfig(config, "outputfileGriddedData", fileNameGriddedData, Config::MUSTSET,  "", "");
    readConfig(config, "type",                  types,               Config::MUSTSET,  "", "");
    readConfig(config, "typeTransmitter",       typeTransmitter,     Config::OPTIONAL, "", "choose transmitter view, e.g. '***G18'");
    readConfig(config, "inputfileResiduals",    fileNameResiduals,   Config::MUSTSET,  "", "GNSS receiver residuals");
    readConfig(config, "R",                     a,                   Config::DEFAULT,  STRING_DEFAULT_GRS80_a, "reference radius for ellipsoidal coordinates");
    readConfig(config, "inverseFlattening",     f,                   Config::DEFAULT,  STRING_DEFAULT_GRS80_f, "reference flattening for ellipsoidal coordinates");
    if(isCreateSchema(config)) return;

    // ============================

    Ellipsoid                        ellipsoid(a, f);
    std::vector<Vector3d>            points;
    std::vector<std::vector<Double>> values(types.size());
    for(const FileName &fileName : fileNameResiduals)
    {
      logStatus<<"read GNSS residuals <"<<fileName<<">"<<Log::endl;
      if(!System::exists(fileName))
      {
        logWarning<<"file not exist -> continue"<<Log::endl;
        continue;
      }

      for(auto &epoch : GnssReceiverArc(InstrumentFile::read(fileName)))
      {
        UInt idObs = 0;
        for(GnssType satType : epoch.satellite)
        {
          Bool                  found = FALSE;
          Double                azimuth=NAN_EXPR, elevation=NAN_EXPR;
          std::vector<GnssType> typesTmp = types;
          std::vector<Double>   valuesPerPoint(types.size(), NAN_EXPR);

          // find type for the satellite system, loop over all obs for this satellite
          UInt idType = std::distance(epoch.obsType.begin(), std::find(epoch.obsType.begin(), epoch.obsType.end(), satType));
          for(; (idType<epoch.obsType.size()) && (epoch.obsType.at(idType)==satType); idType++, idObs++)
          {
            const GnssType type  = epoch.obsType.at(idType) + satType;
            const Double   value = epoch.observation.at(idObs);

            if(!typeTransmitter.hasWildcard(GnssType::PRN)) // isTransmitter
            {
              if(type == (GnssType::AZIMUT    + GnssType::L2)) azimuth   = value;
              if(type == (GnssType::ELEVATION + GnssType::L2)) elevation = value;
            }
            else
            {
              if(type == (GnssType::AZIMUT    + GnssType::L1)) azimuth   = value;
              if(type == (GnssType::ELEVATION + GnssType::L1)) elevation = value;
            }

            UInt idx;
            if(type.isInList(typesTmp, idx) && (type == typeTransmitter) && value && !std::isnan(value))
            {
              valuesPerPoint.at(idx) = value;
              typesTmp.at(idx) = GnssType(static_cast<UInt64>(-1)); // disable types already found
              found = TRUE;
            }
          } // for(idObs)

          if(found)
          {
            if(std::isnan(azimuth) || std::isnan(elevation))
              throw(Exception("file must contain azimuth and elevation."));
            points.push_back(ellipsoid(Angle(azimuth), Angle(elevation), 0));
            for(UInt i=0; i<values.size(); i++)
              values.at(i).push_back(valuesPerPoint.at(i));
          }
        } // for(satType)
      } // for(epoch)
    } // for(idFile)

    // ============================

    logStatus<<"save values to file <"<<fileNameGriddedData<<">"<<Log::endl;
    GriddedData griddedData(ellipsoid, points, std::vector<Double>(points.size(), 1), values);
    writeFileGriddedData(fileNameGriddedData, griddedData);
    MiscGriddedData::printStatistics(griddedData);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/
