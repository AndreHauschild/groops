/***********************************************/
/**
* @file troposphereGlobalMapping.h
*
* @brief Global Mapping Function (GMF).
* @see Troposphere
*
* @author Andreas Kvas
* @date 2024-06-23
*/
/***********************************************/

#ifndef __GROOPS_TROPOSPHEREGLOBALMAPPING__
#define __GROOPS_TROPOSPHEREGLOBALMAPPING__

// Latex documentation
#ifdef DOCSTRING_Troposphere
static const char *docstringTroposphereGlobalMapping = R"(
\subsection{GlobalMapping}\label{troposphereType:globalMapping}

Global Mapping Function (GMF) tropospheric mapping function.
Reference: Boehm, J., A.E. Niell, P. Tregoning, H. Schuh (2006),
Global Mapping Functions (GMF): A new empirical mapping function based on numerical weather model data,
Geoph. Res. Letters, Vol. 33, L07304, doi:10.1029/2005GL025545.

Apriori hydrostatic and non-hydrostatic delays are computed from approximate meteorological
data provided via \configFile{inputfileGpt}{griddedData}.

)";
#endif

/***********************************************/

#include "troposphere.h"

/***** CLASS ***********************************/

/** @brief Global Mapping Function (GMF).
* @ingroup troposphereGroup
* @see Troposphere */
class TroposphereGlobalMapping : public Troposphere
{
  FileName              fileNameGpt;
  Double                a_ht, b_ht, c_ht;
  Vector                ahm, aha, awm, awa;
  Vector                phh, c11h, c10h;
  Vector                cos_lat, dhgt;

public:
  TroposphereGlobalMapping(Config &config);
 ~TroposphereGlobalMapping() {}

  void   init(const std::vector<Vector3d> &stationPositions) override;
  Double slantDelay(const Time &time, UInt stationId, Angle azimuth, Angle elevation) const override;
  Double mappingFunctionHydrostatic(const Time &time, UInt stationId, Angle azimuth, Angle elevation) const override;
  Double mappingFunctionWet(const Time &time, UInt stationId, Angle azimuth, Angle elevation) const override;
  void   mappingFunctionGradient(const Time &time, UInt stationId, Angle azimuth, Angle elevation, Double &dx, Double &dy) const override;
  void   getAprioriValues(const Time&time, UInt stationId, Double &zenithDryDelay, Double &zenithWetDelay, Double &gradientDryNorth,
                          Double &gradientWetNorth, Double &gradientDryEast, Double &gradientWetEast, Double &aDry, Double &aWet) const override;
};

/***********************************************/

#endif
