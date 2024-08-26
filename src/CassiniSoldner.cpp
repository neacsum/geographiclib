/**
 * \file CassiniSoldner.cpp
 * \brief Implementation for GeographicLib::CassiniSoldner class
 *
 * Copyright (c) Charles Karney (2009-2022) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#include <GeographicLib/CassiniSoldner.hpp>

namespace GeographicLib {

  using namespace std;

  CassiniSoldner::CassiniSoldner(const Geodesic& earth)
    : earth_(earth) {}

  CassiniSoldner::CassiniSoldner(real lat0, real lon0, const Geodesic& earth)
    : earth_(earth)
  { Reset(lat0, lon0); }

  void CassiniSoldner::Reset(real lat0, real lon0) {
    meridian_ = earth_.Line(lat0, lon0, real(0),
                            Geodesic::LATITUDE | Geodesic::LONGITUDE |
                            Geodesic::DISTANCE | Geodesic::DISTANCE_IN |
                            Geodesic::AZIMUTH);
    real f = earth_.Flattening();
    Math::sincosd(LatitudeOrigin(), sbet0_, cbet0_);
    sbet0_ *= (1 - f);
    Math::norm(sbet0_, cbet0_);
  }

  void CassiniSoldner::Forward(real lat, real lon, real& x, real& y,
                               real& azi, real& rk) const {
    if (!Init())
      return;
    real dlon = Math::AngDiff(LongitudeOrigin(), lon);
    real sig12, s12, azi1, azi2;
    sig12 = earth_.Inverse(lat, -fabs(dlon), lat, fabs(dlon), s12, azi1, azi2);
    sig12 *= real(0.5);
    s12 *= real(0.5);
    if (s12 == 0) {
      real da = Math::AngDiff(azi1, azi2)/2;
      if (fabs(dlon) <= 90) {
        azi1 = 90 - da;
        azi2 = 90 + da;
      } else {
        azi1 = -90 - da;
        azi2 = -90 + da;
      }
    }
    if (signbit(dlon)) {
      azi2 = azi1;
      s12 = -s12;
      sig12 = -sig12;
    }
    x = s12;
    azi = Math::AngNormalize(azi2);
    GeodesicLine perp(earth_.Line(lat, dlon, azi, Geodesic::GEODESICSCALE));
    real t;
    perp.GenPosition(true, -sig12,
                     Geodesic::GEODESICSCALE,
                     t, t, t, t, t, t, rk, t);

    real salp0, calp0;
    Math::sincosd(perp.EquatorialAzimuth(), salp0, calp0);
    real
      sbet1 = lat >=0 ? calp0 : -calp0,
      cbet1 = fabs(dlon) <= 90 ? fabs(salp0) : -fabs(salp0),
      sbet01 = sbet1 * cbet0_ - cbet1 * sbet0_,
      cbet01 = cbet1 * cbet0_ + sbet1 * sbet0_,
      sig01 = atan2(sbet01, cbet01) / Math::degree();
    meridian_.GenPosition(true, sig01,
                          Geodesic::DISTANCE,
                          t, t, t, y, t, t, t, t);
  }

  void CassiniSoldner::Reverse(real x, real y, real& lat, real& lon,
                               real& azi, real& rk) const {
    if (!Init())
      return;
    real lat1, lon1;
    real azi0, t;
    meridian_.Position(y, lat1, lon1, azi0);
    earth_.Direct(lat1, lon1, azi0 + 90, x, lat, lon, azi, rk, t);
  }

} // namespace GeographicLib
