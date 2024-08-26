/**
 * \file PolygonArea.cpp
 * \brief Implementation for GeographicLib::PolygonAreaT class
 *
 * Copyright (c) Charles Karney (2010-2023) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#include <GeographicLib/PolygonArea.hpp>

namespace GeographicLib {

  using namespace std;

  template<class GeodType>
  int PolygonAreaT<GeodType>::transit(real lon1, real lon2) {
    // Return 1 or -1 if crossing prime meridian in east or west direction.
    // Otherwise return zero.  longitude = +/-0 considered to be positive.
    // This is (should be?) compatible with transitdirect which computes
    // exactly the parity of
    //   int(floor((lon1 + lon12) / 360)) - int(floor(lon1 / 360)))
    real lon12 = Math::AngDiff(lon1, lon2);
    lon1 = Math::AngNormalize(lon1);
    lon2 = Math::AngNormalize(lon2);
    // N.B. lon12 == 0 gives cross = 0
    return
      // edge case lon1 = 180, lon2 = 360->0, lon12 = 180 to give 1
      lon12 > 0 && ((lon1 < 0 && lon2 >= 0) ||
                    // lon12 > 0 && lon1 > 0 && lon2 == 0 implies lon1 == 180
                    (lon1 > 0 && lon2 == 0)) ? 1 :
      // non edge case lon1 = -180, lon2 = -360->-0, lon12 = -180
      (lon12 < 0 && lon1 >= 0 && lon2 < 0 ? -1 : 0);
    // This was the old method (treating +/- 0 as negative).  However, with the
    // new scheme for handling longitude differences this fails on:
    // lon1 = -180, lon2 = -360->-0, lon12 = -180 gives 0 not -1.
    //    return
    //      lon1 <= 0 && lon2 > 0 && lon12 > 0 ? 1 :
    //      (lon2 <= 0 && lon1 > 0 && lon12 < 0 ? -1 : 0);
  }

  // an alternate version of transit to deal with longitudes in the direct
  // problem.
  template<class GeodType>
  int PolygonAreaT<GeodType>::transitdirect(real lon1, real lon2) {
    // Compute exactly the parity of
    //   int(floor(lon2 / 360)) - int(floor(lon1 / 360))
    // C++ C remainder -> [-360, 360]
    // Java % -> (-720, 720) switch to IEEEremainder -> [-360, 360]
    // JS % -> (-720, 720)
    // Python fmod -> (-720, 720) swith to Math.remainder
    // Fortran, Octave skip
    // If mod function gives result in [-360, 360]
    // [0, 360) -> 0; [-360, 0) or 360 -> 1
    // If mod function gives result in (-720, 720)
    // [0, 360) or [-inf, -360) -> 0; [-360, 0) or [360, inf) -> 1
    lon1 = remainder(lon1, real(2 * 360));
    lon2 = remainder(lon2, real(2 * 360));
    return ( (lon2 >= 0 && lon2 < 360 ? 0 : 1) -
             (lon1 >= 0 && lon1 < 360 ? 0 : 1) );
  }

  template<class GeodType>
  void PolygonAreaT<GeodType>::AddPoint(real lat, real lon) {
    if (num_ == 0) {
      lat0_ = lat1_ = lat;
      lon0_ = lon1_ = lon;
    } else {
      real s12, S12, t;
      earth_.GenInverse(lat1_, lon1_, lat, lon, mask_,
                        s12, t, t, t, t, t, S12);
      perimetersum_ += s12;
      if (!polyline_) {
        areasum_ += S12;
        crossings_ += transit(lon1_, lon);
      }
      lat1_ = lat; lon1_ = lon;
    }
    ++num_;
  }

  template<class GeodType>
  void PolygonAreaT<GeodType>::AddEdge(real azi, real s) {
    if (num_) {                 // Do nothing if num_ is zero
      real lat, lon, S12, t;
      earth_.GenDirect(lat1_, lon1_, azi, false, s, mask_,
                       lat, lon, t, t, t, t, t, S12);
      perimetersum_ += s;
      if (!polyline_) {
        areasum_ += S12;
        crossings_ += transitdirect(lon1_, lon);
      }
      lat1_ = lat; lon1_ = lon;
      ++num_;
    }
  }

  template<class GeodType>
  unsigned PolygonAreaT<GeodType>::Compute(bool reverse, bool sign,
                                           real& perimeter, real& area) const
  {
    real s12, S12, t;
    if (num_ < 2) {
      perimeter = 0;
      if (!polyline_)
        area = 0;
      return num_;
    }
    if (polyline_) {
      perimeter = perimetersum_();
      return num_;
    }
    earth_.GenInverse(lat1_, lon1_, lat0_, lon0_, mask_,
                      s12, t, t, t, t, t, S12);
    perimeter = perimetersum_(s12);
    Accumulator<> tempsum(areasum_);
    tempsum += S12;
    int crossings = crossings_ + transit(lon1_, lon0_);
    AreaReduce(tempsum, crossings, reverse, sign);
    area = real(0) + tempsum();
    return num_;
  }

  template<class GeodType>
  unsigned PolygonAreaT<GeodType>::TestPoint(real lat, real lon,
                                             bool reverse, bool sign,
                                             real& perimeter, real& area) const
  {
    if (num_ == 0) {
      perimeter = 0;
      if (!polyline_)
        area = 0;
      return 1;
    }
    perimeter = perimetersum_();
    real tempsum = polyline_ ? 0 : areasum_();
    int crossings = crossings_;
    unsigned num = num_ + 1;
    for (int i = 0; i < (polyline_ ? 1 : 2); ++i) {
      real s12, S12, t;
      earth_.GenInverse(i == 0 ? lat1_ : lat, i == 0 ? lon1_ : lon,
                        i != 0 ? lat0_ : lat, i != 0 ? lon0_ : lon,
                        mask_, s12, t, t, t, t, t, S12);
      perimeter += s12;
      if (!polyline_) {
        tempsum += S12;
        crossings += transit(i == 0 ? lon1_ : lon,
                             i != 0 ? lon0_ : lon);
      }
    }

    if (polyline_)
      return num;

    AreaReduce(tempsum, crossings, reverse, sign);
    area = real(0) + tempsum;
    return num;
  }

  template<class GeodType>
  unsigned PolygonAreaT<GeodType>::TestEdge(real azi, real s,
                                            bool reverse, bool sign,
                                            real& perimeter, real& area) const
  {
    if (num_ == 0) {            // we don't have a starting point!
      perimeter = Math::NaN();
      if (!polyline_)
        area = Math::NaN();
      return 0;
    }
    unsigned num = num_ + 1;
    perimeter = perimetersum_() + s;
    if (polyline_)
      return num;

    real tempsum =  areasum_();
    int crossings = crossings_;
    {
      real lat, lon, s12, S12, t;
      earth_.GenDirect(lat1_, lon1_, azi, false, s, mask_,
                       lat, lon, t, t, t, t, t, S12);
      tempsum += S12;
      crossings += transitdirect(lon1_, lon);
      earth_.GenInverse(lat, lon, lat0_, lon0_, mask_,
                        s12, t, t, t, t, t, S12);
      perimeter += s12;
      tempsum += S12;
      crossings += transit(lon, lon0_);
    }

    AreaReduce(tempsum, crossings, reverse, sign);
    area = real(0) + tempsum;
    return num;
  }

  template<class GeodType>
  template<typename T>
  void PolygonAreaT<GeodType>::AreaReduce(T& area, int crossings,
                                          bool reverse, bool sign) const {
    Remainder(area);
    if (crossings & 1) area += (area < 0 ? 1 : -1) * area0_/2;
    // area is with the clockwise sense.  If !reverse convert to
    // counter-clockwise convention.
    if (!reverse) area *= -1;
    // If sign put area in (-area0_/2, area0_/2], else put area in [0, area0_)
    if (sign) {
      if (area > area0_/2)
        area -= area0_;
      else if (area <= -area0_/2)
        area += area0_;
    } else {
      if (area >= area0_)
        area -= area0_;
      else if (area < 0)
        area += area0_;
    }
  }

  template class GEOGRAPHICLIB_EXPORT PolygonAreaT<Geodesic>;
  template class GEOGRAPHICLIB_EXPORT PolygonAreaT<GeodesicExact>;
  template class GEOGRAPHICLIB_EXPORT PolygonAreaT<Rhumb>;

} // namespace GeographicLib
