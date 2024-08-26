/**
 * \file Georef.hpp
 * \brief Header for GeographicLib::Georef namespace
 *
 * Copyright (c) Charles Karney (2015-2024) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_GEOREF_HPP)
#define GEOGRAPHICLIB_GEOREF_HPP 1

#include <GeographicLib/Constants.hpp>

namespace GeographicLib {

  /**
   * \brief Conversions for the World Geographic Reference System (georef)
   *
   * The World Geographic Reference System is described in
   * - https://en.wikipedia.org/wiki/Georef
   * - https://web.archive.org/web/20161214054445/http://earth-info.nga.mil/GandG/coordsys/grids/georef.pdf
   * .
   * It provides a compact string representation of a geographic area
   * (expressed as latitude and longitude).  The classes GARS and Geohash
   * implement similar compact representations.
   *
   * Example of use:
   * \include examples/Georef.cpp
   **********************************************************************/

  namespace Georef {
    constexpr int maxprec_ = 11;        // approximately equivalent to MGRS class
    constexpr int base_ = 10;           // Base for minutes

    /**
     * Convert from geographic coordinates to georef.
     *
     * @param[in] lat latitude of point (degrees).
     * @param[in] lon longitude of point (degrees).
     * @param[in] prec the precision of the resulting georef.
     * @param[out] georef the georef string.
     * @exception GeographicErr if \e lat is not in [&minus;90&deg;,
     *   90&deg;].
     * @exception std::bad_alloc if memory for \e georef can't be allocated.
     *
     * \e prec specifies the precision of \e georef as follows:
     * - \e prec = &minus;1 (min), 15&deg;
     * - \e prec = 0, 1&deg;
     * - \e prec = 1, converted to \e prec = 2
     * - \e prec = 2, 1'
     * - \e prec = 3, 0.1'
     * - \e prec = 4, 0.01'
     * - \e prec = 5, 0.001'
     * - &hellip;
     * - \e prec = 11 (max), 10<sup>&minus;9</sup>'
     *
     * If \e lat or \e lon is NaN, then \e georef is set to "INVALID".
     **********************************************************************/
    void Forward(real lat, real lon, int prec, std::string& georef);

    /**
     * Convert from Georef to geographic coordinates.
     *
     * @param[in] georef the Georef.
     * @param[out] lat latitude of point (degrees).
     * @param[out] lon longitude of point (degrees).
     * @param[out] prec the precision of \e georef.
     * @param[in] centerp if true (the default) return the center
     *   \e georef, otherwise return the south-west corner.
     * @exception GeographicErr if \e georef is illegal.
     *
     * The case of the letters in \e georef is ignored.  \e prec is in the
     * range [&minus;1, 11] and gives the precision of \e georef as follows:
     * - \e prec = &minus;1 (min), 15&deg;
     * - \e prec = 0, 1&deg;
     * - \e prec = 1, not returned
     * - \e prec = 2, 1'
     * - \e prec = 3, 0.1'
     * - \e prec = 4, 0.01'
     * - \e prec = 5, 0.001'
     * - &hellip;
     * - \e prec = 11 (max), 10<sup>&minus;9</sup>'
     *
     * If the first 3 characters of \e georef are "INV", then \e lat and \e lon
     * are set to NaN and \e prec is unchanged.
     **********************************************************************/
    void Reverse(const std::string& georef, real& lat, real& lon,
                        int& prec, bool centerp = true);

    /**
     * The angular resolution of a Georef.
     *
     * @param[in] prec the precision of the Georef.
     * @return the latitude-longitude resolution (degrees).
     *
     * Internally, \e prec is first put in the range [&minus;1, 11].
     **********************************************************************/
    inline real Resolution(int prec) {
      if (prec < 1)
        return real(prec < 0 ? 15 : 1);
      else {
        // Treat prec = 1 as 2.
        prec = std::max(2, std::min(maxprec_, prec));
        // Need extra real because, since C++11, pow(float, int) returns double
        return 1/(60 * real(std::pow(real(base_), prec - 2)));
      }
    }

    /**
     * The Georef precision required to meet a given geographic resolution.
     *
     * @param[in] res the minimum of resolution in latitude and longitude
     *   (degrees).
     * @return Georef precision.
     *
     * The returned length is in the range [0, 11].
     **********************************************************************/
    inline int Precision(real res) {
      using std::fabs; res = fabs(res);
      for (int prec = 0; prec < maxprec_; ++prec) {
        if (prec == 1)
          continue;
        if (Resolution(prec) <= res)
          return prec;
      }
      return maxprec_;
    }

  } //namespace Georef

} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_GEOREF_HPP
