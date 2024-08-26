/**
 * \file AlbersEqualArea.hpp
 * \brief Header for GeographicLib::AlbersEqualArea class
 *
 * Copyright (c) Charles Karney (2010-2023) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_ALBERSEQUALAREA_HPP)
#define GEOGRAPHICLIB_ALBERSEQUALAREA_HPP 1

#include <GeographicLib/Constants.hpp>

namespace GeographicLib {

  /**
   * \brief Albers equal area conic projection
   *
   * Implementation taken from the report,
   * - J. P. Snyder,
   *   <a href="https://pubs.usgs.gov/publication/pp1395"> Map Projections: A
   *   Working Manual</a>, USGS Professional Paper 1395 (1987),
   *   pp. 101--102.
   *
   * This is a implementation of the equations in Snyder except that divided
   * differences will be [have been] used to transform the expressions into
   * ones which may be evaluated accurately.  [In this implementation, the
   * projection correctly becomes the cylindrical equal area or the azimuthal
   * equal area projection when the standard latitude is the equator or a
   * pole.]
   *
   * The ellipsoid parameters, the standard parallels, and the scale on the
   * standard parallels are set in the constructor.  Internally, the case with
   * two standard parallels is converted into a single standard parallel, the
   * latitude of minimum azimuthal scale, with an azimuthal scale specified on
   * this parallel.  This latitude is also used as the latitude of origin which
   * is returned by AlbersEqualArea::OriginLatitude.  The azimuthal scale on
   * the latitude of origin is given by AlbersEqualArea::CentralScale.  The
   * case with two standard parallels at opposite poles is singular and is
   * disallowed.  The central meridian (which is a trivial shift of the
   * longitude) is specified as the \e lon0 argument of the
   * AlbersEqualArea::Forward and AlbersEqualArea::Reverse functions.
   * AlbersEqualArea::Forward and AlbersEqualArea::Reverse also return the
   * meridian convergence, &gamma;, and azimuthal scale, \e k.  A small square
   * aligned with the cardinal directions is projected to a rectangle with
   * dimensions \e k (in the E-W direction) and 1/\e k (in the N-S direction).
   * The E-W sides of the rectangle are oriented &gamma; degrees
   * counter-clockwise from the \e x axis.  There is no provision in this class
   * for specifying a false easting or false northing or a different latitude
   * of origin.
   *
   * Example of use:
   * \include examples/AlbersEqualArea.cpp
   *
   * <a href="ConicProj.1.html">ConicProj</a> is a command-line utility
   * providing access to the functionality of LambertConformalConic and
   * AlbersEqualArea.
   **********************************************************************/
  class GEOGRAPHICLIB_EXPORT AlbersEqualArea {
  private:
    real a_, f_, fm_, e2_, e_, e2m_, qZ_, qx_;
    real sign_, lat0_, k0_;
    real n0_, m02_, nrho0_, k2_, txi0_, scxi0_, sxi0_;

    real atanhee (real x) const;
    real Datanhee (real x, real y) const;

    // DDatanhee(x,y) = (Datanhee(1,y) - Datanhee(1,x))/(y-x)
    real DDatanhee(real x, real y) const;
    real DDatanhee0(real x, real y) const;
    real DDatanhee1(real x, real y) const;
    real DDatanhee2(real x, real y) const;
    void Init(real sphi1, real cphi1, real sphi2, real cphi2, real k1);
    real txif(real tphi) const;
    real tphif(real txi) const;
  public:

    /**
     * Constructor with a single standard parallel.
     *
     * @param[in] a equatorial radius of ellipsoid (meters).
     * @param[in] f flattening of ellipsoid.  Setting \e f = 0 gives a sphere.
     *   Negative \e f gives a prolate ellipsoid.
     * @param[in] stdlat standard parallel (degrees), the circle of tangency.
     * @param[in] k0 azimuthal scale on the standard parallel.
     * @exception GeographicErr if \e a, (1 &minus; \e f) \e a, or \e k0 is
     *   not positive.
     * @exception GeographicErr if \e stdlat is not in [&minus;90&deg;,
     *   90&deg;].
     **********************************************************************/
    AlbersEqualArea(real a, real f, real stdlat, real k0);

    /**
     * Constructor with two standard parallels.
     *
     * @param[in] a equatorial radius of ellipsoid (meters).
     * @param[in] f flattening of ellipsoid.  Setting \e f = 0 gives a sphere.
     *   Negative \e f gives a prolate ellipsoid.
     * @param[in] stdlat1 first standard parallel (degrees).
     * @param[in] stdlat2 second standard parallel (degrees).
     * @param[in] k1 azimuthal scale on the standard parallels.
     * @exception GeographicErr if \e a, (1 &minus; \e f) \e a, or \e k1 is
     *   not positive.
     * @exception GeographicErr if \e stdlat1 or \e stdlat2 is not in
     *   [&minus;90&deg;, 90&deg;], or if \e stdlat1 and \e stdlat2 are
     *   opposite poles.
     **********************************************************************/
    AlbersEqualArea(real a, real f, real stdlat1, real stdlat2, real k1);

    /**
     * Constructor with two standard parallels specified by sines and cosines.
     *
     * @param[in] a equatorial radius of ellipsoid (meters).
     * @param[in] f flattening of ellipsoid.  Setting \e f = 0 gives a sphere.
     *   Negative \e f gives a prolate ellipsoid.
     * @param[in] sinlat1 sine of first standard parallel.
     * @param[in] coslat1 cosine of first standard parallel.
     * @param[in] sinlat2 sine of second standard parallel.
     * @param[in] coslat2 cosine of second standard parallel.
     * @param[in] k1 azimuthal scale on the standard parallels.
     * @exception GeographicErr if \e a, (1 &minus; \e f) \e a, or \e k1 is
     *   not positive.
     * @exception GeographicErr if \e stdlat1 or \e stdlat2 is not in
     *   [&minus;90&deg;, 90&deg;], or if \e stdlat1 and \e stdlat2 are
     *   opposite poles.
     *
     * This allows parallels close to the poles to be specified accurately.
     * This routine computes the latitude of origin and the azimuthal scale at
     * this latitude.  If \e dlat = abs(\e lat2 &minus; \e lat1) &le; 160&deg;,
     * then the error in the latitude of origin is less than 4.5 &times;
     * 10<sup>&minus;14</sup>d;.
     **********************************************************************/
    AlbersEqualArea(real a, real f,
                    real sinlat1, real coslat1,
                    real sinlat2, real coslat2,
                    real k1);

    /**
     * Set the azimuthal scale for the projection.
     *
     * @param[in] lat (degrees).
     * @param[in] k azimuthal scale at latitude \e lat (default 1).
     * @exception GeographicErr \e k is not positive.
     * @exception GeographicErr if \e lat is not in (&minus;90&deg;,
     *   90&deg;).
     *
     * This allows a "latitude of conformality" to be specified.
     **********************************************************************/
    void SetScale(real lat, real k = real(1));

    /**
     * Forward projection, from geographic to Lambert conformal conic.
     *
     * @param[in] lon0 central meridian longitude (degrees).
     * @param[in] lat latitude of point (degrees).
     * @param[in] lon longitude of point (degrees).
     * @param[out] x easting of point (meters).
     * @param[out] y northing of point (meters).
     * @param[out] gamma meridian convergence at point (degrees).
     * @param[out] k azimuthal scale of projection at point; the radial
     *   scale is the 1/\e k.
     *
     * The latitude origin is given by AlbersEqualArea::LatitudeOrigin().  No
     * false easting or northing is added and \e lat should be in the range
     * [&minus;90&deg;, 90&deg;].  The values of \e x and \e y returned for
     * points which project to infinity (i.e., one or both of the poles) will
     * be large but finite.
     **********************************************************************/
    void Forward(real lon0, real lat, real lon,
                 real& x, real& y, real& gamma, real& k) const;

    /**
     * Reverse projection, from Lambert conformal conic to geographic.
     *
     * @param[in] lon0 central meridian longitude (degrees).
     * @param[in] x easting of point (meters).
     * @param[in] y northing of point (meters).
     * @param[out] lat latitude of point (degrees).
     * @param[out] lon longitude of point (degrees).
     * @param[out] gamma meridian convergence at point (degrees).
     * @param[out] k azimuthal scale of projection at point; the radial
     *   scale is the 1/\e k.
     *
     * The latitude origin is given by AlbersEqualArea::LatitudeOrigin().  No
     * false easting or northing is added.  The value of \e lon returned is in
     * the range [&minus;180&deg;, 180&deg;].  The value of \e lat returned is
     * in the range [&minus;90&deg;, 90&deg;].  If the input point is outside
     * the legal projected space the nearest pole is returned.
     **********************************************************************/
    void Reverse(real lon0, real x, real y,
                 real& lat, real& lon, real& gamma, real& k) const;

    /**
     * AlbersEqualArea::Forward without returning the convergence and
     * scale.
     **********************************************************************/
    void Forward(real lon0, real lat, real lon,
                 real& x, real& y) const {
      real gamma, k;
      Forward(lon0, lat, lon, x, y, gamma, k);
    }

    /**
     * AlbersEqualArea::Reverse without returning the convergence and
     * scale.
     **********************************************************************/
    void Reverse(real lon0, real x, real y,
                 real& lat, real& lon) const {
      real gamma, k;
      Reverse(lon0, x, y, lat, lon, gamma, k);
    }

    /** \name Inspector functions
     **********************************************************************/
    ///@{
    /**
     * @return \e a the equatorial radius of the ellipsoid (meters).  This is
     *   the value used in the constructor.
     **********************************************************************/
    real EquatorialRadius() const { return a_; }

    /**
     * @return \e f the flattening of the ellipsoid.  This is the value used in
     *   the constructor.
     **********************************************************************/
    real Flattening() const { return f_; }

    /**
     * @return latitude of the origin for the projection (degrees).
     *
     * This is the latitude of minimum azimuthal scale and equals the \e stdlat
     * in the 1-parallel constructor and lies between \e stdlat1 and \e stdlat2
     * in the 2-parallel constructors.
     **********************************************************************/
    real OriginLatitude() const { return lat0_; }

    /**
     * @return central scale for the projection.  This is the azimuthal scale
     *   on the latitude of origin.
     **********************************************************************/
    real CentralScale() const { return k0_; }
    ///@}

    /**
     * A global instantiation of AlbersEqualArea with the WGS84 ellipsoid, \e
     * stdlat = 0, and \e k0 = 1.  This degenerates to the cylindrical equal
     * area projection.
     **********************************************************************/
    static const AlbersEqualArea& CylindricalEqualArea();

    /**
     * A global instantiation of AlbersEqualArea with the WGS84 ellipsoid, \e
     * stdlat = 90&deg;, and \e k0 = 1.  This degenerates to the
     * Lambert azimuthal equal area projection.
     **********************************************************************/
    static const AlbersEqualArea& AzimuthalEqualAreaNorth();

    /**
     * A global instantiation of AlbersEqualArea with the WGS84 ellipsoid, \e
     * stdlat = &minus;90&deg;, and \e k0 = 1.  This degenerates to the
     * Lambert azimuthal equal area projection.
     **********************************************************************/
    static const AlbersEqualArea& AzimuthalEqualAreaSouth();
  };

} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_ALBERSEQUALAREA_HPP
