/**
 * \file MagneticCircle.hpp
 * \brief Header for GeographicLib::MagneticCircle class
 *
 * Copyright (c) Charles Karney (2011-2022) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_MAGNETICCIRCLE_HPP)
#define GEOGRAPHICLIB_MAGNETICCIRCLE_HPP 1

#include <vector>
#include <GeographicLib/Constants.hpp>
#include <GeographicLib/CircularEngine.hpp>

namespace GeographicLib {

  /**
   * \brief Geomagnetic field on a circle of latitude
   *
   * Evaluate the earth's magnetic field on a circle of constant height and
   * latitude.  This uses a CircularEngine to pre-evaluate the inner sum of the
   * spherical harmonic sum, allowing the values of the field at several
   * different longitudes to be evaluated rapidly.
   *
   * Use MagneticModel::Circle to create a MagneticCircle object.  (The
   * constructor for this class is private.)
   *
   * Example of use:
   * \include examples/MagneticCircle.cpp
   *
   * <a href="MagneticField.1.html">MagneticField</a> is a command-line utility
   * providing access to the functionality of MagneticModel and MagneticCircle.
   **********************************************************************/

  class GEOGRAPHICLIB_EXPORT MagneticCircle {
  private:

    real a_, f_, lat_, h_, t_, cphi_, sphi_, t1_, dt0_;
    bool interpolate_, constterm_;
    CircularEngine circ0_, circ1_, circ2_;

    MagneticCircle(real a, real f, real lat, real h, real t,
                   real cphi, real sphi, real t1, real dt0,
                   bool interpolate,
                   const CircularEngine& circ0, const CircularEngine& circ1)
      : a_(a)
      , f_(f)
      , lat_(Math::LatFix(lat))
      , h_(h)
      , t_(t)
      , cphi_(cphi)
      , sphi_(sphi)
      , t1_(t1)
      , dt0_(dt0)
      , interpolate_(interpolate)
      , constterm_(false)
      , circ0_(circ0)
      , circ1_(circ1)
    {}

    MagneticCircle(real a, real f, real lat, real h, real t,
                   real cphi, real sphi, real t1, real dt0,
                   bool interpolate,
                   const CircularEngine& circ0, const CircularEngine& circ1,
                   const CircularEngine& circ2)
      : a_(a)
      , f_(f)
      , lat_(lat)
      , h_(h)
      , t_(t)
      , cphi_(cphi)
      , sphi_(sphi)
      , t1_(t1)
      , dt0_(dt0)
      , interpolate_(interpolate)
      , constterm_(true)
      , circ0_(circ0)
      , circ1_(circ1)
      , circ2_(circ2)
    {}

    void Field(real lon, bool diffp,
               real& Bx, real& By, real& Bz,
               real& Bxt, real& Byt, real& Bzt) const;

    void FieldGeocentric(real slam, real clam,
                         real& BX, real& BY, real& BZ,
                         real& BXt, real& BYt, real& BZt) const;

    friend class MagneticModel; // MagneticModel calls the private constructor

  public:

    /**
     * A default constructor for the normal gravity.  This sets up an
     * uninitialized object which can be later replaced by the
     * MagneticModel::Circle.
     **********************************************************************/
    MagneticCircle() : a_(-1) {}

    /** \name Compute the magnetic field
     **********************************************************************/
    ///@{
    /**
     * Evaluate the components of the geomagnetic field at a particular
     * longitude.
     *
     * @param[in] lon longitude of the point (degrees).
     * @param[out] Bx the easterly component of the magnetic field (nanotesla).
     * @param[out] By the northerly component of the magnetic field
     *   (nanotesla).
     * @param[out] Bz the vertical (up) component of the magnetic field
     *   (nanotesla).
     **********************************************************************/
    void operator()(real lon, real& Bx, real& By, real& Bz) const {
      real dummy;
      Field(lon, false, Bx, By, Bz, dummy, dummy, dummy);
    }

    /**
     * Evaluate the components of the geomagnetic field and their time
     * derivatives at a particular longitude.
     *
     * @param[in] lon longitude of the point (degrees).
     * @param[out] Bx the easterly component of the magnetic field (nanotesla).
     * @param[out] By the northerly component of the magnetic field
     *   (nanotesla).
     * @param[out] Bz the vertical (up) component of the magnetic field
     *   (nanotesla).
     * @param[out] Bxt the rate of change of \e Bx (nT/yr).
     * @param[out] Byt the rate of change of \e By (nT/yr).
     * @param[out] Bzt the rate of change of \e Bz (nT/yr).
     **********************************************************************/
    void operator()(real lon, real& Bx, real& By, real& Bz,
                    real& Bxt, real& Byt, real& Bzt) const {
      Field(lon, true, Bx, By, Bz, Bxt, Byt, Bzt);
    }

    /**
     * Evaluate the components of the geomagnetic field and their time
     * derivatives at a particular longitude.
     *
     * @param[in] lon longitude of the point (degrees).
     * @param[out] BX the \e X component of the magnetic field (nT).
     * @param[out] BY the \e Y component of the magnetic field (nT).
     * @param[out] BZ the \e Z component of the magnetic field (nT).
     * @param[out] BXt the rate of change of \e BX (nT/yr).
     * @param[out] BYt the rate of change of \e BY (nT/yr).
     * @param[out] BZt the rate of change of \e BZ (nT/yr).
     **********************************************************************/
    void FieldGeocentric(real lon, real& BX, real& BY, real& BZ,
                         real& BXt, real& BYt, real& BZt) const;
    ///@}

    /** \name Inspector functions
     **********************************************************************/
    ///@{
    /**
     * @return true if the object has been initialized.
     **********************************************************************/
    bool Init() const { return a_ > 0; }
    /**
     * @return \e a the equatorial radius of the ellipsoid (meters).  This is
     *   the value inherited from the MagneticModel object used in the
     *   constructor.
     **********************************************************************/
    real EquatorialRadius() const
    { return Init() ? a_ : Math::NaN(); }
    /**
     * @return \e f the flattening of the ellipsoid.  This is the value
     *   inherited from the MagneticModel object used in the constructor.
     **********************************************************************/
    real Flattening() const
    { return Init() ? f_ : Math::NaN(); }
    /**
     * @return the latitude of the circle (degrees).
     **********************************************************************/
    real Latitude() const
    { return Init() ? lat_ : Math::NaN(); }
    /**
     * @return the height of the circle (meters).
     **********************************************************************/
    real Height() const
    { return Init() ? h_ : Math::NaN(); }
    /**
     * @return the time (fractional years).
     **********************************************************************/
    real Time() const
    { return Init() ? t_ : Math::NaN(); }
    ///@}
  };

} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_MAGNETICCIRCLE_HPP
