/**
 * \file CircularEngine.hpp
 * \brief Header for GeographicLib::CircularEngine class
 *
 * Copyright (c) Charles Karney (2011-2015) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_CIRCULARENGINE_HPP)
#define GEOGRAPHICLIB_CIRCULARENGINE_HPP 1

#include <vector>
#include <GeographicLib/Constants.hpp>
#include <GeographicLib/SphericalEngine.hpp>

#if defined(_MSC_VER)
// Squelch warnings about dll vs vector
#  pragma warning (push)
#  pragma warning (disable: 4251)
#endif

namespace GeographicLib {

  /**
   * \brief Spherical harmonic sums for a circle
   *
   * The class is a companion to SphericalEngine.  If the results of a
   * spherical harmonic sum are needed for several points on a circle of
   * constant latitude \e lat and height \e h, then SphericalEngine::Circle can
   * compute the inner sum, which is independent of longitude \e lon, and
   * produce a CircularEngine object.  CircularEngine::operator()() can
   * then be used to perform the outer sum for particular vales of \e lon.
   * This can lead to substantial improvements in computational speed for high
   * degree sum (approximately by a factor of \e N / 2 where \e N is the
   * maximum degree).
   *
   * CircularEngine is tightly linked to the internals of SphericalEngine.  For
   * that reason, the constructor for this class is private.  Use
   * SphericalHarmonic::Circle, SphericalHarmonic1::Circle, and
   * SphericalHarmonic2::Circle to create instances of this class.
   *
   * CircularEngine stores the coefficients needed to allow the summation over
   * order to be performed in 2 or 6 vectors of length \e M + 1 (depending on
   * whether gradients are to be calculated).  For this reason the constructor
   * may throw a std::bad_alloc exception.
   *
   * Example of use:
   * \include examples/CircularEngine.cpp
   **********************************************************************/

  class GEOGRAPHICLIB_EXPORT CircularEngine {
  private:
    enum normalization {
      FULL = SphericalEngine::FULL,
      SCHMIDT = SphericalEngine::SCHMIDT,
    };
    int mM_;
    bool gradp_;
    unsigned norm_;
    real _a, _r, _u, _t;
    std::vector<real> wc_, ws_, wrc_, wrs_, wtc_, wts_;
    real q_, uq_, uq2_;

    real Value(bool gradp, real sl, real cl,
                     real& gradx, real& grady, real& gradz) const;

    friend class SphericalEngine;
    CircularEngine(int M, bool gradp, unsigned norm,
                   real a, real r, real u, real t)
      : mM_(M)
      , gradp_(gradp)
      , norm_(norm)
      , _a(a)
      , _r(r)
      , _u(u)
      , _t(t)
      , wc_(std::vector<real>(mM_ + 1, 0))
      , ws_(std::vector<real>(mM_ + 1, 0))
      , wrc_(std::vector<real>(gradp_ ? mM_ + 1 : 0, 0))
      , wrs_(std::vector<real>(gradp_ ? mM_ + 1 : 0, 0))
      , wtc_(std::vector<real>(gradp_ ? mM_ + 1 : 0, 0))
      , wts_(std::vector<real>(gradp_ ? mM_ + 1 : 0, 0))
      {
        q_ = _a / _r;
        uq_ = _u * q_;
        uq2_ = Math::sq(uq_);
      }

    void SetCoeff(int m, real wc, real ws)
    { wc_[m] = wc; ws_[m] = ws; }

    void SetCoeff(int m, real wc, real ws,
                  real wrc, real wrs, real wtc, real wts) {
      wc_[m] = wc; ws_[m] = ws;
      if (gradp_) {
        wrc_[m] = wrc; wrs_[m] = wrs;
        wtc_[m] = wtc; wts_[m] = wts;
      }
    }

  public:

    /**
     * A default constructor.  CircularEngine::operator()() on the resulting
     * object returns zero.  The resulting object can be assigned to the result
     * of SphericalHarmonic::Circle.
     **********************************************************************/
    CircularEngine()
      : mM_(-1)
      , gradp_(true)
      , _u(0)
      , _t(1)
      {}

    /**
     * Evaluate the sum for a particular longitude given in terms of its
     * sine and cosine.
     *
     * @param[in] sinlon the sine of the longitude.
     * @param[in] coslon the cosine of the longitude.
     * @return \e V the value of the sum.
     *
     * The arguments must satisfy <i>sinlon</i><sup>2</sup> +
     * <i>coslon</i><sup>2</sup> = 1.
     **********************************************************************/
    real operator()(real sinlon, real coslon) const {
      real dummy;
      return Value(false, sinlon, coslon, dummy, dummy, dummy);
    }

    /**
     * Evaluate the sum for a particular longitude.
     *
     * @param[in] lon the longitude (degrees).
     * @return \e V the value of the sum.
     **********************************************************************/
    real operator()(real lon) const {
      real sinlon, coslon;
      Math::sincosd(lon, sinlon, coslon);
      return (*this)(sinlon, coslon);
    }

    /**
     * Evaluate the sum and its gradient for a particular longitude given in
     * terms of its sine and cosine.
     *
     * @param[in] sinlon the sine of the longitude.
     * @param[in] coslon the cosine of the longitude.
     * @param[out] gradx \e x component of the gradient.
     * @param[out] grady \e y component of the gradient.
     * @param[out] gradz \e z component of the gradient.
     * @return \e V the value of the sum.
     *
     * The gradients will only be computed if the CircularEngine object was
     * created with this capability (e.g., via \e gradp = true in
     * SphericalHarmonic::Circle).  If not, \e gradx, etc., will not be
     * touched.  The arguments must satisfy <i>sinlon</i><sup>2</sup> +
     * <i>coslon</i><sup>2</sup> = 1.
     **********************************************************************/
    real operator()(real sinlon, real coslon,
                          real& gradx, real& grady, real& gradz) const {
      return Value(true, sinlon, coslon, gradx, grady, gradz);
    }

    /**
     * Evaluate the sum and its gradient for a particular longitude.
     *
     * @param[in] lon the longitude (degrees).
     * @param[out] gradx \e x component of the gradient.
     * @param[out] grady \e y component of the gradient.
     * @param[out] gradz \e z component of the gradient.
     * @return \e V the value of the sum.
     *
     * The gradients will only be computed if the CircularEngine object was
     * created with this capability (e.g., via \e gradp = true in
     * SphericalHarmonic::Circle).  If not, \e gradx, etc., will not be
     * touched.
     **********************************************************************/
    real operator()(real lon,
                          real& gradx, real& grady, real& gradz) const {
      real sinlon, coslon;
      Math::sincosd(lon, sinlon, coslon);
      return (*this)(sinlon, coslon, gradx, grady, gradz);
    }
  };

} // namespace GeographicLib

#if defined(_MSC_VER)
#  pragma warning (pop)
#endif

#endif  // GEOGRAPHICLIB_CIRCULARENGINE_HPP
