/**
 * \file GravityCircle.cpp
 * \brief Implementation for GeographicLib::GravityCircle class
 *
 * Copyright (c) Charles Karney (2011-2020) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#include <GeographicLib/GravityCircle.hpp>
#include <fstream>
#include <sstream>
#include <GeographicLib/Geocentric.hpp>

namespace GeographicLib {

  using namespace std;

  GravityCircle::GravityCircle(mask caps, real a, real f, real lat, real h,
                               real Z, real P, real cphi, real sphi,
                               real amodel, real GMmodel,
                               real dzonal0, real corrmult,
                               real gamma0, real gamma, real frot,
                               const CircularEngine& gravitational,
                               const CircularEngine& disturbing,
                               const CircularEngine& correction)
    : caps_(caps)
    , a_(a)
    , f_(f)
    , lat_(Math::LatFix(lat))
    , h_(h)
    , zZ_(Z)
    , pPx_(P)
    , invR_(1 / hypot(pPx_, zZ_))
    , cpsi_(pPx_ * invR_)
    , spsi_(zZ_ * invR_)
    , cphi_(cphi)
    , sphi_(sphi)
    , amodel_(amodel)
    , gGMmodel_(GMmodel)
    , dzonal0_(dzonal0)
    , corrmult_(corrmult)
    , gamma0_(gamma0)
    , gamma_(gamma)
    , frot_(frot)
    , gravitational_(gravitational)
    , disturbing_(disturbing)
    , correction_(correction)
    {}

  real GravityCircle::Gravity(real lon,
                                    real& gx, real& gy, real& gz) const {
    real slam, clam, M[Geocentric::dim2_];
    Math::sincosd(lon, slam, clam);
    real Wres = W(slam, clam, gx, gy, gz);
    Geocentric::Rotation(sphi_, cphi_, slam, clam, M);
    Geocentric::Unrotate(M, gx, gy, gz, gx, gy, gz);
    return Wres;
  }

  real GravityCircle::Disturbance(real lon, real& deltax, real& deltay,
                                        real& deltaz) const {
    real slam, clam, M[Geocentric::dim2_];
    Math::sincosd(lon, slam, clam);
    real Tres = InternalT(slam, clam, deltax, deltay, deltaz, true, true);
    Geocentric::Rotation(sphi_, cphi_, slam, clam, M);
    Geocentric::Unrotate(M, deltax, deltay, deltaz, deltax, deltay, deltaz);
    return Tres;
  }

  real GravityCircle::GeoidHeight(real lon) const {
    if ((caps_ & GEOID_HEIGHT) != GEOID_HEIGHT)
      return Math::NaN();
    real slam, clam, dummy;
    Math::sincosd(lon, slam, clam);
    real T = InternalT(slam, clam, dummy, dummy, dummy, false, false);
    real correction = corrmult_ * correction_(slam, clam);
    return T/gamma0_ + correction;
  }

  void GravityCircle::SphericalAnomaly(real lon,
                                       real& Dg01, real& xi, real& eta) const {
    if ((caps_ & SPHERICAL_ANOMALY) != SPHERICAL_ANOMALY) {
      Dg01 = xi = eta = Math::NaN();
      return;
    }
    real slam, clam;
    Math::sincosd(lon, slam, clam);
    real
      deltax, deltay, deltaz,
      T = InternalT(slam, clam, deltax, deltay, deltaz, true, false);
    // Rotate cartesian into spherical coordinates
    real MC[Geocentric::dim2_];
    Geocentric::Rotation(spsi_, cpsi_, slam, clam, MC);
    Geocentric::Unrotate(MC, deltax, deltay, deltaz, deltax, deltay, deltaz);
    // H+M, Eq 2-151c
    Dg01 = - deltaz - 2 * T * invR_;
    xi  = -(deltay/gamma_) / Math::degree();
    eta = -(deltax/gamma_) / Math::degree();
  }

  real GravityCircle::W(real slam, real clam,
                              real& gX, real& gY, real& gZ) const {
    real Wres = V(slam, clam, gX, gY, gZ) + frot_ * pPx_ / 2;
    gX += frot_ * clam;
    gY += frot_ * slam;
    return Wres;
  }

  real GravityCircle::V(real slam, real clam,
                              real& GX, real& GY, real& GZ) const {
    if ((caps_ & GRAVITY) != GRAVITY) {
      GX = GY = GZ = Math::NaN();
      return Math::NaN();
    }
    real
      Vres = gravitational_(slam, clam, GX, GY, GZ),
      f = gGMmodel_ / amodel_;
    Vres *= f;
    GX *= f;
    GY *= f;
    GZ *= f;
    return Vres;
  }

  real GravityCircle::InternalT(real slam, real clam,
                                      real& deltaX, real& deltaY, real& deltaZ,
                                      bool gradp, bool correct) const {
    if (gradp) {
      if ((caps_ & DISTURBANCE) != DISTURBANCE) {
        deltaX = deltaY = deltaZ = Math::NaN();
        return Math::NaN();
      }
    } else {
      if ((caps_ & DISTURBING_POTENTIAL) != DISTURBING_POTENTIAL)
        return Math::NaN();
    }
    if (dzonal0_ == 0)
      correct = false;
    real T = (gradp
              ? disturbing_(slam, clam, deltaX, deltaY, deltaZ)
              : disturbing_(slam, clam));
    T = (T / amodel_ - (correct ? dzonal0_ : 0) * invR_) * gGMmodel_;
    if (gradp) {
      real f = gGMmodel_ / amodel_;
      deltaX *= f;
      deltaY *= f;
      deltaZ *= f;
      if (correct) {
        real r3 = gGMmodel_ * dzonal0_ * invR_ * invR_ * invR_;
        deltaX += pPx_ * clam * r3;
        deltaY += pPx_ * slam * r3;
        deltaZ += zZ_ * r3;
      }
    }
    return T;
  }

} // namespace GeographicLib
