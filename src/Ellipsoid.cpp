/**
 * \file Ellipsoid.cpp
 * \brief Implementation for GeographicLib::Ellipsoid class
 *
 * Copyright (c) Charles Karney (2012-2023) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#include <GeographicLib/Ellipsoid.hpp>

namespace GeographicLib {

  using namespace std;

  /// \cond SKIP
  Ellipsoid::Ellipsoid(real a, real f)
    : stol_(real(0.01) * sqrt(numeric_limits<real>::epsilon()))
    , a_(a)
    , f_(f)
    , b_(a_ * (1 - f_))
    , e2_(f_ * (2 - f_))
    , e12_(e2_ / (1 - e2_))
    , n_(f_ / (2  - f_))
    , aux_(a_, f_)
    , rm_(aux_.RectifyingRadius(true))
    , c2_(aux_.AuthalicRadiusSquared(true))
  {}
  /// \endcond

  const Ellipsoid& Ellipsoid::WGS84() {
    static const Ellipsoid wgs84(Constants::WGS84_a(), Constants::WGS84_f());
    return wgs84;
  }

  real Ellipsoid::QuarterMeridian() const
  { return Math::pi()/2 * rm_; }

  real Ellipsoid::Area() const
  { return 4 * Math::pi() * c2_; }

  real Ellipsoid::ParametricLatitude(real phi) const {
    return aux_.Convert(AuxLatitude::PHI, AuxLatitude::BETA,
                        Math::LatFix(phi), true);
  }

  real Ellipsoid::InverseParametricLatitude(real beta) const {
    return aux_.Convert(AuxLatitude::BETA, AuxLatitude::PHI,
                        Math::LatFix(beta), true);
  }

  real Ellipsoid::GeocentricLatitude(real phi) const {
    return aux_.Convert(AuxLatitude::PHI, AuxLatitude::THETA,
                        Math::LatFix(phi), true);
  }

  real Ellipsoid::InverseGeocentricLatitude(real theta) const {
    return aux_.Convert(AuxLatitude::THETA, AuxLatitude::PHI,
                        Math::LatFix(theta), true);
  }

  real Ellipsoid::RectifyingLatitude(real phi) const {
    return aux_.Convert(AuxLatitude::PHI, AuxLatitude::MU,
                        Math::LatFix(phi), true);
  }

  real Ellipsoid::InverseRectifyingLatitude(real mu) const {
    return aux_.Convert(AuxLatitude::MU, AuxLatitude::PHI,
                        Math::LatFix(mu), true);
  }

  real Ellipsoid::AuthalicLatitude(real phi) const {
    return aux_.Convert(AuxLatitude::PHI, AuxLatitude::XI,
                        Math::LatFix(phi), true);
  }

  real Ellipsoid::InverseAuthalicLatitude(real xi) const {
    return aux_.Convert(AuxLatitude::XI, AuxLatitude::PHI,
                        Math::LatFix(xi), true);
  }

  real Ellipsoid::ConformalLatitude(real phi) const {
    return aux_.Convert(AuxLatitude::PHI, AuxLatitude::CHI,
                        Math::LatFix(phi), true);
  }

  real Ellipsoid::InverseConformalLatitude(real chi) const {
    return aux_.Convert(AuxLatitude::CHI, AuxLatitude::PHI,
                        Math::LatFix(chi), true);
  }

  real Ellipsoid::IsometricLatitude(real phi) const {
    return aux_.Convert(AuxLatitude::PHI, AuxLatitude::CHI,
                        AuxAngle::degrees(Math::LatFix(phi)), true).lamd();
  }

  real Ellipsoid::InverseIsometricLatitude(real psi) const {
    return  aux_.Convert(AuxLatitude::CHI, AuxLatitude::PHI,
                         AuxAngle::lamd(psi), true).degrees();
  }

  real Ellipsoid::CircleRadius(real phi) const {
    // a * cos(beta)
    AuxAngle beta(aux_.Convert(AuxLatitude::PHI, AuxLatitude::BETA,
                               AuxAngle::degrees(Math::LatFix(phi)),
                               true).normalized());
    return a_ * beta.x();
  }

  real Ellipsoid::CircleHeight(real phi) const {
    // b * sin(beta)
    AuxAngle beta(aux_.Convert(AuxLatitude::PHI, AuxLatitude::BETA,
                               AuxAngle::degrees(Math::LatFix(phi)),
                               true).normalized());
    return b_ * beta.y();
  }

  real Ellipsoid::MeridianDistance(real phi) const {
    return rm_ * aux_.Convert(AuxLatitude::PHI, AuxLatitude::MU,
                              AuxAngle::degrees(Math::LatFix(phi)),
                              true).radians();
  }

  real Ellipsoid::MeridionalCurvatureRadius(real phi) const {
    real v = 1 - e2_ * Math::sq(Math::sind(Math::LatFix(phi)));
    return a_ * (1 - e2_) / (v * sqrt(v));
  }

  real Ellipsoid::TransverseCurvatureRadius(real phi) const {
    real v = 1 - e2_ * Math::sq(Math::sind(Math::LatFix(phi)));
    return a_ / sqrt(v);
  }

  real Ellipsoid::NormalCurvatureRadius(real phi, real azi) const {
    real calp, salp,
      v = 1 - e2_ * Math::sq(Math::sind(Math::LatFix(phi)));
    Math::sincosd(azi, salp, calp);
    return a_ / (sqrt(v) * (Math::sq(calp) * v / (1 - e2_) + Math::sq(salp)));
  }

} // namespace GeographicLib
