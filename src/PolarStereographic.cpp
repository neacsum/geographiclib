/**
 * \file PolarStereographic.cpp
 * \brief Implementation for GeographicLib::PolarStereographic class
 *
 * Copyright (c) Charles Karney (2008-2022) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#include <GeographicLib/PolarStereographic.hpp>

namespace GeographicLib {

  using namespace std;

  PolarStereographic::PolarStereographic(real a, real f, real k0)
    : a_(a)
    , f_(f)
    , e2_(f_ * (2 - f_))
    , es_((f_ < 0 ? -1 : 1) * sqrt(fabs(e2_)))
    , e2m_(1 - e2_)
    , c_( (1 - f_) * exp(Math::eatanhe(real(1), es_)) )
    , k0_(k0)
  {
    if (!(isfinite(a_) && a_ > 0))
      throw GeographicErr("Equatorial radius is not positive");
    if (!(isfinite(f_) && f_ < 1))
      throw GeographicErr("Polar semi-axis is not positive");
    if (!(isfinite(k0_) && k0_ > 0))
      throw GeographicErr("Scale is not positive");
  }

  const PolarStereographic& PolarStereographic::UPS() {
    static const PolarStereographic ups(Constants::WGS84_a(),
                                        Constants::WGS84_f(),
                                        Constants::UPS_k0());
    return ups;
  }

  // This formulation converts to conformal coordinates by tau = tan(phi) and
  // tau' = tan(phi') where phi' is the conformal latitude.  The formulas are:
  //    tau = tan(phi)
  //    secphi = hypot(1, tau)
  //    sig = sinh(e * atanh(e * tau / secphi))
  //    taup = tan(phip) = tau * hypot(1, sig) - sig * hypot(1, tau)
  //    c = (1 - f) * exp(e * atanh(e))
  //
  // Forward:
  //   rho = (2*k0*a/c) / (hypot(1, taup) + taup)  (taup >= 0)
  //       = (2*k0*a/c) * (hypot(1, taup) - taup)  (taup <  0)
  //
  // Reverse:
  //   taup = ((2*k0*a/c) / rho - rho / (2*k0*a/c))/2
  //
  // Scale:
  //   k = (rho/a) * secphi * sqrt((1-e2) + e2 / secphi^2)
  //
  // In limit rho -> 0, tau -> inf, taup -> inf, secphi -> inf, secphip -> inf
  //   secphip = taup = exp(-e * atanh(e)) * tau = exp(-e * atanh(e)) * secphi

  void PolarStereographic::Forward(bool northp, real lat, real lon,
                                   real& x, real& y,
                                   real& gamma, real& k) const {
    lat = Math::LatFix(lat);
    lat *= northp ? 1 : -1;
    real
      tau = Math::tand(lat),
      secphi = hypot(real(1), tau),
      taup = Math::taupf(tau, es_),
      rho = hypot(real(1), taup) + fabs(taup);
    rho = taup >= 0 ? (lat != 90 ? 1/rho : 0) : rho;
    rho *= 2 * k0_ * a_ / c_;
    k = lat != 90 ?
      (rho / a_) * secphi * sqrt(e2m_ + e2_ / Math::sq(secphi)) : k0_;
    Math::sincosd(lon, x, y);
    x *= rho;
    y *= (northp ? -rho : rho);
    gamma = Math::AngNormalize(northp ? lon : -lon);
  }

  void PolarStereographic::Reverse(bool northp, real x, real y,
                                   real& lat, real& lon,
                                   real& gamma, real& k) const {
    real
      rho = hypot(x, y),
      t = rho != 0 ? rho / (2 * k0_ * a_ / c_) :
      Math::sq(numeric_limits<real>::epsilon()),
      taup = (1 / t - t) / 2,
      tau = Math::tauf(taup, es_),
      secphi = hypot(real(1), tau);
    k = rho != 0 ? (rho / a_) * secphi * sqrt(e2m_ + e2_ / Math::sq(secphi)) :
      k0_;
    lat = (northp ? 1 : -1) * Math::atand(tau);
    lon = Math::atan2d(x, northp ? -y : y );
    gamma = Math::AngNormalize(northp ? lon : -lon);
  }

  void PolarStereographic::SetScale(real lat, real k) {
    if (!(isfinite(k) && k > 0))
      throw GeographicErr("Scale is not positive");
    if (!(-90 < lat && lat <= 90))
      throw GeographicErr ("Latitude must be in (-90d, +90d)");
    real x, y, gamma, kold;
    k0_ = 1;
    Forward(true, lat, 0, x, y, gamma, kold);
    k0_ *= k/kold;
  }

} // namespace GeographicLib
