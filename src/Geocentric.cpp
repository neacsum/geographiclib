/**
 * \file Geocentric.cpp
 * \brief Implementation for GeographicLib::Geocentric class
 *
 * Copyright (c) Charles Karney (2008-2022) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#include <GeographicLib/Geocentric.hpp>

namespace GeographicLib {

  using namespace std;

  Geocentric::Geocentric(real a, real f)
    : a_(a)
    , f_(f)
    , e2_(f_ * (2 - f_))
    , e2m_(Math::sq(1 - f_))    // 1 - e2_
    , e2a_(fabs(e2_))
    , e4a_(Math::sq(e2_))
    , maxrad_(2 * a_ / numeric_limits<real>::epsilon())
  {
    if (!(isfinite(a_) && a_ > 0))
      throw GeographicErr("Equatorial radius is not positive");
    if (!(isfinite(f_) && f_ < 1))
      throw GeographicErr("Polar semi-axis is not positive");
  }

  const Geocentric& Geocentric::WGS84() {
    static const Geocentric wgs84(Constants::WGS84_a(), Constants::WGS84_f());
    return wgs84;
  }

  void Geocentric::IntForward(real lat, real lon, real h,
                              real& X, real& Y, real& Z,
                              real M[dim2_]) const {
    real sphi, cphi, slam, clam;
    Math::sincosd(Math::LatFix(lat), sphi, cphi);
    Math::sincosd(lon, slam, clam);
    real n = a_/sqrt(1 - e2_ * Math::sq(sphi));
    Z = (e2m_ * n + h) * sphi;
    X = (n + h) * cphi;
    Y = X * slam;
    X *= clam;
    if (M)
      Rotation(sphi, cphi, slam, clam, M);
  }

  void Geocentric::IntReverse(real X, real Y, real Z,
                              real& lat, real& lon, real& h,
                              real M[dim2_]) const {
    real
      R = hypot(X, Y),
      slam = R != 0 ? Y / R : 0,
      clam = R != 0 ? X / R : 1;
    h = hypot(R, Z);      // Distance to center of earth
    real sphi, cphi;
    if (h > maxrad_) {
      // We really far away (> 12 million light years); treat the earth as a
      // point and h, above, is an acceptable approximation to the height.
      // This avoids overflow, e.g., in the computation of disc below.  It's
      // possible that h has overflowed to inf; but that's OK.
      //
      // Treat the case X, Y finite, but R overflows to +inf by scaling by 2.
      R = hypot(X/2, Y/2);
      slam = R != 0 ? (Y/2) / R : 0;
      clam = R != 0 ? (X/2) / R : 1;
      real H = hypot(Z/2, R);
      sphi = (Z/2) / H;
      cphi = R / H;
    } else if (e4a_ == 0) {
      // Treat the spherical case.  Dealing with underflow in the general case
      // with e2_ = 0 is difficult.  Origin maps to N pole same as with
      // ellipsoid.
      real H = hypot(h == 0 ? 1 : Z, R);
      sphi = (h == 0 ? 1 : Z) / H;
      cphi = R / H;
      h -= a_;
    } else {
      // Treat prolate spheroids by swapping R and Z here and by switching
      // the arguments to phi = atan2(...) at the end.
      real
        p = Math::sq(R / a_),
        q = e2m_ * Math::sq(Z / a_),
        r = (p + q - e4a_) / 6;
      if (f_ < 0) swap(p, q);
      if ( !(e4a_ * q == 0 && r <= 0) ) {
        real
          // Avoid possible division by zero when r = 0 by multiplying
          // equations for s and t by r^3 and r, resp.
          S = e4a_ * p * q / 4, // S = r^3 * s
          r2 = Math::sq(r),
          r3 = r * r2,
          disc = S * (2 * r3 + S);
        real u = r;
        if (disc >= 0) {
          real T3 = S + r3;
          // Pick the sign on the sqrt to maximize abs(T3).  This minimizes
          // loss of precision due to cancellation.  The result is unchanged
          // because of the way the T is used in definition of u.
          T3 += T3 < 0 ? -sqrt(disc) : sqrt(disc); // T3 = (r * t)^3
          // N.B. cbrt always returns the real root.  cbrt(-8) = -2.
          real T = cbrt(T3); // T = r * t
          // T can be zero; but then r2 / T -> 0.
          u += T + (T != 0 ? r2 / T : 0);
        } else {
          // T is complex, but the way u is defined the result is real.
          real ang = atan2(sqrt(-disc), -(S + r3));
          // There are three possible cube roots.  We choose the root which
          // avoids cancellation.  Note that disc < 0 implies that r < 0.
          u += 2 * r * cos(ang / 3);
        }
        real
          v = sqrt(Math::sq(u) + e4a_ * q), // guaranteed positive
          // Avoid loss of accuracy when u < 0.  Underflow doesn't occur in
          // e4 * q / (v - u) because u ~ e^4 when q is small and u < 0.
          uv = u < 0 ? e4a_ * q / (v - u) : u + v, // u+v, guaranteed positive
          // Need to guard against w going negative due to roundoff in uv - q.
          w = fmax(real(0), e2a_ * (uv - q) / (2 * v)),
          // Rearrange expression for k to avoid loss of accuracy due to
          // subtraction.  Division by 0 not possible because uv > 0, w >= 0.
          k = uv / (sqrt(uv + Math::sq(w)) + w),
          k1 = f_ >= 0 ? k : k - e2_,
          k2 = f_ >= 0 ? k + e2_ : k,
          d = k1 * R / k2,
          H = hypot(Z/k1, R/k2);
        sphi = (Z/k1) / H;
        cphi = (R/k2) / H;
        h = (1 - e2m_/k1) * hypot(d, Z);
      } else {                  // e4 * q == 0 && r <= 0
        // This leads to k = 0 (oblate, equatorial plane) and k + e^2 = 0
        // (prolate, rotation axis) and the generation of 0/0 in the general
        // formulas for phi and h.  using the general formula and division by 0
        // in formula for h.  So handle this case by taking the limits:
        // f > 0: z -> 0, k      ->   e2 * sqrt(q)/sqrt(e4 - p)
        // f < 0: R -> 0, k + e2 -> - e2 * sqrt(q)/sqrt(e4 - p)
        real
          zz = sqrt((f_ >= 0 ? e4a_ - p : p) / e2m_),
          xx = sqrt( f_ <  0 ? e4a_ - p : p        ),
          H = hypot(zz, xx);
        sphi = zz / H;
        cphi = xx / H;
        if (Z < 0) sphi = -sphi; // for tiny negative Z (not for prolate)
        h = - a_ * (f_ >= 0 ? e2m_ : 1) * H / e2a_;
      }
    }
    lat = Math::atan2d(sphi, cphi);
    lon = Math::atan2d(slam, clam);
    if (M)
      Rotation(sphi, cphi, slam, clam, M);
  }

  void Geocentric::Rotation(real sphi, real cphi, real slam, real clam,
                            real M[dim2_]) {
    // This rotation matrix is given by the following quaternion operations
    // qrot(lam, [0,0,1]) * qrot(phi, [0,-1,0]) * [1,1,1,1]/2
    // or
    // qrot(pi/2 + lam, [0,0,1]) * qrot(-pi/2 + phi , [-1,0,0])
    // where
    // qrot(t,v) = [cos(t/2), sin(t/2)*v[1], sin(t/2)*v[2], sin(t/2)*v[3]]

    // Local X axis (east) in geocentric coords
    M[0] = -slam;        M[3] =  clam;        M[6] = 0;
    // Local Y axis (north) in geocentric coords
    M[1] = -clam * sphi; M[4] = -slam * sphi; M[7] = cphi;
    // Local Z axis (up) in geocentric coords
    M[2] =  clam * cphi; M[5] =  slam * cphi; M[8] = sphi;
  }

} // namespace GeographicLib
