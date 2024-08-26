/**
 * \file AlbersEqualArea.cpp
 * \brief Implementation for GeographicLib::AlbersEqualArea class
 *
 * Copyright (c) Charles Karney (2010-2022) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#include <GeographicLib/AlbersEqualArea.hpp>

using namespace std;
using real = GeographicLib::real;

namespace GeographicLib {

static constexpr real eps_ = numeric_limits<real>::epsilon (),
  epsx_ = Math::sq (eps_),
  epsx2_ = Math::sq (epsx_);

static const real tol_ = sqrt (eps_), 
  tol0_ = tol_ * sqrt (sqrt (eps_));

// Divided differences
// Definition: Df(x,y) = (f(x)-f(y))/(x-y)
// See:
//   W. M. Kahan and R. J. Fateman,
//   Symbolic computation of divided differences,
//   SIGSAM Bull. 33(2), 7-28 (1999)
//   https://doi.org/10.1145/334714.334716
//   http://www.cs.berkeley.edu/~fateman/papers/divdiff.pdf
//
// General rules
// h(x) = f(g(x)): Dh(x,y) = Df(g(x),g(y))*Dg(x,y)
// h(x) = f(x)*g(x):
//        Dh(x,y) = Df(x,y)*g(x) + Dg(x,y)*f(y)
//                = Df(x,y)*g(y) + Dg(x,y)*f(x)
//                = Df(x,y)*(g(x)+g(y))/2 + Dg(x,y)*(f(x)+f(y))/2
//
// sn(x) = x/sqrt(1+x^2): Dsn(x,y) = (x+y)/((sn(x)+sn(y))*(1+x^2)*(1+y^2))
static
real Dsn (real x, real y, real sx, real sy) {
  // sx = x/hyp(x)
  real t = x * y;
  return t > 0 ? (x + y) * Math::sq ((sx * sy) / t) / (sx + sy) :
    (x - y != 0 ? (sx - sy) / (x - y) : 1);
}

static const int numit_ = 5;   // Newton iterations in Reverse
static const int numit0_ = 20; // Newton iterations in Init
static real hyp (real x) {
  using std::hypot;
  return hypot (real (1), x);
}

// return atanh(sqrt(x))/sqrt(x) - 1, accurate for small x
static real atanhxm1 (real x);

AlbersEqualArea::AlbersEqualArea(real a, real f, real stdlat, real k0)
  : a_(a)
  , f_(f)
  , fm_(1 - f_)
  , e2_(f_ * (2 - f_))
  , e_(sqrt(fabs(e2_)))
  , e2m_(1 - e2_)
  , qZ_(1 + e2m_ * atanhee(real(1)))
  , qx_(qZ_ / ( 2 * e2m_ ))
{
  if (!(isfinite(a_) && a_ > 0))
    throw GeographicErr("Equatorial radius is not positive");
  if (!(isfinite(f_) && f_ < 1))
    throw GeographicErr("Polar semi-axis is not positive");
  if (!(isfinite(k0) && k0 > 0))
    throw GeographicErr("Scale is not positive");
  if (!(fabs(stdlat) <= 90))
    throw GeographicErr("Standard latitude not in [-90d, +90d]");
  real sphi, cphi;
  Math::sincosd(stdlat, sphi, cphi);
  Init(sphi, cphi, sphi, cphi, k0);
}

AlbersEqualArea::AlbersEqualArea(real a, real f, real stdlat1, real stdlat2,
                                  real k1)
  : a_(a)
  , f_(f)
  , fm_(1 - f_)
  , e2_(f_ * (2 - f_))
  , e_(sqrt(fabs(e2_)))
  , e2m_(1 - e2_)
  , qZ_(1 + e2m_ * atanhee(real(1)))
  , qx_(qZ_ / ( 2 * e2m_ ))
{
  if (!(isfinite(a_) && a_ > 0))
    throw GeographicErr("Equatorial radius is not positive");
  if (!(isfinite(f_) && f_ < 1))
    throw GeographicErr("Polar semi-axis is not positive");
  if (!(isfinite(k1) && k1 > 0))
    throw GeographicErr("Scale is not positive");
  if (!(fabs(stdlat1) <= 90))
    throw GeographicErr("Standard latitude 1 not in [-90d, +90d]");
  if (!(fabs (stdlat2) <= 90))
    throw GeographicErr ("Standard latitude 2 not in [-90d, +90d]");
  real sphi1, cphi1, sphi2, cphi2;
  Math::sincosd(stdlat1, sphi1, cphi1);
  Math::sincosd(stdlat2, sphi2, cphi2);
  Init(sphi1, cphi1, sphi2, cphi2, k1);
}

AlbersEqualArea::AlbersEqualArea(real a, real f,
                                  real sinlat1, real coslat1,
                                  real sinlat2, real coslat2,
                                  real k1)
  : a_(a)
  , f_(f)
  , fm_(1 - f_)
  , e2_(f_ * (2 - f_))
  , e_(sqrt(fabs(e2_)))
  , e2m_(1 - e2_)
  , qZ_(1 + e2m_ * atanhee(real(1)))
  , qx_(qZ_ / ( 2 * e2m_ ))
{
  if (!(isfinite(a_) && a_ > 0))
    throw GeographicErr("Equatorial radius is not positive");
  if (!(isfinite(f_) && f_ < 1))
    throw GeographicErr("Polar semi-axis is not positive");
  if (!(isfinite(k1) && k1 > 0))
    throw GeographicErr("Scale is not positive");
  if (signbit (coslat1))
    throw GeographicErr ("Standard latitude 1 not in [-90, +90d]");
  if (signbit (coslat2))
    throw GeographicErr ("Standard latitude 2 not in [-90d, +90d]");
  if (!(fabs(sinlat1) <= 1 && coslat1 <= 1) || (coslat1 == 0 && sinlat1 == 0))
    throw GeographicErr("Bad sine/cosine of standard latitude 1");
  if (!(fabs(sinlat2) <= 1 && coslat2 <= 1) || (coslat2 == 0 && sinlat2 == 0))
    throw GeographicErr("Bad sine/cosine of standard latitude 2");
  if (coslat1 == 0 && coslat2 == 0 && sinlat1 * sinlat2 <= 0)
    throw GeographicErr
      ("Standard latitudes cannot be opposite poles");
  Init(sinlat1, coslat1, sinlat2, coslat2, k1);
}

void AlbersEqualArea::Init(real sphi1, real cphi1,
                            real sphi2, real cphi2, real k1) {
  {
    real r;
    r = hypot(sphi1, cphi1);
    sphi1 /= r; cphi1 /= r;
    r = hypot(sphi2, cphi2);
    sphi2 /= r; cphi2 /= r;
  }
  bool polar = (cphi1 == 0);
  cphi1 = fmax(epsx_, cphi1);   // Avoid singularities at poles
  cphi2 = fmax(epsx_, cphi2);
  // Determine hemisphere of tangent latitude
  sign_ = sphi1 + sphi2 >= 0 ? 1 : -1;
  // Internally work with tangent latitude positive
  sphi1 *= sign_; sphi2 *= sign_;
  if (sphi1 > sphi2) {
    swap(sphi1, sphi2); swap(cphi1, cphi2); // Make phi1 < phi2
  }
  real
    tphi1 = sphi1/cphi1, tphi2 = sphi2/cphi2;

  // q = (1-e^2)*(sphi/(1-e^2*sphi^2) - atanhee(sphi))
  // qZ = q(pi/2) = (1 + (1-e^2)*atanhee(1))
  // atanhee(x) = atanh(e*x)/e
  // q = sxi * qZ
  // dq/dphi = 2*(1-e^2)*cphi/(1-e^2*sphi^2)^2
  //
  // n = (m1^2-m2^2)/(q2-q1) -> sin(phi0) for phi1, phi2 -> phi0
  // C = m1^2 + n*q1 = (m1^2*q2-m2^2*q1)/(q2-q1)
  // let
  //   rho(pi/2)/rho(-pi/2) = (1-s)/(1+s)
  //   s = n*qZ/C
  //     = qZ * (m1^2-m2^2)/(m1^2*q2-m2^2*q1)
  //     = qZ * (scbet2^2 - scbet1^2)/(scbet2^2*q2 - scbet1^2*q1)
  //     = (scbet2^2 - scbet1^2)/(scbet2^2*sxi2 - scbet1^2*sxi1)
  //     = (tbet2^2 - tbet1^2)/(scbet2^2*sxi2 - scbet1^2*sxi1)
  // 1-s = -((1-sxi2)*scbet2^2 - (1-sxi1)*scbet1^2)/
  //         (scbet2^2*sxi2 - scbet1^2*sxi1)
  //
  // Define phi0 to give same value of s, i.e.,
  //  s = sphi0 * qZ / (m0^2 + sphi0*q0)
  //    = sphi0 * scbet0^2 / (1/qZ + sphi0 * scbet0^2 * sxi0)

  real tphi0, C;
  if (polar || tphi1 == tphi2) {
    tphi0 = tphi2;
    C = 1;                    // ignored
  } else {
    real
      tbet1 = fm_ * tphi1, scbet12 = 1 + Math::sq(tbet1),
      tbet2 = fm_ * tphi2, scbet22 = 1 + Math::sq(tbet2),
      txi1 = txif(tphi1), cxi1 = 1/hyp(txi1), sxi1 = txi1 * cxi1,
      txi2 = txif(tphi2), cxi2 = 1/hyp(txi2), sxi2 = txi2 * cxi2,
      dtbet2 = fm_ * (tbet1 + tbet2),
      es1 = 1 - e2_ * Math::sq(sphi1), es2 = 1 - e2_ * Math::sq(sphi2),
      /*
      dsxi = ( (e2_ * sq(sphi2 + sphi1) + es2 + es1) / (2 * es2 * es1) +
                Datanhee(sphi2, sphi1) ) * Dsn(tphi2, tphi1, sphi2, sphi1) /
      ( 2 * qx_ ),
      */
      dsxi = ( (1 + e2_ * sphi1 * sphi2) / (es2 * es1) +
                Datanhee(sphi2, sphi1) ) * Dsn(tphi2, tphi1, sphi2, sphi1) /
      ( 2 * qx_ ),
      den = (sxi2 + sxi1) * dtbet2 + (scbet22 + scbet12) * dsxi,
      // s = (sq(tbet2) - sq(tbet1)) / (scbet22*sxi2 - scbet12*sxi1)
      s = 2 * dtbet2 / den,
      // 1-s = -(sq(scbet2)*(1-sxi2) - sq(scbet1)*(1-sxi1)) /
      //        (scbet22*sxi2 - scbet12*sxi1)
      // Write
      //   sq(scbet)*(1-sxi) = sq(scbet)*(1-sphi) * (1-sxi)/(1-sphi)
      sm1 = -Dsn(tphi2, tphi1, sphi2, sphi1) *
      ( -( ((sphi2 <= 0 ? (1 - sxi2) / (1 - sphi2) :
              Math::sq(cxi2/cphi2) * (1 + sphi2) / (1 + sxi2)) +
            (sphi1 <= 0 ? (1 - sxi1) / (1 - sphi1) :
              Math::sq(cxi1/cphi1) * (1 + sphi1) / (1 + sxi1))) ) *
        (1 + e2_ * (sphi1 + sphi2 + sphi1 * sphi2)) /
        (1 +       (sphi1 + sphi2 + sphi1 * sphi2)) +
        (scbet22 * (sphi2 <= 0 ? 1 - sphi2 :
                    Math::sq(cphi2) / ( 1 + sphi2)) +
          scbet12 * (sphi1 <= 0 ? 1 - sphi1 : Math::sq(cphi1) / ( 1 + sphi1)))
        * (e2_ * (1 + sphi1 + sphi2 + e2_ * sphi1 * sphi2)/(es1 * es2)
        +e2m_ * DDatanhee(sphi1, sphi2) ) / qZ_ ) / den;
    // C = (scbet22*sxi2 - scbet12*sxi1) / (scbet22 * scbet12 * (sx2 - sx1))
    C = den / (2 * scbet12 * scbet22 * dsxi);
    tphi0 = (tphi2 + tphi1)/2;
    real stol = tol0_ * fmax(real(1), fabs(tphi0));
    for (int i = 0;
          i < 2*numit0_ ||
            GEOGRAPHICLIB_PANIC("Convergence failure in AlbersEqualArea");
          ++i) {
      // Solve (scbet0^2 * sphi0) / (1/qZ + scbet0^2 * sphi0 * sxi0) = s
      // for tphi0 by Newton's method on
      // v(tphi0) = (scbet0^2 * sphi0) - s * (1/qZ + scbet0^2 * sphi0 * sxi0)
      //          = 0
      // Alt:
      // (scbet0^2 * sphi0) / (1/qZ - scbet0^2 * sphi0 * (1-sxi0))
      //          = s / (1-s)
      // w(tphi0) = (1-s) * (scbet0^2 * sphi0)
      //             - s  * (1/qZ - scbet0^2 * sphi0 * (1-sxi0))
      //          = (1-s) * (scbet0^2 * sphi0)
      //             - S/qZ  * (1 - scbet0^2 * sphi0 * (qZ-q0))
      // Now
      // qZ-q0 = (1+e2*sphi0)*(1-sphi0)/(1-e2*sphi0^2) +
      //         (1-e2)*atanhee((1-sphi0)/(1-e2*sphi0))
      // In limit sphi0 -> 1, qZ-q0 -> 2*(1-sphi0)/(1-e2), so wrte
      // qZ-q0 = 2*(1-sphi0)/(1-e2) + A + B
      // A = (1-sphi0)*( (1+e2*sphi0)/(1-e2*sphi0^2) - (1+e2)/(1-e2) )
      //   = -e2 *(1-sphi0)^2 * (2+(1+e2)*sphi0) / ((1-e2)*(1-e2*sphi0^2))
      // B = (1-e2)*atanhee((1-sphi0)/(1-e2*sphi0)) - (1-sphi0)
      //   = (1-sphi0)*(1-e2)/(1-e2*sphi0)*
      //     ((atanhee(x)/x-1) - e2*(1-sphi0)/(1-e2))
      // x = (1-sphi0)/(1-e2*sphi0), atanhee(x)/x = atanh(e*x)/(e*x)
      //
      // 1 - scbet0^2 * sphi0 * (qZ-q0)
      //   = 1 - scbet0^2 * sphi0 * (2*(1-sphi0)/(1-e2) + A + B)
      //   = D - scbet0^2 * sphi0 * (A + B)
      // D = 1 - scbet0^2 * sphi0 * 2*(1-sphi0)/(1-e2)
      //   = (1-sphi0)*(1-e2*(1+2*sphi0*(1+sphi0)))/((1-e2)*(1+sphi0))
      // dD/dsphi0 = -2*(1-e2*sphi0^2*(2*sphi0+3))/((1-e2)*(1+sphi0)^2)
      // d(A+B)/dsphi0 = 2*(1-sphi0^2)*e2*(2-e2*(1+sphi0^2))/
      //                 ((1-e2)*(1-e2*sphi0^2)^2)

      real
        scphi02 = 1 + Math::sq(tphi0), scphi0 = sqrt(scphi02),
        // sphi0m = 1-sin(phi0) = 1/( sec(phi0) * (tan(phi0) + sec(phi0)) )
        sphi0 = tphi0 / scphi0, sphi0m = 1/(scphi0 * (tphi0 + scphi0)),
        // scbet0^2 * sphi0
        g = (1 + Math::sq( fm_ * tphi0 )) * sphi0,
        // dg/dsphi0 = dg/dtphi0 * scphi0^3
        dg = e2m_ * scphi02 * (1 + 2 * Math::sq(tphi0)) + e2_,
        D = sphi0m * (1 - e2_*(1 + 2*sphi0*(1+sphi0))) / (e2m_ * (1+sphi0)),
        // dD/dsphi0
        dD = -2 * (1 - e2_*Math::sq(sphi0) * (2*sphi0+3)) /
              (e2m_ * Math::sq(1+sphi0)),
        A = -e2_ * Math::sq(sphi0m) * (2+(1+e2_)*sphi0) /
            (e2m_*(1-e2_*Math::sq(sphi0))),
        B = (sphi0m * e2m_ / (1 - e2_*sphi0) *
              (atanhxm1(e2_ *
                        Math::sq(sphi0m / (1-e2_*sphi0))) - e2_*sphi0m/e2m_)),
        // d(A+B)/dsphi0
        dAB = (2 * e2_ * (2 - e2_ * (1 + Math::sq(sphi0))) /
                (e2m_ * Math::sq(1 - e2_*Math::sq(sphi0)) * scphi02)),
        u = sm1 * g - s/qZ_ * ( D - g * (A + B) ),
        // du/dsphi0
        du = sm1 * dg - s/qZ_ * (dD - dg * (A + B) - g * dAB),
        dtu = -u/du * (scphi0 * scphi02);
      tphi0 += dtu;
      if (!(fabs(dtu) >= stol))
        break;
    }
  }
  txi0_ = txif(tphi0); scxi0_ = hyp(txi0_); sxi0_ = txi0_ / scxi0_;
  n0_ = tphi0/hyp(tphi0);
  m02_ = 1 / (1 + Math::sq(fm_ * tphi0));
  nrho0_ = polar ? 0 : a_ * sqrt(m02_);
  k0_ = sqrt(tphi1 == tphi2 ? 1 : C / (m02_ + n0_ * qZ_ * sxi0_)) * k1;
  k2_ = Math::sq(k0_);
  lat0_ = sign_ * atan(tphi0)/Math::degree();
}

const AlbersEqualArea& AlbersEqualArea::CylindricalEqualArea() {
  static const AlbersEqualArea
    cylindricalequalarea(Constants::WGS84_a(), Constants::WGS84_f(),
                          real(0), real(1), real(0), real(1), real(1));
  return cylindricalequalarea;
}

const AlbersEqualArea& AlbersEqualArea::AzimuthalEqualAreaNorth() {
  static const AlbersEqualArea
    azimuthalequalareanorth(Constants::WGS84_a(), Constants::WGS84_f(),
                            real(1), real(0), real(1), real(0), real(1));
  return azimuthalequalareanorth;
}

const AlbersEqualArea& AlbersEqualArea::AzimuthalEqualAreaSouth() {
  static const AlbersEqualArea
    azimuthalequalareasouth(Constants::WGS84_a(), Constants::WGS84_f(),
                            real(-1), real(0), real(-1), real(0), real(1));
  return azimuthalequalareasouth;
}

real AlbersEqualArea::txif(real tphi) const {
  // sxi = ( sphi/(1-e2*sphi^2) + atanhee(sphi) ) /
  //       ( 1/(1-e2) + atanhee(1) )
  //
  // txi = ( sphi/(1-e2*sphi^2) + atanhee(sphi) ) /
  //       sqrt( ( (1+e2*sphi)*(1-sphi)/( (1-e2*sphi^2) * (1-e2) ) +
  //               atanhee((1-sphi)/(1-e2*sphi)) ) *
  //             ( (1-e2*sphi)*(1+sphi)/( (1-e2*sphi^2) * (1-e2) ) +
  //               atanhee((1+sphi)/(1+e2*sphi)) ) )
  //     = ( tphi/(1-e2*sphi^2) + atanhee(sphi, e2)/cphi ) /
  //       sqrt(
  //       ( (1+e2*sphi)/( (1-e2*sphi^2) * (1-e2) ) + Datanhee(1,  sphi)  ) *
  //       ( (1-e2*sphi)/( (1-e2*sphi^2) * (1-e2) ) + Datanhee(1, -sphi)  ) )
  //
  // This function maintains odd parity
  real
    cphi = 1 / sqrt(1 + Math::sq(tphi)),
    sphi = tphi * cphi,
    es1 = e2_ * sphi,
    es2m1 = 1 - es1 * sphi,   // 1 - e2 * sphi^2
    es2m1a = e2m_ * es2m1;    // (1 - e2 * sphi^2) * (1 - e2)
  return ( tphi / es2m1 + atanhee(sphi) / cphi ) /
    sqrt( ( (1 + es1) / es2m1a + Datanhee(1,  sphi) ) *
          ( (1 - es1) / es2m1a + Datanhee(1, -sphi) ) );
}

real AlbersEqualArea::tphif(real txi) const {
  real
    tphi = txi,
    stol = tol_ * fmax(real(1), fabs(txi));
  // CHECK: min iterations = 1, max iterations = 2; mean = 1.99
  for (int i = 0;
        i < numit_ ||
          GEOGRAPHICLIB_PANIC("Convergence failure in AlbersEqualArea");
        ++i) {
    // dtxi/dtphi = (scxi/scphi)^3 * 2*(1-e^2)/(qZ*(1-e^2*sphi^2)^2)
    real
      txia = txif(tphi),
      tphi2 = Math::sq(tphi),
      scphi2 = 1 + tphi2,
      scterm = scphi2/(1 + Math::sq(txia)),
      dtphi = (txi - txia) * scterm * sqrt(scterm) *
      qx_ * Math::sq(1 - e2_ * tphi2 / scphi2);
    tphi += dtphi;
    if (!(fabs(dtphi) >= stol))
      break;
  }
  return tphi;
}

// return atanh(sqrt(x))/sqrt(x) - 1 = x/3 + x^2/5 + x^3/7 + ...
// typical x < e^2 = 2*f
real atanhxm1(real x) {
  real s = 0;
  if (fabs(x) < real(0.5)) {
    static const real lg2eps_ = -log2(numeric_limits<real>::epsilon() / 2);
    int e;
    (void) frexp(x, &e);
    e = -e;
    // x = [0.5,1) * 2^(-e)
    // estimate n s.t. x^n/(2*n+1) < x/3 * epsilon/2
    // a stronger condition is x^(n-1) < epsilon/2
    // taking log2 of both sides, a stronger condition is
    // (n-1)*(-e) < -lg2eps or (n-1)*e > lg2eps or n > ceiling(lg2eps/e)+1
    int n = x == 0 ? 1 : int(ceil(lg2eps_ / e)) + 1;
    while (n--)               // iterating from n-1 down to 0
      s = x * s + (n ? 1 : 0)/real(2*n + 1);
  } else {
    real xs = sqrt(fabs(x));
    s = (x > 0 ? atanh(xs) : atan(xs)) / xs - 1;
  }
  return s;
}

// return (Datanhee(1,y) - Datanhee(1,x))/(y-x)
real AlbersEqualArea::DDatanhee(real x, real y) const {
  // This function is called with x = sphi1, y = sphi2, phi1 <= phi2, sphi2
  // >= 0, abs(sphi1) <= phi2.  However for safety's sake we enforce x <= y.
  if (y < x) swap(x, y);      // ensure that x <= y
  real q1 = fabs(e2_),
    q2 = fabs(2 * e_ / e2m_ * (1 - x));
  return
    x <= 0 || !(fmin(q1, q2) < real(0.75)) ? DDatanhee0(x, y) :
    (q1 < q2 ? DDatanhee1(x, y) : DDatanhee2(x, y));
}

// Rearrange difference so that 1 - x is in the denominator, then do a
// straight divided difference.
real AlbersEqualArea::DDatanhee0(real x, real y) const {
  return (Datanhee(1, y) - Datanhee(x, y))/(1 - x);
}

// The expansion for e2 small
real AlbersEqualArea::DDatanhee1(real x, real y) const {
  // The series in e2 is
  //   sum( c[l] * e2^l, l, 1, N)
  // where
  //   c[l] = sum( x^i * y^j; i >= 0, j >= 0, i+j < 2*l) / (2*l + 1)
  //        = ( (x-y) - (1-y) * x^(2*l+1) + (1-x) * y^(2*l+1) ) /
  //          ( (2*l+1) * (x-y) * (1-y) * (1-x) )
  // For x = y = 1, c[l] = l
  //
  // In the limit x,y -> 1,
  //
  //   DDatanhee -> e2/(1-e2)^2 = sum(l * e2^l, l, 1, inf)
  //
  // Use if e2 is sufficiently small.
  real s = 0;
  real  z = 1, k = 1, t = 0, c = 0, en = 1;
  while (true) {
    t = y * t + z; c += t; z *= x;
    t = y * t + z; c += t; z *= x;
    k += 2; en *= e2_;
    // Here en[l] = e2^l, k[l] = 2*l + 1,
    // c[l] = sum( x^i * y^j; i >= 0, j >= 0, i+j < 2*l) / (2*l + 1)
    // Taylor expansion is
    // s = sum( c[l] * e2^l, l, 1, N)
    real ds = en * c / k;
    s += ds;
    if (!(fabs(ds) > fabs(s) * eps_/2))
      break;            // Iterate until the added term is sufficiently small
  }
  return s;
}

// The expansion for x (and y) close to 1
real AlbersEqualArea::DDatanhee2(real x, real y) const {
  // If x and y are both close to 1, expand in Taylor series in dx = 1-x and
  // dy = 1-y:
  //
  // DDatanhee = sum(C_m * (dx^(m+1) - dy^(m+1)) / (dx - dy), m, 0, inf)
  //
  // where
  //
  // C_m = sum( (m+2)!! / (m+2-2*k)!! *
  //            ((m+1)/2)! / ((m+1)/2-k)! /
  //            (k! * (2*k-1)!!) *
  //            e2^((m+1)/2+k),
  //           k, 0, (m+1)/2) * (-1)^m / ((m+2) * (1-e2)^(m+2))
  // for m odd, and
  //
  // C_m = sum( 2 * (m+1)!! / (m+1-2*k)!! *
  //            (m/2+1)! / (m/2-k)! /
  //            (k! * (2*k+1)!!) *
  //            e2^(m/2+1+k),
  //           k, 0, m/2)) * (-1)^m / ((m+2) * (1-e2)^(m+2))
  // for m even.
  //
  // Here i!! is the double factorial extended to negative i with
  // i!! = (i+2)!!/(i+2).
  //
  // Note that
  //   (dx^(m+1) - dy^(m+1)) / (dx - dy) =
  //     dx^m + dx^(m-1)*dy ... + dx*dy^(m-1) + dy^m
  //
  // Leading (m = 0) term is e2 / (1 - e2)^2
  //
  // Magnitude of mth term relative to the leading term scales as
  //
  //   2*(2*e/(1-e2)*dx)^m
  //
  // So use series if (2*e/(1-e2)*dx) is sufficiently small
  real s, dx = 1 - x, dy = 1 - y, xy = 1, yy = 1, ee = e2_ / Math::sq(e2m_);
  s = ee;
  for (int m = 1; ; ++m) {
      real c = m + 2, t = c;
    yy *= dy;               // yy = dy^m
    xy = dx * xy + yy;
    // Now xy = dx^m + dx^(m-1)*dy ... + dx*dy^(m-1) + dy^m
    //        = (dx^(m+1) - dy^(m+1)) / (dx - dy)
    // max value = (m+1) * max(dx,dy)^m
    ee /= -e2m_;
    if (m % 2 == 0) ee *= e2_;
    // Now ee = (-1)^m * e2^(floor(m/2)+1) / (1-e2)^(m+2)
    int kmax = (m+1)/2;
    for (int k = kmax - 1; k >= 0; --k) {
      // max coeff is less than 2^(m+1)
      c *= (k + 1) * (2 * (k + m - 2*kmax) + 3);
      c /= (kmax - k) * (2 * (kmax - k) + 1);
      // Horner sum for inner e2_ series
      t = e2_ * t + c;
    }
    // Straight sum for outer m series
    real ds = t * ee * xy / (m + 2);
    s = s + ds;
    if (!(fabs(ds) > fabs(s) * eps_/2))
      break;            // Iterate until the added term is sufficiently small
  }
  return s;
}

void AlbersEqualArea::Forward(real lon0, real lat, real lon,
                              real& x, real& y, real& gamma, real& k) const {
  lon = Math::AngDiff(lon0, lon);
  lat *= sign_;
  real sphi, cphi;
  Math::sincosd(Math::LatFix(lat) * sign_, sphi, cphi);
  cphi = fmax(epsx_, cphi);
  real
    lam = lon * Math::degree(),
    tphi = sphi/cphi, txi = txif(tphi), sxi = txi/hyp(txi),
    dq = qZ_ * Dsn(txi, txi0_, sxi, sxi0_) * (txi - txi0_),
    drho = - a_ * dq / (sqrt(m02_ - n0_ * dq) + nrho0_ / a_),
    theta = k2_ * n0_ * lam, stheta = sin(theta), ctheta = cos(theta),
    t = nrho0_ + n0_ * drho;
  x = t * (n0_ != 0 ? stheta / n0_ : k2_ * lam) / k0_;
  y = (nrho0_ *
        (n0_ != 0 ?
        (ctheta < 0 ? 1 - ctheta : Math::sq(stheta)/(1 + ctheta)) / n0_ :
        0)
        - drho * ctheta) / k0_;
  k = k0_ * (t != 0 ? t * hyp(fm_ * tphi) / a_ : 1);
  y *= sign_;
  gamma = sign_ * theta / Math::degree();
}

void AlbersEqualArea::Reverse(real lon0, real x, real y,
                              real& lat, real& lon,
                              real& gamma, real& k) const {
  y *= sign_;
  real
    nx = k0_ * n0_ * x, ny = k0_ * n0_ * y, y1 =  nrho0_ - ny,
    den = hypot(nx, y1) + nrho0_, // 0 implies origin with polar aspect
    drho = den != 0 ? (k0_*x*nx - 2*k0_*y*nrho0_ + k0_*y*ny) / den : 0,
    // dsxia = scxi0 * dsxi
    dsxia = - scxi0_ * (2 * nrho0_ + n0_ * drho) * drho /
            (Math::sq(a_) * qZ_),
    txi = (txi0_ + dsxia) / sqrt(fmax(1 - dsxia * (2*txi0_ + dsxia), epsx2_)),
    tphi = tphif(txi),
    theta = atan2(nx, y1),
    lam = n0_ != 0 ? theta / (k2_ * n0_) : x / (y1 * k0_);
  gamma = sign_ * theta / Math::degree();
  lat = Math::atand(sign_ * tphi);
  lon = lam / Math::degree();
  lon = Math::AngNormalize(lon + Math::AngNormalize(lon0));
  k = k0_ * (den != 0 ? (nrho0_ + n0_ * drho) * hyp(fm_ * tphi) / a_ : 1);
}

void AlbersEqualArea::SetScale(real lat, real k) {
  if (!(isfinite(k) && k > 0))
    throw GeographicErr("Scale is not positive");
  if (!(fabs (lat) < 90))
    throw GeographicErr ("Latitude for SetScale not in (-90d, +90d)");
  real x, y, gamma, kold;
  Forward(0, lat, 0, x, y, gamma, kold);
  k /= kold;
  k0_ *= k;
  k2_ = Math::sq(k0_);
}

// atanh(      e   * x)/      e   if f > 0
// atan (sqrt(-e2) * x)/sqrt(-e2) if f < 0
// x                              if f = 0
real AlbersEqualArea::atanhee (real x) const {
  using std::atan; using std::atanh;
  return f_ > 0 ? atanh (e_ * x) / e_ : (f_ < 0 ? (atan (e_ * x) / e_) : x);
}

// Datanhee(x,y) = (atanee(x)-atanee(y))/(x-y)
//               = atanhee((x-y)/(1-e^2*x*y))/(x-y)
real AlbersEqualArea::Datanhee (real x, real y) const {
  real t = x - y, d = 1 - e2_ * x * y;
  return t == 0 ? 1 / d :
    (x * y < 0 ? atanhee (x) - atanhee (y) : atanhee (t / d)) / t;
}

} // namespace GeographicLib
