/**
 * \file Math.hpp
 * \brief Header for GeographicLib::Math namespace
 *
 * Copyright (c) Charles Karney (2008-2024) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

// Constants.hpp includes Math.hpp.  Place this include outside Math.hpp's
// include guard to enforce this ordering.
#if !defined (GEOGRAPHICLIB_CONSTANTS_HPP)
#include <GeographicLib/Constants.hpp>
#endif

#if !defined(GEOGRAPHICLIB_MATH_HPP)
#define GEOGRAPHICLIB_MATH_HPP 1

#if !defined(GEOGRAPHICLIB_WORDS_BIGENDIAN)
#  define GEOGRAPHICLIB_WORDS_BIGENDIAN 0
#endif

#if !defined(GEOGRAPHICLIB_HAVE_LONG_DOUBLE)
#  define GEOGRAPHICLIB_HAVE_LONG_DOUBLE 0
#endif

#if !defined(GEOGRAPHICLIB_PRECISION)
/**
 * The precision of floating point numbers used in %GeographicLib.  1 means
 * float (single precision); 2 (the default) means double; 3 means long double;
 * 4 is reserved for quadruple precision.  Nearly all the testing has been
 * carried out with doubles and that's the recommended configuration.  In order
 * for long double to be used, GEOGRAPHICLIB_HAVE_LONG_DOUBLE needs to be
 * defined.  Note that with Microsoft Visual Studio, long double is the same as
 * double.
 **********************************************************************/
#  define GEOGRAPHICLIB_PRECISION 2
#endif

#include <cmath>
#include <limits>

//make sure we can use std::min and std::max 
#undef min
#undef max
#include <algorithm>

#if GEOGRAPHICLIB_PRECISION == 4
#include <memory>
#include <boost/version.hpp>
#include <boost/multiprecision/float128.hpp>
#include <boost/math/special_functions.hpp>
#elif GEOGRAPHICLIB_PRECISION == 5
#include <mpreal.h>
#endif

#if GEOGRAPHICLIB_PRECISION > 3
// volatile keyword makes no sense for multiprec types
#define GEOGRAPHICLIB_VOLATILE
// Signal a convergence failure with multiprec types by throwing an exception
// at loop exit.
#define GEOGRAPHICLIB_PANIC(msg) \
  (throw GeographicLib::GeographicErr(msg), false)
#else
#define GEOGRAPHICLIB_VOLATILE volatile
// Ignore convergence failures with standard floating points types by allowing
// loop to exit cleanly.
#define GEOGRAPHICLIB_PANIC(msg) false
#endif

namespace GeographicLib
{

#if GEOGRAPHICLIB_HAVE_LONG_DOUBLE
  /**
    * The extended precision type for real numbers, used for some testing.
    * This is long double on computers with this type; otherwise it is double.
    **********************************************************************/
  typedef long double extended;
#else
  typedef double extended;
#endif

#if GEOGRAPHICLIB_PRECISION == 2
  /**
    * The real type for %GeographicLib. Nearly all the testing has been done
    * with \e real = double.  However, the algorithms should also work with
    * float and long double (where available).  (<b>CAUTION</b>: reasonable
    * accuracy typically cannot be obtained using floats.)
    **********************************************************************/
  typedef double real;
#elif GEOGRAPHICLIB_PRECISION == 1
  typedef float real;
#elif GEOGRAPHICLIB_PRECISION == 3
  typedef extended real;
#elif GEOGRAPHICLIB_PRECISION == 4
  typedef boost::multiprecision::float128 real;
#elif GEOGRAPHICLIB_PRECISION == 5
  typedef mpfr::mpreal real;
#else
  typedef double real;
#endif

/**
  * \brief Mathematical functions needed by %GeographicLib
  *
  * Define mathematical functions in order to localize system dependencies and
  * to provide generic versions of the functions.  In addition define a real
  * type to be used by %GeographicLib.
  *
  * Example of use:
  * \include examples/Math.cpp
  **********************************************************************/
namespace Math
{

  /**
   * The NaN (not a number)
   *
   * @tparam T the type of the returned value.
   * @return NaN if available, otherwise return the max real of type T.
   **********************************************************************/
  template<typename T = real> constexpr T NaN ();

  /**
   * Infinity
   *
   * @tparam T the type of the returned value.
   * @return infinity if available, otherwise return the max real.
   **********************************************************************/
  template<typename T = real> constexpr T infinity ();

  /**
   * @return the number of bits of precision in a real number.
   **********************************************************************/
  constexpr int digits();

  /**
   * Set the binary precision of a real number.
   *
   * @param[in] ndigits the number of bits of precision.
   * @return the resulting number of bits of precision.
   *
   * This only has an effect when GEOGRAPHICLIB_PRECISION = 5.  See also
   * Utility::set_digits for caveats about when this routine should be
   * called.
   **********************************************************************/
  inline int set_digits (int ndigits) {
#if GEOGRAPHICLIB_PRECISION != 5
    (void)ndigits;
#else
    mpfr::mpreal::set_default_prec (ndigits >= 2 ? ndigits : 2);
#endif
    return digits ();
  }

  /**
   * @return the number of decimal digits of precision in a real number.
   **********************************************************************/
  constexpr int digits10();

  /**
   * Number of additional decimal digits of precision for real relative to
   * double (0 for float).
   **********************************************************************/
  inline constexpr int extra_digits () {
    return
      digits10 () > std::numeric_limits<double>::digits10 ?
      digits10 () - std::numeric_limits<double>::digits10 : 0;
  }

  /**
   * true if the machine is big-endian.
   **********************************************************************/
  const bool bigendian = GEOGRAPHICLIB_WORDS_BIGENDIAN;

  /**
   * @tparam T the type of the returned value.
   * @return &pi;.
   **********************************************************************/
  template<typename T = real> static T pi() {
    using std::atan2;
    static const T pi = atan2(T(0), T(-1));
    return pi;
  }

  /**
   * @tparam T the type of the returned value.
   * @return the number of radians in a degree.
   **********************************************************************/
  template<typename T = real> T degree() {
    static const T degree = pi<T>() / T(180);
    return degree;
  }

  /**
   * Square a number.
   *
   * @tparam T the type of the argument and the returned value.
   * @param[in] x
   * @return <i>x</i><sup>2</sup>.
   **********************************************************************/
  template<typename T> T constexpr sq(T x)
  { return x * x; }

  /**
   * Normalize a two-vector.
   *
   * @tparam T the type of the argument and the returned value.
   * @param[in,out] x on output set to <i>x</i>/hypot(<i>x</i>, <i>y</i>).
   * @param[in,out] y on output set to <i>y</i>/hypot(<i>x</i>, <i>y</i>).
   **********************************************************************/
  template<typename T> void norm(T& x, T& y) {
#if defined(_MSC_VER) && defined(_M_IX86)
    // hypot for Visual Studio (A=win32) fails monotonicity, e.g., with
    //   x  = 0.6102683302836215
    //   y1 = 0.7906090004346522
    //   y2 = y1 + 1e-16
    // the test
    //   hypot(x, y2) >= hypot(x, y1)
    // fails.  Reported 2021-03-14:
    //   https://developercommunity.visualstudio.com/t/1369259
    // See also:
    //   https://bugs.python.org/issue43088
    using std::sqrt; T h = sqrt(x * x + y * y);
#else
    using std::hypot; T h = hypot(x, y);
#endif
    x /= h; y /= h;
  }

 /**
  * The error-free sum of two numbers.
  *
  * @tparam T the type of the argument and the returned value.
  * @param[in] u
  * @param[in] v
  * @param[out] t the exact error given by (\e u + \e v) - \e s.
  * @return \e s = round(\e u + \e v).
  *
  * See D. E. Knuth, TAOCP, Vol 2, 4.2.2, Theorem B.
  *
  * \note \e t can be the same as one of the first two arguments.
  **********************************************************************/
  template<typename T> T sum (T u, T v, T& t) {
    GEOGRAPHICLIB_VOLATILE T s = u + v;
    GEOGRAPHICLIB_VOLATILE T up = s - v;
    GEOGRAPHICLIB_VOLATILE T vpp = s - up;
    up -= u;
    vpp -= v;
    // if s = 0, then t = 0 and give t the same sign as s
    // mpreal needs T(0) here
    t = s != 0 ? T (0) - (up + vpp) : s;
    // u + v =       s      + t
    //       = round(u + v) + t
    return s;
  }

  /**
   * Evaluate a polynomial.
   *
   * @tparam T the type of the arguments and returned value.
   * @param[in] N the order of the polynomial.
   * @param[in] p the coefficient array (of size \e N + 1) with
   *   <i>p</i><sub>0</sub> being coefficient of <i>x</i><sup><i>N</i></sup>.
   * @param[in] x the variable.
   * @return the value of the polynomial.
   *
   * Evaluate &sum;<sub><i>n</i>=0..<i>N</i></sub>
   * <i>p</i><sub><i>n</i></sub> <i>x</i><sup><i>N</i>&minus;<i>n</i></sup>.
   * Return 0 if \e N &lt; 0.  Return <i>p</i><sub>0</sub>, if \e N = 0 (even
   * if \e x is infinite or a nan).  The evaluation uses Horner's method.
   **********************************************************************/
  template<typename T> T polyval(int N, const T p[], T x) {
    // This used to employ Math::fma; but that's too slow and it seemed not
    // to improve the accuracy noticeably.  This might change when there's
    // direct hardware support for fma.
    T y = N < 0 ? 0 : *p++;
    while (--N >= 0) y = y * x + *p++;
    return y;
  }

  /**
   * Normalize an angle.
   *
   * @tparam T the type of the argument and returned value.
   * @param[in] x the angle in degrees.
   * @return the angle reduced to the range [&minus;180&deg;, 180&deg;].
   *
   * The range of \e x is unrestricted.  If the result is &plusmn;0&deg; or
   * &plusmn;180&deg; then the sign is the sign of \e x.
   **********************************************************************/
  template<typename T> T AngNormalize(T x) {
    T y = remainder (x, T (360));
#if GEOGRAPHICLIB_PRECISION == 4
    // boost-quadmath doesn't set the sign of 0 correctly, see
    // https://github.com/boostorg/multiprecision/issues/426
    // Fixed by https://github.com/boostorg/multiprecision/pull/428
    if (y == 0) y = copysign (y, x);
#endif
    return fabs (y) == T (180) ? copysign (T (180), x) : y;
  }

  /**
   * Normalize a latitude.
   *
   * @tparam T the type of the argument and returned value.
   * @param[in] x the angle in degrees.
   * @return x if it is in the range [&minus;90&deg;, 90&deg;], otherwise
   *   return NaN.
   **********************************************************************/
  template<typename T> T LatFix(T x)
  { using std::fabs; return fabs(x) > T(90) ? NaN<T>() : x; }

  /**
   * The exact difference of two angles reduced to
   * [&minus;180&deg;, 180&deg;].
   *
   * @tparam T the type of the arguments and returned value.
   * @param[in] x the first angle in degrees.
   * @param[in] y the second angle in degrees.
   * @param[out] e the error term in degrees.
   * @return \e d, the truncated value of \e y &minus; \e x.
   *
   * This computes \e z = \e y &minus; \e x exactly, reduced to
   * [&minus;180&deg;, 180&deg;]; and then sets \e z = \e d + \e e where \e d
   * is the nearest representable number to \e z and \e e is the truncation
   * error.  If \e z = &plusmn;0&deg; or &plusmn;180&deg;, then the sign of
   * \e d is given by the sign of \e y &minus; \e x.  The maximum absolute
   * value of \e e is 2<sup>&minus;26</sup> (for doubles).
   **********************************************************************/
  template<typename T> T AngDiff(T x, T y, T& e) {
    // Use remainder instead of AngNormalize, since we treat boundary cases
    // later taking account of the error
    T d = sum (remainder (-x, T (360)), remainder (y, T (360)), e);
    // This second sum can only change d if abs(d) < 128, so don't need to
    // apply remainder yet again.
    d = sum (remainder (d, T (360)), e, e);
    // Fix the sign if d = -180, 0, 180.
    if (d == 0 || fabs (d) == 180)
      // If e == 0, take sign from y - x
      // else (e != 0, implies d = +/-180), d and e must have opposite signs
      d = copysign (d, e == 0 ? y - x : -e);
    return d;
  }

  /**
   * Difference of two angles reduced to [&minus;180&deg;, 180&deg;]
   *
   * @tparam T the type of the arguments and returned value.
   * @param[in] x the first angle in degrees.
   * @param[in] y the second angle in degrees.
   * @return \e y &minus; \e x, reduced to the range [&minus;180&deg;,
   *   180&deg;].
   *
   * The result is equivalent to computing the difference exactly, reducing
   * it to [&minus;180&deg;, 180&deg;] and rounding the result.
   **********************************************************************/
  template<typename T> T AngDiff(T x, T y)
  { T e; return AngDiff(x, y, e); }

  /**
   * Coarsen a value close to zero.
   *
   * @tparam T the type of the argument and returned value.
   * @param[in] x
   * @return the coarsened value.
   *
   * The makes the smallest gap in \e x = 1/16 &minus; nextafter(1/16, 0) =
   * 1/2<sup>57</sup> for doubles = 0.8 pm on the earth if \e x is an angle
   * in degrees.  (This is about 2000 times more resolution than we get with
   * angles around 90&deg;.)  We use this to avoid having to deal with near
   * singular cases when \e x is non-zero but tiny (e.g.,
   * 10<sup>&minus;200</sup>).  This sign of &plusmn;0 is preserved.
   **********************************************************************/
  template<typename T> T AngRound(T x) {
    static const T z = T (1) / T (16);
    GEOGRAPHICLIB_VOLATILE T y = fabs (x);
    GEOGRAPHICLIB_VOLATILE T w = z - y;
    // The compiler mustn't "simplify" z - (z - y) to y
    y = w > 0 ? z - w : y;
    return copysign (y, x);
  }

  /**
   * Evaluate the sine and cosine function with the argument in degrees
   *
   * @tparam T the type of the arguments.
   * @param[in] x in degrees.
   * @param[out] sinx sin(<i>x</i>).
   * @param[out] cosx cos(<i>x</i>).
   *
   * The results obey exactly the elementary properties of the trigonometric
   * functions, e.g., sin 9&deg; = cos 81&deg; = &minus; sin 123456789&deg;.
   * If x = &minus;0 or a negative multiple of 180&deg;, then \e sinx =
   * &minus;0; this is the only case where &minus;0 is returned.
   **********************************************************************/
  template<typename T> void sincosd(T x, T& sinx, T& cosx) {
    // In order to minimize round-off errors, this function exactly reduces
    // the argument to the range [-45, 45] before converting it to radians.
    T d, r; int q = 0;
    d = remquo (x, T (90), &q);   // now abs(r) <= 45
    r = d * degree<T> ();
    // g++ -O turns these two function calls into a call to sincos
    T s = sin (r), c = cos (r);
    if (2 * fabs (d) == 90) {
      c = sqrt (1 / T (2));
      s = copysign (c, r);
    }
    else if (3 * fabs (d) == 90) {
      c = sqrt (T (3)) / 2;
      s = copysign (1 / T (2), r);
    }
    switch (unsigned (q) & 3U) {
    case 0U: sinx = s; cosx = c; break;
    case 1U: sinx = c; cosx = -s; break;
    case 2U: sinx = -s; cosx = -c; break;
    default: sinx = -c; cosx = s; break; // case 3U
    }
    // http://www.open-std.org/jtc1/sc22/wg14/www/docs/n1950.pdf
    // mpreal needs T(0) here
    cosx += T (0);                            // special values from F.10.1.12
    if (sinx == 0) sinx = copysign (sinx, x); // special values from F.10.1.13
  }


  /**
   * Evaluate the sine and cosine with reduced argument plus correction
   *
   * @tparam T the type of the arguments.
   * @param[in] x reduced angle in degrees.
   * @param[in] t correction in degrees.
   * @param[out] sinx sin(<i>x</i> + <i>t</i>).
   * @param[out] cosx cos(<i>x</i> + <i>t</i>).
   *
   * This is a variant of Math::sincosd allowing a correction to the angle to
   * be supplied.  \e x must be in [&minus;180&deg;, 180&deg;] and \e t is
   * assumed to be a <i>small</i> correction.  Math::AngRound is applied to
   * the reduced angle to prevent problems with \e x + \e t being extremely
   * close but not exactly equal to one of the four cardinal directions.
   **********************************************************************/
  template<typename T> void sincosde(T x, T t, T& sinx, T& cosx) {
    // In order to minimize round-off errors, this function exactly reduces
    // the argument to the range [-45, 45] before converting it to radians.
    // This implementation allows x outside [-180, 180], but implementations in
    // other languages may not.
    int q = 0;
    T d = AngRound (remquo (x, T (90), &q) + t), // now abs(r) <= 45
      r = d * degree<T> ();
    // g++ -O turns these two function calls into a call to sincos
    T s = sin (r), c = cos (r);
    if (2 * fabs (d) == 90) {
      c = sqrt (1 / T (2));
      s = copysign (c, r);
    }
    else if (3 * fabs (d) == 90) {
      c = sqrt (T (3)) / 2;
      s = copysign (1 / T (2), r);
    }
    switch (unsigned (q) & 3U) {
    case 0U: sinx = s; cosx = c; break;
    case 1U: sinx = c; cosx = -s; break;
    case 2U: sinx = -s; cosx = -c; break;
    default: sinx = -c; cosx = s; break; // case 3U
    }
    // http://www.open-std.org/jtc1/sc22/wg14/www/docs/n1950.pdf
    // mpreal needs T(0) here
    cosx += T (0);                            // special values from F.10.1.12
    if (sinx == 0) sinx = copysign (sinx, x); // special values from F.10.1.13
  }


  /**
   * Evaluate the sine function with the argument in degrees
   *
   * @tparam T the type of the argument and the returned value.
   * @param[in] x in degrees.
   * @return sin(<i>x</i>).
   *
   * The result is +0 for \e x = +0 and positive multiples of 180&deg;.  The
   * result is &minus;0 for \e x = -0 and negative multiples of 180&deg;.
   **********************************************************************/
  template<typename T> T sind(T x) {
    // See sincosd
    int q = 0;
    T d = remquo (x, T (90), &q), // now abs(r) <= 45
      r = d * degree<T> ();
    unsigned p = unsigned (q);
    // r = p & 1U ? cos(r) : sin(r); replaced by ...
    r = p & 1U ? (2 * fabs (d) == 90 ? sqrt (1 / T (2)) :
                  (3 * fabs (d) == 90 ? sqrt (T (3)) / 2 : cos (r))) :
      copysign (2 * fabs (d) == 90 ? sqrt (1 / T (2)) :
               (3 * fabs (d) == 90 ? 1 / T (2) : sin (r)), r);
    if (p & 2U) r = -r;
    if (r == 0) r = copysign (r, x);
    return r;
  }

  /**
   * Evaluate the cosine function with the argument in degrees
   *
   * @tparam T the type of the argument and the returned value.
   * @param[in] x in degrees.
   * @return cos(<i>x</i>).
   *
   * The result is +0 for \e x an odd multiple of 90&deg;.
   **********************************************************************/
  template<typename T> T cosd(T x) {
    // See sincosd
    int q = 0;
    T d = remquo (x, T (90), &q), // now abs(r) <= 45
      r = d * degree<T> ();
    unsigned p = unsigned (q + 1);
    r = p & 1U ? (2 * fabs (d) == 90 ? sqrt (1 / T (2)) :
                  (3 * fabs (d) == 90 ? sqrt (T (3)) / 2 : cos (r))) :
      copysign (2 * fabs (d) == 90 ? sqrt (1 / T (2)) :
               (3 * fabs (d) == 90 ? 1 / T (2) : sin (r)), r);
    if (p & 2U) r = -r;
    // mpreal needs T(0) here
    return T (0) + r;
  }

  /**
   * Evaluate the tangent function with the argument in degrees
   *
   * @tparam T the type of the argument and the returned value.
   * @param[in] x in degrees.
   * @return tan(<i>x</i>).
   *
   * If \e x is an odd multiple of 90&deg;, then a suitably large (but
   * finite) value is returned.
   **********************************************************************/
  template<typename T> T tand(T x) {
    static const T overflow = 1 / sq (std::numeric_limits<T>::epsilon ());
    T s, c;
    sincosd (x, s, c);
    // http://www.open-std.org/jtc1/sc22/wg14/www/docs/n1950.pdf
    T r = s / c;  // special values from F.10.1.14
    // With C++17 this becomes clamp(s / c, -overflow, overflow);
    // Use max/min here (instead of fmax/fmin) to preserve NaN
    return std::min (std::max (r, -overflow), overflow);
  }

  /**
   * Evaluate the atan2 function with the result in degrees
   *
   * @tparam T the type of the arguments and the returned value.
   * @param[in] y
   * @param[in] x
   * @return atan2(<i>y</i>, <i>x</i>) in degrees.
   *
   * The result is in the range [&minus;180&deg; 180&deg;].  N.B.,
   * atan2d(&plusmn;0, &minus;1) = &plusmn;180&deg;.
   **********************************************************************/
  template<typename T> T atan2d(T y, T x) {
    // In order to minimize round-off errors, this function rearranges the
    // arguments so that result of atan2 is in the range [-pi/4, pi/4] before
    // converting it to degrees and mapping the result to the correct
    // quadrant.
    int q = 0;
    if (fabs (y) > fabs (x)) { std::swap (x, y); q = 2; }
    if (std::signbit (x)) { x = -x; ++q; }
    // here x >= 0 and x >= abs(y), so angle is in [-pi/4, pi/4]
    T ang = atan2 (y, x) / degree<T> ();
    switch (q) {
    case 1: ang = copysign (T (180), y) - ang; break;
    case 2: ang = 90 - ang; break;
    case 3: ang = -90 + ang; break;
    default: break;
    }
    return ang;
  }

  /**
   * Evaluate the atan function with the result in degrees
   *
   * @tparam T the type of the argument and the returned value.
   * @param[in] x
   * @return atan(<i>x</i>) in degrees.
   **********************************************************************/
  template<typename T> T atand(T x) { 
    return atan2d (x, T (1)); 
  }

  /**
   * Evaluate <i>e</i> atanh(<i>e x</i>)
   *
   * @tparam T the type of the argument and the returned value.
   * @param[in] x
   * @param[in] es the signed eccentricity =  sign(<i>e</i><sup>2</sup>)
   *    sqrt(|<i>e</i><sup>2</sup>|)
   * @return <i>e</i> atanh(<i>e x</i>)
   *
   * If <i>e</i><sup>2</sup> is negative (<i>e</i> is imaginary), the
   * expression is evaluated in terms of atan.
   **********************************************************************/
  template<typename T> T eatanhe(T x, T es) {
    return es > 0 ? es * atanh (es * x) : -es * atan (es * x);
  }

  /**
   * tan&chi; in terms of tan&phi;
   *
   * @tparam T the type of the argument and the returned value.
   * @param[in] tau &tau; = tan&phi;
   * @param[in] es the signed eccentricity = sign(<i>e</i><sup>2</sup>)
   *   sqrt(|<i>e</i><sup>2</sup>|)
   * @return &tau;&prime; = tan&chi;
   *
   * See Eqs. (7--9) of
   * C. F. F. Karney,
   * <a href="https://doi.org/10.1007/s00190-011-0445-3">
   * Transverse Mercator with an accuracy of a few nanometers,</a>
   * J. Geodesy 85(8), 475--485 (Aug. 2011)
   * (preprint
   * <a href="https://arxiv.org/abs/1002.1417">arXiv:1002.1417</a>).
   **********************************************************************/
  template<typename T> T taupf(T tau, T es) {
    // Need this test, otherwise tau = +/-inf gives taup = nan.
    if (std::isfinite (tau)) {
      T tau1 = hypot (T (1), tau),
        sig = sinh (eatanhe (tau / tau1, es));
      return hypot (T (1), sig) * tau - sig * tau1;
    }
    else
      return tau;
  }

  /**
   * tan&phi; in terms of tan&chi;
   *
   * @tparam T the type of the argument and the returned value.
   * @param[in] taup &tau;&prime; = tan&chi;
   * @param[in] es the signed eccentricity = sign(<i>e</i><sup>2</sup>)
   *   sqrt(|<i>e</i><sup>2</sup>|)
   * @return &tau; = tan&phi;
   *
   * See Eqs. (19--21) of
   * C. F. F. Karney,
   * <a href="https://doi.org/10.1007/s00190-011-0445-3">
   * Transverse Mercator with an accuracy of a few nanometers,</a>
   * J. Geodesy 85(8), 475--485 (Aug. 2011)
   * (preprint
   * <a href="https://arxiv.org/abs/1002.1417">arXiv:1002.1417</a>).
   **********************************************************************/
  template<typename T> T tauf(T taup, T es) {
    static const int numit = 5;
    // min iterations = 1, max iterations = 2; mean = 1.95
    static const T tol = sqrt (std::numeric_limits<T>::epsilon ()) / 10;
    static const T taumax = 2 / sqrt (std::numeric_limits<T>::epsilon ());
    T e2m = 1 - sq (es),
      // To lowest order in e^2, taup = (1 - e^2) * tau = e2m_ * tau; so use
      // tau = taup/e2m as a starting guess. Only 1 iteration is needed for
      // |lat| < 3.35 deg, otherwise 2 iterations are needed.  If, instead, tau
      // = taup is used the mean number of iterations increases to 1.999 (2
      // iterations are needed except near tau = 0).
      //
      // For large tau, taup = exp(-es*atanh(es)) * tau.  Use this as for the
      // initial guess for |taup| > 70 (approx |phi| > 89deg).  Then for
      // sufficiently large tau (such that sqrt(1+tau^2) = |tau|), we can exit
      // with the intial guess and avoid overflow problems.  This also reduces
      // the mean number of iterations slightly from 1.963 to 1.954.
      tau = fabs (taup) > 70 ? taup * exp (eatanhe (T (1), es)) : taup / e2m,
      stol = tol * fmax (T (1), fabs (taup));
    if (!(fabs (tau) < taumax)) return tau; // handles +/-inf and nan
    for (int i = 0;
         i < numit ||
           GEOGRAPHICLIB_PANIC ("Convergence failure in Math::tauf");
         ++i) {
      T taupa = taupf (tau, es),
        dtau = (taup - taupa) * (1 + e2m * sq (tau)) /
        (e2m * hypot (T (1), tau) * hypot (T (1), taupa));
      tau += dtau;
      if (!(fabs (dtau) >= stol))
        break;
    }
    return tau;
  }

  /**
   * Implement hypot with 3 parameters
   *
   * @tparam T the type of the argument and the returned value.
   * @param[in] x
   * @param[in] y
   * @param[in] z
   * @return sqrt(<i>x</i><sup>2</sup> + <i>y</i><sup>2</sup> +
   *   <i>z</i><sup>2</sup>).
   **********************************************************************/
  template<typename T> T hypot3(T x, T y, T z) {
#if __cplusplus < 201703L || GEOGRAPHICLIB_PRECISION == 4
    return sqrt (x * x + y * y + z * z);
#else
    return hypot (x, y, z);
#endif
  }

  /**
   * Swap the bytes of a quantity
   *
   * @tparam T the type of the argument and the returned value.
   * @param[in] x
   * @return x with its bytes swapped.
   **********************************************************************/
  template<typename T> T swab(T x) {
    union {
      T r;
      unsigned char c[sizeof(T)];
    } b;
    b.r = x;
    for (int i = sizeof(T)/2; i--; )
      std::swap(b.c[i], b.c[sizeof(T) - 1 - i]);
    return b.r;
  }

  inline constexpr int digits () {
#if GEOGRAPHICLIB_PRECISION != 5
    return std::numeric_limits<real>::digits;
#else
    return std::numeric_limits<real>::digits ();
#endif
  }

  inline constexpr int digits10 () {
#if GEOGRAPHICLIB_PRECISION != 5
    return std::numeric_limits<real>::digits10;
#else
    return std::numeric_limits<real>::digits10 ();
#endif
  }

  template<typename T> inline constexpr T NaN () {
#if defined(_MSC_VER)
    return std::numeric_limits<T>::has_quiet_NaN ?
      std::numeric_limits<T>::quiet_NaN () :
      (std::numeric_limits<T>::max)();
#else
    return std::numeric_limits<T>::has_quiet_NaN ?
      std::numeric_limits<T>::quiet_NaN () :
      std::numeric_limits<T>::max ();
#endif
  }

  template<typename T> inline constexpr T infinity () {
#if defined(_MSC_VER)
    return std::numeric_limits<T>::has_infinity ?
      std::numeric_limits<T>::infinity () :
      (std::numeric_limits<T>::max)();
#else
    return std::numeric_limits<T>::has_infinity ?
      std::numeric_limits<T>::infinity () :
      std::numeric_limits<T>::max ();
#endif
  }

} // namespace Math

} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_MATH_HPP
