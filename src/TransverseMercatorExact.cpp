/**
 * \file TransverseMercatorExact.cpp
 * \brief Implementation for GeographicLib::TransverseMercatorExact class
 *
 * Copyright (c) Charles Karney (2008-2022) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 *
 * The relevant section of Lee's paper is part V, pp 67--101,
 * <a href="https://doi.org/10.3138/X687-1574-4325-WM62">Conformal
 * Projections Based On Jacobian Elliptic Functions</a>;
 * <a href="https://archive.org/details/conformalproject0000leel/page/92">
 * borrow from archive.org</a>.
 *
 * The method entails using the Thompson Transverse Mercator as an
 * intermediate projection.  The projections from the intermediate
 * coordinates to [\e phi, \e lam] and [\e x, \e y] are given by elliptic
 * functions.  The inverse of these projections are found by Newton's method
 * with a suitable starting guess.
 *
 * This implementation and notation closely follows Lee, with the following
 * exceptions:
 * <center><table>
 * <tr><th>Lee    <th>here    <th>Description
 * <tr><td>x/a    <td>xi      <td>Northing (unit Earth)
 * <tr><td>y/a    <td>eta     <td>Easting (unit Earth)
 * <tr><td>s/a    <td>sigma   <td>xi + i * eta
 * <tr><td>y      <td>x       <td>Easting
 * <tr><td>x      <td>y       <td>Northing
 * <tr><td>k      <td>e       <td>eccentricity
 * <tr><td>k^2    <td>mu      <td>elliptic function parameter
 * <tr><td>k'^2   <td>mv      <td>elliptic function complementary parameter
 * <tr><td>m      <td>k       <td>scale
 * <tr><td>zeta   <td>zeta    <td>complex longitude = Mercator = chi in paper
 * <tr><td>s      <td>sigma   <td>complex GK = zeta in paper
 * </table></center>
 *
 * Minor alterations have been made in some of Lee's expressions in an
 * attempt to control round-off.  For example atanh(sin(phi)) is replaced by
 * asinh(tan(phi)) which maintains accuracy near phi = pi/2.  Such changes
 * are noted in the code.
 **********************************************************************/

#include <GeographicLib/TransverseMercatorExact.hpp>

namespace GeographicLib {

  using namespace std;

  TransverseMercatorExact::TransverseMercatorExact(real a, real f, real k0,
                                                   bool extendp)
    : tol_(numeric_limits<real>::epsilon())
    , tol2_(real(0.1) * tol_)
    , taytol_(pow(tol_, real(0.6)))
    , a_(a)
    , f_(f)
    , k0_(k0)
    , mu_(f_ * (2 - f_))        // e^2
    , mv_(1 - mu_)              // 1 - e^2
    , e_(sqrt(mu_))
    , extendp_(extendp)
    , eEu_(mu_)
    , eEv_(mv_)
  {
    if (!(isfinite(a_) && a_ > 0))
      throw GeographicErr("Equatorial radius is not positive");
    if (!(f_ > 0))
      throw GeographicErr("Flattening is not positive");
    if (!(f_ < 1))
      throw GeographicErr("Polar semi-axis is not positive");
    if (!(isfinite(k0_) && k0_ > 0))
      throw GeographicErr("Scale is not positive");
  }

  const TransverseMercatorExact& TransverseMercatorExact::UTM() {
    static const TransverseMercatorExact utm(Constants::WGS84_a(),
                                             Constants::WGS84_f(),
                                             Constants::UTM_k0());
    return utm;
  }

  void TransverseMercatorExact::zeta(real /*u*/, real snu, real cnu, real dnu,
                                     real /*v*/, real snv, real cnv, real dnv,
                                     real& taup, real& lam) const {
    // Lee 54.17 but write
    // atanh(snu * dnv) = asinh(snu * dnv / sqrt(cnu^2 + mv_ * snu^2 * snv^2))
    // atanh(e_ * snu / dnv) =
    //         asinh(e_ * snu / sqrt(mu_ * cnu^2 + mv_ * cnv^2))
    // Overflow value s.t. atan(overflow) = pi/2
    static const real
      overflow = 1 / Math::sq(numeric_limits<real>::epsilon());
    real
      d1 = sqrt(Math::sq(cnu) + mv_ * Math::sq(snu * snv)),
      d2 = sqrt(mu_ * Math::sq(cnu) + mv_ * Math::sq(cnv)),
      t1 = (d1 != 0 ? snu * dnv / d1 : (signbit(snu) ? -overflow : overflow)),
      t2 = (d2 != 0 ? sinh( e_ * asinh(e_ * snu / d2) ) :
            (signbit(snu) ? -overflow : overflow));
    // psi = asinh(t1) - asinh(t2)
    // taup = sinh(psi)
    taup = t1 * hypot(real(1), t2) - t2 * hypot(real(1), t1);
    lam = (d1 != 0 && d2 != 0) ?
      atan2(dnu * snv, cnu * cnv) - e_ * atan2(e_ * cnu * snv, dnu * cnv) :
      0;
  }

  void TransverseMercatorExact::dwdzeta(real /*u*/,
                                        real snu, real cnu, real dnu,
                                        real /*v*/,
                                        real snv, real cnv, real dnv,
                                        real& du, real& dv) const {
    // Lee 54.21 but write (1 - dnu^2 * snv^2) = (cnv^2 + mu_ * snu^2 * snv^2)
    // (see A+S 16.21.4)
    real d = mv_ * Math::sq(Math::sq(cnv) + mu_ * Math::sq(snu * snv));
    du =  cnu * dnu * dnv * (Math::sq(cnv) - mu_ * Math::sq(snu * snv)) / d;
    dv = -snu * snv * cnv * (Math::sq(dnu * dnv) + mu_ * Math::sq(cnu)) / d;
  }

  // Starting point for zetainv
  bool TransverseMercatorExact::zetainv0(real psi, real lam,
                                         real& u, real& v) const {
    bool retval = false;
    if (psi < -e_ * Math::pi()/4 &&
        lam > (1 - 2 * e_) * Math::pi()/2 &&
        psi < lam - (1 - e_) * Math::pi()/2) {
      // N.B. this branch is normally not taken because psi < 0 is converted
      // psi > 0 by Forward.
      //
      // There's a log singularity at w = w0 = Eu.K() + i * Ev.K(),
      // corresponding to the south pole, where we have, approximately
      //
      //   psi = e_ + i * pi/2 - e_ * atanh(cos(i * (w - w0)/(1 + mu_/2)))
      //
      // Inverting this gives:
      real
        psix = 1 - psi / e_,
        lamx = (Math::pi()/2 - lam) / e_;
      u = asinh(sin(lamx) / hypot(cos(lamx), sinh(psix))) *
        (1 + mu_/2);
      v = atan2(cos(lamx), sinh(psix)) * (1 + mu_/2);
      u = eEu_.K() - u;
      v = eEv_.K() - v;
    } else if (psi < e_ * Math::pi()/2 &&
               lam > (1 - 2 * e_) * Math::pi()/2) {
      // At w = w0 = i * Ev.K(), we have
      //
      //     zeta = zeta0 = i * (1 - e_) * pi/2
      //     zeta' = zeta'' = 0
      //
      // including the next term in the Taylor series gives:
      //
      // zeta = zeta0 - (mv_ * e_) / 3 * (w - w0)^3
      //
      // When inverting this, we map arg(w - w0) = [-90, 0] to
      // arg(zeta - zeta0) = [-90, 180]
      real
        dlam = lam - (1 - e_) * Math::pi()/2,
        rad = hypot(psi, dlam),
        // atan2(dlam-psi, psi+dlam) + 45d gives arg(zeta - zeta0) in range
        // [-135, 225).  Subtracting 180 (since multiplier is negative) makes
        // range [-315, 45).  Multiplying by 1/3 (for cube root) gives range
        // [-105, 15).  In particular the range [-90, 180] in zeta space maps
        // to [-90, 0] in w space as required.
        ang = atan2(dlam-psi, psi+dlam) - real(0.75) * Math::pi();
      // Error using this guess is about 0.21 * (rad/e)^(5/3)
      retval = rad < e_ * taytol_;
      rad = cbrt(3 / (mv_ * e_) * rad);
      ang /= 3;
      u = rad * cos(ang);
      v = rad * sin(ang) + eEv_.K();
    } else {
      // Use spherical TM, Lee 12.6 -- writing atanh(sin(lam) / cosh(psi)) =
      // asinh(sin(lam) / hypot(cos(lam), sinh(psi))).  This takes care of the
      // log singularity at zeta = Eu.K() (corresponding to the north pole)
      v = asinh(sin(lam) / hypot(cos(lam), sinh(psi)));
      u = atan2(sinh(psi), cos(lam));
      // But scale to put 90,0 on the right place
      u *= eEu_.K() / (Math::pi()/2);
      v *= eEu_.K() / (Math::pi()/2);
    }
    return retval;
  }

  // Invert zeta using Newton's method
  void TransverseMercatorExact::zetainv(real taup, real lam,
                                        real& u, real& v) const  {
    real
      psi = asinh(taup),
      scal = 1/hypot(real(1), taup);
    if (zetainv0(psi, lam, u, v))
      return;
    real stol2 = tol2_ / Math::sq(fmax(psi, real(1)));
    // min iterations = 2, max iterations = 6; mean = 4.0
    for (int i = 0, trip = 0;
         i < numit_ ||
           GEOGRAPHICLIB_PANIC
           ("Convergence failure in TransverseMercatorExact");
         ++i) {
      real snu, cnu, dnu, snv, cnv, dnv;
      eEu_.am(u, snu, cnu, dnu);
      eEv_.am(v, snv, cnv, dnv);
      real tau1, lam1, du1, dv1;
      zeta(u, snu, cnu, dnu, v, snv, cnv, dnv, tau1, lam1);
      dwdzeta(u, snu, cnu, dnu, v, snv, cnv, dnv, du1, dv1);
      tau1 -= taup;
      lam1 -= lam;
      tau1 *= scal;
      real
        delu = tau1 * du1 - lam1 * dv1,
        delv = tau1 * dv1 + lam1 * du1;
      u -= delu;
      v -= delv;
      if (trip)
        break;
      real delw2 = Math::sq(delu) + Math::sq(delv);
      if (!(delw2 >= stol2))
        ++trip;
    }
  }

  void TransverseMercatorExact::sigma(real /*u*/, real snu, real cnu, real dnu,
                                      real v, real snv, real cnv, real dnv,
                                      real& xi, real& eta) const {
    // Lee 55.4 writing
    // dnu^2 + dnv^2 - 1 = mu_ * cnu^2 + mv_ * cnv^2
    real d = mu_ * Math::sq(cnu) + mv_ * Math::sq(cnv);
    xi = eEu_.E(snu, cnu, dnu) - mu_ * snu * cnu * dnu / d;
    eta = v - eEv_.E(snv, cnv, dnv) + mv_ * snv * cnv * dnv / d;
  }

  void TransverseMercatorExact::dwdsigma(real /*u*/,
                                         real snu, real cnu, real dnu,
                                         real /*v*/,
                                         real snv, real cnv, real dnv,
                                         real& du, real& dv) const {
    // Reciprocal of 55.9: dw/ds = dn(w)^2/mv_, expanding complex dn(w) using
    // A+S 16.21.4
    real d = mv_ * Math::sq(Math::sq(cnv) + mu_ * Math::sq(snu * snv));
    real
      dnr = dnu * cnv * dnv,
      dni = - mu_ * snu * cnu * snv;
    du = (Math::sq(dnr) - Math::sq(dni)) / d;
    dv = 2 * dnr * dni / d;
  }

  // Starting point for sigmainv
  bool TransverseMercatorExact::sigmainv0(real xi, real eta,
                                          real& u, real& v) const {
    bool retval = false;
    if (eta > real(1.25) * eEv_.KE() ||
        (xi < -real(0.25) * eEu_.E() && xi < eta - eEv_.KE())) {
      // sigma as a simple pole at w = w0 = Eu.K() + i * Ev.K() and sigma is
      // approximated by
      //
      // sigma = (Eu.E() + i * Ev.KE()) + 1/(w - w0)
      real
        x = xi - eEu_.E(),
        y = eta - eEv_.KE(),
        r2 = Math::sq(x) + Math::sq(y);
      u = eEu_.K() + x/r2;
      v = eEv_.K() - y/r2;
    } else if ((eta > real(0.75) * eEv_.KE() && xi < real(0.25) * eEu_.E())
               || eta > eEv_.KE()) {
      // At w = w0 = i * Ev.K(), we have
      //
      //     sigma = sigma0 = i * Ev.KE()
      //     sigma' = sigma'' = 0
      //
      // including the next term in the Taylor series gives:
      //
      // sigma = sigma0 - mv_ / 3 * (w - w0)^3
      //
      // When inverting this, we map arg(w - w0) = [-pi/2, -pi/6] to
      // arg(sigma - sigma0) = [-pi/2, pi/2]
      // mapping arg = [-pi/2, -pi/6] to [-pi/2, pi/2]
      real
        deta = eta - eEv_.KE(),
        rad = hypot(xi, deta),
        // Map the range [-90, 180] in sigma space to [-90, 0] in w space.  See
        // discussion in zetainv0 on the cut for ang.
        ang = atan2(deta-xi, xi+deta) - real(0.75) * Math::pi();
      // Error using this guess is about 0.068 * rad^(5/3)
      retval = rad < 2 * taytol_;
      rad = cbrt(3 / mv_ * rad);
      ang /= 3;
      u = rad * cos(ang);
      v = rad * sin(ang) + eEv_.K();
    } else {
      // Else use w = sigma * Eu.K/Eu.E (which is correct in the limit e_ -> 0)
      u = xi * eEu_.K()/eEu_.E();
      v = eta * eEu_.K()/eEu_.E();
    }
    return retval;
  }

  // Invert sigma using Newton's method
  void TransverseMercatorExact::sigmainv(real xi, real eta,
                                         real& u, real& v) const {
    if (sigmainv0(xi, eta, u, v))
      return;
    // min iterations = 2, max iterations = 7; mean = 3.9
    for (int i = 0, trip = 0;
         i < numit_ ||
           GEOGRAPHICLIB_PANIC
           ("Convergence failure in TransverseMercatorExact");
         ++i) {
      real snu, cnu, dnu, snv, cnv, dnv;
      eEu_.am(u, snu, cnu, dnu);
      eEv_.am(v, snv, cnv, dnv);
      real xi1, eta1, du1, dv1;
      sigma(u, snu, cnu, dnu, v, snv, cnv, dnv, xi1, eta1);
      dwdsigma(u, snu, cnu, dnu, v, snv, cnv, dnv, du1, dv1);
      xi1 -= xi;
      eta1 -= eta;
      real
        delu = xi1 * du1 - eta1 * dv1,
        delv = xi1 * dv1 + eta1 * du1;
      u -= delu;
      v -= delv;
      if (trip)
        break;
      real delw2 = Math::sq(delu) + Math::sq(delv);
      if (!(delw2 >= tol2_))
        ++trip;
    }
  }

  void TransverseMercatorExact::Scale(real tau, real /*lam*/,
                                      real snu, real cnu, real dnu,
                                      real snv, real cnv, real dnv,
                                      real& gamma, real& k) const {
    real sec2 = 1 + Math::sq(tau);    // sec(phi)^2
    // Lee 55.12 -- negated for our sign convention.  gamma gives the bearing
    // (clockwise from true north) of grid north
    gamma = atan2(mv_ * snu * snv * cnv, cnu * dnu * dnv);
    // Lee 55.13 with nu given by Lee 9.1 -- in sqrt change the numerator
    // from
    //
    //    (1 - snu^2 * dnv^2) to (mv_ * snv^2 + cnu^2 * dnv^2)
    //
    // to maintain accuracy near phi = 90 and change the denomintor from
    //
    //    (dnu^2 + dnv^2 - 1) to (mu_ * cnu^2 + mv_ * cnv^2)
    //
    // to maintain accuracy near phi = 0, lam = 90 * (1 - e).  Similarly
    // rewrite sqrt term in 9.1 as
    //
    //    mv_ + mu_ * c^2 instead of 1 - mu_ * sin(phi)^2
    k = sqrt(mv_ + mu_ / sec2) * sqrt(sec2) *
      sqrt( (mv_ * Math::sq(snv) + Math::sq(cnu * dnv)) /
            (mu_ * Math::sq(cnu) + mv_ * Math::sq(cnv)) );
  }

  void TransverseMercatorExact::Forward(real lon0, real lat, real lon,
                                        real& x, real& y,
                                        real& gamma, real& k) const {
    lat = Math::LatFix(lat);
    lon = Math::AngDiff(lon0, lon);
    // Explicitly enforce the parity
    int
      latsign = (!extendp_ && signbit(lat)) ? -1 : 1,
      lonsign = (!extendp_ && signbit(lon)) ? -1 : 1;
    lon *= lonsign;
    lat *= latsign;
    bool backside = !extendp_ && lon > 90;
    if (backside) {
      if (lat == 0)
        latsign = -1;
      lon = 180 - lon;
    }
    real
      lam = lon * Math::degree(),
      tau = Math::tand(lat);

    // u,v = coordinates for the Thompson TM, Lee 54
    real u, v;
    if (lat == 90) {
      u = eEu_.K();
      v = 0;
    } else if (lat == 0 && lon == 90 * (1 - e_)) {
      u = 0;
      v = eEv_.K();
    } else
      // tau = tan(phi), taup = sinh(psi)
      zetainv(Math::taupf(tau, e_), lam, u, v);

    real snu, cnu, dnu, snv, cnv, dnv;
    eEu_.am(u, snu, cnu, dnu);
    eEv_.am(v, snv, cnv, dnv);

    real xi, eta;
    sigma(u, snu, cnu, dnu, v, snv, cnv, dnv, xi, eta);
    if (backside)
      xi = 2 * eEu_.E() - xi;
    y = xi * a_ * k0_ * latsign;
    x = eta * a_ * k0_ * lonsign;

    if (lat == 90) {
      gamma = lon;
      k = 1;
    } else {
      // Recompute (tau, lam) from (u, v) to improve accuracy of Scale
      zeta(u, snu, cnu, dnu, v, snv, cnv, dnv, tau, lam);
      tau = Math::tauf(tau, e_);
      Scale(tau, lam, snu, cnu, dnu, snv, cnv, dnv, gamma, k);
      gamma /= Math::degree();
    }
    if (backside)
      gamma = 180 - gamma;
    gamma *= latsign * lonsign;
    k *= k0_;
  }

  void TransverseMercatorExact::Reverse(real lon0, real x, real y,
                                        real& lat, real& lon,
                                        real& gamma, real& k) const {
    // This undoes the steps in Forward.
    real
      xi = y / (a_ * k0_),
      eta = x / (a_ * k0_);
    // Explicitly enforce the parity
    int
      xisign = (!extendp_ && signbit(xi)) ? -1 : 1,
      etasign = (!extendp_ && signbit(eta)) ? -1 : 1;
    xi *= xisign;
    eta *= etasign;
    bool backside = !extendp_ && xi > eEu_.E();
    if (backside)
      xi = 2 * eEu_.E()- xi;

    // u,v = coordinates for the Thompson TM, Lee 54
    real u, v;
    if (xi == 0 && eta == eEv_.KE()) {
      u = 0;
      v = eEv_.K();
    } else
      sigmainv(xi, eta, u, v);

    real snu, cnu, dnu, snv, cnv, dnv;
    eEu_.am(u, snu, cnu, dnu);
    eEv_.am(v, snv, cnv, dnv);
    real phi, lam, tau;
    if (v != 0 || u != eEu_.K()) {
      zeta(u, snu, cnu, dnu, v, snv, cnv, dnv, tau, lam);
      tau = Math::tauf(tau, e_);
      phi = atan(tau);
      lat = phi / Math::degree();
      lon = lam / Math::degree();
      Scale(tau, lam, snu, cnu, dnu, snv, cnv, dnv, gamma, k);
      gamma /= Math::degree();
    } else {
      lat = 90;
      lon = lam = gamma = 0;
      k = 1;
    }

    if (backside)
      lon = 180 - lon;
    lon *= etasign;
    lon = Math::AngNormalize(lon + Math::AngNormalize(lon0));
    lat *= xisign;
    if (backside)
      gamma = 180 - gamma;
    gamma *= xisign * etasign;
    k *= k0_;
  }

} // namespace GeographicLib
