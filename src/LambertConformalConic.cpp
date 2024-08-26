/**
 * \file LambertConformalConic.cpp
 * \brief Implementation for GeographicLib::LambertConformalConic class
 *
 * Copyright (c) Charles Karney (2010-2022) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#include <GeographicLib/LambertConformalConic.hpp>

namespace GeographicLib {

  using namespace std;

  LambertConformalConic::LambertConformalConic(real a, real f,
                                               real stdlat, real k0)
    : eps_(numeric_limits<real>::epsilon())
    , epsx_(Math::sq(eps_))
    , ahypover_(Math::digits() * log(real(numeric_limits<real>::radix)) + 2)
    , a_(a)
    , f_(f)
    , fm_(1 - f_)
    , e2_(f_ * (2 - f_))
    , es_((f_ < 0 ? -1 : 1) * sqrt(fabs(e2_)))
  {
    if (!(isfinite(a_) && a_ > 0))
      throw GeographicErr("Equatorial radius is not positive");
    if (!(isfinite(f_) && f_ < 1))
      throw GeographicErr("Polar semi-axis is not positive");
    if (!(isfinite(k0) && k0 > 0))
      throw GeographicErr("Scale is not positive");
    if (!(fabs (stdlat) <= 90))
      throw GeographicErr ("Standard latitude not in [-90d, +90d]");
    real sphi, cphi;
    Math::sincosd(stdlat, sphi, cphi);
    Init(sphi, cphi, sphi, cphi, k0);
  }

  LambertConformalConic::LambertConformalConic(real a, real f,
                                               real stdlat1, real stdlat2,
                                               real k1)
    : eps_(numeric_limits<real>::epsilon())
    , epsx_(Math::sq(eps_))
    , ahypover_(Math::digits() * log(real(numeric_limits<real>::radix)) + 2)
    , a_(a)
    , f_(f)
    , fm_(1 - f_)
    , e2_(f_ * (2 - f_))
    , es_((f_ < 0 ? -1 : 1) * sqrt(fabs(e2_)))
  {
    if (!(isfinite(a_) && a_ > 0))
      throw GeographicErr("Equatorial radius is not positive");
    if (!(isfinite(f_) && f_ < 1))
      throw GeographicErr("Polar semi-axis is not positive");
    if (!(isfinite(k1) && k1 > 0))
      throw GeographicErr("Scale is not positive");
    if (!(fabs (stdlat1) <= 90))
      throw GeographicErr ("Standard latitude 1 not in [-90d, +90d]");
    if (!(fabs (stdlat2) <= 90))
      throw GeographicErr ("Standard latitude 2 not in [-90d, +90d]");
    real sphi1, cphi1, sphi2, cphi2;
    Math::sincosd(stdlat1, sphi1, cphi1);
    Math::sincosd(stdlat2, sphi2, cphi2);
    Init(sphi1, cphi1, sphi2, cphi2, k1);
  }

  LambertConformalConic::LambertConformalConic(real a, real f,
                                               real sinlat1, real coslat1,
                                               real sinlat2, real coslat2,
                                               real k1)
    : eps_(numeric_limits<real>::epsilon())
    , epsx_(Math::sq(eps_))
    , ahypover_(Math::digits() * log(real(numeric_limits<real>::radix)) + 2)
    , a_(a)
    , f_(f)
    , fm_(1 - f_)
    , e2_(f_ * (2 - f_))
    , es_((f_ < 0 ? -1 : 1) * sqrt(fabs(e2_)))
  {
    if (!(isfinite(a_) && a_ > 0))
      throw GeographicErr("Equatorial radius is not positive");
    if (!(isfinite(f_) && f_ < 1))
      throw GeographicErr("Polar semi-axis is not positive");
    if (!(isfinite(k1) && k1 > 0))
      throw GeographicErr("Scale is not positive");
    if (signbit (coslat1))
      throw GeographicErr ("Standard latitude 1 not in [-90d, +90d]");
    if (signbit (coslat2))
      throw GeographicErr ("Standard latitude 2 not in [-90d, +90d]");
    if (!(fabs(sinlat1) <= 1 && coslat1 <= 1) || (coslat1 == 0 && sinlat1 == 0))
      throw GeographicErr("Bad sine/cosine of standard latitude 1");
    if (!(fabs(sinlat2) <= 1 && coslat2 <= 1) || (coslat2 == 0 && sinlat2 == 0))
      throw GeographicErr("Bad sine/cosine of standard latitude 2");
    if (coslat1 == 0 || coslat2 == 0)
      if (!(coslat1 == coslat2 && sinlat1 == sinlat2))
        throw GeographicErr
          ("Standard latitudes must be equal is either is a pole");
    Init(sinlat1, coslat1, sinlat2, coslat2, k1);
  }

  void LambertConformalConic::Init(real sphi1, real cphi1,
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
      tphi1 = sphi1/cphi1, tphi2 = sphi2/cphi2, tphi0;
    //
    // Snyder: 15-8: n = (log(m1) - log(m2))/(log(t1)-log(t2))
    //
    // m = cos(bet) = 1/sec(bet) = 1/sqrt(1+tan(bet)^2)
    // bet = parametric lat, tan(bet) = (1-f)*tan(phi)
    //
    // t = tan(pi/4-chi/2) = 1/(sec(chi) + tan(chi)) = sec(chi) - tan(chi)
    // log(t) = -asinh(tan(chi)) = -psi
    // chi = conformal lat
    // tan(chi) = tan(phi)*cosh(xi) - sinh(xi)*sec(phi)
    // xi = eatanhe(sin(phi)), eatanhe(x) = e * atanh(e*x)
    //
    // n = (log(sec(bet2))-log(sec(bet1)))/(asinh(tan(chi2))-asinh(tan(chi1)))
    //
    // Let log(sec(bet)) = b(tphi), asinh(tan(chi)) = c(tphi)
    // Then n = Db(tphi2, tphi1)/Dc(tphi2, tphi1)
    // In limit tphi2 -> tphi1, n -> sphi1
    //
    real
      tbet1 = fm_ * tphi1, scbet1 = hyp(tbet1),
      tbet2 = fm_ * tphi2, scbet2 = hyp(tbet2);
    real
      scphi1 = 1/cphi1,
      xi1 = Math::eatanhe(sphi1, es_), shxi1 = sinh(xi1), chxi1 = hyp(shxi1),
      tchi1 = chxi1 * tphi1 - shxi1 * scphi1, scchi1 = hyp(tchi1),
      scphi2 = 1/cphi2,
      xi2 = Math::eatanhe(sphi2, es_), shxi2 = sinh(xi2), chxi2 = hyp(shxi2),
      tchi2 = chxi2 * tphi2 - shxi2 * scphi2, scchi2 = hyp(tchi2),
      psi1 = asinh(tchi1);
    if (tphi2 - tphi1 != 0) {
      // Db(tphi2, tphi1)
      real num = Dlog1p(Math::sq(tbet2)/(1 + scbet2),
                        Math::sq(tbet1)/(1 + scbet1))
        * Dhyp(tbet2, tbet1, scbet2, scbet1) * fm_;
      // Dc(tphi2, tphi1)
      real den = Dasinh(tphi2, tphi1, scphi2, scphi1)
        - Deatanhe(sphi2, sphi1) * Dsn(tphi2, tphi1, sphi2, sphi1);
      n_ = num/den;

      if (n_ < 1/real(4))
        nc_ = sqrt((1 - n_) * (1 + n_));
      else {
        // Compute nc = cos(phi0) = sqrt((1 - n) * (1 + n)), evaluating 1 - n
        // carefully.  First write
        //
        // Dc(tphi2, tphi1) * (tphi2 - tphi1)
        //   = log(tchi2 + scchi2) - log(tchi1 + scchi1)
        //
        // then den * (1 - n) =
        // (log((tchi2 + scchi2)/(2*scbet2)) -
        //  log((tchi1 + scchi1)/(2*scbet1))) / (tphi2 - tphi1)
        // = Dlog1p(a2, a1) * (tchi2+scchi2 + tchi1+scchi1)/(4*scbet1*scbet2)
        //   * fm * Q
        //
        // where
        // a1 = ( (tchi1 - scbet1) + (scchi1 - scbet1) ) / (2 * scbet1)
        // Q = ((scbet2 + scbet1)/fm)/((scchi2 + scchi1)/D(tchi2, tchi1))
        //     - (tbet2 + tbet1)/(scbet2 + scbet1)
        real t;
        {
          real
            // s1 = (scbet1 - scchi1) * (scbet1 + scchi1)
            s1 = (tphi1 * (2 * shxi1 * chxi1 * scphi1 - e2_ * tphi1) -
                  Math::sq(shxi1) * (1 + 2 * Math::sq(tphi1))),
            s2 = (tphi2 * (2 * shxi2 * chxi2 * scphi2 - e2_ * tphi2) -
                  Math::sq(shxi2) * (1 + 2 * Math::sq(tphi2))),
            // t1 = scbet1 - tchi1
            t1 = tchi1 < 0 ? scbet1 - tchi1 : (s1 + 1)/(scbet1 + tchi1),
            t2 = tchi2 < 0 ? scbet2 - tchi2 : (s2 + 1)/(scbet2 + tchi2),
            a2 = -(s2 / (scbet2 + scchi2) + t2) / (2 * scbet2),
            a1 = -(s1 / (scbet1 + scchi1) + t1) / (2 * scbet1);
          t = Dlog1p(a2, a1) / den;
        }
        // multiply by (tchi2 + scchi2 + tchi1 + scchi1)/(4*scbet1*scbet2) * fm
        t *= ( ( (tchi2 >= 0 ? scchi2 + tchi2 : 1/(scchi2 - tchi2)) +
                 (tchi1 >= 0 ? scchi1 + tchi1 : 1/(scchi1 - tchi1)) ) /
               (4 * scbet1 * scbet2) ) * fm_;

        // Rewrite
        // Q = (1 - (tbet2 + tbet1)/(scbet2 + scbet1)) -
        //     (1 - ((scbet2 + scbet1)/fm)/((scchi2 + scchi1)/D(tchi2, tchi1)))
        //   = tbm - tam
        // where
        real tbm = ( ((tbet1 > 0 ? 1/(scbet1+tbet1) : scbet1 - tbet1) +
                      (tbet2 > 0 ? 1/(scbet2+tbet2) : scbet2 - tbet2)) /
                     (scbet1+scbet2) );

        // tam = (1 - ((scbet2+scbet1)/fm)/((scchi2+scchi1)/D(tchi2, tchi1)))
        //
        // Let
        //   (scbet2 + scbet1)/fm = scphi2 + scphi1 + dbet
        //   (scchi2 + scchi1)/D(tchi2, tchi1) = scphi2 + scphi1 + dchi
        // then
        //   tam = D(tchi2, tchi1) * (dchi - dbet) / (scchi1 + scchi2)
        real
          // D(tchi2, tchi1)
          dtchi = den / Dasinh(tchi2, tchi1, scchi2, scchi1),
          // (scbet2 + scbet1)/fm - (scphi2 + scphi1)
          dbet = (e2_/fm_) * ( 1 / (scbet2 + fm_ * scphi2) +
                               1 / (scbet1 + fm_ * scphi1) );

        // dchi = (scchi2 + scchi1)/D(tchi2, tchi1) - (scphi2 + scphi1)
        // Let
        //    tzet = chxiZ * tphi - shxiZ * scphi
        //    tchi = tzet + nu
        //    scchi = sczet + mu
        // where
        //    xiZ = eatanhe(1), shxiZ = sinh(xiZ), chxiZ = cosh(xiZ)
        //    nu =   scphi * (shxiZ - shxi) - tphi * (chxiZ - chxi)
        //    mu = - scphi * (chxiZ - chxi) + tphi * (shxiZ - shxi)
        // then
        // dchi = ((mu2 + mu1) - D(nu2, nu1) * (scphi2 + scphi1)) /
        //         D(tchi2, tchi1)
        real
          xiZ = Math::eatanhe(real(1), es_),
          shxiZ = sinh(xiZ), chxiZ = hyp(shxiZ),
          // These are differences not divided differences
          // dxiZ1 = xiZ - xi1; dshxiZ1 = shxiZ - shxi; dchxiZ1 = chxiZ - chxi
          dxiZ1 = Deatanhe(real(1), sphi1)/(scphi1*(tphi1+scphi1)),
          dxiZ2 = Deatanhe(real(1), sphi2)/(scphi2*(tphi2+scphi2)),
          dshxiZ1 = Dsinh(xiZ, xi1, shxiZ, shxi1, chxiZ, chxi1) * dxiZ1,
          dshxiZ2 = Dsinh(xiZ, xi2, shxiZ, shxi2, chxiZ, chxi2) * dxiZ2,
          dchxiZ1 = Dhyp(shxiZ, shxi1, chxiZ, chxi1) * dshxiZ1,
          dchxiZ2 = Dhyp(shxiZ, shxi2, chxiZ, chxi2) * dshxiZ2,
          // mu1 + mu2
          amu12 = (- scphi1 * dchxiZ1 + tphi1 * dshxiZ1
                   - scphi2 * dchxiZ2 + tphi2 * dshxiZ2),
          // D(xi2, xi1)
          dxi = Deatanhe(sphi1, sphi2) * Dsn(tphi2, tphi1, sphi2, sphi1),
          // D(nu2, nu1)
          dnu12 =
          ( (f_ * 4 * scphi2 * dshxiZ2 > f_ * scphi1 * dshxiZ1 ?
             // Use divided differences
             (dshxiZ1 + dshxiZ2)/2 * Dhyp(tphi1, tphi2, scphi1, scphi2)
             - ( (scphi1 + scphi2)/2
                 * Dsinh(xi1, xi2, shxi1, shxi2, chxi1, chxi2) * dxi ) :
             // Use ratio of differences
             (scphi2 * dshxiZ2 - scphi1 * dshxiZ1)/(tphi2 - tphi1))
            + ( (tphi1 + tphi2)/2 * Dhyp(shxi1, shxi2, chxi1, chxi2)
                * Dsinh(xi1, xi2, shxi1, shxi2, chxi1, chxi2) * dxi )
            - (dchxiZ1 + dchxiZ2)/2 ),
          // dtchi * dchi
          dchia = (amu12 - dnu12 * (scphi2 + scphi1)),
          tam = (dchia - dtchi * dbet) / (scchi1 + scchi2);
        t *= tbm - tam;
        nc_ = sqrt(fmax(real(0), t) * (1 + n_));
      }
      {
        real r = hypot(n_, nc_);
        n_ /= r;
        nc_ /= r;
      }
      tphi0 = n_ / nc_;
    } else {
      tphi0 = tphi1;
      nc_ = 1/hyp(tphi0);
      n_ = tphi0 * nc_;
      if (polar)
        nc_ = 0;
    }

    scbet0_ = hyp(fm_ * tphi0);
    real shxi0 = sinh(Math::eatanhe(n_, es_));
    tchi0_ = tphi0 * hyp(shxi0) - shxi0 * hyp(tphi0); scchi0_ = hyp(tchi0_);
    psi0_ = asinh(tchi0_);

    lat0_ = atan(sign_ * tphi0) / Math::degree();
    t0nm1_ = expm1(- n_ * psi0_); // Snyder's t0^n - 1
    // a * k1 * m1/t1^n = a * k1 * m2/t2^n = a * k1 * n * (Snyder's F)
    // = a * k1 / (scbet1 * exp(-n * psi1))
    scale_ = a_ * k1 / scbet1 *
      // exp(n * psi1) = exp(- (1 - n) * psi1) * exp(psi1)
      // with (1-n) = nc^2/(1+n) and exp(-psi1) = scchi1 + tchi1
      exp( - (Math::sq(nc_)/(1 + n_)) * psi1 )
      * (tchi1 >= 0 ? scchi1 + tchi1 : 1 / (scchi1 - tchi1));
    // Scale at phi0 = k0 = k1 * (scbet0*exp(-n*psi0))/(scbet1*exp(-n*psi1))
    //                    = k1 * scbet0/scbet1 * exp(n * (psi1 - psi0))
    // psi1 - psi0 = Dasinh(tchi1, tchi0) * (tchi1 - tchi0)
    k0_ = k1 * (scbet0_/scbet1) *
      exp( - (Math::sq(nc_)/(1 + n_)) *
           Dasinh(tchi1, tchi0_, scchi1, scchi0_) * (tchi1 - tchi0_))
      * (tchi1 >= 0 ? scchi1 + tchi1 : 1 / (scchi1 - tchi1)) /
      (scchi0_ + tchi0_);
    nrho0_ = polar ? 0 : a_ * k0_ / scbet0_;
    {
      // Figure drhomax_ using code at beginning of Forward with lat = -90
      real
        sphi = -1, cphi =  epsx_,
        tphi = sphi/cphi,
        scphi = 1/cphi, shxi = sinh(Math::eatanhe(sphi, es_)),
        tchi = hyp(shxi) * tphi - shxi * scphi, scchi = hyp(tchi),
        psi = asinh(tchi),
        dpsi = Dasinh(tchi, tchi0_, scchi, scchi0_) * (tchi - tchi0_);
      drhomax_ = - scale_ * (2 * nc_ < 1 && dpsi != 0 ?
                             (exp(Math::sq(nc_)/(1 + n_) * psi ) *
                              (tchi > 0 ? 1/(scchi + tchi) : (scchi - tchi))
                              - (t0nm1_ + 1))/(-n_) :
                             Dexp(-n_ * psi, -n_ * psi0_) * dpsi);
    }
  }

  const LambertConformalConic& LambertConformalConic::Mercator() {
    static const LambertConformalConic mercator(Constants::WGS84_a(),
                                                Constants::WGS84_f(),
                                                real(0), real(1));
    return mercator;
  }

  void LambertConformalConic::Forward(real lon0, real lat, real lon,
                                      real& x, real& y,
                                      real& gamma, real& k) const {
    lon = Math::AngDiff(lon0, lon);
    // From Snyder, we have
    //
    // theta = n * lambda
    // x = rho * sin(theta)
    //   = (nrho0 + n * drho) * sin(theta)/n
    // y = rho0 - rho * cos(theta)
    //   = nrho0 * (1-cos(theta))/n - drho * cos(theta)
    //
    // where nrho0 = n * rho0, drho = rho - rho0
    // and drho is evaluated with divided differences
    real sphi, cphi;
    Math::sincosd(Math::LatFix(lat) * sign_, sphi, cphi);
    cphi = fmax(epsx_, cphi);
    real
      lam = lon * Math::degree(),
      tphi = sphi/cphi, scbet = hyp(fm_ * tphi),
      scphi = 1/cphi, shxi = sinh(Math::eatanhe(sphi, es_)),
      tchi = hyp(shxi) * tphi - shxi * scphi, scchi = hyp(tchi),
      psi = asinh(tchi),
      theta = n_ * lam, stheta = sin(theta), ctheta = cos(theta),
      dpsi = Dasinh(tchi, tchi0_, scchi, scchi0_) * (tchi - tchi0_),
      drho = - scale_ * (2 * nc_ < 1 && dpsi != 0 ?
                         (exp(Math::sq(nc_)/(1 + n_) * psi ) *
                          (tchi > 0 ? 1/(scchi + tchi) : (scchi - tchi))
                          - (t0nm1_ + 1))/(-n_) :
                         Dexp(-n_ * psi, -n_ * psi0_) * dpsi);
    x = (nrho0_ + n_ * drho) * (n_ != 0 ? stheta / n_ : lam);
    y = nrho0_ *
      (n_ != 0 ?
       (ctheta < 0 ? 1 - ctheta : Math::sq(stheta)/(1 + ctheta)) / n_ : 0)
      - drho * ctheta;
    k = k0_ * (scbet/scbet0_) /
      (exp( - (Math::sq(nc_)/(1 + n_)) * dpsi )
       * (tchi >= 0 ? scchi + tchi : 1 / (scchi - tchi)) / (scchi0_ + tchi0_));
    y *= sign_;
    gamma = sign_ * theta / Math::degree();
  }

  void LambertConformalConic::Reverse(real lon0, real x, real y,
                                      real& lat, real& lon,
                                      real& gamma, real& k) const {
    // From Snyder, we have
    //
    //        x = rho * sin(theta)
    // rho0 - y = rho * cos(theta)
    //
    // rho = hypot(x, rho0 - y)
    // drho = (n*x^2 - 2*y*nrho0 + n*y^2)/(hypot(n*x, nrho0-n*y) + nrho0)
    // theta = atan2(n*x, nrho0-n*y)
    //
    // From drho, obtain t^n-1
    // psi = -log(t), so
    // dpsi = - Dlog1p(t^n-1, t0^n-1) * drho / scale
    y *= sign_;
    real
      // Guard against 0 * inf in computation of ny
      nx = n_ * x, ny = n_ != 0 ? n_ * y : 0, y1 = nrho0_ - ny,
      den = hypot(nx, y1) + nrho0_, // 0 implies origin with polar aspect
      // isfinite test is to avoid inf/inf
      drho = ((den != 0 && isfinite(den))
              ? (x*nx + y * (ny - 2*nrho0_)) / den
              : den);
    drho = fmin(drho, drhomax_);
    if (n_ == 0)
      drho = fmax(drho, -drhomax_);
    real
      tnm1 = t0nm1_ + n_ * drho/scale_,
      dpsi = (den == 0 ? 0 :
              (tnm1 + 1 != 0 ? - Dlog1p(tnm1, t0nm1_) * drho / scale_ :
               ahypover_));
    real tchi;
    if (2 * n_ <= 1) {
      // tchi = sinh(psi)
      real
        psi = psi0_ + dpsi, tchia = sinh(psi), scchi = hyp(tchia),
        dtchi = Dsinh(psi, psi0_, tchia, tchi0_, scchi, scchi0_) * dpsi;
      tchi = tchi0_ + dtchi;    // Update tchi using divided difference
    } else {
      // tchi = sinh(-1/n * log(tn))
      //      = sinh((1-1/n) * log(tn) - log(tn))
      //      = + sinh((1-1/n) * log(tn)) * cosh(log(tn))
      //        - cosh((1-1/n) * log(tn)) * sinh(log(tn))
      // (1-1/n) = - nc^2/(n*(1+n))
      // cosh(log(tn)) = (tn + 1/tn)/2; sinh(log(tn)) = (tn - 1/tn)/2
      real
        tn = tnm1 + 1 == 0 ? epsx_ : tnm1 + 1,
        sh = sinh( -Math::sq(nc_)/(n_ * (1 + n_)) *
                   (2 * tn > 1 ? log1p(tnm1) : log(tn)) );
      tchi = sh * (tn + 1/tn)/2 - hyp(sh) * (tnm1 * (tn + 1)/tn)/2;
    }

    // log(t) = -asinh(tan(chi)) = -psi
    gamma = atan2(nx, y1);
    real
      tphi = Math::tauf(tchi, es_),
      scbet = hyp(fm_ * tphi), scchi = hyp(tchi),
      lam = n_ != 0 ? gamma / n_ : x / y1;
    lat = Math::atand(sign_ * tphi);
    lon = lam / Math::degree();
    lon = Math::AngNormalize(lon + Math::AngNormalize(lon0));
    k = k0_ * (scbet/scbet0_) /
      (exp(nc_ != 0 ? - (Math::sq(nc_)/(1 + n_)) * dpsi : 0)
       * (tchi >= 0 ? scchi + tchi : 1 / (scchi - tchi)) / (scchi0_ + tchi0_));
    gamma /= sign_ * Math::degree();
  }

  void LambertConformalConic::SetScale(real lat, real k) {
    if (!(isfinite(k) && k > 0))
      throw GeographicErr("Scale is not positive");
    if (!(fabs (lat) <= 90))
      throw GeographicErr ("Latitude for SetScale not in [-90d, +90d]");
    if (fabs(lat) == 90 && !(nc_ == 0 && lat * n_ > 0))
      throw GeographicErr("Incompatible polar latitude in SetScale");
    real x, y, gamma, kold;
    Forward(0, lat, 0, x, y, gamma, kold);
    k /= kold;
    scale_ *= k;
    k0_ *= k;
  }

} // namespace GeographicLib
