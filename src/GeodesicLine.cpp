/**
 * \file GeodesicLine.cpp
 * \brief Implementation for GeographicLib::GeodesicLine class
 *
 * Copyright (c) Charles Karney (2009-2023) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 *
 * This is a reformulation of the geodesic problem.  The notation is as
 * follows:
 * - at a general point (no suffix or 1 or 2 as suffix)
 *   - phi = latitude
 *   - beta = latitude on auxiliary sphere
 *   - omega = longitude on auxiliary sphere
 *   - lambda = longitude
 *   - alpha = azimuth of great circle
 *   - sigma = arc length along great circle
 *   - s = distance
 *   - tau = scaled distance (= sigma at multiples of pi/2)
 * - at northwards equator crossing
 *   - beta = phi = 0
 *   - omega = lambda = 0
 *   - alpha = alpha0
 *   - sigma = s = 0
 * - a 12 suffix means a difference, e.g., s12 = s2 - s1.
 * - s and c prefixes mean sin and cos
 **********************************************************************/

#include <GeographicLib/GeodesicLine.hpp>

namespace GeographicLib {

  using namespace std;

  void GeodesicLine::LineInit(const Geodesic& g,
                              real lat1, real lon1,
                              real azi1, real salp1, real calp1,
                              unsigned caps) {
    tiny_ = g.tiny_;
    lat1_ = Math::LatFix(lat1);
    lon1_ = lon1;
    azi1_ = azi1;
    salp1_ = salp1;
    calp1_ = calp1;
    a_ = g.a_;
    f_ = g.f_;
    b_ = g.b_;
    c2_ = g.c2_;
    f1_ = g.f1_;
    // Always allow latitude and azimuth and unrolling of longitude
    caps_ = caps | LATITUDE | AZIMUTH | LONG_UNROLL;

    real cbet1, sbet1;
    Math::sincosd(Math::AngRound(lat1_), sbet1, cbet1); sbet1 *= f1_;
    // Ensure cbet1 = +epsilon at poles
    Math::norm(sbet1, cbet1); cbet1 = fmax(tiny_, cbet1);
    dn1_ = sqrt(1 + g.ep2_ * Math::sq(sbet1));

    // Evaluate alp0 from sin(alp1) * cos(bet1) = sin(alp0),
    salp0_ = salp1_ * cbet1; // alp0 in [0, pi/2 - |bet1|]
    // Alt: calp0 = hypot(sbet1, calp1 * cbet1).  The following
    // is slightly better (consider the case salp1 = 0).
    calp0_ = hypot(calp1_, salp1_ * sbet1);
    // Evaluate sig with tan(bet1) = tan(sig1) * cos(alp1).
    // sig = 0 is nearest northward crossing of equator.
    // With bet1 = 0, alp1 = pi/2, we have sig1 = 0 (equatorial line).
    // With bet1 =  pi/2, alp1 = -pi, sig1 =  pi/2
    // With bet1 = -pi/2, alp1 =  0 , sig1 = -pi/2
    // Evaluate omg1 with tan(omg1) = sin(alp0) * tan(sig1).
    // With alp0 in (0, pi/2], quadrants for sig and omg coincide.
    // No atan2(0,0) ambiguity at poles since cbet1 = +epsilon.
    // With alp0 = 0, omg1 = 0 for alp1 = 0, omg1 = pi for alp1 = pi.
    ssig1_ = sbet1; somg1_ = salp0_ * sbet1;
    csig1_ = comg1_ = sbet1 != 0 || calp1_ != 0 ? cbet1 * calp1_ : 1;
    Math::norm(ssig1_, csig1_); // sig1 in (-pi, pi]
    // Math::norm(somg1_, comg1_); -- don't need to normalize!

    a13_ = s13_ = Math::NaN();
    exact_ = g.exact_;
    if (exact_) {
      lineexact_.LineInit(g.geodexact_, lat1, lon1, azi1, salp1, calp1, caps);
      return;
    }

    k2_ = Math::sq(calp0_) * g.ep2_;
    real eps = k2_ / (2 * (1 + sqrt(1 + k2_)) + k2_);

    if (caps_ & CAP_C1) {
      aA1m1_ = Geodesic::A1m1f(eps);
      Geodesic::C1f(eps, cC1a_);
      bB11_ = Geodesic::SinCosSeries(true, ssig1_, csig1_, cC1a_, nC1_);
      real s = sin(bB11_), c = cos(bB11_);
      // tau1 = sig1 + B11
      stau1_ = ssig1_ * c + csig1_ * s;
      ctau1_ = csig1_ * c - ssig1_ * s;
      // Not necessary because C1pa reverts C1a
      //    bB11_ = -SinCosSeries(true, stau1_, ctau1_, cC1pa_, nC1p_);
    }

    if (caps_ & CAP_C1p)
      Geodesic::C1pf(eps, cC1pa_);

    if (caps_ & CAP_C2) {
      aA2m1_ = Geodesic::A2m1f(eps);
      Geodesic::C2f(eps, cC2a_);
      bB21_ = Geodesic::SinCosSeries(true, ssig1_, csig1_, cC2a_, nC2_);
    }

    if (caps_ & CAP_C3) {
      g.C3f(eps, cC3a_);
      aA3c_ = -f_ * salp0_ * g.A3f(eps);
      bB31_ = Geodesic::SinCosSeries(true, ssig1_, csig1_, cC3a_, nC3_-1);
    }

    if (caps_ & CAP_C4) {
      g.C4f(eps, _cC4a);
      // Multiplier = a^2 * e^2 * cos(alpha0) * sin(alpha0)
      aA4_ = Math::sq(a_) * calp0_ * salp0_ * g.e2_;
      bB41_ = Geodesic::SinCosSeries(false, ssig1_, csig1_, _cC4a, nC4_);
    }

  }

  GeodesicLine::GeodesicLine(const Geodesic& g,
                             real lat1, real lon1, real azi1,
                             unsigned caps) {
    azi1 = Math::AngNormalize(azi1);
    real salp1, calp1;
    // Guard against underflow in salp0.  Also -0 is converted to +0.
    Math::sincosd(Math::AngRound(azi1), salp1, calp1);
    LineInit(g, lat1, lon1, azi1, salp1, calp1, caps);
  }

  GeodesicLine::GeodesicLine(const Geodesic& g,
                             real lat1, real lon1,
                             real azi1, real salp1, real calp1,
                             unsigned caps, bool arcmode, real s13_a13) {
    LineInit(g, lat1, lon1, azi1, salp1, calp1, caps);
    GenSetDistance(arcmode, s13_a13);
  }

  real GeodesicLine::GenPosition(bool arcmode, real s12_a12,
                                       unsigned outmask,
                                       real& lat2, real& lon2, real& azi2,
                                       real& s12, real& m12,
                                       real& M12, real& M21,
                                       real& S12) const {
    if (exact_)
      return lineexact_.GenPosition(arcmode, s12_a12, outmask,
                                    lat2, lon2, azi2,
                                    s12, m12, M12, M21, S12);
    outmask &= caps_ & OUT_MASK;
    if (!( Init() && (arcmode || (caps_ & (OUT_MASK & DISTANCE_IN))) ))
      // Uninitialized or impossible distance calculation requested
      return Math::NaN();

    // Avoid warning about uninitialized B12.
    real sig12, ssig12, csig12, B12 = 0, AB1 = 0;
    if (arcmode) {
      // Interpret s12_a12 as spherical arc length
      sig12 = s12_a12 * Math::degree();
      Math::sincosd(s12_a12, ssig12, csig12);
    } else {
      // Interpret s12_a12 as distance
      real
        tau12 = s12_a12 / (b_ * (1 + aA1m1_)),
        s = sin(tau12),
        c = cos(tau12);
      // tau2 = tau1 + tau12
      B12 = - Geodesic::SinCosSeries(true,
                                     stau1_ * c + ctau1_ * s,
                                     ctau1_ * c - stau1_ * s,
                                     cC1pa_, nC1p_);
      sig12 = tau12 - (B12 - bB11_);
      ssig12 = sin(sig12); csig12 = cos(sig12);
      if (fabs(f_) > 0.01) {
        // Reverted distance series is inaccurate for |f| > 1/100, so correct
        // sig12 with 1 Newton iteration.  The following table shows the
        // approximate maximum error for a = WGS_a() and various f relative to
        // GeodesicExact.
        //     erri = the error in the inverse solution (nm)
        //     errd = the error in the direct solution (series only) (nm)
        //     errda = the error in the direct solution
        //             (series + 1 Newton) (nm)
        //
        //       f     erri  errd errda
        //     -1/5    12e6 1.2e9  69e6
        //     -1/10  123e3  12e6 765e3
        //     -1/20   1110 108e3  7155
        //     -1/50  18.63 200.9 27.12
        //     -1/100 18.63 23.78 23.37
        //     -1/150 18.63 21.05 20.26
        //      1/150 22.35 24.73 25.83
        //      1/100 22.35 25.03 25.31
        //      1/50  29.80 231.9 30.44
        //      1/20   5376 146e3  10e3
        //      1/10  829e3  22e6 1.5e6
        //      1/5   157e6 3.8e9 280e6
        real
          ssig2 = ssig1_ * csig12 + csig1_ * ssig12,
          csig2 = csig1_ * csig12 - ssig1_ * ssig12;
        B12 = Geodesic::SinCosSeries(true, ssig2, csig2, cC1a_, nC1_);
        real serr = (1 + aA1m1_) * (sig12 + (B12 - bB11_)) - s12_a12 / b_;
        sig12 = sig12 - serr / sqrt(1 + k2_ * Math::sq(ssig2));
        ssig12 = sin(sig12); csig12 = cos(sig12);
        // Update B12 below
      }
    }

    real ssig2, csig2, sbet2, cbet2, salp2, calp2;
    // sig2 = sig1 + sig12
    ssig2 = ssig1_ * csig12 + csig1_ * ssig12;
    csig2 = csig1_ * csig12 - ssig1_ * ssig12;
    real dn2 = sqrt(1 + k2_ * Math::sq(ssig2));
    if (outmask & (DISTANCE | REDUCEDLENGTH | GEODESICSCALE)) {
      if (arcmode || fabs(f_) > 0.01)
        B12 = Geodesic::SinCosSeries(true, ssig2, csig2, cC1a_, nC1_);
      AB1 = (1 + aA1m1_) * (B12 - bB11_);
    }
    // sin(bet2) = cos(alp0) * sin(sig2)
    sbet2 = calp0_ * ssig2;
    // Alt: cbet2 = hypot(csig2, salp0 * ssig2);
    cbet2 = hypot(salp0_, calp0_ * csig2);
    if (cbet2 == 0)
      // I.e., salp0 = 0, csig2 = 0.  Break the degeneracy in this case
      cbet2 = csig2 = tiny_;
    // tan(alp0) = cos(sig2)*tan(alp2)
    salp2 = salp0_; calp2 = calp0_ * csig2; // No need to normalize

    if (outmask & DISTANCE)
      s12 = arcmode ? b_ * ((1 + aA1m1_) * sig12 + AB1) : s12_a12;

    if (outmask & LONGITUDE) {
      // tan(omg2) = sin(alp0) * tan(sig2)
      real somg2 = salp0_ * ssig2, comg2 = csig2,  // No need to normalize
        E = copysign(real(1), salp0_);       // east-going?
      // omg12 = omg2 - omg1
      real omg12 = outmask & LONG_UNROLL
        ? E * (sig12
               - (atan2(    ssig2, csig2) - atan2(    ssig1_, csig1_))
               + (atan2(E * somg2, comg2) - atan2(E * somg1_, comg1_)))
        : atan2(somg2 * comg1_ - comg2 * somg1_,
                comg2 * comg1_ + somg2 * somg1_);
      real lam12 = omg12 + aA3c_ *
        ( sig12 + (Geodesic::SinCosSeries(true, ssig2, csig2, cC3a_, nC3_-1)
                   - bB31_));
      real lon12 = lam12 / Math::degree();
      lon2 = outmask & LONG_UNROLL ? lon1_ + lon12 :
        Math::AngNormalize(Math::AngNormalize(lon1_) +
                           Math::AngNormalize(lon12));
    }

    if (outmask & LATITUDE)
      lat2 = Math::atan2d(sbet2, f1_ * cbet2);

    if (outmask & AZIMUTH)
      azi2 = Math::atan2d(salp2, calp2);

    if (outmask & (REDUCEDLENGTH | GEODESICSCALE)) {
      real
        B22 = Geodesic::SinCosSeries(true, ssig2, csig2, cC2a_, nC2_),
        AB2 = (1 + aA2m1_) * (B22 - bB21_),
        J12 = (aA1m1_ - aA2m1_) * sig12 + (AB1 - AB2);
      if (outmask & REDUCEDLENGTH)
        // Add parens around (csig1_ * ssig2) and (ssig1_ * csig2) to ensure
        // accurate cancellation in the case of coincident points.
        m12 = b_ * ((dn2 * (csig1_ * ssig2) - dn1_ * (ssig1_ * csig2))
                    - csig1_ * csig2 * J12);
      if (outmask & GEODESICSCALE) {
        real t = k2_ * (ssig2 - ssig1_) * (ssig2 + ssig1_) / (dn1_ + dn2);
        M12 = csig12 + (t *  ssig2 -  csig2 * J12) * ssig1_ / dn1_;
        M21 = csig12 - (t * ssig1_ - csig1_ * J12) *  ssig2 /  dn2;
      }
    }

    if (outmask & AREA) {
      real
        B42 = Geodesic::SinCosSeries(false, ssig2, csig2, _cC4a, nC4_);
      real salp12, calp12;
      if (calp0_ == 0 || salp0_ == 0) {
        // alp12 = alp2 - alp1, used in atan2 so no need to normalize
        salp12 = salp2 * calp1_ - calp2 * salp1_;
        calp12 = calp2 * calp1_ + salp2 * salp1_;
        // We used to include here some patch up code that purported to deal
        // with nearly meridional geodesics properly.  However, this turned out
        // to be wrong once salp1_ = -0 was allowed (via
        // Geodesic::InverseLine).  In fact, the calculation of {s,c}alp12
        // was already correct (following the IEEE rules for handling signed
        // zeros).  So the patch up code was unnecessary (as well as
        // dangerous).
      } else {
        // tan(alp) = tan(alp0) * sec(sig)
        // tan(alp2-alp1) = (tan(alp2) -tan(alp1)) / (tan(alp2)*tan(alp1)+1)
        // = calp0 * salp0 * (csig1-csig2) / (salp0^2 + calp0^2 * csig1*csig2)
        // If csig12 > 0, write
        //   csig1 - csig2 = ssig12 * (csig1 * ssig12 / (1 + csig12) + ssig1)
        // else
        //   csig1 - csig2 = csig1 * (1 - csig12) + ssig12 * ssig1
        // No need to normalize
        salp12 = calp0_ * salp0_ *
          (csig12 <= 0 ? csig1_ * (1 - csig12) + ssig12 * ssig1_ :
           ssig12 * (csig1_ * ssig12 / (1 + csig12) + ssig1_));
        calp12 = Math::sq(salp0_) + Math::sq(calp0_) * csig1_ * csig2;
      }
      S12 = c2_ * atan2(salp12, calp12) + aA4_ * (B42 - bB41_);
    }

    return arcmode ? s12_a12 : sig12 / Math::degree();
  }

  void GeodesicLine::SetDistance(real s13) {
    s13_ = s13;
    real t;
    // This will set a13_ to NaN if the GeodesicLine doesn't have the
    // DISTANCE_IN capability.
    a13_ = GenPosition(false, s13_, 0u, t, t, t, t, t, t, t, t);
  }

  void GeodesicLine::SetArc(real a13) {
    a13_ = a13;
    // In case the GeodesicLine doesn't have the DISTANCE capability.
    s13_ = Math::NaN();
    real t;
    GenPosition(true, a13_, DISTANCE, t, t, t, s13_, t, t, t, t);
  }

  void GeodesicLine::GenSetDistance(bool arcmode, real s13_a13) {
    arcmode ? SetArc(s13_a13) : SetDistance(s13_a13);
  }

} // namespace GeographicLib
