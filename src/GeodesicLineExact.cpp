/**
 * \file GeodesicLineExact.cpp
 * \brief Implementation for GeographicLib::GeodesicLineExact class
 *
 * Copyright (c) Charles Karney (2012-2022) <karney@alum.mit.edu> and licensed
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

#include <GeographicLib/GeodesicLineExact.hpp>

namespace GeographicLib {

  using namespace std;

  void GeodesicLineExact::LineInit(const GeodesicExact& g,
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
    e2_ = g.e2_;
    nC4_ = g.nC4_;
    // Always allow latitude and azimuth and unrolling of longitude
    caps_ = caps | LATITUDE | AZIMUTH | LONG_UNROLL;

    real cbet1, sbet1;
    Math::sincosd(Math::AngRound(lat1_), sbet1, cbet1); sbet1 *= f1_;
    // Ensure cbet1 = +epsilon at poles
    Math::norm(sbet1, cbet1); cbet1 = fmax(tiny_, cbet1);
    dn1_ = (f_ >= 0 ? sqrt(1 + g.ep2_ * Math::sq(sbet1)) :
            sqrt(1 - e2_ * Math::sq(cbet1)) / f1_);

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
    // Without normalization we have schi1 = somg1.
    cchi1_ = f1_ * dn1_ * comg1_;
    Math::norm(ssig1_, csig1_); // sig1 in (-pi, pi]
    // Math::norm(somg1_, comg1_); -- don't need to normalize!
    // Math::norm(_schi1, cchi1_); -- don't need to normalize!

    k2_ = Math::sq(calp0_) * g.ep2_;
    eE_.Reset(-k2_, -g.ep2_, 1 + k2_, 1 + g.ep2_);

    if (caps_ & CAP_E) {
      eE0_ = eE_.E() / (Math::pi() / 2);
      eE1_ = eE_.deltaE(ssig1_, csig1_, dn1_);
      real s = sin(eE1_), c = cos(eE1_);
      // tau1 = sig1 + B11
      stau1_ = ssig1_ * c + csig1_ * s;
      ctau1_ = csig1_ * c - ssig1_ * s;
      // Not necessary because Einv inverts E
      //    eE1_ = -eE_.deltaEinv(stau1_, ctau1_);
    }

    if (caps_ & CAP_D) {
      dD0_ = eE_.D() / (Math::pi() / 2);
      dD1_ = eE_.deltaD(ssig1_, csig1_, dn1_);
    }

    if (caps_ & CAP_H) {
      hH0_ = eE_.H() / (Math::pi() / 2);
      hH1_ = eE_.deltaH(ssig1_, csig1_, dn1_);
    }

    if (caps_ & CAP_C4) {
      // Multiplier = a^2 * e^2 * cos(alpha0) * sin(alpha0)
      aA4_ = Math::sq(a_) * calp0_ * salp0_ * e2_;
      if (aA4_ == 0)
        bB41_ = 0;
      else {
        GeodesicExact::I4Integrand i4(g.ep2_, k2_);
        cC4a_.resize(nC4_);
        g.fft_.transform(i4, cC4a_.data());
        bB41_ = DST::integral(ssig1_, csig1_, cC4a_.data(), nC4_);
      }
    }

    a13_ = s13_ = Math::NaN();
  }

  GeodesicLineExact::GeodesicLineExact(const GeodesicExact& g,
                                       real lat1, real lon1, real azi1,
                                       unsigned caps) {
    azi1 = Math::AngNormalize(azi1);
    real salp1, calp1;
    // Guard against underflow in salp0.  Also -0 is converted to +0.
    Math::sincosd(Math::AngRound(azi1), salp1, calp1);
    LineInit(g, lat1, lon1, azi1, salp1, calp1, caps);
  }

  GeodesicLineExact::GeodesicLineExact(const GeodesicExact& g,
                                       real lat1, real lon1,
                                       real azi1, real salp1, real calp1,
                                       unsigned caps,
                                       bool arcmode, real s13_a13) {
    LineInit(g, lat1, lon1, azi1, salp1, calp1, caps);
    GenSetDistance(arcmode, s13_a13);
  }

  real GeodesicLineExact::GenPosition(bool arcmode, real s12_a12,
                                            unsigned outmask,
                                            real& lat2, real& lon2, real& azi2,
                                            real& s12, real& m12,
                                            real& M12, real& M21,
                                            real& S12) const {
    outmask &= caps_ & OUT_MASK;
    if (!( Init() && (arcmode || (caps_ & (OUT_MASK & DISTANCE_IN))) ))
      // Uninitialized or impossible distance calculation requested
      return Math::NaN();

    // Avoid warning about uninitialized B12.
    real sig12, ssig12, csig12, E2 = 0, AB1 = 0;
    if (arcmode) {
      // Interpret s12_a12 as spherical arc length
      sig12 = s12_a12 * Math::degree();
      Math::sincosd(s12_a12, ssig12, csig12);
    } else {
      // Interpret s12_a12 as distance
      real
        tau12 = s12_a12 / (b_ * eE0_),
        s = sin(tau12),
        c = cos(tau12);
      // tau2 = tau1 + tau12
      E2 = - eE_.deltaEinv(stau1_ * c + ctau1_ * s, ctau1_ * c - stau1_ * s);
      sig12 = tau12 - (E2 - eE1_);
      ssig12 = sin(sig12);
      csig12 = cos(sig12);
    }

    real ssig2, csig2, sbet2, cbet2, salp2, calp2;
    // sig2 = sig1 + sig12
    ssig2 = ssig1_ * csig12 + csig1_ * ssig12;
    csig2 = csig1_ * csig12 - ssig1_ * ssig12;
    real dn2 = eE_.Delta(ssig2, csig2);
    if (outmask & (DISTANCE | REDUCEDLENGTH | GEODESICSCALE)) {
      if (arcmode) {
        E2 = eE_.deltaE(ssig2, csig2, dn2);
      }
      AB1 = eE0_ * (E2 - eE1_);
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
      s12 = arcmode ? b_ * (eE0_ * sig12 + AB1) : s12_a12;

    if (outmask & LONGITUDE) {
      real somg2 = salp0_ * ssig2, comg2 = csig2,  // No need to normalize
        E = copysign(real(1), salp0_);       // east-going?
      // Without normalization we have schi2 = somg2.
      real cchi2 =  f1_ * dn2 *  comg2;
      real chi12 = outmask & LONG_UNROLL
        ? E * (sig12
               - (atan2(    ssig2, csig2) - atan2(    ssig1_, csig1_))
               + (atan2(E * somg2, cchi2) - atan2(E * somg1_, cchi1_)))
        : atan2(somg2 * cchi1_ - cchi2 * somg1_,
                cchi2 * cchi1_ + somg2 * somg1_);
      real lam12 = chi12 -
        e2_/f1_ * salp0_ * hH0_ *
        (sig12 + (eE_.deltaH(ssig2, csig2, dn2) - hH1_));
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
      real J12 = k2_ * dD0_ * (sig12 + (eE_.deltaD(ssig2, csig2, dn2) - dD1_));
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
      real B42 = aA4_ == 0 ? 0 :
        DST::integral(ssig2, csig2, cC4a_.data(), nC4_);
      real salp12, calp12;
      if (calp0_ == 0 || salp0_ == 0) {
        // alp12 = alp2 - alp1, used in atan2 so no need to normalize
        salp12 = salp2 * calp1_ - calp2 * salp1_;
        calp12 = calp2 * calp1_ + salp2 * salp1_;
        // We used to include here some patch up code that purported to deal
        // with nearly meridional geodesics properly.  However, this turned out
        // to be wrong once salp1_ = -0 was allowed (via
        // GeodesicExact::InverseLine).  In fact, the calculation of {s,c}alp12
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

  void GeodesicLineExact::SetDistance(real s13) {
    s13_ = s13;
    real t;
    // This will set a13_ to NaN if the GeodesicLineExact doesn't have the
    // DISTANCE_IN capability.
    a13_ = GenPosition(false, s13_, 0u, t, t, t, t, t, t, t, t);
  }

  void GeodesicLineExact::SetArc(real a13) {
    a13_ = a13;
    // In case the GeodesicLineExact doesn't have the DISTANCE capability.
    s13_ = Math::NaN();
    real t;
    GenPosition(true, a13_, DISTANCE, t, t, t, s13_, t, t, t, t);
  }

  void GeodesicLineExact::GenSetDistance(bool arcmode, real s13_a13) {
    arcmode ? SetArc(s13_a13) : SetDistance(s13_a13);
  }

} // namespace GeographicLib
