/**
 * \file CircularEngine.cpp
 * \brief Implementation for GeographicLib::CircularEngine class
 *
 * Copyright (c) Charles Karney (2011) <karney@alum.mit.edu> and licensed under
 * the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#include <GeographicLib/CircularEngine.hpp>

namespace GeographicLib {

  using namespace std;

  real CircularEngine::Value(bool gradp, real sl, real cl,
                                   real& gradx, real& grady, real& gradz) const
  {
    gradp = gradp_ && gradp;
    const vector<real>& root( SphericalEngine::sqrttable() );

    // Initialize outer sum
    real vc  = 0, vc2  = 0, vs  = 0, vs2  = 0;   // v [N + 1], v [N + 2]
    // vr, vt, vl and similar w variable accumulate the sums for the
    // derivatives wrt r, theta, and lambda, respectively.
    real vrc = 0, vrc2 = 0, vrs = 0, vrs2 = 0;   // vr[N + 1], vr[N + 2]
    real vtc = 0, vtc2 = 0, vts = 0, vts2 = 0;   // vt[N + 1], vt[N + 2]
    real vlc = 0, vlc2 = 0, vls = 0, vls2 = 0;   // vl[N + 1], vl[N + 2]
    for (int m = mM_; m >= 0; --m) {             // m = M .. 0
      // Now Sc[m] = wc, Ss[m] = ws
      // Sc'[m] = wtc, Ss'[m] = wtc
      if (m) {
        real v, A, B;           // alpha[m], beta[m + 1]
        switch (norm_) {
        case FULL:
          v = root[2] * root[2 * m + 3] / root[m + 1];
          A = cl * v * uq_;
          B = - v * root[2 * m + 5] / (root[8] * root[m + 2]) * uq2_;
          break;
        case SCHMIDT:
          v = root[2] * root[2 * m + 1] / root[m + 1];
          A = cl * v * uq_;
          B = - v * root[2 * m + 3] / (root[8] * root[m + 2]) * uq2_;
          break;
        default:
          A = B = 0;
        }
        v = A * vc  + B * vc2  +  wc_[m] ; vc2  = vc ; vc  = v;
        v = A * vs  + B * vs2  +  ws_[m] ; vs2  = vs ; vs  = v;
        if (gradp) {
          v = A * vrc + B * vrc2 +  wrc_[m]; vrc2 = vrc; vrc = v;
          v = A * vrs + B * vrs2 +  wrs_[m]; vrs2 = vrs; vrs = v;
          v = A * vtc + B * vtc2 +  wtc_[m]; vtc2 = vtc; vtc = v;
          v = A * vts + B * vts2 +  wts_[m]; vts2 = vts; vts = v;
          v = A * vlc + B * vlc2 + m*ws_[m]; vlc2 = vlc; vlc = v;
          v = A * vls + B * vls2 - m*wc_[m]; vls2 = vls; vls = v;
        }
      } else {
        real A, B, qs;
        switch (norm_) {
        case FULL:
          A = root[3] * uq_;       // F[1]/(q*cl) or F[1]/(q*sl)
          B = - root[15]/2 * uq2_; // beta[1]/q
          break;
        case SCHMIDT:
          A = uq_;
          B = - root[3]/2 * uq2_;
          break;
        default:
          A = B = 0;
        }
        qs = q_ / SphericalEngine::scale();
        vc = qs * (wc_[m] + A * (cl * vc + sl * vs ) + B * vc2);
        if (gradp) {
          qs /= _r;
          // The components of the gradient in circular coordinates are
          // r: dV/dr
          // theta: 1/r * dV/dtheta
          // lambda: 1/(r*u) * dV/dlambda
          vrc =    - qs * (wrc_[m] + A * (cl * vrc + sl * vrs) + B * vrc2);
          vtc =      qs * (wtc_[m] + A * (cl * vtc + sl * vts) + B * vtc2);
          vlc = qs / _u * (          A * (cl * vlc + sl * vls) + B * vlc2);
        }
      }
    }

    if (gradp) {
      // Rotate into cartesian (geocentric) coordinates
      gradx = cl * (_u * vrc + _t * vtc) - sl * vlc;
      grady = sl * (_u * vrc + _t * vtc) + cl * vlc;
      gradz =           _t * vrc - _u * vtc                ;
    }
    return vc;
  }

} // namespace GeographicLib
