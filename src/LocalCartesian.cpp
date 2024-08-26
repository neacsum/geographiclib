/**
 * \file LocalCartesian.cpp
 * \brief Implementation for GeographicLib::LocalCartesian class
 *
 * Copyright (c) Charles Karney (2008-2015) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#include <GeographicLib/LocalCartesian.hpp>

namespace GeographicLib {

  using namespace std;

  void LocalCartesian::Reset(real lat0, real lon0, real h0) {
    lat0_ = Math::LatFix(lat0);
    lon0_ = Math::AngNormalize(lon0);
    h0_ = h0;
    earth_.Forward(lat0_, lon0_, h0_, x0_, y0_, z0_);
    real sphi, cphi, slam, clam;
    Math::sincosd(lat0_, sphi, cphi);
    Math::sincosd(lon0_, slam, clam);
    Geocentric::Rotation(sphi, cphi, slam, clam, r_);
  }

  void LocalCartesian::MatrixMultiply(real M[dim2_]) const {
    // M = r' . M
    real t[dim2_];
    copy(M, M + dim2_, t);
    for (size_t i = 0; i < dim2_; ++i) {
      size_t row = i / dim_, col = i % dim_;
      M[i] = r_[row] * t[col] + r_[row+3] * t[col+3] + r_[row+6] * t[col+6];
    }
  }

  void LocalCartesian::IntForward(real lat, real lon, real h,
                                  real& x, real& y, real& z,
                                  real M[dim2_]) const {
    real xc, yc, zc;
    earth_.IntForward(lat, lon, h, xc, yc, zc, M);
    xc -= x0_; yc -= y0_; zc -= z0_;
    x = r_[0] * xc + r_[3] * yc + r_[6] * zc;
    y = r_[1] * xc + r_[4] * yc + r_[7] * zc;
    z = r_[2] * xc + r_[5] * yc + r_[8] * zc;
    if (M)
      MatrixMultiply(M);
  }

  void LocalCartesian::IntReverse(real x, real y, real z,
                                  real& lat, real& lon, real& h,
                                  real M[dim2_]) const {
    real
      xc = x0_ + r_[0] * x + r_[1] * y + r_[2] * z,
      yc = y0_ + r_[3] * x + r_[4] * y + r_[5] * z,
      zc = z0_ + r_[6] * x + r_[7] * y + r_[8] * z;
    earth_.IntReverse(xc, yc, zc, lat, lon, h, M);
    if (M)
      MatrixMultiply(M);
  }

} // namespace GeographicLib
