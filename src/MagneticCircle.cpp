/**
 * \file MagneticCircle.cpp
 * \brief Implementation for GeographicLib::MagneticCircle class
 *
 * Copyright (c) Charles Karney (2011-2021) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#include <GeographicLib/MagneticCircle.hpp>
#include <fstream>
#include <sstream>
#include <GeographicLib/Geocentric.hpp>

namespace GeographicLib {

  using namespace std;

  void MagneticCircle::FieldGeocentric(real slam, real clam,
                                       real& BX, real& BY, real& BZ,
                                       real& BXt, real& BYt, real& BZt) const {
    real BXc = 0, BYc = 0, BZc = 0;
    circ0_(slam, clam, BX, BY, BZ);
    circ1_(slam, clam, BXt, BYt, BZt);
    if (constterm_)
      circ2_(slam, clam, BXc, BYc, BZc);
    if (interpolate_) {
      BXt = (BXt - BX) / dt0_;
      BYt = (BYt - BY) / dt0_;
      BZt = (BZt - BZ) / dt0_;
    }
    BX += t1_ * BXt + BXc;
    BY += t1_ * BYt + BYc;
    BZ += t1_ * BZt + BZc;

    BXt *= - a_;
    BYt *= - a_;
    BZt *= - a_;

    BX *= - a_;
    BY *= - a_;
    BZ *= - a_;
  }

  void MagneticCircle::FieldGeocentric(real lon,
                                       real& BX, real& BY, real& BZ,
                                       real& BXt, real& BYt, real& BZt) const {
    real slam, clam;
    Math::sincosd(lon, slam, clam);
    FieldGeocentric(slam, clam, BX, BY, BZ, BXt, BYt, BZt);
  }

  void MagneticCircle::Field(real lon, bool diffp,
                             real& Bx, real& By, real& Bz,
                             real& Bxt, real& Byt, real& Bzt) const {
    real slam, clam;
    Math::sincosd(lon, slam, clam);
    real M[Geocentric::dim2_];
    Geocentric::Rotation(sphi_, cphi_, slam, clam, M);
    real BX, BY, BZ, BXt, BYt, BZt; // Components in geocentric basis
    FieldGeocentric(slam, clam, BX, BY, BZ, BXt, BYt, BZt);
    if (diffp)
      Geocentric::Unrotate(M, BXt, BYt, BZt, Bxt, Byt, Bzt);
    Geocentric::Unrotate(M, BX, BY, BZ, Bx, By, Bz);
  }

} // namespace GeographicLib
