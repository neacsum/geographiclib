/**
 * \file MagneticModel.cpp
 * \brief Implementation for GeographicLib::MagneticModel class
 *
 * Copyright (c) Charles Karney (2011-2021) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#include <GeographicLib/MagneticModel.hpp>
#include <fstream>
#include <GeographicLib/SphericalEngine.hpp>
#include <GeographicLib/MagneticCircle.hpp>
#include <GeographicLib/Utility.hpp>

#if !defined(GEOGRAPHICLIB_DATA)
#  if defined(_WIN32)
#    define GEOGRAPHICLIB_DATA "C:/ProgramData/GeographicLib"
#  else
#    define GEOGRAPHICLIB_DATA "/usr/local/share/GeographicLib"
#  endif
#endif

#if !defined(GEOGRAPHICLIB_MAGNETIC_DEFAULT_NAME)
#  define GEOGRAPHICLIB_MAGNETIC_DEFAULT_NAME "wmm2020"
#endif

#if defined(_MSC_VER)
// Squelch warnings about unsafe use of getenv
#  pragma warning (disable: 4996)
#endif

namespace GeographicLib {

  using namespace std;

  MagneticModel::MagneticModel(const std::string& name, const std::string& path,
                               const Geocentric& earth, int Nmax, int Mmax)
    : name_(name)
    , dir_(path)
    , description_("NONE")
    , date_("UNKNOWN")
    , t0_(Math::NaN())
    , dt0_(1)
    , tmin_(Math::NaN())
    , tmax_(Math::NaN())
    , a_(Math::NaN())
    , hmin_(Math::NaN())
    , hmax_(Math::NaN())
    , nNmodels_(1)
    , nNconstants_(0)
    , nmx_(-1)
    , mmx_(-1)
    , norm_(SphericalHarmonic::SCHMIDT)
    , earth_(earth)
  {
    if (dir_.empty())
      dir_ = DefaultMagneticPath();
    bool truncate = Nmax >= 0 || Mmax >= 0;
    if (truncate) {
      if (Nmax >= 0 && Mmax < 0) Mmax = Nmax;
      if (Nmax < 0) Nmax = numeric_limits<int>::max();
      if (Mmax < 0) Mmax = numeric_limits<int>::max();
    }
    ReadMetadata(name_);
    gG_.resize(nNmodels_ + 1 + nNconstants_);
    hH_.resize(nNmodels_ + 1 + nNconstants_);
    {
      string coeff = filename_ + ".cof";
      ifstream coeffstr(coeff.c_str(), ios::binary);
      if (!coeffstr.good())
        throw GeographicErr("Error opening " + coeff);
      char id[idlength_ + 1];
      coeffstr.read(id, idlength_);
      if (!coeffstr.good())
        throw GeographicErr("No header in " + coeff);
      id[idlength_] = '\0';
      if (id_ != string(id))
        throw GeographicErr("ID mismatch: " + id_ + " vs " + id);
      for (int i = 0; i < nNmodels_ + 1 + nNconstants_; ++i) {
        int N, M;
        if (truncate) { N = Nmax; M = Mmax; }
        SphericalEngine::coeff::readcoeffs(coeffstr, N, M, gG_[i], hH_[i],
                                           truncate);
        if (!(M < 0 || gG_[i][0] == 0))
          throw GeographicErr("A degree 0 term is not permitted");
        harm_.push_back(SphericalHarmonic(gG_[i], hH_[i], N, N, M, a_, norm_));
        nmx_ = max(nmx_, harm_.back().Coefficients().nmx());
        mmx_ = max(mmx_, harm_.back().Coefficients().mmx());
      }
      int pos = int(coeffstr.tellg());
      coeffstr.seekg(0, ios::end);
      if (pos != coeffstr.tellg())
        throw GeographicErr("Extra data in " + coeff);
    }
  }

  void MagneticModel::ReadMetadata(const string& name) {
    const char* spaces = " \t\n\v\f\r";
    filename_ = dir_ + "/" + name + ".wmm";
    ifstream metastr(filename_.c_str());
    if (!metastr.good())
      throw GeographicErr("Cannot open " + filename_);
    string line;
    getline(metastr, line);
    if (!(line.size() >= 6 && line.substr(0,5) == "WMMF-"))
      throw GeographicErr(filename_ + " does not contain WMMF-n signature");
    string::size_type n = line.find_first_of(spaces, 5);
    if (n != string::npos)
      n -= 5;
    string version(line, 5, n);
    if (!(version == "1" || version == "2"))
      throw GeographicErr("Unknown version in " + filename_ + ": " + version);
    string key, val;
    while (getline(metastr, line)) {
      if (!Utility::ParseLine(line, key, val))
        continue;
      // Process key words
      if (key == "Name")
        name_ = val;
      else if (key == "Description")
        description_ = val;
      else if (key == "ReleaseDate")
        date_ = val;
      else if (key == "Radius")
        a_ = Utility::val<real>(val);
      else if (key == "Type") {
        if (!(val == "Linear" || val == "linear"))
          throw GeographicErr("Only linear models are supported");
      } else if (key == "Epoch")
        t0_ = Utility::val<real>(val);
      else if (key == "DeltaEpoch")
        dt0_ = Utility::val<real>(val);
      else if (key == "NumModels")
        nNmodels_ = Utility::val<int>(val);
      else if (key == "NumConstants")
        nNconstants_ = Utility::val<int>(val);
      else if (key == "MinTime")
        tmin_ = Utility::val<real>(val);
      else if (key == "MaxTime")
        tmax_ = Utility::val<real>(val);
      else if (key == "MinHeight")
        hmin_ = Utility::val<real>(val);
      else if (key == "MaxHeight")
        hmax_ = Utility::val<real>(val);
      else if (key == "Normalization") {
        if (val == "FULL" || val == "Full" || val == "full")
          norm_ = SphericalHarmonic::FULL;
        else if (val == "SCHMIDT" || val == "Schmidt" || val == "schmidt")
          norm_ = SphericalHarmonic::SCHMIDT;
        else
          throw GeographicErr("Unknown normalization " + val);
      } else if (key == "ByteOrder") {
        if (val == "Big" || val == "big")
          throw GeographicErr("Only little-endian ordering is supported");
        else if (!(val == "Little" || val == "little"))
          throw GeographicErr("Unknown byte ordering " + val);
      } else if (key == "ID")
        id_ = val;
      // else unrecognized keywords are skipped
    }
    // Check values
    if (!(isfinite(a_) && a_ > 0))
      throw GeographicErr("Reference radius must be positive");
    if (!(t0_ > 0))
      throw GeographicErr("Epoch time not defined");
    if (tmin_ >= tmax_)
      throw GeographicErr("Min time exceeds max time");
    if (hmin_ >= hmax_)
      throw GeographicErr("Min height exceeds max height");
    if (int(id_.size()) != idlength_)
      throw GeographicErr("Invalid ID");
    if (nNmodels_ < 1)
      throw GeographicErr("NumModels must be positive");
    if (!(nNconstants_ == 0 || nNconstants_ == 1))
      throw GeographicErr("NumConstants must be 0 or 1");
    if (!(dt0_ > 0)) {
      if (nNmodels_ > 1)
        throw GeographicErr("DeltaEpoch must be positive");
      else
        dt0_ = 1;
    }
  }

  void MagneticModel::FieldGeocentric(real t, real X, real Y, real Z,
                                      real& BX, real& BY, real& BZ,
                                      real& BXt, real& BYt, real& BZt) const {
    t -= t0_;
    int n = max(min(int(floor(t / dt0_)), nNmodels_ - 1), 0);
    bool interpolate = n + 1 < nNmodels_;
    t -= n * dt0_;
    // Components in geocentric basis
    // initial values to suppress warning
    real BXc = 0, BYc = 0, BZc = 0;
    harm_[n](X, Y, Z, BX, BY, BZ);
    harm_[n + 1](X, Y, Z, BXt, BYt, BZt);
    if (nNconstants_)
      harm_[nNmodels_ + 1](X, Y, Z, BXc, BYc, BZc);
    if (interpolate) {
      // Convert to a time derivative
      BXt = (BXt - BX) / dt0_;
      BYt = (BYt - BY) / dt0_;
      BZt = (BZt - BZ) / dt0_;
    }
    BX += t * BXt + BXc;
    BY += t * BYt + BYc;
    BZ += t * BZt + BZc;

    BXt = BXt * - a_;
    BYt = BYt * - a_;
    BZt = BZt * - a_;

    BX *= - a_;
    BY *= - a_;
    BZ *= - a_;
  }

  void MagneticModel::Field(real t, real lat, real lon, real h, bool diffp,
                            real& Bx, real& By, real& Bz,
                            real& Bxt, real& Byt, real& Bzt) const {
    real X, Y, Z;
    real M[Geocentric::dim2_];
    earth_.IntForward(lat, lon, h, X, Y, Z, M);
    // Components in geocentric basis
    // initial values to suppress warning
    real BX = 0, BY = 0, BZ = 0, BXt = 0, BYt = 0, BZt = 0;
    FieldGeocentric(t, X, Y, Z, BX, BY, BZ, BXt, BYt, BZt);
    if (diffp)
      Geocentric::Unrotate(M, BXt, BYt, BZt, Bxt, Byt, Bzt);
    Geocentric::Unrotate(M, BX, BY, BZ, Bx, By, Bz);
  }

  MagneticCircle MagneticModel::Circle(real t, real lat, real h) const {
    real t1 = t - t0_;
    int n = max(min(int(floor(t1 / dt0_)), nNmodels_ - 1), 0);
    bool interpolate = n + 1 < nNmodels_;
    t1 -= n * dt0_;
    real X, Y, Z, M[Geocentric::dim2_];
    earth_.IntForward(lat, 0, h, X, Y, Z, M);
    // Y = 0, cphi = M[7], sphi = M[8];

    return (nNconstants_ == 0 ?
            MagneticCircle(a_, earth_.f_, lat, h, t,
                           M[7], M[8], t1, dt0_, interpolate,
                           harm_[n].Circle(X, Z, true),
                           harm_[n + 1].Circle(X, Z, true)) :
            MagneticCircle(a_, earth_.f_, lat, h, t,
                           M[7], M[8], t1, dt0_, interpolate,
                           harm_[n].Circle(X, Z, true),
                           harm_[n + 1].Circle(X, Z, true),
                           harm_[nNmodels_ + 1].Circle(X, Z, true)));
  }

  void MagneticModel::FieldComponents(real Bx, real By, real Bz,
                                      real Bxt, real Byt, real Bzt,
                                      real& H, real& F, real& D, real& I,
                                      real& Ht, real& Ft,
                                      real& Dt, real& It) {
    H = hypot(Bx, By);
    Ht = H != 0 ? (Bx * Bxt + By * Byt) / H : hypot(Bxt, Byt);
    D = H != 0 ? Math::atan2d(Bx, By) : Math::atan2d(Bxt, Byt);
    Dt = (H != 0 ? (By * Bxt - Bx * Byt) / Math::sq(H) : 0) / Math::degree();
    F = hypot(H, Bz);
    Ft = F != 0 ? (H * Ht + Bz * Bzt) / F : hypot(Ht, Bzt);
    I = F != 0 ? Math::atan2d(-Bz, H) : Math::atan2d(-Bzt, Ht);
    It = (F != 0 ? (Bz * Ht - H * Bzt) / Math::sq(F) : 0) / Math::degree();
  }

  string MagneticModel::DefaultMagneticPath() {
    string path;
    char* magneticpath = getenv("GEOGRAPHICLIB_MAGNETIC_PATH");
    if (magneticpath)
      path = string(magneticpath);
    if (!path.empty())
      return path;
    char* datapath = getenv("GEOGRAPHICLIB_DATA");
    if (datapath)
      path = string(datapath);
    return (!path.empty() ? path : string(GEOGRAPHICLIB_DATA)) + "/magnetic";
  }

  string MagneticModel::DefaultMagneticName() {
    string name;
    char* magneticname = getenv("GEOGRAPHICLIB_MAGNETIC_NAME");
    if (magneticname)
      name = string(magneticname);
    return !name.empty() ? name : string(GEOGRAPHICLIB_MAGNETIC_DEFAULT_NAME);
  }

} // namespace GeographicLib
