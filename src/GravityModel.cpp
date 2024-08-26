/**
 * \file GravityModel.cpp
 * \brief Implementation for GeographicLib::GravityModel class
 *
 * Copyright (c) Charles Karney (2011-2023) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#include <GeographicLib/GravityModel.hpp>
#include <fstream>
#include <limits>
#include <GeographicLib/SphericalEngine.hpp>
#include <GeographicLib/GravityCircle.hpp>
#include <GeographicLib/Utility.hpp>

#if !defined(GEOGRAPHICLIB_DATA)
#  if defined(_WIN32)
#    define GEOGRAPHICLIB_DATA "C:/ProgramData/GeographicLib"
#  else
#    define GEOGRAPHICLIB_DATA "/usr/local/share/GeographicLib"
#  endif
#endif

#if !defined(GEOGRAPHICLIB_GRAVITY_DEFAULT_NAME)
#  define GEOGRAPHICLIB_GRAVITY_DEFAULT_NAME "egm96"
#endif

#if defined(_MSC_VER)
// Squelch warnings about unsafe use of getenv
#  pragma warning (disable: 4996)
#endif

namespace GeographicLib {

  using namespace std;

  GravityModel::GravityModel(const std::string& name, const std::string& path,
                             int Nmax, int Mmax)
    : name_(name)
    , dir_(path)
    , description_("NONE")
    , date_("UNKNOWN")
    , amodel_(Math::NaN())
    , gGMmodel_(Math::NaN())
    , zeta0_(0)
    , corrmult_(1)
    , nmx_(-1)
    , mmx_(-1)
    , norm_(SphericalHarmonic::FULL)
  {
    if (dir_.empty())
      dir_ = DefaultGravityPath();
    bool truncate = Nmax >= 0 || Mmax >= 0;
    if (truncate) {
      if (Nmax >= 0 && Mmax < 0) Mmax = Nmax;
      if (Nmax < 0) Nmax = numeric_limits<int>::max();
      if (Mmax < 0) Mmax = numeric_limits<int>::max();
    }
    ReadMetadata(name_);
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
      int N, M;
      if (truncate) { N = Nmax; M = Mmax; }
      SphericalEngine::coeff::readcoeffs(coeffstr, N, M, cCx_, sSx_, truncate);
      if (!(N >= 0 && M >= 0))
        throw GeographicErr("Degree and order must be at least 0");
      if (cCx_[0] != 0)
        throw GeographicErr("The degree 0 term should be zero");
      cCx_[0] = 1;              // Include the 1/r term in the sum
      gravitational_ = SphericalHarmonic(cCx_, sSx_, N, N, M, amodel_, norm_);
      if (truncate) { N = Nmax; M = Mmax; }
      SphericalEngine::coeff::readcoeffs(coeffstr, N, M, cCC_, cCS_, truncate);
      if (N < 0) {
        N = M = 0;
        cCC_.resize(1, real(0));
      }
      cCC_[0] += zeta0_ / corrmult_;
      correction_ = SphericalHarmonic(cCC_, cCS_, N, N, M, real(1), norm_);
      int pos = int(coeffstr.tellg());
      coeffstr.seekg(0, ios::end);
      if (pos != coeffstr.tellg())
        throw GeographicErr("Extra data in " + coeff);
    }
    int nmx = gravitational_.Coefficients().nmx();
    nmx_ = max(nmx, correction_.Coefficients().nmx());
    mmx_ = max(gravitational_.Coefficients().mmx(),
               correction_.Coefficients().mmx());
    // Adjust the normalization of the normal potential to match the model.
    real mult = earth_.gGM_ / gGMmodel_;
    real amult = Math::sq(earth_.a_ / amodel_);
    // The 0th term in zonal_ should be is 1 + dzonal0_.  Instead set it to 1
    // to give exact cancellation with the (0,0) term in the model and account
    // for dzonal0_ separately.
    zonal_.clear(); zonal_.push_back(1);
    dzonal0_ = (earth_.MassConstant() - gGMmodel_) / gGMmodel_;
    for (int n = 2; n <= nmx; n += 2) {
      // Only include as many normal zonal terms as matter.  Figuring the limit
      // in this way works because the coefficients of the normal potential
      // (which is smooth) decay much more rapidly that the corresponding
      // coefficient of the model potential (which is bumpy).  Typically this
      // goes out to n = 18.
      mult *= amult;
      real
        r = cCx_[n],                                       // the model term
        s = - mult * earth_.Jn(n) / sqrt(real(2 * n + 1)), // the normal term
        t = r - s;                                         // the difference
      if (t == r)               // the normal term is negligible
        break;
      zonal_.push_back(0);      // index = n - 1; the odd terms are 0
      zonal_.push_back(s);
    }
    int nmx1 = int(zonal_.size()) - 1;
    disturbing_ = SphericalHarmonic1(cCx_, sSx_,
                                     gravitational_.Coefficients().N(),
                                     nmx, gravitational_.Coefficients().mmx(),
                                     zonal_,
                                     zonal_, // This is not accessed!
                                     nmx1, nmx1, 0,
                                     amodel_,
                                     SphericalHarmonic1::normalization(norm_));
  }

  void GravityModel::ReadMetadata(const string& name) {
    const char* spaces = " \t\n\v\f\r";
    filename_ = dir_ + "/" + name + ".egm";
    ifstream metastr(filename_.c_str());
    if (!metastr.good())
      throw GeographicErr("Cannot open " + filename_);
    string line;
    getline(metastr, line);
    if (!(line.size() >= 6 && line.substr(0,5) == "EGMF-"))
      throw GeographicErr(filename_ + " does not contain EGMF-n signature");
    string::size_type n = line.find_first_of(spaces, 5);
    if (n != string::npos)
      n -= 5;
    string version(line, 5, n);
    if (version != "1")
      throw GeographicErr("Unknown version in " + filename_ + ": " + version);
    string key, val;
    real a = Math::NaN(), GM = a, omega = a, f = a, J2 = a;
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
      else if (key == "ModelRadius")
        amodel_ = Utility::val<real>(val);
      else if (key == "ModelMass")
        gGMmodel_ = Utility::val<real>(val);
      else if (key == "AngularVelocity")
        omega = Utility::val<real>(val);
      else if (key == "ReferenceRadius")
        a = Utility::val<real>(val);
      else if (key == "ReferenceMass")
        GM = Utility::val<real>(val);
      else if (key == "Flattening")
        f = Utility::fract<real>(val);
      else if (key == "DynamicalFormFactor")
        J2 = Utility::fract<real>(val);
      else if (key == "HeightOffset")
        zeta0_ = Utility::fract<real>(val);
      else if (key == "CorrectionMultiplier")
        corrmult_ = Utility::fract<real>(val);
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
    if (!(isfinite(amodel_) && amodel_ > 0))
      throw GeographicErr("Model radius must be positive");
    if (!(isfinite(gGMmodel_) && gGMmodel_ > 0))
      throw GeographicErr("Model mass constant must be positive");
    if (!(isfinite(corrmult_) && corrmult_ > 0))
      throw GeographicErr("Correction multiplier must be positive");
    if (!(isfinite(zeta0_)))
      throw GeographicErr("Height offset must be finite");
    if (int(id_.size()) != idlength_)
      throw GeographicErr("Invalid ID");
    if (isfinite(f) && isfinite(J2))
      throw GeographicErr("Cannot specify both f and J2");
    earth_ = NormalGravity(a, GM, omega,
                           isfinite(f) ? f : J2, isfinite(f));
  }

  real GravityModel::InternalT(real X, real Y, real Z,
                                     real& deltaX, real& deltaY, real& deltaZ,
                                     bool gradp, bool correct) const {
    // If correct, then produce the correct T = W - U.  Otherwise, neglect the
    // n = 0 term (which is proportial to the difference in the model and
    // reference values of GM).
    if (dzonal0_ == 0)
      // No need to do the correction
      correct = false;
    real T, invR = correct ? 1 / hypot(hypot(X, Y), Z) : 1;
    if (gradp) {
      // initial values to suppress warnings
      deltaX = deltaY = deltaZ = 0;
      T = disturbing_(-1, X, Y, Z, deltaX, deltaY, deltaZ);
      real f = gGMmodel_ / amodel_;
      deltaX *= f;
      deltaY *= f;
      deltaZ *= f;
      if (correct) {
        invR = gGMmodel_ * dzonal0_ * invR * invR * invR;
        deltaX += X * invR;
        deltaY += Y * invR;
        deltaZ += Z * invR;
      }
    } else
      T = disturbing_(-1, X, Y, Z);
    T = (T / amodel_ - (correct ? dzonal0_ : 0) * invR) * gGMmodel_;
    return T;
  }

  real GravityModel::V(real X, real Y, real Z,
                             real& GX, real& GY, real& GZ) const {
    real
      Vres = gravitational_(X, Y, Z, GX, GY, GZ),
      f = gGMmodel_ / amodel_;
    Vres *= f;
    GX *= f;
    GY *= f;
    GZ *= f;
    return Vres;
  }

  real GravityModel::W(real X, real Y, real Z,
                             real& gX, real& gY, real& gZ) const {
    real fX, fY,
      Wres = V(X, Y, Z, gX, gY, gZ) + earth_.Phi(X, Y, fX, fY);
    gX += fX;
    gY += fY;
    return Wres;
  }

  void GravityModel::SphericalAnomaly(real lat, real lon, real h,
                                      real& Dg01, real& xi, real& eta) const {
    real X, Y, Z, M[Geocentric::dim2_];
    earth_.Earth().IntForward(lat, lon, h, X, Y, Z, M);
    real
      deltax, deltay, deltaz,
      T = InternalT(X, Y, Z, deltax, deltay, deltaz, true, false),
      clam = M[3], slam = -M[0],
      P = hypot(X, Y),
      R = hypot(P, Z),
      // theta is geocentric latitude
      ctheta = R != 0 ? P / R : M[7],
      stheta = R != 0 ? Z / R : M[8];
    // Rotate cartesian into spherical coordinates
    real MC[Geocentric::dim2_];
    Geocentric::Rotation(stheta, ctheta, slam, clam, MC);
    Geocentric::Unrotate(MC, deltax, deltay, deltaz, deltax, deltay, deltaz);
    // H+M, Eq 2-151c
    Dg01 = - deltaz - 2 * T / R;
    real gammaX, gammaY, gammaZ;
    earth_.U(X, Y, Z, gammaX, gammaY, gammaZ);
    real gamma = hypot( hypot(gammaX, gammaY), gammaZ);
    xi  = -(deltay/gamma) / Math::degree();
    eta = -(deltax/gamma) / Math::degree();
  }

  real GravityModel::GeoidHeight(real lat, real lon) const
  {
    real X, Y, Z;
    earth_.Earth().IntForward(lat, lon, 0, X, Y, Z, NULL);
    real
      gamma0 = earth_.SurfaceGravity(lat),
      dummy,
      T = InternalT(X, Y, Z, dummy, dummy, dummy, false, false),
      invR = 1 / hypot(hypot(X, Y), Z),
      correction = corrmult_ * correction_(invR * X, invR * Y, invR * Z);
    // zeta0_ has been included in correction_
    return T/gamma0 + correction;
  }

  real GravityModel::Gravity(real lat, real lon, real h,
                                   real& gx, real& gy, real& gz) const {
    real X, Y, Z, M[Geocentric::dim2_];
    earth_.Earth().IntForward(lat, lon, h, X, Y, Z, M);
    real Wres = W(X, Y, Z, gx, gy, gz);
    Geocentric::Unrotate(M, gx, gy, gz, gx, gy, gz);
    return Wres;
  }
  real GravityModel::Disturbance(real lat, real lon, real h,
                                       real& deltax, real& deltay,
                                       real& deltaz) const {
    real X, Y, Z, M[Geocentric::dim2_];
    earth_.Earth().IntForward(lat, lon, h, X, Y, Z, M);
    real Tres = InternalT(X, Y, Z, deltax, deltay, deltaz, true, true);
    Geocentric::Unrotate(M, deltax, deltay, deltaz, deltax, deltay, deltaz);
    return Tres;
  }

  GravityCircle GravityModel::Circle(real lat, real h, unsigned caps) const {
    if (h != 0)
      // Disallow invoking GeoidHeight unless h is zero.
      caps &= ~(CAP_GAMMA0 | CAP_C);
    real X, Y, Z, M[Geocentric::dim2_];
    earth_.Earth().IntForward(lat, 0, h, X, Y, Z, M);
    // Y = 0, cphi = M[7], sphi = M[8];
    real
      invR = 1 / hypot(X, Z),
      gamma0 = (caps & CAP_GAMMA0 ?earth_.SurfaceGravity(lat)
                : Math::NaN()),
      fx, fy, fz, gamma;
    if (caps & CAP_GAMMA) {
      earth_.U(X, Y, Z, fx, fy, fz); // fy = 0
      gamma = hypot(fx, fz);
    } else
      gamma = Math::NaN();
    earth_.Phi(X, Y, fx, fy);
    return GravityCircle(GravityCircle::mask(caps),
                         earth_.a_, earth_.f_, lat, h, Z, X, M[7], M[8],
                         amodel_, gGMmodel_, dzonal0_, corrmult_,
                         gamma0, gamma, fx,
                         caps & CAP_G ?
                         gravitational_.Circle(X, Z, true) :
                         CircularEngine(),
                         // N.B. If CAP_DELTA is set then CAP_T should be too.
                         caps & CAP_T ?
                         disturbing_.Circle(-1, X, Z, (caps&CAP_DELTA) != 0) :
                         CircularEngine(),
                         caps & CAP_C ?
                         correction_.Circle(invR * X, invR * Z, false) :
                         CircularEngine());
  }

  string GravityModel::DefaultGravityPath() {
    string path;
    char* gravitypath = getenv("GEOGRAPHICLIB_GRAVITY_PATH");
    if (gravitypath)
      path = string(gravitypath);
    if (!path.empty())
      return path;
    char* datapath = getenv("GEOGRAPHICLIB_DATA");
    if (datapath)
      path = string(datapath);
    return (!path.empty() ? path : string(GEOGRAPHICLIB_DATA)) + "/gravity";
  }

  string GravityModel::DefaultGravityName() {
    string name;
    char* gravityname = getenv("GEOGRAPHICLIB_GRAVITY_NAME");
    if (gravityname)
      name = string(gravityname);
    return !name.empty() ? name : string(GEOGRAPHICLIB_GRAVITY_DEFAULT_NAME);
  }

} // namespace GeographicLib
