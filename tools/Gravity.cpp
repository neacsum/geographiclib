/**
 * \file Gravity.cpp
 * \brief Command line utility for evaluating gravity fields
 *
 * Copyright (c) Charles Karney (2011-2022) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 *
 * See the <a href="Gravity.1.html">man page</a> for usage information.
 **********************************************************************/

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <GeographicLib/GravityModel.hpp>
#include <GeographicLib/GravityCircle.hpp>
#include <GeographicLib/DMS.hpp>
#include <GeographicLib/Utility.hpp>

int usage (int retval, bool brief);

int main(int argc, const char* const argv[]) {
  try {
    using namespace GeographicLib;
    Utility::set_digits();
    bool verbose = false, longfirst = false;
    std::string dir;
    std::string model = GravityModel::DefaultGravityName();
    std::string istring, ifile, ofile, cdelim;
    char lsep = ';';
    real lat = 0, h = 0;
    bool circle = false;
    int prec = -1, Nmax = -1, Mmax = -1;
    enum {
      GRAVITY = 0,
      DISTURBANCE = 1,
      ANOMALY = 2,
      UNDULATION = 3,
    };
    unsigned mode = GRAVITY;
    for (int m = 1; m < argc; ++m) {
      std::string arg(argv[m]);
      if (arg == "-n") {
        if (++m == argc) return usage(1, true);
        model = argv[m];
      } else if (arg == "-d") {
        if (++m == argc) return usage(1, true);
        dir = argv[m];
      } else if (arg == "-N") {
        if (++m == argc) return usage(1, true);
        try {
          Nmax = Utility::val<int>(std::string(argv[m]));
          if (Nmax < 0) {
            std::cerr << "Maximum degree " << argv[m] << " is negative\n";
            return 1;
          }
        }
        catch (const std::exception&) {
          std::cerr << "Precision " << argv[m] << " is not a number\n";
          return 1;
        }
      } else if (arg == "-M") {
        if (++m == argc) return usage(1, true);
        try {
          Mmax = Utility::val<int>(std::string(argv[m]));
          if (Mmax < 0) {
            std::cerr << "Maximum order " << argv[m] << " is negative\n";
            return 1;
          }
        }
        catch (const std::exception&) {
          std::cerr << "Precision " << argv[m] << " is not a number\n";
          return 1;
        }
      } else if (arg == "-G")
        mode = GRAVITY;
      else if (arg == "-D")
        mode = DISTURBANCE;
      else if (arg == "-A")
        mode = ANOMALY;
      else if (arg == "-H")
        mode = UNDULATION;
      else if (arg == "-c") {
        if (m + 2 >= argc) return usage(1, true);
        try {
          using std::fabs;
          DMS::flag ind;
          lat = DMS::Decode(std::string(argv[++m]), ind);
          if (ind == DMS::LONGITUDE)
            throw GeographicErr("Bad hemisphere letter on latitude");
          if (!(fabs(lat) <= 90))
            throw GeographicErr("Latitude not in [-90d, +90d]");
          h = Utility::val<real>(std::string(argv[++m]));
          circle = true;
        }
        catch (const std::exception& e) {
          std::cerr << "Error decoding argument of " << arg << ": "
                    << e.what() << "\n";
          return 1;
        }
      } else if (arg == "-w")
        longfirst = !longfirst;
      else if (arg == "-p") {
        if (++m == argc) return usage(1, true);
        try {
          prec = Utility::val<int>(std::string(argv[m]));
        }
        catch (const std::exception&) {
          std::cerr << "Precision " << argv[m] << " is not a number\n";
          return 1;
        }
      } else if (arg == "-v")
        verbose = true;
      else if (arg == "--input-string") {
        if (++m == argc) return usage(1, true);
        istring = argv[m];
      } else if (arg == "--input-file") {
        if (++m == argc) return usage(1, true);
        ifile = argv[m];
      } else if (arg == "--output-file") {
        if (++m == argc) return usage(1, true);
        ofile = argv[m];
      } else if (arg == "--line-separator") {
        if (++m == argc) return usage(1, true);
        if (std::string(argv[m]).size() != 1) {
          std::cerr << "Line separator must be a single character\n";
          return 1;
        }
        lsep = argv[m][0];
      } else if (arg == "--comment-delimiter") {
        if (++m == argc) return usage(1, true);
        cdelim = argv[m];
      } else if (arg == "--version") {
        std::cout << argv[0] << ": GeographicLib version "
                  << GEOGRAPHICLIB_VERSION_STRING << "\n";
        return 0;
      } else {
        int retval = usage(!(arg == "-h" || arg == "--help"), arg != "--help");
        if (arg == "-h")
          std::cout<< "\nDefault gravity path = \""
                   << GravityModel::DefaultGravityPath()
                   << "\"\nDefault gravity name = \""
                   << GravityModel::DefaultGravityName()
                   << "\"\n";
        return retval;
      }
    }

    if (!ifile.empty() && !istring.empty()) {
      std::cerr << "Cannot specify --input-string and --input-file together\n";
      return 1;
    }
    if (ifile == "-") ifile.clear();
    std::ifstream infile;
    std::istringstream instring;
    if (!ifile.empty()) {
      infile.open(ifile.c_str());
      if (!infile.is_open()) {
        std::cerr << "Cannot open " << ifile << " for reading\n";
        return 1;
      }
    } else if (!istring.empty()) {
      std::string::size_type m = 0;
      while (true) {
        m = istring.find(lsep, m);
        if (m == std::string::npos)
          break;
        istring[m] = '\n';
      }
      instring.str(istring);
    }
    std::istream* input = !ifile.empty() ? &infile :
      (!istring.empty() ? &instring : &std::cin);

    std::ofstream outfile;
    if (ofile == "-") ofile.clear();
    if (!ofile.empty()) {
      outfile.open(ofile.c_str());
      if (!outfile.is_open()) {
        std::cerr << "Cannot open " << ofile << " for writing\n";
        return 1;
      }
    }
    std::ostream* output = !ofile.empty() ? &outfile : &std::cout;

    switch (mode) {
    case GRAVITY:
      prec = std::min(16 + Math::extra_digits(), prec < 0 ? 5 : prec);
      break;
    case DISTURBANCE:
    case ANOMALY:
      prec = std::min(14 + Math::extra_digits(), prec < 0 ? 3 : prec);
      break;
    case UNDULATION:
    default:
      prec = std::min(12 + Math::extra_digits(), prec < 0 ? 4 : prec);
      break;
    }
    int retval = 0;
    try {
      using std::isfinite;
      const GravityModel g(model, dir, Nmax, Mmax);
      if (circle) {
        if (!isfinite(h))
          throw GeographicErr("Bad height");
        else if (mode == UNDULATION && h != 0)
          throw GeographicErr("Height should be zero for geoid undulations");
      }
      if (verbose) {
        std::cerr << "Gravity file: " << g.GravityFile()      << "\n"
                  << "Name: "         << g.GravityModelName() << "\n"
                  << "Description: "  << g.Description()      << "\n"
                  << "Date & Time: "  << g.DateTime()         << "\n";
      }
      unsigned mask = (mode == GRAVITY ? GravityModel::GRAVITY :
                       (mode == DISTURBANCE ? GravityModel::DISTURBANCE :
                        (mode == ANOMALY ? GravityModel::SPHERICAL_ANOMALY :
                         GravityModel::GEOID_HEIGHT))); // mode == UNDULATION
      const GravityCircle c(circle ? g.Circle(lat, h, mask) : GravityCircle());
      std::string s, eol, stra, strb;
      std::istringstream str;
      while (std::getline(*input, s)) {
        try {
          eol = "\n";
          if (!cdelim.empty()) {
            std::string::size_type m = s.find(cdelim);
            if (m != std::string::npos) {
              eol = " " + s.substr(m) + "\n";
              s = s.substr(0, m);
            }
          }
          str.clear(); str.str(s);
          real lon;
          if (circle) {
            if (!(str >> strb))
              throw GeographicErr("Incomplete input: " + s);
            DMS::flag ind;
            lon = DMS::Decode(strb, ind);
            if (ind == DMS::LATITUDE)
              throw GeographicErr("Bad hemisphere letter on " + strb);
          } else {
            if (!(str >> stra >> strb))
              throw GeographicErr("Incomplete input: " + s);
            DMS::DecodeLatLon(stra, strb, lat, lon, longfirst);
            h = 0;
            if (!(str >> h))    // h is optional
              str.clear();
            if (mode == UNDULATION && h != 0)
                throw GeographicErr("Height must be zero for geoid heights");
          }
          if (str >> stra)
            throw GeographicErr("Extra junk in input: " + s);
          switch (mode) {
          case GRAVITY:
            {
              real gx, gy, gz;
              if (circle) {
                c.Gravity(lon, gx, gy, gz);
              } else {
                g.Gravity(lat, lon, h, gx, gy, gz);
              }
              *output << Utility::str(gx, prec) << " "
                      << Utility::str(gy, prec) << " "
                      << Utility::str(gz, prec) << eol;
            }
            break;
          case DISTURBANCE:
            {
              real deltax, deltay, deltaz;
              if (circle) {
                c.Disturbance(lon, deltax, deltay, deltaz);
              } else {
                g.Disturbance(lat, lon, h, deltax, deltay, deltaz);
              }
              // Convert to mGals
              *output << Utility::str(deltax * 100000, prec) << " "
                      << Utility::str(deltay * 100000, prec) << " "
                      << Utility::str(deltaz * 100000, prec)
                      << eol;
            }
            break;
          case ANOMALY:
            {
              real Dg01, xi, eta;
              if (circle)
                c.SphericalAnomaly(lon, Dg01, xi, eta);
              else
                g.SphericalAnomaly(lat, lon, h, Dg01, xi, eta);
              Dg01 *= 100000;   // Convert to mGals
              xi *= 3600;       // Convert to arcsecs
              eta *= 3600;
              *output << Utility::str(Dg01, prec) << " "
                      << Utility::str(xi, prec) << " "
                      << Utility::str(eta, prec) << eol;
            }
            break;
          case UNDULATION:
          default:
            {
              real N = circle ? c.GeoidHeight(lon) : g.GeoidHeight(lat, lon);
              *output << Utility::str(N, prec) << eol;
            }
            break;
          }
        }
        catch (const std::exception& e) {
          *output << "ERROR: " << e.what() << "\n";
          retval = 1;
        }
      }
    }
    catch (const std::exception& e) {
      std::cerr << "Error reading " << model << ": " << e.what() << "\n";
      retval = 1;
    }
    return retval;
  }
  catch (const std::exception& e) {
    std::cerr << "Caught exception: " << e.what() << "\n";
    return 1;
  }
  catch (...) {
    std::cerr << "Caught unknown exception\n";
    return 1;
  }
}

int usage (int retval, bool brief) {
  if (brief)
    (retval ? std::cerr : std::cout) << "Usage:\n"
    "    Gravity [ -n name ] [ -d dir ] [ -N Nmax ] [ -M Mmax ] [ -G | -D | -A |\n"
    "    -H ] [ -c lat h ] [ -w ] [ -p prec ] [ -v ] [ --comment-delimiter\n"
    "    commentdelim ] [ --version | -h | --help ] [ --input-file infile |\n"
    "    --input-string instring ] [ --line-separator linesep ] [ --output-file\n"
    "    outfile ]\n"
    "\n"
    "For full documentation type:\n"
    "    Gravity --help\n"
    "or visit:\n"
    "    https://geographiclib.sourceforge.io/C++/2.4/Gravity.1.html\n";
  else
    (retval ? std::cerr : std::cout) << "Man page:\n"
    "NAME\n"
    "       Gravity -- compute the earth's gravity field\n"
    "\n"
    "SYNOPSIS\n"
    "       Gravity [ -n name ] [ -d dir ] [ -N Nmax ] [ -M Mmax ] [ -G | -D | -A |\n"
    "       -H ] [ -c lat h ] [ -w ] [ -p prec ] [ -v ] [ --comment-delimiter\n"
    "       commentdelim ] [ --version | -h | --help ] [ --input-file infile |\n"
    "       --input-string instring ] [ --line-separator linesep ] [ --output-file\n"
    "       outfile ]\n"
    "\n"
    "DESCRIPTION\n"
    "       Gravity reads in positions on standard input and prints out the\n"
    "       gravitational field on standard output.\n"
    "\n"
    "       The input line is of the form lat lon h.  lat and lon are the latitude\n"
    "       and longitude expressed as decimal degrees or degrees, minutes, and\n"
    "       seconds; for details on the allowed formats for latitude and longitude,\n"
    "       see the \"GEOGRAPHIC COORDINATES\" section of GeoConvert(1).  h is the\n"
    "       height above the ellipsoid in meters; this quantity is optional and\n"
    "       defaults to 0.  Alternatively, the gravity field can be computed at\n"
    "       various points on a circle of latitude (constant lat and h) via the -c\n"
    "       option; in this case only the longitude should be given on the input\n"
    "       lines.  The quantities printed out are governed by the -G (default),\n"
    "       -D, -A, or -H options.\n"
    "\n"
    "       All the supported gravity models, except for grs80, use WGS84 as the\n"
    "       reference ellipsoid a = 6378137 m, f = 1/298.257223563, omega =\n"
    "       7292115e-11 rad/s, and GM = 3986004.418e8 m^3/s^2.\n"
    "\n"
    "OPTIONS\n"
    "       -n name\n"
    "           use gravity field model name instead of the default \"egm96\".  See\n"
    "           \"MODELS\".\n"
    "\n"
    "       -d dir\n"
    "           read gravity models from dir instead of the default.  See \"MODELS\".\n"
    "\n"
    "       -N Nmax\n"
    "           limit the degree of the model to Nmax.\n"
    "\n"
    "       -M Mmax\n"
    "           limit the order of the model to Mmax.\n"
    "\n"
    "       -G  compute the acceleration due to gravity (including the centrifugal\n"
    "           acceleration due the the earth's rotation) g.  The output consists\n"
    "           of gx gy gz (all in m/s^2), where the x, y, and z components are in\n"
    "           easterly, northerly, and up directions, respectively.  Usually gz\n"
    "           is negative.\n"
    "\n"
    "       -D  compute the gravity disturbance delta = g - gamma, where gamma is\n"
    "           the \"normal\" gravity due to the reference ellipsoid .  The output\n"
    "           consists of deltax deltay deltaz (all in mGal, 1 mGal = 10^-5\n"
    "           m/s^2), where the x, y, and z components are in easterly,\n"
    "           northerly, and up directions, respectively.  Note that deltax = gx,\n"
    "           because gammax = 0.\n"
    "\n"
    "       -A  computes the gravitational anomaly.  The output consists of 3 items\n"
    "           Dg01 xi eta, where Dg01 is in mGal (1 mGal = 10^-5 m/s^2) and xi\n"
    "           and eta are in arcseconds.  The gravitational anomaly compares the\n"
    "           gravitational field g at P with the normal gravity gamma at Q where\n"
    "           the P is vertically above Q and the gravitational potential at P\n"
    "           equals the normal potential at Q.  Dg01 gives the difference in the\n"
    "           magnitudes of these two vectors and xi and eta give the difference\n"
    "           in their directions (as northerly and easterly components).  The\n"
    "           calculation uses a spherical approximation to match the results of\n"
    "           the NGA's synthesis programs.\n"
    "\n"
    "       -H  compute the height of the geoid above the reference ellipsoid (in\n"
    "           meters).  In this case, h should be zero.  The results accurately\n"
    "           match the results of the NGA's synthesis programs.  GeoidEval(1)\n"
    "           can compute geoid heights much more quickly by interpolating on a\n"
    "           grid of precomputed results; however the results from GeoidEval(1)\n"
    "           are only accurate to a few millimeters.\n"
    "\n"
    "       -c lat h\n"
    "           evaluate the field on a circle of latitude given by lat and h\n"
    "           instead of reading these quantities from the input lines.  In this\n"
    "           case, Gravity can calculate the field considerably more quickly.\n"
    "           If geoid heights are being computed (the -H option), then h must be\n"
    "           zero.\n"
    "\n"
    "       -w  toggle the longitude first flag (it starts off); if the flag is on,\n"
    "           then on input and output, longitude precedes latitude (except that,\n"
    "           on input, this can be overridden by a hemisphere designator, N, S,\n"
    "           E, W).\n"
    "\n"
    "       -p prec\n"
    "           set the output precision to prec.  By default prec is 5 for\n"
    "           acceleration due to gravity, 3 for the gravity disturbance and\n"
    "           anomaly, and 4 for the geoid height.\n"
    "\n"
    "       -v  print information about the gravity model on standard error before\n"
    "           processing the input.\n"
    "\n"
    "       --comment-delimiter commentdelim\n"
    "           set the comment delimiter to commentdelim (e.g., \"#\" or \"//\").  If\n"
    "           set, the input lines will be scanned for this delimiter and, if\n"
    "           found, the delimiter and the rest of the line will be removed prior\n"
    "           to processing and subsequently appended to the output line\n"
    "           (separated by a space).\n"
    "\n"
    "       --version\n"
    "           print version and exit.\n"
    "\n"
    "       -h  print usage, the default gravity path and name, and exit.\n"
    "\n"
    "       --help\n"
    "           print full documentation and exit.\n"
    "\n"
    "       --input-file infile\n"
    "           read input from the file infile instead of from standard input; a\n"
    "           file name of \"-\" stands for standard input.\n"
    "\n"
    "       --input-string instring\n"
    "           read input from the string instring instead of from standard input.\n"
    "           All occurrences of the line separator character (default is a\n"
    "           semicolon) in instring are converted to newlines before the reading\n"
    "           begins.\n"
    "\n"
    "       --line-separator linesep\n"
    "           set the line separator character to linesep.  By default this is a\n"
    "           semicolon.\n"
    "\n"
    "       --output-file outfile\n"
    "           write output to the file outfile instead of to standard output; a\n"
    "           file name of \"-\" stands for standard output.\n"
    "\n"
    "MODELS\n"
    "       Gravity computes the gravity field using one of the following models\n"
    "\n"
    "           egm84, earth gravity model 1984.  See\n"
    "             https://earth-info.nga.mil/index.php?dir=wgs84&action=wgs84#tab_egm84\n"
    "           egm96, earth gravity model 1996.  See\n"
    "             https://earth-info.nga.mil/index.php?dir=wgs84&action=wgs84#tab_egm96\n"
    "           egm2008, earth gravity model 2008.  See\n"
    "             https://earth-info.nga.mil/index.php?dir=wgs84&action=wgs84#tab_egm2008\n"
    "           wgs84, world geodetic system 1984.  This returns the normal\n"
    "             gravity for the WGS84 ellipsoid.\n"
    "           grs80, geodetic reference system 1980.  This returns the normal\n"
    "             gravity for the GRS80 ellipsoid.\n"
    "\n"
    "       These models approximate the gravitation field above the surface of the\n"
    "       earth.  By default, the \"egm96\" gravity model is used.  This may\n"
    "       changed by setting the environment variable\n"
    "       \"GEOGRAPHICLIB_GRAVITY_NAME\" or with the -n option.\n"
    "\n"
    "       The gravity models will be loaded from a directory specified at compile\n"
    "       time.  This may changed by setting the environment variables\n"
    "       \"GEOGRAPHICLIB_GRAVITY_PATH\" or \"GEOGRAPHICLIB_DATA\", or with the -d\n"
    "       option.  The -h option prints the default gravity path and name.  Use\n"
    "       the -v option to ascertain the full path name of the data file.\n"
    "\n"
    "       Instructions for downloading and installing gravity models are\n"
    "       available at\n"
    "       <https://geographiclib.sourceforge.io/C++/doc/gravity.html#gravityinst>.\n"
    "\n"
    "ENVIRONMENT\n"
    "       GEOGRAPHICLIB_GRAVITY_NAME\n"
    "           Override the compile-time default gravity name of \"egm96\".  The -h\n"
    "           option reports the value of GEOGRAPHICLIB_GRAVITY_NAME, if defined,\n"
    "           otherwise it reports the compile-time value.  If the -n name option\n"
    "           is used, then name takes precedence.\n"
    "\n"
    "       GEOGRAPHICLIB_GRAVITY_PATH\n"
    "           Override the compile-time default gravity path.  This is typically\n"
    "           \"/usr/local/share/GeographicLib/gravity\" on Unix-like systems and\n"
    "           \"C:/ProgramData/GeographicLib/gravity\" on Windows systems.  The -h\n"
    "           option reports the value of GEOGRAPHICLIB_GRAVITY_PATH, if defined,\n"
    "           otherwise it reports the compile-time value.  If the -d dir option\n"
    "           is used, then dir takes precedence.\n"
    "\n"
    "       GEOGRAPHICLIB_DATA\n"
    "           Another way of overriding the compile-time default gravity path.\n"
    "           If it is set (and if GEOGRAPHICLIB_GRAVITY_PATH is not set), then\n"
    "           $GEOGRAPHICLIB_DATA/gravity is used.\n"
    "\n"
    "ERRORS\n"
    "       An illegal line of input will print an error message to standard output\n"
    "       beginning with \"ERROR:\" and causes Gravity to return an exit code of 1.\n"
    "       However, an error does not cause Gravity to terminate; following lines\n"
    "       will be converted.\n"
    "\n"
    "EXAMPLES\n"
    "       The gravity field from EGM2008 at the top of Mount Everest\n"
    "\n"
    "           echo 27:59:17N 86:55:32E 8820 | Gravity -n egm2008\n"
    "           => -0.00001 0.00103 -9.76782\n"
    "\n"
    "SEE ALSO\n"
    "       GeoConvert(1), GeoidEval(1), geographiclib-get-gravity(8).\n"
    "\n"
    "AUTHOR\n"
    "       Gravity was written by Charles Karney.\n"
    "\n"
    "HISTORY\n"
    "       Gravity was added to GeographicLib,\n"
    "       <https://geographiclib.sourceforge.io>, in version 1.16.\n"
    ;
  return retval;
}
