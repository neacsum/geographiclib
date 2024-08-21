/**
 * \file TransverseMercatorProj.cpp
 * \brief Command line utility for transverse Mercator projections
 *
 * Copyright (c) Charles Karney (2008-2023) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 *
 * See the <a href="TransverseMercatorProj.1.html">man page</a> for usage
 * information.
 **********************************************************************/

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <GeographicLib/TransverseMercator.hpp>
#include <GeographicLib/DMS.hpp>
#include <GeographicLib/Utility.hpp>

#if defined(_MSC_VER)
// Squelch warnings about potentially uninitialized local variables
#  pragma warning (disable: 4701)
#endif

int usage (int retval, bool brief);

int main(int argc, const char* const argv[]) {
  try {
    using namespace GeographicLib;
    typedef Math::real real;
    Utility::set_digits();
    bool exact = true, extended = false, reverse = false, longfirst = false;
    real
      a = Constants::WGS84_a(),
      f = Constants::WGS84_f(),
      k0 = Constants::UTM_k0(),
      lon0 = 0;
    int prec = 6;
    std::string istring, ifile, ofile, cdelim;
    char lsep = ';';

    for (int m = 1; m < argc; ++m) {
      std::string arg(argv[m]);
      if (arg == "-r")
        reverse = true;
      else if (arg == "-t") {
        exact = true;
        extended = true;
      } else if (arg == "-s") {
        exact = false;
        extended = false;
      } else if (arg == "-l") {
        if (++m >= argc) return usage(1, true);
        try {
          DMS::flag ind;
          lon0 = DMS::Decode(std::string(argv[m]), ind);
          if (ind == DMS::LATITUDE)
            throw GeographicErr("Bad hemisphere");
          lon0 = Math::AngNormalize(lon0);
        }
        catch (const std::exception& e) {
          std::cerr << "Error decoding argument of " << arg << ": "
                    << e.what() << "\n";
          return 1;
        }
      } else if (arg == "-k") {
        if (++m >= argc) return usage(1, true);
        try {
          k0 = Utility::val<real>(std::string(argv[m]));
        }
        catch (const std::exception& e) {
          std::cerr << "Error decoding argument of " << arg << ": "
                    << e.what() << "\n";
          return 1;
        }
      } else if (arg == "-e") {
        if (m + 2 >= argc) return usage(1, true);
        try {
          a = Utility::val<real>(std::string(argv[m + 1]));
          f = Utility::fract<real>(std::string(argv[m + 2]));
        }
        catch (const std::exception& e) {
          std::cerr << "Error decoding arguments of -e: " << e.what() << "\n";
          return 1;
        }
        m += 2;
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
      } else if (arg == "--input-string") {
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
      } else
        return usage(!(arg == "-h" || arg == "--help"), arg != "--help");
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

    const TransverseMercator TM(a, f, k0, exact, extended);

    // Max precision = 10: 0.1 nm in distance, 10^-15 deg (= 0.11 nm),
    // 10^-11 sec (= 0.3 nm).
    prec = std::min(10 + Math::extra_digits(), std::max(0, prec));
    std::string s, eol, stra, strb, strc;
    std::istringstream str;
    int retval = 0;
    std::cout << std::fixed;
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
        real lat, lon, x, y;
        if (!(str >> stra >> strb))
          throw GeographicErr("Incomplete input: " + s);
        if (reverse) {
          x = Utility::val<real>(stra);
          y = Utility::val<real>(strb);
        } else
          DMS::DecodeLatLon(stra, strb, lat, lon, longfirst);
        if (str >> strc)
          throw GeographicErr("Extraneous input: " + strc);
        real gamma, k;
        if (reverse) {
          TM.Reverse(lon0, x, y, lat, lon, gamma, k);
          *output << Utility::str(longfirst ? lon : lat, prec + 5) << " "
                  << Utility::str(longfirst ? lat : lon, prec + 5) << " "
                  << Utility::str(gamma, prec + 6) << " "
                  << Utility::str(k, prec + 6) << eol;
        } else {
          TM.Forward(lon0, lat, lon, x, y, gamma, k);
          *output << Utility::str(x, prec) << " "
                  << Utility::str(y, prec) << " "
                  << Utility::str(gamma, prec + 6) << " "
                  << Utility::str(k, prec + 6) << eol;
        }
      }
      catch (const std::exception& e) {
        *output << "ERROR: " << e.what() << "\n";
        retval = 1;
      }
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
    "    TransverseMercatorProj [ -s | -t ] [ -l lon0 ] [ -k k0 ] [ -r ] [ -e a\n"
    "    f ] [ -w ] [ -p prec ] [ --comment-delimiter commentdelim ] [ --version\n"
    "    | -h | --help ] [ --input-file infile | --input-string instring ] [\n"
    "    --line-separator linesep ] [ --output-file outfile ]\n"
    "\n"
    "For full documentation type:\n"
    "    TransverseMercatorProj --help\n"
    "or visit:\n"
    "    https://geographiclib.sourceforge.io/C++/2.4/TransverseMercatorProj.1.html\n";
  else
    (retval ? std::cerr : std::cout) << "Man page:\n"
    "NAME\n"
    "       TransverseMercatorProj -- perform transverse Mercator projection\n"
    "\n"
    "SYNOPSIS\n"
    "       TransverseMercatorProj [ -s | -t ] [ -l lon0 ] [ -k k0 ] [ -r ] [ -e a\n"
    "       f ] [ -w ] [ -p prec ] [ --comment-delimiter commentdelim ] [ --version\n"
    "       | -h | --help ] [ --input-file infile | --input-string instring ] [\n"
    "       --line-separator linesep ] [ --output-file outfile ]\n"
    "\n"
    "DESCRIPTION\n"
    "       Perform the transverse Mercator projections.  Convert geodetic\n"
    "       coordinates to transverse Mercator coordinates.  The central meridian\n"
    "       is given by lon0.  The longitude of origin is the equator.  The scale\n"
    "       on the central meridian is k0.  By default an implementation of the\n"
    "       exact transverse Mercator projection is used.\n"
    "\n"
    "       Geodetic coordinates are provided on standard input as a set of lines\n"
    "       containing (blank separated) latitude and longitude (decimal degrees or\n"
    "       degrees, minutes, seconds); for detils on the allowed formats for\n"
    "       latitude and longitude, see the \"GEOGRAPHIC COORDINATES\" section of\n"
    "       GeoConvert(1).  For each set of geodetic coordinates, the corresponding\n"
    "       projected easting, x, and northing, y, (meters) are printed on standard\n"
    "       output together with the meridian convergence gamma (degrees) and scale\n"
    "       k.  The meridian convergence is the bearing of grid north (the y axis)\n"
    "       measured clockwise from true north.\n"
    "\n"
    "OPTIONS\n"
    "       -s  use the sixth-order Krueger series approximation to the transverse\n"
    "           Mercator projection instead of the exact projection.\n"
    "\n"
    "       -t  use the exact algorithm with the \"EXTENDED DOMAIN\".\n"
    "\n"
    "       -l lon0\n"
    "           specify the longitude of origin lon0 (degrees, default 0).\n"
    "\n"
    "       -k k0\n"
    "           specify the scale k0 on the central meridian (default 0.9996).\n"
    "\n"
    "       -r  perform the reverse projection.  x and y are given on standard\n"
    "           input and each line of standard output gives latitude, longitude,\n"
    "           gamma, and k.\n"
    "\n"
    "       -e a f\n"
    "           specify the ellipsoid via the equatorial radius, a and the\n"
    "           flattening, f.  Setting f = 0 results in a sphere.  Specify f < 0\n"
    "           for a prolate ellipsoid.  A simple fraction, e.g., 1/297, is\n"
    "           allowed for f.  By default, the WGS84 ellipsoid is used, a =\n"
    "           6378137 m, f = 1/298.257223563.  If the exact algorithm is used, f\n"
    "           must be positive.\n"
    "\n"
    "       -w  on input and output, longitude precedes latitude (except that on\n"
    "           input this can be overridden by a hemisphere designator, N, S, E,\n"
    "           W).\n"
    "\n"
    "       -p prec\n"
    "           set the output precision to prec (default 6).  prec is the number\n"
    "           of digits after the decimal point for lengths (in meters).  For\n"
    "           latitudes and longitudes (in degrees), the number of digits after\n"
    "           the decimal point is prec + 5.  For the convergence (in degrees)\n"
    "           and scale, the number of digits after the decimal point is prec +\n"
    "           6.\n"
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
    "       -h  print usage and exit.\n"
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
    "EXTENDED DOMAIN\n"
    "       The exact transverse Mercator projection has a branch point on the\n"
    "       equator at longitudes (relative to lon0) of +/- (1 - e) 90 = 82.636...,\n"
    "       where e is the eccentricity of the ellipsoid.  The standard convention\n"
    "       for handling this branch point is to map positive (negative) latitudes\n"
    "       into positive (negative) northings y; i.e., a branch cut is placed on\n"
    "       the equator.  With the extended domain, the northern sheet of the\n"
    "       projection is extended into the south hemisphere by pushing the branch\n"
    "       cut south from the branch points.  See the reference below for details.\n"
    "\n"
    "EXAMPLES\n"
    "          echo 0 90 | TransverseMercatorProj\n"
    "          => 25953592.84 9997964.94 90 18.40\n"
    "          echo 260e5 100e5 | TransverseMercatorProj -r\n"
    "          => -0.02 90.00 90.01 18.48\n"
    "\n"
    "ERRORS\n"
    "       An illegal line of input will print an error message to standard output\n"
    "       beginning with \"ERROR:\" and causes TransverseMercatorProj to return an\n"
    "       exit code of 1.  However, an error does not cause\n"
    "       TransverseMercatorProj to terminate; following lines will be converted.\n"
    "\n"
    "AUTHOR\n"
    "       TransverseMercatorProj was written by Charles Karney.\n"
    "\n"
    "SEE ALSO\n"
    "       The algorithms for the transverse Mercator projection are described in\n"
    "       C. F. F. Karney, Transverse Mercator with an accuracy of a few\n"
    "       nanometers, J. Geodesy 85(8), 475-485 (Aug. 2011); DOI\n"
    "       <https://doi.org/10.1007/s00190-011-0445-3>; preprint\n"
    "       <https://arxiv.org/abs/1002.1417>.  The explanation of the extended\n"
    "       domain of the projection with the -t option is given in Section 5 of\n"
    "       this paper.\n"
    "\n"
    "HISTORY\n"
    "       TransverseMercatorProj was added to GeographicLib,\n"
    "       <https://geographiclib.sourceforge.io>, in 2009-01.  Prior to version\n"
    "       1.9 it was called TransverseMercatorTest (and its interface was\n"
    "       slightly different).\n"
    ;
  return retval;
}
