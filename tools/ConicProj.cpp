/**
 * \file ConicProj.cpp
 * \brief Command line utility for conical projections
 *
 * Copyright (c) Charles Karney (2009-2017) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 *
 * See the <a href="ConicProj.1.html">man page</a> for usage information.
 **********************************************************************/

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <GeographicLib/LambertConformalConic.hpp>
#include <GeographicLib/AlbersEqualArea.hpp>
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
    bool lcc = false, albers = false, reverse = false, longfirst = false;
    real lat1 = 0, lat2 = 0, lon0 = 0, k1 = 1;
    real
      a = Constants::WGS84_a(),
      f = Constants::WGS84_f();
    int prec = 6;
    std::string istring, ifile, ofile, cdelim;
    char lsep = ';';

    for (int m = 1; m < argc; ++m) {
      std::string arg(argv[m]);
      if (arg == "-r")
        reverse = true;
      else if (arg == "-c" || arg == "-a") {
        lcc = arg == "-c";
        albers = arg == "-a";
        if (m + 2 >= argc) return usage(1, true);
        try {
          for (int i = 0; i < 2; ++i) {
            DMS::flag ind;
            (i ? lat2 : lat1) = DMS::Decode(std::string(argv[++m]), ind);
            if (ind == DMS::LONGITUDE)
              throw GeographicErr("Bad hemisphere");
          }
        }
        catch (const std::exception& e) {
          std::cerr << "Error decoding arguments of " << arg << ": "
                    << e.what() << "\n";
          return 1;
        }
      } else if (arg == "-l") {
        if (++m == argc) return usage(1, true);
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
        if (++m == argc) return usage(1, true);
        try {
          k1 = Utility::val<real>(std::string(argv[m]));
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

    if (!(lcc || albers)) {
      std::cerr << "Must specify \"-c lat1 lat2\" or "
                << "\"-a lat1 lat2\"\n";
      return 1;
    }

    const LambertConformalConic lproj =
      lcc ? LambertConformalConic(a, f, lat1, lat2, k1)
      : LambertConformalConic(1, 0, 0, 0, 1);
    const AlbersEqualArea aproj =
      albers ? AlbersEqualArea(a, f, lat1, lat2, k1)
      : AlbersEqualArea(1, 0, 0, 0, 1);

    // Max precision = 10: 0.1 nm in distance, 10^-15 deg (= 0.11 nm),
    // 10^-11 sec (= 0.3 nm).
    prec = std::min(10 + Math::extra_digits(), std::max(0, prec));
    std::string s, eol, stra, strb, strc;
    std::istringstream str;
    int retval = 0;
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
        real lat, lon, x, y, gamma, k;
        if (!(str >> stra >> strb))
          throw GeographicErr("Incomplete input: " + s);
        if (reverse) {
          x = Utility::val<real>(stra);
          y = Utility::val<real>(strb);
        } else
          DMS::DecodeLatLon(stra, strb, lat, lon, longfirst);
        if (str >> strc)
          throw GeographicErr("Extraneous input: " + strc);
        if (reverse) {
          if (lcc)
            lproj.Reverse(lon0, x, y, lat, lon, gamma, k);
          else
            aproj.Reverse(lon0, x, y, lat, lon, gamma, k);
          *output << Utility::str(longfirst ? lon : lat, prec + 5) << " "
                  << Utility::str(longfirst ? lat : lon, prec + 5) << " "
                  << Utility::str(gamma, prec + 6) << " "
                  << Utility::str(k, prec + 6) << eol;
        } else {
          if (lcc)
            lproj.Forward(lon0, lat, lon, x, y, gamma, k);
          else
            aproj.Forward(lon0, lat, lon, x, y, gamma, k);
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
    "    ConicProj ( -c | -a ) lat1 lat2 [ -l lon0 ] [ -k k1 ] [ -r ] [ -e a f ]\n"
    "    [ -w ] [ -p prec ] [ --comment-delimiter commentdelim ] [ --version |\n"
    "    -h | --help ] [ --input-file infile | --input-string instring ] [\n"
    "    --line-separator linesep ] [ --output-file outfile ]\n"
    "\n"
    "For full documentation type:\n"
    "    ConicProj --help\n"
    "or visit:\n"
    "    https://geographiclib.sourceforge.io/C++/2.4/ConicProj.1.html\n";
  else
    (retval ? std::cerr : std::cout) << "Man page:\n"
    "NAME\n"
    "       ConicProj -- perform conic projections\n"
    "\n"
    "SYNOPSIS\n"
    "       ConicProj ( -c | -a ) lat1 lat2 [ -l lon0 ] [ -k k1 ] [ -r ] [ -e a f ]\n"
    "       [ -w ] [ -p prec ] [ --comment-delimiter commentdelim ] [ --version |\n"
    "       -h | --help ] [ --input-file infile | --input-string instring ] [\n"
    "       --line-separator linesep ] [ --output-file outfile ]\n"
    "\n"
    "DESCRIPTION\n"
    "       Perform one of two conic projections geodesics.  Convert geodetic\n"
    "       coordinates to either Lambert conformal conic or Albers equal area\n"
    "       coordinates.  The standard latitudes lat1 and lat2 are specified by\n"
    "       that the -c option (for Lambert conformal conic) or the -a option (for\n"
    "       Albers equal area).  At least one of these options must be given (the\n"
    "       last one given is used).  Specify lat1 = lat2, to obtain the case with\n"
    "       a single standard parallel.  The central meridian is given by lon0.\n"
    "       The longitude of origin is given by the latitude of minimum (azimuthal)\n"
    "       scale for Lambert conformal conic (Albers equal area).  The (azimuthal)\n"
    "       scale on the standard parallels is k1.\n"
    "\n"
    "       Geodetic coordinates are provided on standard input as a set of lines\n"
    "       containing (blank separated) latitude and longitude (decimal degrees or\n"
    "       degrees, minutes, seconds);  for details on the allowed formats for\n"
    "       latitude and longitude, see the \"GEOGRAPHIC COORDINATES\" section of\n"
    "       GeoConvert(1).  For each set of geodetic coordinates, the corresponding\n"
    "       projected easting, x, and northing, y, (meters) are printed on standard\n"
    "       output together with the meridian convergence gamma (degrees) and\n"
    "       (azimuthal) scale k.  For Albers equal area, the radial scale is 1/k.\n"
    "       The meridian convergence is the bearing of the y axis measured\n"
    "       clockwise from true north.\n"
    "\n"
    "       Special cases of the Lambert conformal projection are the Mercator\n"
    "       projection (the standard latitudes equal and opposite) and the polar\n"
    "       stereographic projection (both standard latitudes correspond to the\n"
    "       same pole).  Special cases of the Albers equal area projection are the\n"
    "       cylindrical equal area projection (the standard latitudes equal and\n"
    "       opposite), the Lambert azimuthal equal area projection (both standard\n"
    "       latitude corresponds to the same pole), and the Lambert equal area\n"
    "       conic projection (one standard parallel is at a pole).\n"
    "\n"
    "OPTIONS\n"
    "       -c lat1 lat2\n"
    "           use the Lambert conformal conic projection with standard parallels\n"
    "           lat1 and lat2.\n"
    "\n"
    "       -a lat1 lat2\n"
    "           use the Albers equal area projection with standard parallels lat1\n"
    "           and lat2.\n"
    "\n"
    "       -l lon0\n"
    "           specify the longitude of origin lon0 (degrees, default 0).\n"
    "\n"
    "       -k k1\n"
    "           specify the (azimuthal) scale k1 on the standard parallels (default\n"
    "           1).\n"
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
    "           6378137 m, f = 1/298.257223563.\n"
    "\n"
    "       -w  toggle the longitude first flag (it starts off); if the flag is on,\n"
    "           then on input and output, longitude precedes latitude (except that,\n"
    "           on input, this can be overridden by a hemisphere designator, N, S,\n"
    "           E, W).\n"
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
    "EXAMPLES\n"
    "          echo 39.95N 75.17W | ConicProj -c 40d58 39d56 -l 77d45W\n"
    "          => 220445 -52372 1.67 1.0\n"
    "          echo 220445 -52372 | ConicProj -c 40d58 39d56 -l 77d45W -r\n"
    "          => 39.95 -75.17 1.67 1.0\n"
    "\n"
    "ERRORS\n"
    "       An illegal line of input will print an error message to standard output\n"
    "       beginning with \"ERROR:\" and causes ConicProj to return an exit code of\n"
    "       1.  However, an error does not cause ConicProj to terminate; following\n"
    "       lines will be converted.\n"
    "\n"
    "AUTHOR\n"
    "       ConicProj was written by Charles Karney.\n"
    "\n"
    "HISTORY\n"
    "       ConicProj was added to GeographicLib,\n"
    "       <https://geographiclib.sourceforge.io>, in version 1.9.\n"
    ;
  return retval;
}
