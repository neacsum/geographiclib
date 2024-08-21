/**
 * \file GeodesicProj.cpp
 * \brief Command line utility for geodesic projections
 *
 * Copyright (c) Charles Karney (2009-2017) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 *
 * See the <a href="GeodesicProj.1.html">man page</a> for usage information.
 **********************************************************************/

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <GeographicLib/Geodesic.hpp>
#include <GeographicLib/AzimuthalEquidistant.hpp>
#include <GeographicLib/CassiniSoldner.hpp>
#include <GeographicLib/Gnomonic.hpp>
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
    bool azimuthal = false, cassini = false, gnomonic = false, reverse = false,
      longfirst = false;
    real lat0 = 0, lon0 = 0;
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
      else if (arg == "-c" || arg == "-z" || arg == "-g") {
        cassini = azimuthal =  gnomonic = false;
        cassini = arg == "-c";
        azimuthal = arg == "-z";
        gnomonic = arg == "-g";
        if (m + 2 >= argc) return usage(1, true);
        try {
          DMS::DecodeLatLon(std::string(argv[m + 1]), std::string(argv[m + 2]),
                            lat0, lon0, longfirst);
        }
        catch (const std::exception& e) {
          std::cerr << "Error decoding arguments of " << arg << ": "
                    << e.what() << "\n";
          return 1;
        }
        m += 2;
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

    if (!(azimuthal || cassini || gnomonic)) {
      std::cerr << "Must specify \"-z lat0 lon0\" or "
                << "\"-c lat0 lon0\" or \"-g lat0 lon0\"\n";
      return 1;
    }

    const Geodesic geod(a, f);
    const CassiniSoldner cs = cassini ?
      CassiniSoldner(lat0, lon0, geod) : CassiniSoldner(geod);
    const AzimuthalEquidistant az(geod);
    const Gnomonic gn(geod);

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
        real lat, lon, x, y, azi, rk;
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
          if (cassini)
            cs.Reverse(x, y, lat, lon, azi, rk);
          else if (azimuthal)
            az.Reverse(lat0, lon0, x, y, lat, lon, azi, rk);
          else
            gn.Reverse(lat0, lon0, x, y, lat, lon, azi, rk);
          *output << Utility::str(longfirst ? lon : lat, prec + 5) << " "
                  << Utility::str(longfirst ? lat : lon, prec + 5) << " "
                  << Utility::str(azi, prec + 5) << " "
                  << Utility::str(rk, prec + 6) << eol;
        } else {
          if (cassini)
            cs.Forward(lat, lon, x, y, azi, rk);
          else if (azimuthal)
            az.Forward(lat0, lon0, lat, lon, x, y, azi, rk);
          else
            gn.Forward(lat0, lon0, lat, lon, x, y, azi, rk);
          *output << Utility::str(x, prec) << " "
                  << Utility::str(y, prec) << " "
                  << Utility::str(azi, prec + 5) << " "
                  << Utility::str(rk, prec + 6) << eol;
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
    "    GeodesicProj ( -z | -c | -g ) lat0 lon0 [ -r ] [ -e a f ] [ -w ] [ -p\n"
    "    prec ] [ --comment-delimiter commentdelim ] [ --version | -h | --help ]\n"
    "    [ --input-file infile | --input-string instring ] [ --line-separator\n"
    "    linesep ] [ --output-file outfile ]\n"
    "\n"
    "For full documentation type:\n"
    "    GeodesicProj --help\n"
    "or visit:\n"
    "    https://geographiclib.sourceforge.io/C++/2.4/GeodesicProj.1.html\n";
  else
    (retval ? std::cerr : std::cout) << "Man page:\n"
    "NAME\n"
    "       GeodesicProj -- perform projections based on geodesics\n"
    "\n"
    "SYNOPSIS\n"
    "       GeodesicProj ( -z | -c | -g ) lat0 lon0 [ -r ] [ -e a f ] [ -w ] [ -p\n"
    "       prec ] [ --comment-delimiter commentdelim ] [ --version | -h | --help ]\n"
    "       [ --input-file infile | --input-string instring ] [ --line-separator\n"
    "       linesep ] [ --output-file outfile ]\n"
    "\n"
    "DESCRIPTION\n"
    "       Perform projections based on geodesics.  Convert geodetic coordinates\n"
    "       to either azimuthal equidistant, Cassini-Soldner, or gnomonic\n"
    "       coordinates.  The center of the projection (lat0, lon0) is specified by\n"
    "       either the -c option (for Cassini-Soldner), the -z option (for\n"
    "       azimuthal equidistant), or the -g option (for gnomonic).  At least one\n"
    "       of these options must be given (the last one given is used).\n"
    "\n"
    "       Geodetic coordinates are provided on standard input as a set of lines\n"
    "       containing (blank separated) latitude and longitude (decimal degrees or\n"
    "       degrees, minutes, seconds); for details on the allowed formats for\n"
    "       latitude and longitude, see the \"GEOGRAPHIC COORDINATES\" section of\n"
    "       GeoConvert(1).  For each set of geodetic coordinates, the corresponding\n"
    "       projected coordinates x, y (meters) are printed on standard output\n"
    "       together with the azimuth azi (degrees) and reciprocal scale rk.  For\n"
    "       Cassini-Soldner, azi is the bearing of the easting direction and the\n"
    "       scale in the easting direction is 1 and the scale in the northing\n"
    "       direction is 1/rk.  For azimuthal equidistant and gnomonic, azi is the\n"
    "       bearing of the radial direction and the scale in the azimuthal\n"
    "       direction is 1/rk.  For azimuthal equidistant and gnomonic, the scales\n"
    "       in the radial direction are 1 and 1/rk^2, respectively.\n"
    "\n"
    "OPTIONS\n"
    "       -z lat0 lon0\n"
    "           use the azimuthal equidistant projection centered at latitude =\n"
    "           lat0, longitude = lon0.  The -w flag can be used to swap the\n"
    "           default order of the 2 coordinates, provided that it appears before\n"
    "           -z.\n"
    "\n"
    "       -c lat0 lon0\n"
    "           use the Cassini-Soldner projection centered at latitude = lat0,\n"
    "           longitude = lon0.  The -w flag can be used to swap the default\n"
    "           order of the 2 coordinates, provided that it appears before -c.\n"
    "\n"
    "       -g lat0 lon0\n"
    "           use the ellipsoidal gnomonic projection centered at latitude =\n"
    "           lat0, longitude = lon0.  The -w flag can be used to swap the\n"
    "           default order of the 2 coordinates, provided that it appears before\n"
    "           -g.\n"
    "\n"
    "       -r  perform the reverse projection.  x and y are given on standard\n"
    "           input and each line of standard output gives latitude, longitude,\n"
    "           azi, and rk.\n"
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
    "           latitudes, longitudes, and azimuths (in degrees), the number of\n"
    "           digits after the decimal point is prec + 5.  For the scale, the\n"
    "           number of digits after the decimal point is prec + 6.\n"
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
    "          echo 48.648 -2.007 | GeodesicProj -c 48.836 2.337\n"
    "          => -319919 -11791 86.7 0.999\n"
    "          echo -319919 -11791 | GeodesicProj -c 48.836 2.337 -r\n"
    "          => 48.648 -2.007 86.7 0.999\n"
    "\n"
    "ERRORS\n"
    "       An illegal line of input will print an error message to standard output\n"
    "       beginning with \"ERROR:\" and causes GeodesicProj to return an exit code\n"
    "       of 1.  However, an error does not cause GeodesicProj to terminate;\n"
    "       following lines will be converted.\n"
    "\n"
    "SEE ALSO\n"
    "       The ellipsoidal gnomonic projection is derived in Section 8 of C. F. F.\n"
    "       Karney, Algorithms for geodesics, J. Geodesy 87, 43-55 (2013); DOI\n"
    "       <https://doi.org/10.1007/s00190-012-0578-z>; addenda:\n"
    "       <https://geographiclib.sourceforge.io/geod-addenda.html>.\n"
    "\n"
    "AUTHOR\n"
    "       GeodesicProj was written by Charles Karney.\n"
    "\n"
    "HISTORY\n"
    "       GeodesicProj was added to GeographicLib,\n"
    "       <https://geographiclib.sourceforge.io>, in 2009-08.  Prior to version\n"
    "       1.9 it was called EquidistantTest.\n"
    ;
  return retval;
}
