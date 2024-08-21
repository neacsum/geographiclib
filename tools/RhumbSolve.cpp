/**
 * \file RhumbSolve.cpp
 * \brief Command line utility for rhumb line calculations
 *
 * Copyright (c) Charles Karney (2014-2023) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 *
 * See the <a href="RhumbSolve.1.html">man page</a> for usage information.
 **********************************************************************/

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include <limits>
#include <GeographicLib/Rhumb.hpp>
#include <GeographicLib/DMS.hpp>
#include <GeographicLib/Utility.hpp>

#if defined(_MSC_VER)
// Squelch warnings about potentially uninitialized local variables
#  pragma warning (disable: 4701)
#endif

int usage (int retval, bool brief);

using namespace GeographicLib;
typedef Math::real real;

std::string LatLonString(real lat, real lon, int prec, bool dms, char dmssep,
                         bool longfirst) {
  using namespace GeographicLib;
  std::string
    latstr = dms ? DMS::Encode(lat, prec + 5, DMS::LATITUDE, dmssep) :
    DMS::Encode(lat, prec + 5, DMS::NUMBER),
    lonstr = dms ? DMS::Encode(lon, prec + 5, DMS::LONGITUDE, dmssep) :
    DMS::Encode(lon, prec + 5, DMS::NUMBER);
  return
    (longfirst ? lonstr : latstr) + " " + (longfirst ? latstr : lonstr);
}

std::string AzimuthString(real azi, int prec, bool dms, char dmssep) {
  return dms ? DMS::Encode(azi, prec + 5, DMS::AZIMUTH, dmssep) :
    DMS::Encode(azi, prec + 5, DMS::NUMBER);
}

int main(int argc, const char* const argv[]) {
  try {
    Utility::set_digits();
    bool linecalc = false, inverse = false, dms = false, exact = false,
      unroll = false, longfirst = false;
    real
      a = Constants::WGS84_a(),
      f = Constants::WGS84_f();
    real lat1, lon1, azi12 = Math::NaN(), lat2, lon2, s12, S12;
    int prec = 3;
    std::string istring, ifile, ofile, cdelim;
    char lsep = ';', dmssep = char(0);

    for (int m = 1; m < argc; ++m) {
      std::string arg(argv[m]);
      if (arg == "-i") {
        inverse = true;
        linecalc = false;
      } else if (arg == "-L") {
        inverse = false;
        linecalc = true;
        if (m + 3 >= argc) return usage(1, true);
        try {
          DMS::DecodeLatLon(std::string(argv[m + 1]), std::string(argv[m + 2]),
                            lat1, lon1, longfirst);
          azi12 = DMS::DecodeAzimuth(std::string(argv[m + 3]));
        }
        catch (const std::exception& e) {
          std::cerr << "Error decoding arguments of -L: " << e.what() << "\n";
          return 1;
        }
        m += 3;
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
      } else if (arg == "-u")
        unroll = true;
      else if (arg == "-d") {
        dms = true;
        dmssep = '\0';
      } else if (arg == "-:") {
        dms = true;
        dmssep = ':';
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
      } else if (arg == "-E")
        exact = true;
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

    const Rhumb rh(a, f, exact);
    const RhumbLine rhl(linecalc ? rh.Line(lat1, lon1, azi12) :
                        rh.Line(0, 0, Math::qd));
    // Max precision = 10: 0.1 nm in distance, 10^-15 deg (= 0.11 nm),
    // 10^-11 sec (= 0.3 nm).
    prec = std::min(10 + Math::extra_digits(), std::max(0, prec));
    std::string s, eol, slat1, slon1, slat2, slon2, sazi, ss12, strc;
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
        if (linecalc) {
          if (!(str >> ss12))
            throw GeographicErr("Incomplete input: " + s);
          if (str >> strc)
            throw GeographicErr("Extraneous input: " + strc);
          s12 = Utility::val<real>(ss12);
          rhl.GenPosition(s12, Rhumb::ALL | (unroll ? Rhumb::LONG_UNROLL : 0),
                          lat2, lon2, S12);
          *output << LatLonString(lat2, lon2, prec, dms, dmssep, longfirst)
                  << " " << Utility::str(S12, std::max(prec-7, 0)) << eol;
        } else if (inverse) {
          if (!(str >> slat1 >> slon1 >> slat2 >> slon2))
            throw GeographicErr("Incomplete input: " + s);
          if (str >> strc)
            throw GeographicErr("Extraneous input: " + strc);
          DMS::DecodeLatLon(slat1, slon1, lat1, lon1, longfirst);
          DMS::DecodeLatLon(slat2, slon2, lat2, lon2, longfirst);
          rh.Inverse(lat1, lon1, lat2, lon2, s12, azi12, S12);
          *output << AzimuthString(azi12, prec, dms, dmssep) << " "
                  << Utility::str(s12, prec) << " "
                  << Utility::str(S12, std::max(prec-7, 0)) << eol;
        } else {                // direct
          if (!(str >> slat1 >> slon1 >> sazi >> ss12))
            throw GeographicErr("Incomplete input: " + s);
          if (str >> strc)
            throw GeographicErr("Extraneous input: " + strc);
          DMS::DecodeLatLon(slat1, slon1, lat1, lon1, longfirst);
          azi12 = DMS::DecodeAzimuth(sazi);
          s12 = Utility::val<real>(ss12);
          rh.GenDirect(lat1, lon1, azi12, s12,
                       Rhumb::ALL | (unroll ? Rhumb::LONG_UNROLL : 0),
                       lat2, lon2, S12);
          *output << LatLonString(lat2, lon2, prec, dms, dmssep, longfirst)
                  << " " << Utility::str(S12, std::max(prec-7, 0)) << eol;
        }
      }
      catch (const std::exception& e) {
        // Write error message cout so output lines match input lines
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
    "    RhumbSolve [ -i | -L lat1 lon1 azi12 ] [ -e a f ] [ -u ] [ -d | -: ] [\n"
    "    -w ] [ -p prec ] [ -E ] [ --comment-delimiter commentdelim ] [\n"
    "    --version | -h | --help ] [ --input-file infile | --input-string\n"
    "    instring ] [ --line-separator linesep ] [ --output-file outfile ]\n"
    "\n"
    "For full documentation type:\n"
    "    RhumbSolve --help\n"
    "or visit:\n"
    "    https://geographiclib.sourceforge.io/C++/2.4/RhumbSolve.1.html\n";
  else
    (retval ? std::cerr : std::cout) << "Man page:\n"
    "NAME\n"
    "       RhumbSolve -- perform rhumb line calculations\n"
    "\n"
    "SYNOPSIS\n"
    "       RhumbSolve [ -i | -L lat1 lon1 azi12 ] [ -e a f ] [ -u ] [ -d | -: ] [\n"
    "       -w ] [ -p prec ] [ -E ] [ --comment-delimiter commentdelim ] [\n"
    "       --version | -h | --help ] [ --input-file infile | --input-string\n"
    "       instring ] [ --line-separator linesep ] [ --output-file outfile ]\n"
    "\n"
    "DESCRIPTION\n"
    "       The path with constant heading between two points on the ellipsoid at\n"
    "       (lat1, lon1) and (lat2, lon2) is called the rhumb line or loxodrome.\n"
    "       Its length is s12 and the rhumb line has a forward azimuth azi12 along\n"
    "       its length.  The quantity S12 is the area between the rhumb line from\n"
    "       point 1 to point 2 and the equator; i.e., it is the area, measured\n"
    "       counter-clockwise, of the geodesic quadrilateral with corners\n"
    "       (lat1,lon1), (0,lon1), (0,lon2), and (lat2,lon2).  The longitude\n"
    "       becomes indeterminate when a rhumb line passes through a pole, and\n"
    "       RhumbSolve reports NaNs for the longitude and the area in this case.\n"
    "\n"
    "       NOTE: the rhumb line is not the shortest path between two points; that\n"
    "       is the geodesic and it is calculated by GeodSolve(1).\n"
    "\n"
    "       RhumbSolve operates in one of three modes:\n"
    "\n"
    "       1.  By default, RhumbSolve accepts lines on the standard input\n"
    "           containing lat1 lon1 azi12 s12 and prints lat2 lon2 S12 on standard\n"
    "           output.  This is the direct calculation.\n"
    "\n"
    "       2.  With the -i option, RhumbSolve performs the inverse calculation.\n"
    "           It reads lines containing lat1 lon1 lat2 lon2 and prints the values\n"
    "           of azi12 s12 S12 for the corresponding shortest rhumb lines.\n"
    "\n"
    "       3.  Command line arguments -L lat1 lon1 azi12 specify a rhumb line.\n"
    "           RhumbSolve then accepts a sequence of s12 values (one per line) on\n"
    "           standard input and prints lat2 lon2 S12 for each.  This generates a\n"
    "           sequence of points on a rhumb line.\n"
    "\n"
    "OPTIONS\n"
    "       -i  perform an inverse calculation (see 2 above).\n"
    "\n"
    "       -L lat1 lon1 azi12\n"
    "           line mode (see 3 above); generate a sequence of points along the\n"
    "           rhumb line specified by lat1 lon1 azi12.  The -w flag can be used\n"
    "           to swap the default order of the 2 geographic coordinates, provided\n"
    "           that it appears before -L.\n"
    "\n"
    "       -e a f\n"
    "           specify the ellipsoid via the equatorial radius, a and the\n"
    "           flattening, f.  Setting f = 0 results in a sphere.  Specify f < 0\n"
    "           for a prolate ellipsoid.  A simple fraction, e.g., 1/297, is\n"
    "           allowed for f.  By default, the WGS84 ellipsoid is used, a =\n"
    "           6378137 m, f = 1/298.257223563.\n"
    "\n"
    "       -u  unroll the longitude.  Normally, on output longitudes are reduced\n"
    "           to lie in [-180deg,180deg).  However with this option, the returned\n"
    "           longitude lon2 is \"unrolled\" so that lon2 - lon1 indicates how\n"
    "           often and in what sense the geodesic has encircled the earth.\n"
    "\n"
    "       -d  output angles as degrees, minutes, seconds instead of decimal\n"
    "           degrees.\n"
    "\n"
    "       -:  like -d, except use : as a separator instead of the d, ', and \"\n"
    "           delimiters.\n"
    "\n"
    "       -w  on input and output, longitude precedes latitude (except that on\n"
    "           input this can be overridden by a hemisphere designator, N, S, E,\n"
    "           W).\n"
    "\n"
    "       -p prec\n"
    "           set the output precision to prec (default 3); prec is the precision\n"
    "           relative to 1 m.  See \"PRECISION\".\n"
    "\n"
    "       -E  By default, the rhumb line calculations are carried out using\n"
    "           series expansions valid for |f| < 0.01.  If -E is supplied, exact\n"
    "           equations for the rhumb line are used and the area integral is\n"
    "           computed with an accurate fit based on this exact equations; these\n"
    "           are valid for arbitrary eccentricities.\n"
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
    "INPUT\n"
    "       RhumbSolve measures all angles in degrees, all lengths (s12) in meters,\n"
    "       and all areas (S12) in meters^2.  On input angles (latitude, longitude,\n"
    "       azimuth, arc length) can be as decimal degrees or degrees, minutes,\n"
    "       seconds.  For example, \"40d30\", \"40d30'\", \"40:30\", \"40.5d\", and 40.5\n"
    "       are all equivalent.  By default, latitude precedes longitude for each\n"
    "       point (the -w flag switches this convention); however on input either\n"
    "       may be given first by appending (or prepending) N or S to the latitude\n"
    "       and E or W to the longitude.  Azimuths are measured clockwise from\n"
    "       north; however this may be overridden with E or W.\n"
    "\n"
    "       For details on the allowed formats for angles, see the \"GEOGRAPHIC\n"
    "       COORDINATES\" section of GeoConvert(1).\n"
    "\n"
    "PRECISION\n"
    "       prec gives precision of the output with prec = 0 giving 1 m precision,\n"
    "       prec = 3 giving 1 mm precision, etc.  prec is the number of digits\n"
    "       after the decimal point for lengths.  For decimal degrees, the number\n"
    "       of digits after the decimal point is prec + 5.  For DMS (degree,\n"
    "       minute, seconds) output, the number of digits after the decimal point\n"
    "       in the seconds component is prec + 1.  The minimum value of prec is 0\n"
    "       and the maximum is 10.\n"
    "\n"
    "ERRORS\n"
    "       An illegal line of input will print an error message to standard output\n"
    "       beginning with \"ERROR:\" and causes RhumbSolve to return an exit code of\n"
    "       1.  However, an error does not cause RhumbSolve to terminate; following\n"
    "       lines will be converted.\n"
    "\n"
    "ACCURACY\n"
    "       The algorithm used by RhumbSolve uses either series expansions or (if\n"
    "       -E is specified) exact formulas for computing the rhumb line and the\n"
    "       area.  These series are formulas are accurate for |f| < 0.01 and the\n"
    "       exact formulas apply for any value of the flattening.  The computation\n"
    "       of rhumb lines and the area involves the ratio of differences and, for\n"
    "       nearly east- or west-going rhumb lines, this might result in a large\n"
    "       loss of accuracy.  However, this problem is avoided by the use of\n"
    "       divided differences. For the WGS84 ellipsoid, the error is about 10\n"
    "       nanometers using either method.\n"
    "\n"
    "EXAMPLES\n"
    "       Route from JFK Airport to Singapore Changi Airport:\n"
    "\n"
    "          echo 40:38:23N 073:46:44W 01:21:33N 103:59:22E |\n"
    "          RhumbSolve -i -: -p 0\n"
    "\n"
    "          103:34:58.2 18523563 45921660958919\n"
    "\n"
    "       N.B. This is not the route typically taken by aircraft because it's\n"
    "       considerably longer than the geodesic given by GeodSolve(1).\n"
    "\n"
    "       Waypoints on the route at intervals of 2000km:\n"
    "\n"
    "          for ((i = 0; i <= 20; i += 2)); do echo ${i}000000;done |\n"
    "          RhumbSolve -L 40:38:23N 073:46:44W 103:34:58.2 -: -p 0\n"
    "\n"
    "          40:38:23.0N 073:46:44.0W 0\n"
    "          36:24:30.3N 051:28:26.4W 9817078307821\n"
    "          32:10:26.8N 030:20:57.3W 18224745682005\n"
    "          27:56:13.2N 010:10:54.2W 25358020327741\n"
    "          23:41:50.1N 009:12:45.5E 31321269267102\n"
    "          19:27:18.7N 027:59:22.1E 36195163180159\n"
    "          15:12:40.2N 046:17:01.1E 40041499143669\n"
    "          10:57:55.9N 064:12:52.8E 42906570007050\n"
    "          06:43:07.3N 081:53:28.8E 44823504180200\n"
    "          02:28:16.2N 099:24:54.5E 45813843358737\n"
    "          01:46:36.0S 116:52:59.7E 45888525219677\n"
    "\n"
    "SEE ALSO\n"
    "       GeoConvert(1), GeodSolve(1).\n"
    "\n"
    "       An online version of this utility is availbable at\n"
    "       <https://geographiclib.sourceforge.io/cgi-bin/RhumbSolve>.\n"
    "\n"
    "       An online version of this utility is availbable at\n"
    "       <https://geographiclib.sourceforge.io/cgi-bin/RhumbSolve>.\n"
    "\n"
    "       This solution for rhumb line is described in C. F. F. Karney, The area\n"
    "       of rhumb polygons, Technical Report, SRI International (2023); URL:\n"
    "       <https://arxiv.org/abs/2303.03219>.\n"
    "\n"
    "       The Wikipedia page, Rhumb line,\n"
    "       <https://en.wikipedia.org/wiki/Rhumb_line>.\n"
    "\n"
    "AUTHOR\n"
    "       RhumbSolve was written by Charles Karney.\n"
    "\n"
    "HISTORY\n"
    "       RhumbSolve was added to GeographicLib,\n"
    "       <https://geographiclib.sourceforge.io>, in version 1.37 and\n"
    "       substantially rewritten in version 2.2.\n"
    ;
  return retval;
}
