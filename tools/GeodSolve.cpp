/**
 * \file GeodSolve.cpp
 * \brief Command line utility for geodesic calculations
 *
 * Copyright (c) Charles Karney (2009-2023) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 *
 * See the <a href="GeodSolve.1.html">man page</a> for usage information.
 **********************************************************************/

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <GeographicLib/Geodesic.hpp>
#include <GeographicLib/GeodesicLine.hpp>
#include <GeographicLib/DMS.hpp>
#include <GeographicLib/Utility.hpp>

#if defined(_MSC_VER)
// Squelch warnings about potentially uninitialized local variables
#  pragma warning (disable: 4701)
#endif

int usage (int retval, bool brief);

typedef GeographicLib::Math::real real;

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
  using namespace GeographicLib;
  return dms ? DMS::Encode(azi, prec + 5, DMS::AZIMUTH, dmssep) :
    DMS::Encode(azi, prec + 5, DMS::NUMBER);
}

std::string DistanceStrings(real s12, real a12,
                            bool full, bool arcmode, int prec, bool dms) {
  using namespace GeographicLib;
  std::string s;
  if (full || !arcmode)
    s += Utility::str(s12, prec);
  if (full)
    s += " ";
  if (full || arcmode)
    s += DMS::Encode(a12, prec + 5, dms ? DMS::NONE : DMS::NUMBER);
  return s;
}

real ReadDistance(const std::string& s, bool arcmode, bool fraction = false) {
  using namespace GeographicLib;
  return fraction ? Utility::fract<real>(s) :
    (arcmode ? DMS::DecodeAngle(s) : Utility::val<real>(s));
}

int main(int argc, const char* const argv[]) {
  try {
    using namespace GeographicLib;
    enum { NONE = 0, LINE, DIRECT, INVERSE };
    Utility::set_digits();
    bool inverse = false, arcmode = false,
      dms = false, full = false, exact = false, unroll = false,
      longfirst = false, azi2back = false, fraction = false,
      arcmodeline = false;
    real
      a = Constants::WGS84_a(),
      f = Constants::WGS84_f();
    real lat1, lon1, azi1, lat2, lon2, azi2, s12, m12, a12, M12, M21, S12,
      mult = 1;
    int linecalc = NONE, prec = 3;
    std::string istring, ifile, ofile, cdelim;
    char lsep = ';', dmssep = char(0);

    for (int m = 1; m < argc; ++m) {
      std::string arg(argv[m]);
      if (arg == "-i") {
        inverse = true;
        linecalc = NONE;
      } else if (arg == "-a")
        arcmode = !arcmode;
      else if (arg == "-F")
        fraction = true;
      else if (arg == "-L") {
        inverse = false;
        linecalc = LINE;
        if (m + 3 >= argc) return usage(1, true);
        try {
          DMS::DecodeLatLon(std::string(argv[m + 1]), std::string(argv[m + 2]),
                            lat1, lon1, longfirst);
          azi1 = DMS::DecodeAzimuth(std::string(argv[m + 3]));
        }
        catch (const std::exception& e) {
          std::cerr << "Error decoding arguments of -L: " << e.what() << "\n";
          return 1;
        }
        m += 3;
      } else if (arg == "-D") {
        inverse = false;
        linecalc = DIRECT;
        if (m + 4 >= argc) return usage(1, true);
        try {
          DMS::DecodeLatLon(std::string(argv[m + 1]), std::string(argv[m + 2]),
                            lat1, lon1, longfirst);
          azi1 = DMS::DecodeAzimuth(std::string(argv[m + 3]));
          s12 = ReadDistance(std::string(argv[m + 4]), arcmode);
          arcmodeline = arcmode;
        }
        catch (const std::exception& e) {
          std::cerr << "Error decoding arguments of -D: " << e.what() << "\n";
          return 1;
        }
        m += 4;
      } else if (arg == "-I") {
        inverse = false;
        linecalc = INVERSE;
        if (m + 4 >= argc) return usage(1, true);
        try {
          DMS::DecodeLatLon(std::string(argv[m + 1]), std::string(argv[m + 2]),
                            lat1, lon1, longfirst);
          DMS::DecodeLatLon(std::string(argv[m + 3]), std::string(argv[m + 4]),
                            lat2, lon2, longfirst);
        }
        catch (const std::exception& e) {
          std::cerr << "Error decoding arguments of -I: " << e.what() << "\n";
          return 1;
        }
        m += 4;
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
      else if (arg == "-b")
        azi2back = true;
      else if (arg == "-f")
        full = true;
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

    unsigned outmask = Geodesic::LATITUDE | Geodesic::LONGITUDE |
      Geodesic::AZIMUTH;        // basic output quantities
    outmask |= inverse ? Geodesic::DISTANCE : // distance-related flags
      (arcmode ? Geodesic::NONE : Geodesic::DISTANCE_IN);
    // longitude unrolling
    outmask |= unroll ? Geodesic::LONG_UNROLL : Geodesic::NONE;
    // full output -- don't use Geodesic::ALL since this includes DISTANCE_IN
    outmask |= full ? (Geodesic::DISTANCE | Geodesic::REDUCEDLENGTH |
                       Geodesic::GEODESICSCALE | Geodesic::AREA) :
      Geodesic::NONE;

    const Geodesic geods(a, f, exact);
    GeodesicLine ls;
    if (linecalc) {
      if (linecalc == LINE) fraction = false;
      ls = linecalc == DIRECT ?
        geods.GenDirectLine(lat1, lon1, azi1, arcmodeline, s12, outmask) :
        linecalc == INVERSE ?
        geods.InverseLine(lat1, lon1, lat2, lon2, outmask) :
        // linecalc == LINE
        geods.Line(lat1, lon1, azi1, outmask);
      mult = fraction ? ls.GenDistance(arcmode) : 1;
      if (linecalc == INVERSE) azi1 = ls.Azimuth();
    }

    // Max precision = 10: 0.1 nm in distance, 10^-15 deg (= 0.11 nm),
    // 10^-11 sec (= 0.3 nm).
    prec = std::min(10 + Math::extra_digits(), std::max(0, prec));
    std::string s, eol, slat1, slon1, slat2, slon2, sazi1, ss12, strc;
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
        if (inverse) {
          if (!(str >> slat1 >> slon1 >> slat2 >> slon2))
            throw GeographicErr("Incomplete input: " + s);
          if (str >> strc)
            throw GeographicErr("Extraneous input: " + strc);
          DMS::DecodeLatLon(slat1, slon1, lat1, lon1, longfirst);
          DMS::DecodeLatLon(slat2, slon2, lat2, lon2, longfirst);
          a12 = geods.GenInverse(lat1, lon1, lat2, lon2, outmask,
                                 s12, azi1, azi2, m12, M12, M21, S12);
          if (full) {
            if (unroll) {
              real e;
              lon2 = lon1 + Math::AngDiff(lon1, lon2, e);
              lon2 += e;
            } else {
              lon1 = Math::AngNormalize(lon1);
              lon2 = Math::AngNormalize(lon2);
            }
            *output << LatLonString(lat1, lon1, prec, dms, dmssep, longfirst)
                    << " ";
          }
          *output << AzimuthString(azi1, prec, dms, dmssep) << " ";
          if (full)
            *output << LatLonString(lat2, lon2, prec, dms, dmssep, longfirst)
                    << " ";
          if (azi2back) {
            using std::copysign;
            // map +/-0 -> -/+180; +/-180 -> -/+0
            // this depends on abs(azi2) <= 180
            azi2 = copysign(azi2 + copysign(real(Math::hd), -azi2), -azi2);
          }
          *output << AzimuthString(azi2, prec, dms, dmssep) << " "
                  << DistanceStrings(s12, a12, full, arcmode, prec, dms);
          if (full)
            *output << " " << Utility::str(m12, prec)
                    << " " << Utility::str(M12, prec+7)
                    << " " << Utility::str(M21, prec+7)
                    << " " << Utility::str(S12, std::max(prec-7, 0));
          *output << eol;
        } else {
          if (linecalc) {
            if (!(str >> ss12))
              throw GeographicErr("Incomplete input: " + s);
            if (str >> strc)
              throw GeographicErr("Extraneous input: " + strc);
            // In fraction mode input is read as a distance
            s12 = ReadDistance(ss12, !fraction && arcmode, fraction) * mult;
            a12 = ls.GenPosition(arcmode, s12, outmask,
                                 lat2, lon2, azi2, s12, m12, M12, M21, S12);
          } else {
            if (!(str >> slat1 >> slon1 >> sazi1 >> ss12))
              throw GeographicErr("Incomplete input: " + s);
            if (str >> strc)
              throw GeographicErr("Extraneous input: " + strc);
            DMS::DecodeLatLon(slat1, slon1, lat1, lon1, longfirst);
            azi1 = DMS::DecodeAzimuth(sazi1);
            s12 = ReadDistance(ss12, arcmode);
            a12 = geods.GenDirect(lat1, lon1, azi1, arcmode, s12, outmask,
                                  lat2, lon2, azi2, s12, m12, M12, M21, S12);
          }
          if (full)
            *output
              << LatLonString(lat1, unroll ? lon1 : Math::AngNormalize(lon1),
                              prec, dms, dmssep, longfirst)
              << " " << AzimuthString(azi1, prec, dms, dmssep) << " ";
          if (azi2back) {
            using std::copysign;
            // map +/-0 -> -/+180; +/-180 -> -/+0
            // this depends on abs(azi2) <= 180
            azi2 = copysign(azi2 + copysign(real(Math::hd), -azi2), -azi2);
          }
          *output << LatLonString(lat2, lon2, prec, dms, dmssep, longfirst)
                  << " " << AzimuthString(azi2, prec, dms, dmssep);
          if (full)
            *output << " "
                    << DistanceStrings(s12, a12, full, arcmode, prec, dms)
                    << " " << Utility::str(m12, prec)
                    << " " << Utility::str(M12, prec+7)
                    << " " << Utility::str(M21, prec+7)
                    << " " << Utility::str(S12, std::max(prec-7, 0));
          *output << eol;
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
    "    GeodSolve [ -i | -L lat1 lon1 azi1 | -D lat1 lon1 azi1 s13 | -I lat1\n"
    "    lon1 lat3 lon3 ] [ -a ] [ -e a f ] [ -u ] [ -F ] [ -d | -: ] [ -w ] [\n"
    "    -b ] [ -f ] [ -p prec ] [ -E ] [ --comment-delimiter commentdelim ] [\n"
    "    --version | -h | --help ] [ --input-file infile | --input-string\n"
    "    instring ] [ --line-separator linesep ] [ --output-file outfile ]\n"
    "\n"
    "For full documentation type:\n"
    "    GeodSolve --help\n"
    "or visit:\n"
    "    https://geographiclib.sourceforge.io/C++/2.4/GeodSolve.1.html\n";
  else
    (retval ? std::cerr : std::cout) << "Man page:\n"
    "NAME\n"
    "       GeodSolve -- perform geodesic calculations\n"
    "\n"
    "SYNOPSIS\n"
    "       GeodSolve [ -i | -L lat1 lon1 azi1 | -D lat1 lon1 azi1 s13 | -I lat1\n"
    "       lon1 lat3 lon3 ] [ -a ] [ -e a f ] [ -u ] [ -F ] [ -d | -: ] [ -w ] [\n"
    "       -b ] [ -f ] [ -p prec ] [ -E ] [ --comment-delimiter commentdelim ] [\n"
    "       --version | -h | --help ] [ --input-file infile | --input-string\n"
    "       instring ] [ --line-separator linesep ] [ --output-file outfile ]\n"
    "\n"
    "DESCRIPTION\n"
    "       The shortest path between two points on the ellipsoid at (lat1, lon1)\n"
    "       and (lat2, lon2) is called the geodesic.  Its length is s12 and the\n"
    "       geodesic from point 1 to point 2 has forward azimuths azi1 and azi2 at\n"
    "       the two end points.\n"
    "\n"
    "       GeodSolve operates in one of three modes:\n"
    "\n"
    "       1.  By default, GeodSolve accepts lines on the standard input\n"
    "           containing lat1 lon1 azi1 s12 and prints lat2 lon2 azi2 on standard\n"
    "           output.  This is the direct geodesic calculation.\n"
    "\n"
    "       2.  With the -i option, GeodSolve performs the inverse geodesic\n"
    "           calculation.  It reads lines containing lat1 lon1 lat2 lon2 and\n"
    "           prints the corresponding values of azi1 azi2 s12.\n"
    "\n"
    "       3.  Command line arguments -L lat1 lon1 azi1 specify a geodesic line.\n"
    "           GeodSolve then accepts a sequence of s12 values (one per line) on\n"
    "           standard input and prints lat2 lon2 azi2 for each.  This generates\n"
    "           a sequence of points on a single geodesic.  Command line arguments\n"
    "           -D and -I work similarly with the geodesic line defined in terms of\n"
    "           a direct or inverse geodesic calculation, respectively.\n"
    "\n"
    "OPTIONS\n"
    "       -i  perform an inverse geodesic calculation (see 2 above).\n"
    "\n"
    "       -L lat1 lon1 azi1\n"
    "           line mode (see 3 above); generate a sequence of points along the\n"
    "           geodesic specified by lat1 lon1 azi1.  The -w flag can be used to\n"
    "           swap the default order of the 2 geographic coordinates, provided\n"
    "           that it appears before -L.\n"
    "\n"
    "       -D lat1 lon1 azi1 s13\n"
    "           line mode (see 3 above); generate a sequence of points along the\n"
    "           geodesic specified by lat1 lon1 azi1 s13.  The -w flag can be used\n"
    "           to swap the default order of the 2 geographic coordinates, provided\n"
    "           that it appears before -D.  Similarly, the -a flag can be used to\n"
    "           change the interpretation of s13 to a13, provided that it appears\n"
    "           before -D.\n"
    "\n"
    "       -I lat1 lon1 lat3 lon3\n"
    "           line mode (see 3 above); generate a sequence of points along the\n"
    "           geodesic specified by lat1 lon1 lat3 lon3.  The -w flag can be used\n"
    "           to swap the default order of the 2 geographic coordinates, provided\n"
    "           that it appears before -I.\n"
    "\n"
    "       -a  toggle the arc mode flag (it starts off); if this flag is on, then\n"
    "           on input and output s12 is replaced by a12 the arc length (in\n"
    "           degrees) on the auxiliary sphere.  See \"AUXILIARY SPHERE\".\n"
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
    "           often and in what sense the geodesic has encircled the earth.  Use\n"
    "           the -f option, to get both longitudes printed.\n"
    "\n"
    "       -F  fractional mode.  This only has any effect with the -D and -I\n"
    "           options (and is otherwise ignored).  The values read on standard\n"
    "           input are interpreted as fractional distances to point 3, i.e., as\n"
    "           s12/s13 instead of s12.  If arc mode is in effect, then the values\n"
    "           denote fractional arc length, i.e., a12/a13.  The fractional\n"
    "           distances can be entered as a simple fraction, e.g., 3/4.\n"
    "\n"
    "       -d  output angles as degrees, minutes, seconds instead of decimal\n"
    "           degrees.\n"
    "\n"
    "       -:  like -d, except use : as a separator instead of the d, ', and \"\n"
    "           delimiters.\n"
    "\n"
    "       -w  toggle the longitude first flag (it starts off); if the flag is on,\n"
    "           then on input and output, longitude precedes latitude (except that,\n"
    "           on input, this can be overridden by a hemisphere designator, N, S,\n"
    "           E, W).\n"
    "\n"
    "       -b  report the back azimuth at point 2 instead of the forward azimuth.\n"
    "\n"
    "       -f  full output; each line of output consists of 12 quantities: lat1\n"
    "           lon1 azi1 lat2 lon2 azi2 s12 a12 m12 M12 M21 S12.  a12 is described\n"
    "           in \"AUXILIARY SPHERE\".  The four quantities m12, M12, M21, and S12\n"
    "           are described in \"ADDITIONAL QUANTITIES\".\n"
    "\n"
    "       -p prec\n"
    "           set the output precision to prec (default 3); prec is the precision\n"
    "           relative to 1 m.  See \"PRECISION\".\n"
    "\n"
    "       -E  use \"exact\" algorithms (based on elliptic integrals) for the\n"
    "           geodesic calculations.  These are more accurate than the (default)\n"
    "           series expansions for |f| > 0.02.\n"
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
    "       GeodSolve measures all angles in degrees and all lengths (s12) in\n"
    "       meters, and all areas (S12) in meters^2.  On input angles (latitude,\n"
    "       longitude, azimuth, arc length) can be as decimal degrees or degrees,\n"
    "       minutes, seconds.  For example, \"40d30\", \"40d30'\", \"40:30\", \"40.5d\",\n"
    "       and 40.5 are all equivalent.  By default, latitude precedes longitude\n"
    "       for each point (the -w flag switches this convention); however on input\n"
    "       either may be given first by appending (or prepending) N or S to the\n"
    "       latitude and E or W to the longitude.  Azimuths are measured clockwise\n"
    "       from north; however this may be overridden with E or W.\n"
    "\n"
    "       For details on the allowed formats for angles, see the \"GEOGRAPHIC\n"
    "       COORDINATES\" section of GeoConvert(1).\n"
    "\n"
    "AUXILIARY SPHERE\n"
    "       Geodesics on the ellipsoid can be transferred to the auxiliary sphere\n"
    "       on which the distance is measured in terms of the arc length a12\n"
    "       (measured in degrees) instead of s12.  In terms of a12, 180 degrees is\n"
    "       the distance from one equator crossing to the next or from the minimum\n"
    "       latitude to the maximum latitude.  Geodesics with a12 > 180 degrees do\n"
    "       not correspond to shortest paths.  With the -a flag, s12 (on both input\n"
    "       and output) is replaced by a12.  The -a flag does not affect the full\n"
    "       output given by the -f flag (which always includes both s12 and a12).\n"
    "\n"
    "ADDITIONAL QUANTITIES\n"
    "       The -f flag reports four additional quantities.\n"
    "\n"
    "       The reduced length of the geodesic, m12, is defined such that if the\n"
    "       initial azimuth is perturbed by dazi1 (radians) then the second point\n"
    "       is displaced by m12 dazi1 in the direction perpendicular to the\n"
    "       geodesic.  m12 is given in meters.  On a curved surface the reduced\n"
    "       length obeys a symmetry relation, m12 + m21 = 0.  On a flat surface, we\n"
    "       have m12 = s12.\n"
    "\n"
    "       M12 and M21 are geodesic scales.  If two geodesics are parallel at\n"
    "       point 1 and separated by a small distance dt, then they are separated\n"
    "       by a distance M12 dt at point 2.  M21 is defined similarly (with the\n"
    "       geodesics being parallel to one another at point 2).  M12 and M21 are\n"
    "       dimensionless quantities.  On a flat surface, we have M12 = M21 = 1.\n"
    "\n"
    "       If points 1, 2, and 3 lie on a single geodesic, then the following\n"
    "       addition rules hold:\n"
    "\n"
    "          s13 = s12 + s23,\n"
    "          a13 = a12 + a23,\n"
    "          S13 = S12 + S23,\n"
    "          m13 = m12 M23 + m23 M21,\n"
    "          M13 = M12 M23 - (1 - M12 M21) m23 / m12,\n"
    "          M31 = M32 M21 - (1 - M23 M32) m12 / m23.\n"
    "\n"
    "       Finally, S12 is the area between the geodesic from point 1 to point 2\n"
    "       and the equator; i.e., it is the area, measured counter-clockwise, of\n"
    "       the geodesic quadrilateral with corners (lat1,lon1), (0,lon1),\n"
    "       (0,lon2), and (lat2,lon2).  It is given in meters^2.\n"
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
    "       beginning with \"ERROR:\" and causes GeodSolve to return an exit code of\n"
    "       1.  However, an error does not cause GeodSolve to terminate; following\n"
    "       lines will be converted.\n"
    "\n"
    "ACCURACY\n"
    "       Using the (default) series solution, GeodSolve is accurate to about 15\n"
    "       nm (15 nanometers) for the WGS84 ellipsoid.  The approximate maximum\n"
    "       error (expressed as a distance) for an ellipsoid with the same\n"
    "       equatorial radius as the WGS84 ellipsoid and different values of the\n"
    "       flattening is\n"
    "\n"
    "          |f|     error\n"
    "          0.01    25 nm\n"
    "          0.02    30 nm\n"
    "          0.05    10 um\n"
    "          0.1    1.5 mm\n"
    "          0.2    300 mm\n"
    "\n"
    "       If -E is specified, GeodSolve is accurate to about 40 nm (40\n"
    "       nanometers) for the WGS84 ellipsoid.  The approximate maximum error\n"
    "       (expressed as a distance) for an ellipsoid with a quarter meridian of\n"
    "       10000 km and different values of the a/b = 1 - f is\n"
    "\n"
    "          1-f    error (nm)\n"
    "          1/128   387\n"
    "          1/64    345\n"
    "          1/32    269\n"
    "          1/16    210\n"
    "          1/8     115\n"
    "          1/4      69\n"
    "          1/2      36\n"
    "            1      15\n"
    "            2      25\n"
    "            4      96\n"
    "            8     318\n"
    "           16     985\n"
    "           32    2352\n"
    "           64    6008\n"
    "          128   19024\n"
    "\n"
    "MULTIPLE SOLUTIONS\n"
    "       The shortest distance returned for the inverse problem is (obviously)\n"
    "       uniquely defined.  However, in a few special cases there are multiple\n"
    "       azimuths which yield the same shortest distance.  Here is a catalog of\n"
    "       those cases:\n"
    "\n"
    "       lat1 = -lat2 (with neither point at a pole)\n"
    "           If azi1 = azi2, the geodesic is unique.  Otherwise there are two\n"
    "           geodesics and the second one is obtained by setting [azi1,azi2] =\n"
    "           [azi2,azi1], [M12,M21] = [M21,M12], S12 = -S12.  (This occurs when\n"
    "           the longitude difference is near +/-180 for oblate ellipsoids.)\n"
    "\n"
    "       lon2 = lon1 +/- 180 (with neither point at a pole)\n"
    "           If azi1 = 0 or +/-180, the geodesic is unique.  Otherwise there are\n"
    "           two geodesics and the second one is obtained by setting [azi1,azi2]\n"
    "           = [-azi1,-azi2], S12 = -S12.  (This occurs when lat2 is near -lat1\n"
    "           for prolate ellipsoids.)\n"
    "\n"
    "       Points 1 and 2 at opposite poles\n"
    "           There are infinitely many geodesics which can be generated by\n"
    "           setting [azi1,azi2] = [azi1,azi2] + [d,-d], for arbitrary d.  (For\n"
    "           spheres, this prescription applies when points 1 and 2 are\n"
    "           antipodal.)\n"
    "\n"
    "       s12 = 0 (coincident points)\n"
    "           There are infinitely many geodesics which can be generated by\n"
    "           setting [azi1,azi2] = [azi1,azi2] + [d,d], for arbitrary d.\n"
    "\n"
    "EXAMPLES\n"
    "       Route from JFK Airport to Singapore Changi Airport:\n"
    "\n"
    "          echo 40:38:23N 073:46:44W 01:21:33N 103:59:22E |\n"
    "          GeodSolve -i -: -p 0\n"
    "\n"
    "          003:18:29.9 177:29:09.2 15347628\n"
    "\n"
    "       Equally spaced waypoints on the route:\n"
    "\n"
    "          for ((i = 0; i <= 10; ++i)); do echo $i/10; done |\n"
    "          GeodSolve -I 40:38:23N 073:46:44W 01:21:33N 103:59:22E -F -: -p 0\n"
    "\n"
    "          40:38:23.0N 073:46:44.0W 003:18:29.9\n"
    "          54:24:51.3N 072:25:39.6W 004:18:44.1\n"
    "          68:07:37.7N 069:40:42.9W 006:44:25.4\n"
    "          81:38:00.4N 058:37:53.9W 017:28:52.7\n"
    "          83:43:26.0N 080:37:16.9E 156:26:00.4\n"
    "          70:20:29.2N 097:01:29.4E 172:31:56.4\n"
    "          56:38:36.0N 100:14:47.6E 175:26:10.5\n"
    "          42:52:37.1N 101:43:37.2E 176:34:28.6\n"
    "          29:03:57.0N 102:39:34.8E 177:07:35.2\n"
    "          15:13:18.6N 103:22:08.0E 177:23:44.7\n"
    "          01:21:33.0N 103:59:22.0E 177:29:09.2\n"
    "\n"
    "SEE ALSO\n"
    "       GeoConvert(1).\n"
    "\n"
    "       An online version of this utility is availbable at\n"
    "       <https://geographiclib.sourceforge.io/cgi-bin/GeodSolve>.\n"
    "\n"
    "       The algorithms are described in C. F. F. Karney, Algorithms for\n"
    "       geodesics, J. Geodesy 87, 43-55 (2013); DOI:\n"
    "       <https://doi.org/10.1007/s00190-012-0578-z>; addenda:\n"
    "       <https://geographiclib.sourceforge.io/geod-addenda.html>.\n"
    "\n"
    "       The Wikipedia page, Geodesics on an ellipsoid,\n"
    "       <https://en.wikipedia.org/wiki/Geodesics_on_an_ellipsoid>.\n"
    "\n"
    "AUTHOR\n"
    "       GeodSolve was written by Charles Karney.\n"
    "\n"
    "HISTORY\n"
    "       GeodSolve was added to GeographicLib,\n"
    "       <https://geographiclib.sourceforge.io>, in 2009-03.  Prior to version\n"
    "       1.30, it was called Geod.  (The name was changed to avoid a conflict\n"
    "       with the geod utility in proj.4.)\n"
    ;
  return retval;
}
