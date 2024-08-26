/**
 * \file CartConvert.cpp
 * \brief Command line utility for geodetic to cartesian coordinate conversions
 *
 * Copyright (c) Charles Karney (2009-2017) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 *
 * See the <a href="CartConvert.1.html">man page</a> for usage information.
 **********************************************************************/

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <GeographicLib/Geocentric.hpp>
#include <GeographicLib/LocalCartesian.hpp>
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
    Utility::set_digits();
    bool localcartesian = false, reverse = false, longfirst = false;
    real
      a = Constants::WGS84_a(),
      f = Constants::WGS84_f();
    int prec = 6;
    real lat0 = 0, lon0 = 0, h0 = 0;
    std::string istring, ifile, ofile, cdelim;
    char lsep = ';';

    for (int m = 1; m < argc; ++m) {
      std::string arg(argv[m]);
      if (arg == "-r")
        reverse = true;
      else if (arg == "-l") {
        localcartesian = true;
        if (m + 3 >= argc) return usage(1, true);
        try {
          DMS::DecodeLatLon(std::string(argv[m + 1]), std::string(argv[m + 2]),
                            lat0, lon0, longfirst);
          h0 = Utility::val<real>(std::string(argv[m + 3]));
        }
        catch (const std::exception& e) {
          std::cerr << "Error decoding arguments of -l: " << e.what() << "\n";
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
      }  else if (arg == "--input-string") {
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

    const Geocentric ec(a, f);
    const LocalCartesian lc(lat0, lon0, h0, ec);

    // Max precision = 10: 0.1 nm in distance, 10^-15 deg (= 0.11 nm),
    // 10^-11 sec (= 0.3 nm).
    prec = std::min(10 + Math::extra_digits(), std::max(0, prec));
    std::string s, eol, stra, strb, strc, strd;
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
        // initial values to suppress warnings
        real lat, lon, h, x = 0, y = 0, z = 0;
        if (!(str >> stra >> strb >> strc))
          throw GeographicErr("Incomplete input: " + s);
        if (reverse) {
          x = Utility::val<real>(stra);
          y = Utility::val<real>(strb);
          z = Utility::val<real>(strc);
        } else {
          DMS::DecodeLatLon(stra, strb, lat, lon, longfirst);
          h = Utility::val<real>(strc);
        }
        if (str >> strd)
          throw GeographicErr("Extraneous input: " + strd);
        if (reverse) {
          if (localcartesian)
            lc.Reverse(x, y, z, lat, lon, h);
          else
            ec.Reverse(x, y, z, lat, lon, h);
          *output << Utility::str(longfirst ? lon : lat, prec + 5) << " "
                  << Utility::str(longfirst ? lat : lon, prec + 5) << " "
                  << Utility::str(h, prec) << eol;
        } else {
          if (localcartesian)
            lc.Forward(lat, lon, h, x, y, z);
          else
            ec.Forward(lat, lon, h, x, y, z);
          *output << Utility::str(x, prec) << " "
                  << Utility::str(y, prec) << " "
                  << Utility::str(z, prec) << eol;
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
    "    CartConvert [ -r ] [ -l lat0 lon0 h0 ] [ -e a f ] [ -w ] [ -p prec ] [\n"
    "    --comment-delimiter commentdelim ] [ --version | -h | --help ] [\n"
    "    --input-file infile | --input-string instring ] [ --line-separator\n"
    "    linesep ] [ --output-file outfile ]\n"
    "\n"
    "For full documentation type:\n"
    "    CartConvert --help\n"
    "or visit:\n"
    "    https://geographiclib.sourceforge.io/C++/2.4/CartConvert.1.html\n";
  else
    (retval ? std::cerr : std::cout) << "Man page:\n"
    "NAME\n"
    "       CartConvert -- convert geodetic coordinates to geocentric or local\n"
    "       cartesian\n"
    "\n"
    "SYNOPSIS\n"
    "       CartConvert [ -r ] [ -l lat0 lon0 h0 ] [ -e a f ] [ -w ] [ -p prec ] [\n"
    "       --comment-delimiter commentdelim ] [ --version | -h | --help ] [\n"
    "       --input-file infile | --input-string instring ] [ --line-separator\n"
    "       linesep ] [ --output-file outfile ]\n"
    "\n"
    "DESCRIPTION\n"
    "       Convert geodetic coordinates to either geocentric or local cartesian\n"
    "       coordinates.  Geocentric coordinates have the origin at the center of\n"
    "       the earth, with the z axis going thru the north pole, and the x axis\n"
    "       thru latitude = 0, longitude = 0.  By default, the conversion is to\n"
    "       geocentric coordinates.  Specifying -l lat0 lon0 h0 causes a local\n"
    "       coordinate system to be used with the origin at latitude = lat0,\n"
    "       longitude = lon0, height = h0, z normal to the ellipsoid and y due\n"
    "       north.\n"
    "\n"
    "       Geodetic coordinates are provided on standard input as a set of lines\n"
    "       containing (blank separated) latitude, longitude (decimal degrees or\n"
    "       degrees, minutes and seconds), and height above the ellipsoid (meters);\n"
    "       for details on the allowed formats for latitude and longitude, see the\n"
    "       \"GEOGRAPHIC COORDINATES\" section of GeoConvert(1).  For each set of\n"
    "       geodetic coordinates, the corresponding cartesian coordinates x, y, z\n"
    "       (meters) are printed on standard output.\n"
    "\n"
    "OPTIONS\n"
    "       -r  perform the reverse projection.  x, y, z are given on standard\n"
    "           input and each line of standard output gives latitude, longitude,\n"
    "           height.  In general there are multiple solutions and the result\n"
    "           which minimizes the absolute value of height is returned, i.e.,\n"
    "           (latitude, longitude) corresponds to the closest point on the\n"
    "           ellipsoid.\n"
    "\n"
    "       -l lat0 lon0 h0\n"
    "           specifies conversions to and from a local cartesion coordinate\n"
    "           systems with origin lat0 lon0 h0, instead of a geocentric\n"
    "           coordinate system.  The -w flag can be used to swap the default\n"
    "           order of the 2 geographic coordinates, provided that it appears\n"
    "           before -l.\n"
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
    "           of digits after the decimal point for geocentric and local\n"
    "           cartesion coordinates and for the height (in meters).  For\n"
    "           latitudes and longitudes (in degrees), the number of digits after\n"
    "           the decimal point is prec + 5.\n"
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
    "          echo 33.3 44.4 6000 | CartConvert\n"
    "          => 3816209.60 3737108.55 3485109.57\n"
    "          echo 33.3 44.4 6000 | CartConvert -l 33 44 20\n"
    "          => 37288.97 33374.29 5783.64\n"
    "          echo 30000 30000 0 | CartConvert -r\n"
    "          => 6.483 45 -6335709.73\n"
    "\n"
    "ERRORS\n"
    "       An illegal line of input will print an error message to standard output\n"
    "       beginning with \"ERROR:\" and causes CartConvert to return an exit code\n"
    "       of 1.  However, an error does not cause CartConvert to terminate;\n"
    "       following lines will be converted.\n"
    "\n"
    "SEE ALSO\n"
    "       The algorithm for converting geocentric to geodetic coordinates is\n"
    "       given in Appendix B of C. F. F. Karney, Geodesics on an ellipsoid of\n"
    "       revolution, Feb. 2011; preprint <https://arxiv.org/abs/1102.1215>.\n"
    "\n"
    "AUTHOR\n"
    "       CartConvert was written by Charles Karney.\n"
    "\n"
    "HISTORY\n"
    "       CartConvert was added to GeographicLib,\n"
    "       <https://geographiclib.sourceforge.io>, in 2009-02.  Prior to 2009-03\n"
    "       it was called ECEFConvert.\n"
    ;
  return retval;
}
