/**
 * \file Planimeter.cpp
 * \brief Command line utility for measuring the area of geodesic polygons
 *
 * Copyright (c) Charles Karney (2010-2023) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 *
 * See the <a href="Planimeter.1.html">man page</a> for usage information.
 **********************************************************************/

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <GeographicLib/PolygonArea.hpp>
#include <GeographicLib/DMS.hpp>
#include <GeographicLib/Utility.hpp>
#include <GeographicLib/GeoCoords.hpp>
#include <GeographicLib/AuxLatitude.hpp>

int usage (int retval, bool brief);

int main(int argc, const char* const argv[]) {
  try {
    using namespace GeographicLib;
    typedef Math::real real;
    Utility::set_digits();
    enum { GEODESIC, AUTHALIC, RHUMB };
    real
      a = Constants::WGS84_a(),
      f = Constants::WGS84_f();
    bool reverse = false, sign = true, polyline = false, longfirst = false,
      exact = false, geoconvert_compat = false;
    int linetype = GEODESIC;
    int prec = 6;
    std::string istring, ifile, ofile, cdelim;
    char lsep = ';';

    for (int m = 1; m < argc; ++m) {
      std::string arg(argv[m]);
      if (arg == "-r")
        reverse = !reverse;
      else if (arg == "-s")
        sign = !sign;
      else if (arg == "-l")
        polyline = !polyline;
      else if (arg == "-e") {
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
      } else if (arg == "-G")
        linetype = GEODESIC;
      else if (arg == "-Q")
        linetype = AUTHALIC;
      else if (arg == "-R")
        linetype = RHUMB;
      else if (arg == "-E")
        exact = true;
      else if (arg == "--geoconvert-input")
        geoconvert_compat = true;
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

    // Linetype is one of GEODESIC, AUTHALIC, RHUMB
    const AuxLatitude ellip(a, f);
    if (linetype == AUTHALIC) {
      // adjusting a and f to correspond to the authalic sphere
      using std::sqrt;
      a = sqrt(ellip.AuthalicRadiusSquared(exact));
      f = 0;
    }
    const Geodesic geod(a, linetype != RHUMB ? f : 0,
                        exact && linetype != RHUMB);
    const Rhumb rhumb(a, linetype == RHUMB ? f : 0,
                      exact && linetype == RHUMB);
    PolygonArea poly(geod, polyline);
    PolygonAreaRhumb polyr(rhumb, polyline);
    GeoCoords p;

    // Max precision = 10: 0.1 nm in distance, 10^-15 deg (= 0.11 nm),
    // 10^-11 sec (= 0.3 nm).
    prec = std::min(10 + Math::extra_digits(), std::max(0, prec));
    std::string s, eol("\n");
    real perimeter, area;
    unsigned num;
    std::istringstream str;
    std::string slat, slon, junk;
    real lat = 0, lon = 0;
    while (std::getline(*input, s)) {
      if (!cdelim.empty()) {
        std::string::size_type m = s.find(cdelim);
        if (m != std::string::npos) {
          eol = " " + s.substr(m) + "\n";
          s = s.substr(0, m);
        }
      }
      bool endpoly = s.empty();
      if (!endpoly) {
        try {
          using std::isnan;
          if (geoconvert_compat) {
            p.Reset(s, true, longfirst);
            lat = p.Latitude(); lon = p.Longitude();
          } else {
            str.clear(); str.str(s);
            if (!(str >> slat >> slon))
              throw GeographicErr("incomplete input");
            if (str >> junk)
              throw GeographicErr("extra input");
            DMS::DecodeLatLon(slat, slon, lat, lon, longfirst);
          }
          if (isnan(lat) || isnan(lon))
            endpoly = true;
        }
        catch (const GeographicErr&) {
          endpoly = true;
        }
      }
      if (endpoly) {
        num =
          linetype == RHUMB ? polyr.Compute(reverse, sign, perimeter, area) :
          poly.Compute(reverse, sign, perimeter, area); // geodesic + authalic
        if (num > 0) {
          *output << num << " " << Utility::str(perimeter, prec);
          if (!polyline) {
            *output << " " << Utility::str(area, std::max(0, prec - 5));
          }
          *output << eol;
        }
        linetype == RHUMB ? polyr.Clear() : poly.Clear();
        eol = "\n";
      } else {
        linetype == RHUMB ? polyr.AddPoint(lat, lon) :
          poly.AddPoint
          (linetype == AUTHALIC ?
           ellip.Convert(AuxLatitude::PHI, AuxLatitude::XI, lat, exact) : lat,
           lon);
      }
    }
    num =
      linetype == RHUMB ? polyr.Compute(reverse, sign, perimeter, area) :
      poly.Compute(reverse, sign, perimeter, area);
    if (num > 0) {
      *output << num << " " << Utility::str(perimeter, prec);
      if (!polyline) {
        *output << " " << Utility::str(area, std::max(0, prec - 5));
      }
      *output << eol;
    }
      linetype == RHUMB ? polyr.Clear() : poly.Clear();
    eol = "\n";
    return 0;
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
    "    Planimeter [ -r ] [ -s ] [ -l ] [ -e a f ] [ -w ] [ -p prec ] [ -G | -Q\n"
    "    | -R ] [ -E ] [ --geoconvert-input ] [ --comment-delimiter commentdelim\n"
    "    ] [ --version | -h | --help ] [ --input-file infile | --input-string\n"
    "    instring ] [ --line-separator linesep ] [ --output-file outfile ]\n"
    "\n"
    "For full documentation type:\n"
    "    Planimeter --help\n"
    "or visit:\n"
    "    https://geographiclib.sourceforge.io/C++/2.4/Planimeter.1.html\n";
  else
    (retval ? std::cerr : std::cout) << "Man page:\n"
    "NAME\n"
    "       Planimeter -- compute the area of geodesic polygons\n"
    "\n"
    "SYNOPSIS\n"
    "       Planimeter [ -r ] [ -s ] [ -l ] [ -e a f ] [ -w ] [ -p prec ] [ -G | -Q\n"
    "       | -R ] [ -E ] [ --geoconvert-input ] [ --comment-delimiter commentdelim\n"
    "       ] [ --version | -h | --help ] [ --input-file infile | --input-string\n"
    "       instring ] [ --line-separator linesep ] [ --output-file outfile ]\n"
    "\n"
    "DESCRIPTION\n"
    "       Measure the area of a geodesic polygon.  Reads polygon vertices from\n"
    "       standard input, one per line.  Vertices are be given as latitude and\n"
    "       longitude.  By default latitude precedes longitude, however this\n"
    "       convention is reversed with the -w flag and a hemisphere designator (N,\n"
    "       S, E, W) can be used to disambiguate the coordinates.  The end of\n"
    "       input, a blank line, or a line which can't be interpreted as a vertex\n"
    "       signals the end of one polygon and the start of the next.  For each\n"
    "       polygon print a summary line with the number of points, the perimeter\n"
    "       (in meters), and the area (in meters^2).\n"
    "\n"
    "       The edges of the polygon are given by the shortest geodesic (or rhumb\n"
    "       line) between consecutive vertices.  In certain cases, there may be two\n"
    "       or many such shortest path, and in that case, the polygon is not\n"
    "       uniquely specified by its vertices.  For geodesics, this only happens\n"
    "       with very long edges (for the WGS84 ellipsoid, any edge shorter than\n"
    "       19970 km is uniquely specified by its end points).  In such cases,\n"
    "       insert an additional vertex near the middle of the long edge to define\n"
    "       the boundary of the polygon.\n"
    "\n"
    "       By default, polygons traversed in a counter-clockwise direction return\n"
    "       a positive area and those traversed in a clockwise direction return a\n"
    "       negative area.  This sign convention is reversed if the -r option is\n"
    "       given.\n"
    "\n"
    "       Of course, encircling an area in the clockwise direction is equivalent\n"
    "       to encircling the rest of the ellipsoid in the counter-clockwise\n"
    "       direction.  The default interpretation used by Planimeter is the one\n"
    "       that results in a smaller magnitude of area; i.e., the magnitude of the\n"
    "       area is less than or equal to one half the total area of the ellipsoid.\n"
    "       If the -s option is given, then the interpretation used is the one that\n"
    "       results in a positive area; i.e., the area is positive and less than\n"
    "       the total area of the ellipsoid.\n"
    "\n"
    "       Arbitrarily complex polygons are allowed.  In the case of self-\n"
    "       intersecting polygons the area is accumulated \"algebraically\", e.g.,\n"
    "       the areas of the 2 loops in a figure-8 polygon will partially cancel.\n"
    "       Polygons may include one or both poles.  There is no need to close the\n"
    "       polygon.\n"
    "\n"
    "OPTIONS\n"
    "       -r  toggle whether counter-clockwise traversal of the polygon returns a\n"
    "           positive (the default) or negative result.\n"
    "\n"
    "       -s  toggle whether to return a signed result (the default) or not.\n"
    "\n"
    "       -l  toggle whether the vertices represent a polygon (the default) or a\n"
    "           polyline.  For a polyline, the number of points and the length of\n"
    "           the path joining them is returned; the path is not closed and the\n"
    "           area is not reported.\n"
    "\n"
    "       -e a f\n"
    "           specify the ellipsoid via the equatorial radius, a and the\n"
    "           flattening, f.  Setting f = 0 results in a sphere.  Specify f < 0\n"
    "           for a prolate ellipsoid.  A simple fraction, e.g., 1/297, is\n"
    "           allowed for f.  By default, the WGS84 ellipsoid is used, a =\n"
    "           6378137 m, f = 1/298.257223563.  If entering vertices as UTM/UPS or\n"
    "           MGRS coordinates, use the default ellipsoid, since the conversion\n"
    "           of these coordinates to latitude and longitude always uses the\n"
    "           WGS84 parameters.\n"
    "\n"
    "       -w  toggle the longitude first flag (it starts off); if the flag is on,\n"
    "           then when reading geographic coordinates, longitude precedes\n"
    "           latitude (this can be overridden by a hemisphere designator, N, S,\n"
    "           E, W).\n"
    "\n"
    "       -p prec\n"
    "           set the output precision to prec (default 6); the perimeter is\n"
    "           given (in meters) with prec digits after the decimal point; the\n"
    "           area is given (in meters^2) with (prec - 5) digits after the\n"
    "           decimal point.\n"
    "\n"
    "       -G  the edges joining the vertices are geodesics.  This is the default\n"
    "           option and is recommended for terrestrial applications.  This\n"
    "           option, -G, and the following two options, -Q and -R, are mutually\n"
    "           exclusive.\n"
    "\n"
    "       -Q  map the points to the authalic sphere and compute the area of the\n"
    "           resulting spherical polygon.  The area will be reasonable accurate\n"
    "           provided that the edges are sufficiently short.  The perimeter\n"
    "           calculation is not accurate.\n"
    "\n"
    "       -R  the edges joining the vertices are rhumb lines instead of\n"
    "           geodesics.\n"
    "\n"
    "       -E  use the exact equations for the geodesic -G, authalic -Q, and rhumb\n"
    "           -R calculations instead of series expansions.  For the geodesic and\n"
    "           rhumb methods, the area is computed by applying discrete sine\n"
    "           transforms to the integrand in the expression for the area.  These\n"
    "           are more accurate, albeit slower, than the (default) series\n"
    "           expansions for |f| > 0.02 and will give high accuracy for -99 < f <\n"
    "           0.99.  It is not necessary to specify this option for terrestrial\n"
    "           applications.\n"
    "\n"
    "       --geoconvert-input\n"
    "           The input lines are interpreted in the same way as GeoConvert(1)\n"
    "           allowing the coordinates for the vertices to be given as UTM/UPS or\n"
    "           MGRS coordinates, as well as latitude and longitude.  CAUTION:\n"
    "           GeoConvert assumes the coordinates refer to the WGS84 ellipsoid\n"
    "           (disregarding the -e flag) and MGRS coordinates signify the center\n"
    "           of the corresponding MGRS square.\n"
    "\n"
    "       --comment-delimiter commentdelim\n"
    "           set the comment delimiter to commentdelim (e.g., \"#\" or \"//\").  If\n"
    "           set, the input lines will be scanned for this delimiter and, if\n"
    "           found, the delimiter and the rest of the line will be removed prior\n"
    "           to processing.  For a given polygon, the last such string found\n"
    "           will be appended to the output line (separated by a space).\n"
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
    "       Example (the area of the 100km MGRS square 18SWK)\n"
    "\n"
    "          Planimeter --geoconvert-input <<EOF\n"
    "          18n 500000 4400000\n"
    "          18n 600000 4400000\n"
    "          18n 600000 4500000\n"
    "          18n 500000 4500000\n"
    "          EOF\n"
    "          => 4 400139.532959 10007388597.2\n"
    "\n"
    "       The following code takes the output from gdalinfo and reports the area\n"
    "       covered by the data (assuming the edges of the image are geodesics).\n"
    "\n"
    "          #! /bin/sh\n"
    "          grep -E '^((Upper|Lower) (Left|Right)|Center) ' |\n"
    "          sed -e 's/d /d/g' -e \"s/' /'/g\" | tr -s '(),\\r\\t' ' ' | awk '{\n"
    "              if ($1 $2 == \"UpperLeft\")\n"
    "                  ul = $6 \" \" $5;\n"
    "              else if ($1 $2 == \"LowerLeft\")\n"
    "                  ll = $6 \" \" $5;\n"
    "              else if ($1 $2 == \"UpperRight\")\n"
    "                  ur = $6 \" \" $5;\n"
    "              else if ($1 $2 == \"LowerRight\")\n"
    "                  lr = $6 \" \" $5;\n"
    "              else if ($1 == \"Center\") {\n"
    "                  printf \"%s\\n%s\\n%s\\n%s\\n\\n\", ul, ll, lr, ur;\n"
    "                  ul = ll = ur = lr = \"\";\n"
    "              }\n"
    "          }\n"
    "          ' | Planimeter | cut -f3 -d' '\n"
    "\n"
    "ACCURACY\n"
    "       Using the -G option (the default), the accuracy was estimated by\n"
    "       computing the error in the area for 10^7 approximately regular polygons\n"
    "       on the WGS84 ellipsoid.  The centers and the orientations of the\n"
    "       polygons were uniformly distributed, the number of vertices was log-\n"
    "       uniformly distributed in [3, 300], and the center to vertex distance\n"
    "       log-uniformly distributed in [0.1 m, 9000 km].\n"
    "\n"
    "       The maximum error in the perimeter was 200 nm, and the maximum error in\n"
    "       the area was\n"
    "\n"
    "          0.0013 m^2 for perimeter < 10 km\n"
    "          0.0070 m^2 for perimeter < 100 km\n"
    "          0.070 m^2 for perimeter < 1000 km\n"
    "          0.11 m^2 for all perimeters\n"
    "\n"
    "SEE ALSO\n"
    "       GeoConvert(1), GeodSolve(1).\n"
    "\n"
    "       An online version of this utility is availbable at\n"
    "       <https://geographiclib.sourceforge.io/cgi-bin/Planimeter>.\n"
    "\n"
    "       The algorithm for the area of geodesic polygon is given in Section 6 of\n"
    "       C. F. F. Karney, Algorithms for geodesics, J. Geodesy 87, 43-55 (2013);\n"
    "       DOI <https://doi.org/10.1007/s00190-012-0578-z>; addenda:\n"
    "       <https://geographiclib.sourceforge.io/geod-addenda.html>.\n"
    "\n"
    "       The algorithm for the area of a rhumb polygon is given in Section 3 of\n"
    "       C. F. F. Karney, The area of rhumb polygons, Technical Report, SRI\n"
    "       International (2023); URL: <https://arxiv.org/abs/2303.03219>.\n"
    "\n"
    "AUTHOR\n"
    "       Planimeter was written by Charles Karney.\n"
    "\n"
    "HISTORY\n"
    "       Planimeter was added to GeographicLib,\n"
    "       <https://geographiclib.sourceforge.io>, in version 1.4.\n"
    ;
  return retval;
}
