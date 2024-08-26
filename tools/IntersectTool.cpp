/**
 * \file IntersectTool.cpp
 * \brief Command line utility for geodesic intersections
 *
 * Copyright (c) Charles Karney (2023) <karney@alum.mit.edu> and licensed under
 * the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 *
 * See the <a href="IntersectTool.1.html">man page</a> for usage information.
 **********************************************************************/

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <GeographicLib/Geodesic.hpp>
#include <GeographicLib/DMS.hpp>
#include <GeographicLib/Utility.hpp>
#include <GeographicLib/Intersect.hpp>

int usage (int retval, bool brief);

using namespace GeographicLib;

int main(int argc, const char* const argv[]) {
  try {
    enum { CLOSE = 0, OFFSET, NEXT, SEGMENT };
    Utility::set_digits();
    real
      a = Constants::WGS84_a(),
      f = Constants::WGS84_f(),
      maxdist = -1;
    bool exact = false, check = false, longfirst = false;
    int prec = 3, mode = CLOSE;
    std::string istring, ifile, ofile, cdelim;
    char lsep = ';';

    for (int m = 1; m < argc; ++m) {
      std::string arg(argv[m]);
      if (arg == "-e") {
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
      } else if (arg == "-E")
        exact = true;
      else if (arg == "-p") {
        if (++m == argc) return usage(1, true);
        try {
          prec = Utility::val<int>(std::string(argv[m]));
        }
        catch (const std::exception&) {
          std::cerr << "Precision " << argv[m] << " is not a number\n";
          return 1;
        }
      } else if (arg == "-R") {
        if (++m == argc) return usage(1, true);
        try {
          maxdist = Utility::val<real>(std::string(argv[m]));
        }
        catch (const std::exception&) {
          std::cerr << "Maxdist " << argv[m] << " is not a number\n";
          return 1;
        }
        if (!(maxdist >= 0)) {
          std::cerr << "Maxdist must be nonnegative\n";
          return 1;
        }
      } else if (arg == "-c")
        mode = CLOSE;
      else if (arg == "-o")
        mode = OFFSET;
      else if (arg == "-n")
        mode = NEXT;
      else if (arg == "-i")
        mode = SEGMENT;
      else if (arg == "-C")
        check = true;
      else if (arg == "-w")
        longfirst = true;
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

    Geodesic geod(a, f, exact);
    Intersect intersect(geod);
    real latX1, lonX1, aziX, latY1, lonY1, aziY, latX2, lonX2, latY2, lonY2,
      x0 = 0, y0 = 0, x, y;
    std::string inp[8], s, sc, eol;
    std::istringstream str;
    int retval = 0,
      ninp = mode == CLOSE ? 6 : (mode == NEXT ? 4 :
                                  8); // mode == OFFSET || mode == SEGMENT
    GeodesicLine lineX, lineY;
    unsigned caps = Intersect::LineCaps;
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
        for (int i = 0; i < ninp; ++i) {
          if (!(str >> inp[i]))
          throw GeographicErr("Incomplete input: " + s);
        }
        if (str >> sc)
          throw GeographicErr("Extraneous input: " + sc);
        if (mode == CLOSE || mode == OFFSET) {
          DMS::DecodeLatLon(inp[0], inp[1], latX1, lonX1, longfirst);
          aziX = DMS::DecodeAzimuth(inp[2]);
          DMS::DecodeLatLon(inp[3], inp[4], latY1, lonY1, longfirst);
          aziY = DMS::DecodeAzimuth(inp[5]);
          if (mode == OFFSET) {
            x0 = Utility::val<real>(inp[6]);
            y0 = Utility::val<real>(inp[7]);
          } else
            x0 = y0 = 0;
          lineX = geod.Line(latX1, lonX1, aziX, caps);
          lineY = geod.Line(latY1, lonY1, aziY, caps);
        } else if (mode == NEXT) {
          DMS::DecodeLatLon(inp[0], inp[1], latX1, lonX1, longfirst);
          aziX = DMS::DecodeAzimuth(inp[2]);
          aziY = DMS::DecodeAzimuth(inp[3]);
          lineX = geod.Line(latX1, lonX1, aziX, caps);
          lineY = geod.Line(latX1, lonX1, aziY, caps);
        } else {                // mode == SEGMENT
          DMS::DecodeLatLon(inp[0], inp[1], latX1, lonX1, longfirst);
          DMS::DecodeLatLon(inp[2], inp[3], latX2, lonX2, longfirst);
          DMS::DecodeLatLon(inp[4], inp[5], latY1, lonY1, longfirst);
          DMS::DecodeLatLon(inp[6], inp[7], latY2, lonY2, longfirst);
          lineX = geod.InverseLine(latX1, lonX1, latX2, lonX2,
                                   Intersect::LineCaps);
          lineY = geod.InverseLine(latY1, lonY1, latY2, lonY2,
                                   Intersect::LineCaps);
          x0 = lineX.Distance()/2;
          y0 = lineY.Distance()/2;
        }
        std::pair<real, real> p0(x0, y0);
        if (maxdist < 0) {
          int segmode = 0, c;
          auto p = mode == CLOSE || mode == OFFSET ?
            intersect.Closest(lineX, lineY, p0, &c) :
            mode == NEXT ? intersect.Next(lineX, lineY, &c) :
            intersect.Segment(lineX, lineY, segmode, &c);
          x = p.first; y = p.second;
          *output << Utility::str(x, prec) << " "
                  << Utility::str(y, prec) << " " << c;
          if (mode == SEGMENT)
            *output << " " << segmode;
          *output << eol;
          if (check) {
            lineX.Position(x, latX2, lonX2);
            lineY.Position(y, latY2, lonY2);
            real sXY;
            geod.Inverse(latX2, lonX2, latY2, lonY2, sXY);
            std::cerr << Utility::str(longfirst ? lonX2 : latX2, prec+5) << " "
                      << Utility::str(longfirst ? latX2 : lonX2, prec+5) << " "
                      << Utility::str(longfirst ? lonY2 : latY2, prec+5) << " "
                      << Utility::str(longfirst ? latY2 : lonY2, prec+5) << " "
                      << Utility::str(sXY, prec) << eol;
          }
        } else {
          std::vector<int> c;
          auto v = intersect.All(lineX, lineY, maxdist, c, p0);
          unsigned n = unsigned(v.size());
          for (unsigned i = 0; i < n; ++i) {
            x = v[i].first; y = v[i].second;
            *output << Utility::str(x, prec) << " " << Utility::str(y, prec)
                    << " " << c[i] << " "
                    << Utility::str(Intersect::Dist(v[i], p0), prec)
                    << eol;
            if (check) {
              lineX.Position(x, latX2, lonX2);
              lineY.Position(y, latY2, lonY2);
              real sXY;
              geod.Inverse(latX2, lonX2, latY2, lonY2, sXY);
              std::cerr << Utility::str(longfirst ? lonX2 : latX2, prec+5) << " "
                        << Utility::str(longfirst ? latX2 : lonX2, prec+5) << " "
                        << Utility::str(longfirst ? lonY2 : latY2, prec+5) << " "
                        << Utility::str(longfirst ? latY2 : lonY2, prec+5) << " "
                        << Utility::str(sXY, prec) << eol;
            }
          }
          *output << "nan nan 0 nan" << eol;
          if (check)
            std::cerr << "nan nan nan nan nan" << eol;
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
    "    IntersectTool [ -c | -n | -i | -o | [ -R maxdist ] [ -e a f] [ -E ] [\n"
    "    -w ] [ -p prec ] [ --comment-delimiter commentdelim ] [ --version | -h\n"
    "    | --help ] [ --input-file infile | --input-string instring ] [\n"
    "    --line-separator linesep ] [ --output-file outfile ]\n"
    "\n"
    "For full documentation type:\n"
    "    IntersectTool --help\n"
    "or visit:\n"
    "    https://geographiclib.sourceforge.io/C++/2.4/IntersectTool.1.html\n";
  else
    (retval ? std::cerr : std::cout) << "Man page:\n"
    "NAME\n"
    "       IntersectTool -- perform rhumb line calculations\n"
    "\n"
    "SYNOPSIS\n"
    "       IntersectTool [ -c | -n | -i | -o | [ -R maxdist ] [ -e a f] [ -E ] [\n"
    "       -w ] [ -p prec ] [ --comment-delimiter commentdelim ] [ --version | -h\n"
    "       | --help ] [ --input-file infile | --input-string instring ] [\n"
    "       --line-separator linesep ] [ --output-file outfile ]\n"
    "\n"
    "DESCRIPTION\n"
    "       IntersectTool finds the intersection of two geodesics X and Y.  The\n"
    "       geodesics may either be specified as a location and an azimuth, latX\n"
    "       lonX aziX, or as the shortest geodesic segment between two locations,\n"
    "       latX1 lonX1 and latX2 lonX2.  The intersection is then specified as the\n"
    "       displacements, x and y, along the geodesics X and Y from the starting\n"
    "       points to the intersection.  In the case of the intersection of\n"
    "       geodesic segments, the starting point is first point specified for X or\n"
    "       Y.\n"
    "\n"
    "       Usually this tool returns the closest intersection defined as the one\n"
    "       that minimizes the \"L1\" distance, |x| + |y|.  However, it is possible\n"
    "       to specify an \"origin\" x0 and y0 when determining closeness so that the\n"
    "       intersection which minimizes |x - x0| + |y - y0| is returned.\n"
    "\n"
    "       In the case of intersecting segments the origin is taken to be the\n"
    "       midpoints of the segments; x0 is half the distance from X1 to X2.  In\n"
    "       addition a flag is returned specifying whether the intersection is\n"
    "       \"within\" the segments.\n"
    "\n"
    "       The tool also returns a \"coincidence indicator\" c.  This is typically\n"
    "       0.  However if the geodesics lie on top of one another at the point of\n"
    "       intersection, then c is set to 1, if they are parallel, and -1, if they\n"
    "       are antiparallel.\n"
    "\n"
    "       IntersectTool operates in one of three modes:\n"
    "\n"
    "       1.  With the -c option (the default), IntersectTool accepts lines on\n"
    "           the standard input containing latX lonX aziX latY lonY aziY,\n"
    "           specifying two geodesic lines X and Y, and prints the location of\n"
    "           the closest intersection x y c on standard output.\n"
    "\n"
    "       2.  With the -n option, IntersectTool accepts lines on the standard\n"
    "           input containing latX lonX aziX aziY aziY, specifying a point where\n"
    "           two geodesic lines X and Y intersect, and prints the location of\n"
    "           the next closest intersection x y c on standard output.\n"
    "\n"
    "       3.  With the -i option, IntersectTool accepts lines on the standard\n"
    "           input containing latX1 lonX1 latX2 lonX2 latY1 lonY1 latY2 lonY2,\n"
    "           specifying two geodesic segments X1-X2 and Y1-Y2, and prints x y c\n"
    "           k on standard output.  Here k is a flag in [-4,4] specifying\n"
    "           whether the intersection is within the segments (0) or not (non-\n"
    "           zero).  x and y give the distances from X1 and Y1 respectively.  k\n"
    "           is set to 3 kx + ky where kx = -1 if x < 0, 0 if 0 <= x <= sx, 1 if\n"
    "           sx < x, and similarly for ky; sx is the length of the segment\n"
    "           X1-X2.\n"
    "\n"
    "       4.  With the -o option, IntersectTool accepts lines on the standard\n"
    "           input containing latX lonX aziX latY lonY aziY x0 y0, specifying\n"
    "           two geodesic lines X and Y and two offsets, and prints x y c on\n"
    "           standard output where [x, y] is the intersection closest to [x0,\n"
    "           y0].\n"
    "\n"
    "OPTIONS\n"
    "       -c  find the closest intersection (see 1 above).\n"
    "\n"
    "       -n  find the intersection closest to a given intersection (see 2\n"
    "           above).\n"
    "\n"
    "       -i  find the intersection of two geodesic segments (see 3 above).\n"
    "\n"
    "       -o  find the closest intersection with an offset.\n"
    "\n"
    "       -R maxdist\n"
    "           modifies the four modes to return all the intersections within an\n"
    "           L1 distance, maxdist, of the relevant origin: [0, 0] for -c and -n,\n"
    "           the midpoints of the segments for -i, and [x0, y0] for -o.  For\n"
    "           each intersection, x y c z is printed on standard output.  Here z\n"
    "           is the L1 distance of the intersection from the origin and the\n"
    "           intersections are sorted by the distances.  A line \"nan nan 0 nan\"\n"
    "           is added after the intersections, so that the output can be\n"
    "           associated with the correct lines of the input.  The number of\n"
    "           intersections scales as (maxdist/(pi a))^2.\n"
    "\n"
    "       -C  check the intersections.  For each computed intersection, print on\n"
    "           standard error a line latX lonX latY lonY sXY giving the computed\n"
    "           positions of the intersections points on X and Y and the distance\n"
    "           between them.  If -w is specified, the longitude is given before\n"
    "           the latitude.\n"
    "\n"
    "       -e a f\n"
    "           specify the ellipsoid via the equatorial radius, a and the\n"
    "           flattening, f.  Setting f = 0 results in a sphere.  Specify f < 0\n"
    "           for a prolate ellipsoid.  A simple fraction, e.g., 1/297, is\n"
    "           allowed for f.  By default, the WGS84 ellipsoid is used, a =\n"
    "           6378137 m, f = 1/298.257223563.\n"
    "\n"
    "       -E  use \"exact\" algorithms (based on elliptic integrals) for the\n"
    "           geodesic calculations.  These are more accurate than the (default)\n"
    "           series expansions for |f| > 0.02.\n"
    "\n"
    "       -w  on input, longitude precedes latitude (except that on input this\n"
    "           can be overridden by a hemisphere designator, N, S, E, W).\n"
    "\n"
    "       -p prec\n"
    "           set the output precision to prec (default 3); prec is the precision\n"
    "           relative to 1 m.  See \"PRECISION\".\n"
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
    "       IntersectTool measures all angles in degrees and all lengths in meters.\n"
    "       On input angles (latitude, longitude, azimuth) can be as decimal\n"
    "       degrees or degrees, minutes, seconds.  For example, \"40d30\", \"40d30'\",\n"
    "       \"40:30\", \"40.5d\", and 40.5 are all equivalent.  By default, latitude\n"
    "       precedes longitude for each point (the -w flag switches this\n"
    "       convention); however either may be given first by appending (or\n"
    "       prepending) N or S to the latitude and E or W to the longitude.\n"
    "       Azimuths are measured clockwise from north; however this may be\n"
    "       overridden with E or W.\n"
    "\n"
    "       For details on the allowed formats for angles, see the \"GEOGRAPHIC\n"
    "       COORDINATES\" section of GeoConvert(1).\n"
    "\n"
    "PRECISION\n"
    "       prec gives precision of the output with prec = 0 giving 1 m precision,\n"
    "       prec = 3 giving 1 mm precision, etc.  prec is the number of digits\n"
    "       after the decimal point for lengths.  The latitude and longitude\n"
    "       printed to standard error with the -C option are given in decimal\n"
    "       degrees with prec + 5 digits after the decimal point.  The minimum\n"
    "       value of prec is 0 and the maximum is 10.\n"
    "\n"
    "ERRORS\n"
    "       An illegal line of input will print an error message to standard output\n"
    "       beginning with \"ERROR:\" and causes IntersectTool to return an exit code\n"
    "       of 1.  However, an error does not cause IntersectTool to terminate;\n"
    "       following lines will be converted.\n"
    "\n"
    "ACCURACY\n"
    "       This tool will give nearly full double precision accuracy for |f| <\n"
    "       0.02.  If the -E option is given, full accuracy is achieved for -1/4 <\n"
    "       f < 1/5.  The tool had not been tested outside this range.\n"
    "\n"
    "EXAMPLES\n"
    "       A vessel leaves Plymouth 50N 4W on a geodesic path with initial heading\n"
    "       147.7W.  When will it first cross the equator?\n"
    "\n"
    "          echo 50N 4W 147.7W 0 0 90 | IntersectTool -c -p 0 -C\n"
    "\n"
    "          6058049 -3311253 0\n"
    "          0.00000 -29.74549 -0.00000 -29.74549 0\n"
    "\n"
    "       Answer: after 6058km at longitude 29.7W.  When will it cross the date\n"
    "       line, longitude 180E?  Here we need to use -R because there a closer\n"
    "       intersection on the prime meridian:\n"
    "\n"
    "          echo 50N 4W 147.7W 0 180 0 | IntersectTool -c -p 0 -C -R 2.6e7\n"
    "\n"
    "          -494582 14052230 0 14546812\n"
    "          53.69260 0.00000 53.69260 0.00000 0\n"
    "          19529110 -5932344 0 25461454\n"
    "          -53.51867 180.00000 -53.51867 180.00000 0\n"
    "          nan nan 0 nan\n"
    "          nan nan nan nan nan\n"
    "\n"
    "       We want the second result: after 19529 km at latitude 53.5S.\n"
    "\n"
    "SEE ALSO\n"
    "       GeoConvert(1), GeodSolve(1).\n"
    "\n"
    "       This solution for intersections is described in C. F. F. Karney,\n"
    "       Geodesic intersections, J. Surveying Eng. 150(3), 04024005:1-9 (2024),\n"
    "       DOI: <https://doi.org/10.1061/JSUED2.SUENG-1483>; preprint\n"
    "       <https://arxiv.org/abs/2308.00495>.  It is based on the work of S.\n"
    "       Baseldga and J. C. Martinez-Llario, Intersection and point-to-line\n"
    "       solutions for geodesics on the ellipsoid, Stud. Geophys. Geod. 62,\n"
    "       353-363 (2018); DOI: <https://doi.org/10.1007/s11200-017-1020-z>;\n"
    "\n"
    "AUTHOR\n"
    "       IntersectTool was written by Charles Karney.\n"
    "\n"
    "HISTORY\n"
    "       IntersectTool was added to GeographicLib,\n"
    "       <https://geographiclib.sourceforge.io>, in version 2.3.\n"
    ;
  return retval;
}
