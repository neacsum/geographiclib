/**
 * \file GeoidEval.cpp
 * \brief Command line utility for evaluating geoid heights
 *
 * Copyright (c) Charles Karney (2009-2017) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 *
 * See the <a href="GeoidEval.1.html">man page</a> for usage information.
 **********************************************************************/

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <GeographicLib/Geoid.hpp>
#include <GeographicLib/DMS.hpp>
#include <GeographicLib/Utility.hpp>
#include <GeographicLib/GeoCoords.hpp>

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
    bool cacheall = false, cachearea = false, verbose = false, cubic = true;
    real caches, cachew, cachen, cachee;
    std::string dir;
    std::string geoid = Geoid::DefaultGeoidName();
    Geoid::convertflag heightmult = Geoid::NONE;
    std::string istring, ifile, ofile, cdelim;
    char lsep = ';';
    bool northp = false, longfirst = false;
    int zonenum = UTMUPS::INVALID;

    for (int m = 1; m < argc; ++m) {
      std::string arg(argv[m]);
      if (arg == "-a") {
        cacheall = true;
        cachearea = false;
      }
      else if (arg == "-c") {
        if (m + 4 >= argc) return usage(1, true);
        cacheall = false;
        cachearea = true;
        try {
          DMS::DecodeLatLon(std::string(argv[m + 1]), std::string(argv[m + 2]),
                            caches, cachew, longfirst);
          DMS::DecodeLatLon(std::string(argv[m + 3]), std::string(argv[m + 4]),
                            cachen, cachee, longfirst);
        }
        catch (const std::exception& e) {
          std::cerr << "Error decoding argument of -c: " << e.what() << "\n";
          return 1;
        }
        m += 4;
      } else if (arg == "--msltohae")
        heightmult = Geoid::GEOIDTOELLIPSOID;
      else if (arg == "--haetomsl")
        heightmult = Geoid::ELLIPSOIDTOGEOID;
      else if (arg == "-w")
        longfirst = !longfirst;
      else if (arg == "-z") {
        if (++m == argc) return usage(1, true);
        std::string zone = argv[m];
        try {
          UTMUPS::DecodeZone(zone, zonenum, northp);
        }
        catch (const std::exception& e) {
          std::cerr << "Error decoding zone: " << e.what() << "\n";
          return 1;
        }
        if (!(zonenum >= UTMUPS::MINZONE && zonenum <= UTMUPS::MAXZONE)) {
          std::cerr << "Illegal zone " << zone << "\n";
          return 1;
        }
      } else if (arg == "-n") {
        if (++m == argc) return usage(1, true);
        geoid = argv[m];
      } else if (arg == "-d") {
        if (++m == argc) return usage(1, true);
        dir = argv[m];
      } else if (arg == "-l")
        cubic = false;
      else if (arg == "-v")
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
          std::cout
            << "\nDefault geoid path = \""   << Geoid::DefaultGeoidPath()
            << "\"\nDefault geoid name = \"" << Geoid::DefaultGeoidName()
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

    int retval = 0;
    try {
      const Geoid g(geoid, dir, cubic);
      try {
        if (cacheall)
          g.CacheAll();
        else if (cachearea)
          g.CacheArea(caches, cachew, cachen, cachee);
      }
      catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\nProceeding without a cache\n";
      }
      if (verbose) {
        std::cerr << "Geoid file: "    << g.GeoidFile()     << "\n"
                  << "Description: "   << g.Description()   << "\n"
                  << "Interpolation: " << g.Interpolation() << "\n"
                  << "Date & Time: "   << g.DateTime()      << "\n"
                  << "Offset (m): "    << g.Offset()        << "\n"
                  << "Scale (m): "     << g.Scale()         << "\n"
                  << "Max error (m): " << g.MaxError()      << "\n"
                  << "RMS error (m): " << g.RMSError()      << "\n";
        if (g.Cache())
          std::cerr
            << "Caching:"
            << "\n SW Corner: " << g.CacheSouth() << " " << g.CacheWest()
            << "\n NE Corner: " << g.CacheNorth() << " " << g.CacheEast()
            << "\n";
      }

      GeoCoords p;
      std::string s, eol, suff;
      const char* spaces = " \t\n\v\f\r,"; // Include comma as space
      while (std::getline(*input, s)) {
        try {
          eol = "\n";
          if (!cdelim.empty()) {
            std::string::size_type m = s.find(cdelim);
            if (m != std::string::npos) {
              eol = " " + s.substr(m) + "\n";
              std::string::size_type m1 =
                m > 0 ? s.find_last_not_of(spaces, m - 1) : std::string::npos;
              s = s.substr(0, m1 != std::string::npos ? m1 + 1 : m);
            }
          }
          real height = 0;
          if (zonenum != UTMUPS::INVALID) {
            // Expect "easting northing" if heightmult == 0, or
            // "easting northing height" if heightmult != 0.
            std::string::size_type pa = 0, pb = 0;
            real easting = 0, northing = 0;
            for (int i = 0; i < (heightmult ? 3 : 2); ++i) {
              if (pb == std::string::npos)
                throw GeographicErr("Incomplete input: " + s);
              // Start of i'th token
              pa = s.find_first_not_of(spaces, pb);
              if (pa == std::string::npos)
                throw GeographicErr("Incomplete input: " + s);
              // End of i'th token
              pb = s.find_first_of(spaces, pa);
              (i == 2 ? height : (i == 0 ? easting : northing)) =
                Utility::val<real>(s.substr(pa, (pb == std::string::npos ?
                                                 pb : pb - pa)));
            }
            p.Reset(zonenum, northp, easting, northing);
            if (heightmult) {
              suff = pb == std::string::npos ? "" : s.substr(pb);
              s = s.substr(0, pa);
            }
          } else {
            if (heightmult) {
              // Treat last token as height
              // pb = last char of last token
              // pa = last char preceding white space
              // px = last char of 2nd last token
              std::string::size_type pb = s.find_last_not_of(spaces);
              std::string::size_type pa = s.find_last_of(spaces, pb);
              if (pa == std::string::npos || pb == std::string::npos)
                throw GeographicErr("Incomplete input: " + s);
              height = Utility::val<real>(s.substr(pa + 1, pb - pa));
              s = s.substr(0, pa + 1);
            }
            p.Reset(s, true, longfirst);
          }
          if (heightmult) {
            real h = g(p.Latitude(), p.Longitude());
            *output << s
                    << Utility::str(height + real(heightmult) * h, 4)
                    << suff << eol;
          } else {
            real h = g(p.Latitude(), p.Longitude());
            *output << Utility::str(h, 4) << eol;
          }
        }
        catch (const std::exception& e) {
          *output << "ERROR: " << e.what() << "\n";
          retval = 1;
        }
      }
    }
    catch (const std::exception& e) {
      std::cerr << "Error reading " << geoid << ": " << e.what() << "\n";
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
    "    GeoidEval [ -n name ] [ -d dir ] [ -l ] [ -a | -c south west north east\n"
    "    ] [ -w ] [ -z zone ] [ --msltohae ] [ --haetomsl ] [ -v ] [\n"
    "    --comment-delimiter commentdelim ] [ --version | -h | --help ] [\n"
    "    --input-file infile | --input-string instring ] [ --line-separator\n"
    "    linesep ] [ --output-file outfile ]\n"
    "\n"
    "For full documentation type:\n"
    "    GeoidEval --help\n"
    "or visit:\n"
    "    https://geographiclib.sourceforge.io/C++/2.4/GeoidEval.1.html\n";
  else
    (retval ? std::cerr : std::cout) << "Man page:\n"
    "NAME\n"
    "       GeoidEval -- look up geoid heights\n"
    "\n"
    "SYNOPSIS\n"
    "       GeoidEval [ -n name ] [ -d dir ] [ -l ] [ -a | -c south west north east\n"
    "       ] [ -w ] [ -z zone ] [ --msltohae ] [ --haetomsl ] [ -v ] [\n"
    "       --comment-delimiter commentdelim ] [ --version | -h | --help ] [\n"
    "       --input-file infile | --input-string instring ] [ --line-separator\n"
    "       linesep ] [ --output-file outfile ]\n"
    "\n"
    "DESCRIPTION\n"
    "       GeoidEval reads in positions on standard input and prints out the\n"
    "       corresponding heights of the geoid model above the WGS84 ellipsoid on\n"
    "       standard output.\n"
    "\n"
    "       Positions are given as latitude and longitude, UTM/UPS, or MGRS, in any\n"
    "       of the formats accepted by GeoConvert(1).  (MGRS coordinates signify\n"
    "       the center of the corresponding MGRS square.)  If the -z option is\n"
    "       specified then the specified zone is prepended to each line of input\n"
    "       (which must be in UTM/UPS coordinates).  This allows a file with UTM\n"
    "       eastings and northings in a single zone to be used as standard input.\n"
    "\n"
    "       More accurate results for the geoid height are provided by Gravity(1).\n"
    "       This utility can also compute the direction of gravity accurately.\n"
    "\n"
    "       The height of the geoid above the ellipsoid, N, is sometimes called the\n"
    "       geoid undulation.  It can be used to convert a height above the\n"
    "       ellipsoid, h, to the corresponding height above the geoid (the\n"
    "       orthometric height, roughly the height above mean sea level), H, using\n"
    "       the relations\n"
    "\n"
    "           h = N + H,   H = -N + h.\n"
    "\n"
    "OPTIONS\n"
    "       -n name\n"
    "           use name for the geoid model instead of the default \"egm96-5\".  See\n"
    "           \"GEOIDS\".\n"
    "\n"
    "       -d dir\n"
    "           read geoid data from dir instead of the default.  See \"GEOIDS\".\n"
    "\n"
    "       -l  use bilinear interpolation instead of cubic.  See \"INTERPOLATION\".\n"
    "\n"
    "       -a  cache the entire data set in memory.  See \"CACHE\".\n"
    "\n"
    "       -c south west north east\n"
    "           cache the data bounded by south west north east in memory.  The\n"
    "           first two arguments specify the SW corner of the cache and the last\n"
    "           two arguments specify the NE corner.  The -w flag specifies that\n"
    "           longitude precedes latitude for these corners, provided that it\n"
    "           appears before -c.  See \"CACHE\".\n"
    "\n"
    "       -w  toggle the longitude first flag (it starts off); if the flag is on,\n"
    "           then when reading geographic coordinates, longitude precedes\n"
    "           latitude (this can be overridden by a hemisphere designator, N, S,\n"
    "           E, W).\n"
    "\n"
    "       -z zone\n"
    "           prefix each line of input by zone, e.g., \"38n\".  This should be\n"
    "           used when the input consists of UTM/UPS eastings and northings.\n"
    "\n"
    "       --msltohae\n"
    "           standard input should include a final token on each line which is\n"
    "           treated as a height (in meters) above the geoid and the output\n"
    "           echoes the input line with the height converted to height above\n"
    "           ellipsoid (HAE).  If -z zone is specified then the third token is\n"
    "           treated as the height; this makes it possible to convert LIDAR data\n"
    "           where each line consists of: easting northing height intensity.\n"
    "\n"
    "       --haetomsl\n"
    "           this is similar to --msltohae except that the height token is\n"
    "           treated as a height (in meters) above the ellipsoid and the output\n"
    "           echoes the input line with the height converted to height above the\n"
    "           geoid (MSL).\n"
    "\n"
    "       -v  print information about the geoid model on standard error before\n"
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
    "       -h  print usage, the default path and name for geoid models, and exit.\n"
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
    "GEOIDS\n"
    "       GeoidEval computes geoid heights by interpolating on the data in a\n"
    "       regularly spaced table (see \"INTERPOLATION\").  The following geoid\n"
    "       grids are available (however, some may not be installed):\n"
    "\n"
    "                                         bilinear error    cubic error\n"
    "          name         geoid    grid     max      rms      max      rms\n"
    "          egm84-30     EGM84    30'      1.546 m  70 mm    0.274 m  14 mm\n"
    "          egm84-15     EGM84    15'      0.413 m  18 mm    0.021 m  1.2 mm\n"
    "          egm96-15     EGM96    15'      1.152 m  40 mm    0.169 m  7.0 mm\n"
    "          egm96-5      EGM96     5'      0.140 m  4.6 mm   .0032 m  0.7 mm\n"
    "          egm2008-5    EGM2008   5'      0.478 m  12 mm    0.294 m  4.5 mm\n"
    "          egm2008-2_5  EGM2008   2.5'    0.135 m  3.2 mm   0.031 m  0.8 mm\n"
    "          egm2008-1    EGM2008   1'      0.025 m  0.8 mm   .0022 m  0.7 mm\n"
    "\n"
    "       By default, the \"egm96-5\" geoid model is used.  This may changed by\n"
    "       setting the environment variable \"GEOGRAPHICLIB_GEOID_NAME\" or with the\n"
    "       -n option.  The errors listed here are estimates of the quantization\n"
    "       and interpolation errors in the reported heights compared to the\n"
    "       specified geoid.\n"
    "\n"
    "       The geoid model data will be loaded from a directory specified at\n"
    "       compile time.  This may changed by setting the environment variables\n"
    "       \"GEOGRAPHICLIB_GEOID_PATH\" or \"GEOGRAPHICLIB_DATA\", or with the -d\n"
    "       option.  The -h option prints the default geoid path and name.  Use the\n"
    "       -v option to ascertain the full path name of the data file.\n"
    "\n"
    "       Instructions for downloading and installing geoid data are available at\n"
    "       <https://geographiclib.sourceforge.io/C++/doc/geoid.html#geoidinst>.\n"
    "\n"
    "       NOTE: all the geoids above apply to the WGS84 ellipsoid (a = 6378137 m,\n"
    "       f = 1/298.257223563) only.\n"
    "\n"
    "INTERPOLATION\n"
    "       Cubic interpolation is used to compute the geoid height unless -l is\n"
    "       specified in which case bilinear interpolation is used.  The cubic\n"
    "       interpolation is based on a least-squares fit of a cubic polynomial to\n"
    "       a 12-point stencil\n"
    "\n"
    "          . 1 1 .\n"
    "          1 2 2 1\n"
    "          1 2 2 1\n"
    "          . 1 1 .\n"
    "\n"
    "       The cubic is constrained to be independent of longitude when evaluating\n"
    "       the height at one of the poles.  Cubic interpolation is considerably\n"
    "       more accurate than bilinear; however it results in small\n"
    "       discontinuities in the returned height on cell boundaries.\n"
    "\n"
    "CACHE\n"
    "       By default, the data file is randomly read to compute the geoid heights\n"
    "       at the input positions.  Usually this is sufficient for interactive\n"
    "       use.  If many heights are to be computed, use -c south west north east\n"
    "       to notify GeoidEval to read a rectangle of data into memory; heights\n"
    "       within the this rectangle can then be computed without any disk access.\n"
    "       If -a is specified all the geoid data is read; in the case of\n"
    "       \"egm2008-1\", this requires about 0.5 GB of RAM.  The evaluation of\n"
    "       heights outside the cached area causes the necessary data to be read\n"
    "       from disk.  Use the -v option to verify the size of the cache.\n"
    "\n"
    "       Regardless of whether any cache is requested (with the -a or -c\n"
    "       options), the data for the last grid cell in cached.  This allows the\n"
    "       geoid height along a continuous path to be returned with little disk\n"
    "       overhead.\n"
    "\n"
    "ENVIRONMENT\n"
    "       GEOGRAPHICLIB_GEOID_NAME\n"
    "           Override the compile-time default geoid name of \"egm96-5\".  The -h\n"
    "           option reports the value of GEOGRAPHICLIB_GEOID_NAME, if defined,\n"
    "           otherwise it reports the compile-time value.  If the -n name option\n"
    "           is used, then name takes precedence.\n"
    "\n"
    "       GEOGRAPHICLIB_GEOID_PATH\n"
    "           Override the compile-time default geoid path.  This is typically\n"
    "           \"/usr/local/share/GeographicLib/geoids\" on Unix-like systems and\n"
    "           \"C:/ProgramData/GeographicLib/geoids\" on Windows systems.  The -h\n"
    "           option reports the value of GEOGRAPHICLIB_GEOID_PATH, if defined,\n"
    "           otherwise it reports the compile-time value.  If the -d dir option\n"
    "           is used, then dir takes precedence.\n"
    "\n"
    "       GEOGRAPHICLIB_DATA\n"
    "           Another way of overriding the compile-time default geoid path.  If\n"
    "           it is set (and if GEOGRAPHICLIB_GEOID_PATH is not set), then\n"
    "           $GEOGRAPHICLIB_DATA/geoids is used.\n"
    "\n"
    "ERRORS\n"
    "       An illegal line of input will print an error message to standard output\n"
    "       beginning with \"ERROR:\" and causes GeoidEval to return an exit code of\n"
    "       1.  However, an error does not cause GeoidEval to terminate; following\n"
    "       lines will be converted.\n"
    "\n"
    "ABBREVIATIONS\n"
    "       The geoid is usually approximated by an \"earth gravity model\". The\n"
    "       models published by the NGA are:\n"
    "\n"
    "       EGM84\n"
    "           An earth gravity model published by the NGA in 1984,\n"
    "           <https://earth-info.nga.mil/index.php?dir=wgs84&action=wgs84#tab_egm84>.\n"
    "\n"
    "       EGM96\n"
    "           An earth gravity model published by the NGA in 1996,\n"
    "           <https://earth-info.nga.mil/index.php?dir=wgs84&action=wgs84#tab_egm96>.\n"
    "\n"
    "       EGM2008\n"
    "           An earth gravity model published by the NGA in 2008,\n"
    "           <https://earth-info.nga.mil/index.php?dir=wgs84&action=wgs84#tab_egm2008>.\n"
    "\n"
    "       WGS84\n"
    "           World Geodetic System 1984, <https://en.wikipedia.org/wiki/WGS84>.\n"
    "\n"
    "       HAE Height above the WGS84 ellipsoid.\n"
    "\n"
    "       MSL Mean sea level, used as a convenient short hand for the geoid.\n"
    "           (However, typically, the geoid differs by a few meters from mean\n"
    "           sea level.)\n"
    "\n"
    "EXAMPLES\n"
    "       The height of the EGM96 geoid at Timbuktu\n"
    "\n"
    "           echo 16:46:33N 3:00:34W | GeoidEval\n"
    "           => 28.7068 -0.02e-6 -1.73e-6\n"
    "\n"
    "       The first number returned is the height of the geoid and the 2nd and\n"
    "       3rd are its slopes in the northerly and easterly directions.\n"
    "\n"
    "       Convert a point in UTM zone 18n from MSL to HAE\n"
    "\n"
    "          echo 531595 4468135 23 | GeoidEval --msltohae -z 18n\n"
    "          => 531595 4468135 -10.842\n"
    "\n"
    "SEE ALSO\n"
    "       GeoConvert(1), Gravity(1), geographiclib-get-geoids(8).\n"
    "\n"
    "       An online version of this utility is availbable at\n"
    "       <https://geographiclib.sourceforge.io/cgi-bin/GeoidEval>.\n"
    "\n"
    "AUTHOR\n"
    "       GeoidEval was written by Charles Karney.\n"
    "\n"
    "HISTORY\n"
    "       GeoidEval was added to GeographicLib,\n"
    "       <https://geographiclib.sourceforge.io>, in 2009-09.\n"
    ;
  return retval;
}
