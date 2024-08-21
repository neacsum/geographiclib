/**
 * \file GeoConvert.cpp
 * \brief Command line utility for geographic coordinate conversions
 *
 * Copyright (c) Charles Karney (2008-2017) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 *
 * See the <a href="GeoConvert.1.html">man page</a> for usage information.
 **********************************************************************/

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <GeographicLib/GeoCoords.hpp>
#include <GeographicLib/DMS.hpp>
#include <GeographicLib/Utility.hpp>
#include <GeographicLib/MGRS.hpp>

int usage (int retval, bool brief);

int main(int argc, const char* const argv[]) {
  try {
    using namespace GeographicLib;
    typedef Math::real real;
    Utility::set_digits();
    enum { GEOGRAPHIC, DMS, UTMUPS, MGRS, CONVERGENCE };
    int outputmode = GEOGRAPHIC;
    int prec = 0;
    int zone = UTMUPS::MATCH;
    bool centerp = true, longfirst = false;
    std::string istring, ifile, ofile, cdelim;
    char lsep = ';', dmssep = char(0);
    bool sethemisphere = false, northp = false, abbrev = true, latch = false;

    for (int m = 1; m < argc; ++m) {
      std::string arg(argv[m]);
      if (arg == "-g")
        outputmode = GEOGRAPHIC;
      else if (arg == "-d") {
        outputmode = DMS;
        dmssep = '\0';
      } else if (arg == "-:") {
        outputmode = DMS;
        dmssep = ':';
      } else if (arg == "-u")
        outputmode = UTMUPS;
      else if (arg == "-m")
        outputmode = MGRS;
      else if (arg == "-c")
        outputmode = CONVERGENCE;
      else if (arg == "-n")
        centerp = false;
      else if (arg == "-z") {
        if (++m == argc) return usage(1, true);
        std::string zonestr(argv[m]);
        try {
          UTMUPS::DecodeZone(zonestr, zone, northp);
          sethemisphere = true;
        }
        catch (const std::exception&) {
          std::istringstream str(zonestr);
          char c;
          if (!(str >> zone) || (str >> c)) {
            std::cerr << "Zone " << zonestr
                      << " is not a number or zone+hemisphere\n";
            return 1;
          }
          if (!(zone >= UTMUPS::MINZONE && zone <= UTMUPS::MAXZONE)) {
            std::cerr << "Zone " << zone << " not in [0, 60]\n";
            return 1;
          }
          sethemisphere = false;
        }
        latch = false;
      } else if (arg == "-s") {
        zone = UTMUPS::STANDARD;
        sethemisphere = false;
        latch = false;
      } else if (arg == "-S") {
        zone = UTMUPS::STANDARD;
        sethemisphere = false;
        latch = true;
      } else if (arg == "-t") {
        zone = UTMUPS::UTM;
        sethemisphere = false;
        latch = false;
      } else if (arg == "-T") {
        zone = UTMUPS::UTM;
        sethemisphere = false;
        latch = true;
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
      } else if (arg == "-l")
        abbrev = false;
      else if (arg == "-a")
        abbrev = true;
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
        MGRS::Check();
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

    GeoCoords p;
    std::string s, eol;
    std::string os;
    int retval = 0;

    while (std::getline(*input, s)) {
      eol = "\n";
      try {
        if (!cdelim.empty()) {
          std::string::size_type m = s.find(cdelim);
          if (m != std::string::npos) {
            eol = " " + s.substr(m) + "\n";
            s = s.substr(0, m);
          }
        }
        p.Reset(s, centerp, longfirst);
        p.SetAltZone(zone);
        switch (outputmode) {
        case GEOGRAPHIC:
          os = p.GeoRepresentation(prec, longfirst);
          break;
        case DMS:
          os = p.DMSRepresentation(prec, longfirst, dmssep);
          break;
        case UTMUPS:
          os = (sethemisphere
                ? p.AltUTMUPSRepresentation(northp, prec, abbrev)
                : p.AltUTMUPSRepresentation(prec, abbrev));
          break;
        case MGRS:
          os = p.AltMGRSRepresentation(prec);
          break;
        case CONVERGENCE:
          {
            real
              gamma = p.AltConvergence(),
              k = p.AltScale();
            int prec1 = std::max(-5, std::min(Math::extra_digits() + 8, prec));
            os = Utility::str(gamma, prec1 + 5) + " "
              + Utility::str(k, prec1 + 7);
          }
        }
        if (latch &&
            zone < UTMUPS::MINZONE && p.AltZone() >= UTMUPS::MINZONE) {
          zone = p.AltZone();
          northp = p.Northp();
          sethemisphere = true;
          latch = false;
        }
      }
      catch (const std::exception& e) {
        // Write error message to cout so output lines match input lines
        os = std::string("ERROR: ") + e.what();
        retval = 1;
      }
      *output << os << eol;
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
    "    GeoConvert [ -g | -d | -: | -u | -m | -c ] [ -z zone | -s | -t | -S |\n"
    "    -T ] [ -n ] [ -w ] [ -p prec ] [ -l | -a ] [ --comment-delimiter\n"
    "    commentdelim ] [ --version | -h | --help ] [ --input-file infile |\n"
    "    --input-string instring ] [ --line-separator linesep ] [ --output-file\n"
    "    outfile ]\n"
    "\n"
    "For full documentation type:\n"
    "    GeoConvert --help\n"
    "or visit:\n"
    "    https://geographiclib.sourceforge.io/C++/2.4/GeoConvert.1.html\n";
  else
    (retval ? std::cerr : std::cout) << "Man page:\n"
    "NAME\n"
    "       GeoConvert -- convert geographic coordinates\n"
    "\n"
    "SYNOPSIS\n"
    "       GeoConvert [ -g | -d | -: | -u | -m | -c ] [ -z zone | -s | -t | -S |\n"
    "       -T ] [ -n ] [ -w ] [ -p prec ] [ -l | -a ] [ --comment-delimiter\n"
    "       commentdelim ] [ --version | -h | --help ] [ --input-file infile |\n"
    "       --input-string instring ] [ --line-separator linesep ] [ --output-file\n"
    "       outfile ]\n"
    "\n"
    "DESCRIPTION\n"
    "       GeoConvert reads from standard input interpreting each line as a\n"
    "       geographic coordinate and prints the coordinate in the format specified\n"
    "       by the options on standard output.  The input is interpreted in one of\n"
    "       three different ways depending on how many space or comma delimited\n"
    "       tokens there are on the line.  The options -g, -d, -u, and -m govern\n"
    "       the format of output.  In all cases, the WGS84 model of the earth is\n"
    "       used (a = 6378137 m, f = 1/298.257223563).\n"
    "\n"
    "       geographic\n"
    "           2 tokens (output options -g, -d, or -:) given as latitude longitude\n"
    "           using decimal degrees or degrees, minutes, and seconds.  Latitude\n"
    "           is given first (unless the -w option is given).  See \"GEOGRAPHIC\n"
    "           COORDINATES\" for a description of the format.  For example, the\n"
    "           following are all equivalent\n"
    "\n"
    "               33.3 44.4\n"
    "               E44.4 N33.3\n"
    "               33d18'N 44d24'E\n"
    "               44d24 33d18N\n"
    "               33:18 +44:24\n"
    "\n"
    "       UTM/UPS\n"
    "           3 tokens (output option -u) given as zone+hemisphere easting\n"
    "           northing or easting northing zone+hemisphere, where hemisphere is\n"
    "           either n (or north) or s (or south).  The zone is absent for a UPS\n"
    "           specification.  For example,\n"
    "\n"
    "               38n 444140.54 3684706.36\n"
    "               444140.54 3684706.36 38n\n"
    "               s 2173854.98 2985980.58\n"
    "               2173854.98 2985980.58 s\n"
    "\n"
    "       MRGS\n"
    "           1 token (output option -m) is used to specify the center of an MGRS\n"
    "           grid square.  For example,\n"
    "\n"
    "               38SMB4484\n"
    "               38SMB44140847064\n"
    "\n"
    "OPTIONS\n"
    "       -g  output latitude and longitude using decimal degrees.  Default\n"
    "           output mode.\n"
    "\n"
    "       -d  output latitude and longitude using degrees, minutes, and seconds\n"
    "           (DMS).\n"
    "\n"
    "       -:  like -d, except use : as a separator instead of the d, ', and \"\n"
    "           delimiters.\n"
    "\n"
    "       -u  output UTM or UPS.\n"
    "\n"
    "       -m  output MGRS.\n"
    "\n"
    "       -c  output meridian convergence and scale for the corresponding UTM or\n"
    "           UPS projection.  The meridian convergence is the bearing of grid\n"
    "           north given as degrees clockwise from true north.\n"
    "\n"
    "       -z zone\n"
    "           set the zone to zone for output.  Use either 0 < zone <= 60 for a\n"
    "           UTM zone or zone = 0 for UPS.  Alternatively use a zone+hemisphere\n"
    "           designation, e.g., 38n.  See \"ZONE\".\n"
    "\n"
    "       -s  use the standard UPS and UTM zones.\n"
    "\n"
    "       -t  similar to -s but forces UPS regions to the closest UTM zone.\n"
    "\n"
    "       -S or -T\n"
    "           behave the same as -s and -t, respectively, until the first legal\n"
    "           conversion is performed.  For subsequent points, the zone and\n"
    "           hemisphere of that conversion are used.  This enables a sequence of\n"
    "           points to be converted into UTM or UPS using a consistent\n"
    "           coordinate system.\n"
    "\n"
    "       -n  on input, MGRS coordinates refer to the south-west corner of the\n"
    "           MGRS square instead of the center; see \"MGRS\".\n"
    "\n"
    "       -w  toggle the longitude first flag (it starts off); if the flag is on,\n"
    "           then on input and output, longitude precedes latitude (except that,\n"
    "           on input, this can be overridden by a hemisphere designator, N, S,\n"
    "           E, W).\n"
    "\n"
    "       -p prec\n"
    "           set the output precision to prec (default 0); prec is the precision\n"
    "           relative to 1 m.  See \"PRECISION\".\n"
    "\n"
    "       -l  on output, UTM/UPS uses the long forms north and south to designate\n"
    "           the hemisphere instead of n or s.\n"
    "\n"
    "       -a  on output, UTM/UPS uses the abbreviations n and s to designate the\n"
    "           hemisphere instead of north or south; this is the default\n"
    "           representation.\n"
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
    "PRECISION\n"
    "       prec gives precision of the output with prec = 0 giving 1 m precision,\n"
    "       prec = 3 giving 1 mm precision, etc.  prec is the number of digits\n"
    "       after the decimal point for UTM/UPS.  For MGRS, The number of digits\n"
    "       per coordinate is 5 + prec; prec = -6 results in just the grid zone.\n"
    "       For decimal degrees, the number of digits after the decimal point is 5\n"
    "       + prec.  For DMS (degree, minute, seconds) output, the number of digits\n"
    "       after the decimal point in the seconds components is 1 + prec; if this\n"
    "       is negative then use minutes (prec = -2 or -3) or degrees (prec <= -4)\n"
    "       as the least significant component.  Print convergence, resp. scale,\n"
    "       with 5 + prec, resp. 7 + prec, digits after the decimal point.  The\n"
    "       minimum value of prec is -5 (-6 for MGRS) and the maximum is 9 for\n"
    "       UTM/UPS, 9 for decimal degrees, 10 for DMS, 6 for MGRS, and 8 for\n"
    "       convergence and scale.\n"
    "\n"
    "GEOGRAPHIC COORDINATES\n"
    "       The utility accepts geographic coordinates, latitude and longitude, in\n"
    "       a number of common formats.  Latitude precedes longitude, unless the -w\n"
    "       option is given which switches this convention.  On input, either\n"
    "       coordinate may be given first by appending or prepending N or S to the\n"
    "       latitude and E or W to the longitude.  These hemisphere designators\n"
    "       carry an implied sign, positive for N and E and negative for S and W.\n"
    "       This sign multiplies any +/- sign prefixing the coordinate.  The\n"
    "       coordinates may be given as decimal degree or as degrees, minutes, and\n"
    "       seconds.  d, ', and \" are used to denote degrees, minutes, and seconds,\n"
    "       with the least significant designator optional.  (See \"QUOTING\" for how\n"
    "       to quote the characters ' and \" when entering coordinates on the\n"
    "       command line.)  Alternatively, : (colon) may be used to separate the\n"
    "       various components.  Only the final component of coordinate can include\n"
    "       a decimal point, and the minutes and seconds components must be less\n"
    "       than 60.\n"
    "\n"
    "       It is also possible to carry out addition or subtraction operations in\n"
    "       geographic coordinates.  If the coordinate includes interior signs\n"
    "       (i.e., not at the beginning or immediately after an initial hemisphere\n"
    "       designator), then the coordinate is split before such signs; the pieces\n"
    "       are parsed separately and the results summed.  For example the point\n"
    "       15\" east of 39N 70W is\n"
    "\n"
    "           39N 70W+0:0:15E\n"
    "\n"
    "       WARNING: \"Exponential\" notation is not recognized for geographic\n"
    "       coordinates.  Thus 7.0E1 is illegal, while 7.0E+1 is parsed as (7.0E) +\n"
    "       (+1), yielding the same result as 8.0E.\n"
    "\n"
    "       Various unicode characters (encoded with UTF-8) may also be used to\n"
    "       denote degrees, minutes, and seconds, e.g., the degree, prime, and\n"
    "       double prime symbols; in addition two single quotes can be used to\n"
    "       represent \".\n"
    "\n"
    "       The other GeographicLib utilities use the same rules for interpreting\n"
    "       geographic coordinates; in addition, azimuths and arc lengths are\n"
    "       interpreted the same way.\n"
    "\n"
    "QUOTING\n"
    "       Unfortunately the characters ' and \" have special meanings in many\n"
    "       shells and have to be entered with care.  However note (1) that the\n"
    "       trailing designator is optional and that (2) you can use colons as a\n"
    "       separator character.  Thus 10d20' can be entered as 10d20 or 10:20 and\n"
    "       10d20'30\" can be entered as 10:20:30.\n"
    "\n"
    "       Unix shells (sh, bash, tsch)\n"
    "           The characters ' and \" can be quoted by preceding them with a \\\n"
    "           (backslash); or you can quote a string containing ' with a pair of\n"
    "           \"s.  The two alternatives are illustrated by\n"
    "\n"
    "              echo 10d20\\'30\\\" \"20d30'40\" | GeoConvert -d -p -1\n"
    "              => 10d20'30\"N 020d30'40\"E\n"
    "\n"
    "           Quoting of command line arguments is similar\n"
    "\n"
    "              GeoConvert -d -p -1 --input-string \"10d20'30\\\" 20d30'40\"\n"
    "              => 10d20'30\"N 020d30'40\"E\n"
    "\n"
    "       Windows command shell (cmd)\n"
    "           The ' character needs no quoting; the \" character can either be\n"
    "           quoted by a ^ or can be represented by typing ' twice.  (This\n"
    "           quoting is usually unnecessary because the trailing designator can\n"
    "           be omitted.)  Thus\n"
    "\n"
    "              echo 10d20'30'' 20d30'40 | GeoConvert -d -p -1\n"
    "              => 10d20'30\"N 020d30'40\"E\n"
    "\n"
    "           Use \\ to quote the \" character in a command line argument\n"
    "\n"
    "              GeoConvert -d -p -1 --input-string \"10d20'30\\\" 20d30'40\"\n"
    "              => 10d20'30\"N 020d30'40\"E\n"
    "\n"
    "       Input from a file\n"
    "           No quoting need be done if the input from a file.  Thus each line\n"
    "           of the file \"input.txt\" should just contain the plain coordinates.\n"
    "\n"
    "             GeoConvert -d -p -1 < input.txt\n"
    "\n"
    "MGRS\n"
    "       MGRS coordinates represent a square patch of the earth, thus\n"
    "       \"38SMB4488\" is in zone \"38n\" with 444km <= easting < 445km and 3688km\n"
    "       <= northing < 3689km.  Consistent with this representation, coordinates\n"
    "       are truncated (instead of rounded) to the requested precision.  When an\n"
    "       MGRS coordinate is provided as input, GeoConvert treats this as a\n"
    "       representative point within the square.  By default, this\n"
    "       representative point is the center of the square (\"38n 444500 3688500\"\n"
    "       in the example above).  (This leads to a stable conversion between MGRS\n"
    "       and geographic coordinates.)  However, if the -n option is given then\n"
    "       the south-west corner of the square is returned instead (\"38n 444000\n"
    "       3688000\" in the example above).\n"
    "\n"
    "ZONE\n"
    "       If the input is geographic, GeoConvert uses the standard rules of\n"
    "       selecting UTM vs UPS and for assigning the UTM zone (with the Norway\n"
    "       and Svalbard exceptions).  If the input is UTM/UPS or MGRS, then the\n"
    "       choice between UTM and UPS and the UTM zone mirrors the input.  The -z\n"
    "       zone, -s, and -t options allow these rules to be overridden with zone =\n"
    "       0 being used to indicate UPS.  For example, the point\n"
    "\n"
    "          79.9S 6.1E\n"
    "\n"
    "       corresponds to possible MGRS coordinates\n"
    "\n"
    "          32CMS4324728161 (standard UTM zone = 32)\n"
    "          31CEM6066227959 (neighboring UTM zone = 31)\n"
    "            BBZ1945517770 (neighboring UPS zone)\n"
    "\n"
    "       then\n"
    "\n"
    "          echo 79.9S 6.1E      | GeoConvert -p -3 -m       => 32CMS4328\n"
    "          echo 31CEM6066227959 | GeoConvert -p -3 -m       => 31CEM6027\n"
    "          echo 31CEM6066227959 | GeoConvert -p -3 -m -s    => 32CMS4328\n"
    "          echo 31CEM6066227959 | GeoConvert -p -3 -m -z 0  =>   BBZ1917\n"
    "\n"
    "       Is zone is specified with a hemisphere, then this is honored when\n"
    "       printing UTM coordinates:\n"
    "\n"
    "          echo -1 3 | GeoConvert -u         => 31s 500000 9889470\n"
    "          echo -1 3 | GeoConvert -u -z 31   => 31s 500000 9889470\n"
    "          echo -1 3 | GeoConvert -u -z 31s  => 31s 500000 9889470\n"
    "          echo -1 3 | GeoConvert -u -z 31n  => 31n 500000 -110530\n"
    "\n"
    "       NOTE: the letter in the zone specification for UTM is a hemisphere\n"
    "       designator n or s and not an MGRS latitude band letter.  Convert the\n"
    "       MGRS latitude band letter to a hemisphere as follows: replace C thru M\n"
    "       by s (or south); replace N thru X by n (or north).\n"
    "\n"
    "EXAMPLES\n"
    "          echo 38SMB4488 | GeoConvert         => 33.33424 44.40363\n"
    "          echo 38SMB4488 | GeoConvert -: -p 1 => 33:20:03.25N 044:2413.06E\n"
    "          echo 38SMB4488 | GeoConvert -u      => 38n 444500 3688500\n"
    "          echo E44d24 N33d20 | GeoConvert -m -p -3 => 38SMB4488\n"
    "\n"
    "       GeoConvert can be used to do simple arithmetic using degree, minutes,\n"
    "       and seconds.  For example, sometimes data is tiled in 15 second squares\n"
    "       tagged by the DMS representation of the SW corner.  The tags of the\n"
    "       tile at 38:59:45N 077:02:00W and its 8 neighbors are then given by\n"
    "\n"
    "           t=0:0:15\n"
    "           for y in -$t +0 +$t; do\n"
    "               for x in -$t +0 +$t; do\n"
    "                   echo 38:59:45N$y 077:02:00W$x\n"
    "               done\n"
    "           done | GeoConvert -: -p -1 | tr -d ': '\n"
    "           =>\n"
    "           385930N0770215W\n"
    "           385930N0770200W\n"
    "           385930N0770145W\n"
    "           385945N0770215W\n"
    "           385945N0770200W\n"
    "           385945N0770145W\n"
    "           390000N0770215W\n"
    "           390000N0770200W\n"
    "           390000N0770145W\n"
    "\n"
    "ERRORS\n"
    "       An illegal line of input will print an error message to standard output\n"
    "       beginning with \"ERROR:\" and causes GeoConvert to return an exit code of\n"
    "       1.  However, an error does not cause GeoConvert to terminate; following\n"
    "       lines will be converted.\n"
    "\n"
    "ABBREVIATIONS\n"
    "       UTM Universal Transverse Mercator,\n"
    "           <https://en.wikipedia.org/wiki/Universal_Transverse_Mercator_coordinate_system>.\n"
    "\n"
    "       UPS Universal Polar Stereographic,\n"
    "           <https://en.wikipedia.org/wiki/Universal_Polar_Stereographic>.\n"
    "\n"
    "       MGRS\n"
    "           Military Grid Reference System,\n"
    "           <https://en.wikipedia.org/wiki/Military_grid_reference_system>.\n"
    "\n"
    "       WGS84\n"
    "           World Geodetic System 1984, <https://en.wikipedia.org/wiki/WGS84>.\n"
    "\n"
    "SEE ALSO\n"
    "       An online version of this utility is availbable at\n"
    "       <https://geographiclib.sourceforge.io/cgi-bin/GeoConvert>.\n"
    "\n"
    "       The algorithms for the transverse Mercator projection are described in\n"
    "       C. F. F. Karney, Transverse Mercator with an accuracy of a few\n"
    "       nanometers, J. Geodesy 85(8), 475-485 (Aug. 2011); DOI\n"
    "       <https://doi.org/10.1007/s00190-011-0445-3>; preprint\n"
    "       <https://arxiv.org/abs/1002.1417>.\n"
    "\n"
    "AUTHOR\n"
    "       GeoConvert was written by Charles Karney.\n"
    "\n"
    "HISTORY\n"
    "       GeoConvert was added to GeographicLib,\n"
    "       <https://geographiclib.sourceforge.io>, in 2009-01.\n"
    ;
  return retval;
}
