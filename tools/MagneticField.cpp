/**
 * \file MagneticField.cpp
 * \brief Command line utility for evaluating magnetic fields
 *
 * Copyright (c) Charles Karney (2011-2022) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 *
 * See the <a href="MagneticField.1.html">man page</a> for usage information.
 **********************************************************************/

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <GeographicLib/MagneticModel.hpp>
#include <GeographicLib/MagneticCircle.hpp>
#include <GeographicLib/DMS.hpp>
#include <GeographicLib/Utility.hpp>

int usage (int retval, bool brief);

int main(int argc, const char* const argv[]) {
  try {
    using namespace GeographicLib;
    Utility::set_digits();
    bool verbose = false, longfirst = false;
    std::string dir;
    std::string model = MagneticModel::DefaultMagneticName();
    std::string istring, ifile, ofile, cdelim;
    char lsep = ';';
    real time = 0, lat = 0, h = 0;
    bool timeset = false, circle = false, rate = false;
    real hguard = 500000, tguard = 50;
    int prec = 1, Nmax = -1, Mmax = -1;

    for (int m = 1; m < argc; ++m) {
      std::string arg(argv[m]);
      if (arg == "-n") {
        if (++m == argc) return usage(1, true);
        model = argv[m];
      } else if (arg == "-d") {
        if (++m == argc) return usage(1, true);
        dir = argv[m];
      } else if (arg == "-N") {
        if (++m == argc) return usage(1, true);
        try {
          Nmax = Utility::val<int>(std::string(argv[m]));
          if (Nmax < 0) {
            std::cerr << "Maximum degree " << argv[m] << " is negative\n";
            return 1;
          }
        }
        catch (const std::exception&) {
          std::cerr << "Precision " << argv[m] << " is not a number\n";
          return 1;
        }
      } else if (arg == "-M") {
        if (++m == argc) return usage(1, true);
        try {
          Mmax = Utility::val<int>(std::string(argv[m]));
          if (Mmax < 0) {
            std::cerr << "Maximum order " << argv[m] << " is negative\n";
            return 1;
          }
        }
        catch (const std::exception&) {
          std::cerr << "Precision " << argv[m] << " is not a number\n";
          return 1;
        }
      } else if (arg == "-t") {
        if (++m == argc) return usage(1, true);
        try {
          time = Utility::fractionalyear<real>(std::string(argv[m]));
          timeset = true;
          circle = false;
        }
        catch (const std::exception& e) {
          std::cerr << "Error decoding argument of " << arg << ": "
                    << e.what() << "\n";
          return 1;
        }
      } else if (arg == "-c") {
        if (m + 3 >= argc) return usage(1, true);
        try {
          using std::fabs;
          time = Utility::fractionalyear<real>(std::string(argv[++m]));
          DMS::flag ind;
          lat = DMS::Decode(std::string(argv[++m]), ind);
          if (ind == DMS::LONGITUDE)
            throw GeographicErr("Bad hemisphere letter on latitude");
          if (!(fabs(lat) <= 90))
            throw GeographicErr("Latitude not in [-90d, +90d]");
          h = Utility::val<real>(std::string(argv[++m]));
          timeset = false;
          circle = true;
        }
        catch (const std::exception& e) {
          std::cerr << "Error decoding argument of " << arg << ": "
                    << e.what() << "\n";
          return 1;
        }
      } else if (arg == "-r")
        rate = !rate;
      else if (arg == "-w")
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
      } else if (arg == "-T") {
        if (++m == argc) return usage(1, true);
        try {
          tguard = Utility::val<real>(std::string(argv[m]));
        }
        catch (const std::exception& e) {
          std::cerr << "Error decoding argument of " << arg << ": "
                    << e.what() << "\n";
          return 1;
        }
      } else if (arg == "-H") {
        if (++m == argc) return usage(1, true);
        try {
          hguard = Utility::val<real>(std::string(argv[m]));
        }
        catch (const std::exception& e) {
          std::cerr << "Error decoding argument of " << arg << ": "
                    << e.what() << "\n";
          return 1;
        }
      } else if (arg == "-v")
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
          std::cout<< "\nDefault magnetic path = \""
                   << MagneticModel::DefaultMagneticPath()
                   << "\"\nDefault magnetic name = \""
                   << MagneticModel::DefaultMagneticName()
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

    using std::fmax;
    tguard = fmax(real(0), tguard);
    hguard = fmax(real(0), hguard);
    prec = std::min(10 + Math::extra_digits(), std::max(0, prec));
    int retval = 0;
    try {
      using std::isfinite;
      const MagneticModel m(model, dir, Geocentric::WGS84(), Nmax, Mmax);
      if ((timeset || circle)
          && (!isfinite(time) ||
              time < m.MinTime() - tguard ||
              time > m.MaxTime() + tguard))
        throw GeographicErr("Time " + Utility::str(time) +
                            " too far outside allowed range [" +
                            Utility::str(m.MinTime()) + "," +
                            Utility::str(m.MaxTime()) + "]");
      if (circle
          && (!isfinite(h) ||
              h < m.MinHeight() - hguard ||
              h > m.MaxHeight() + hguard))
        throw GeographicErr("Height " + Utility::str(h/1000) +
                            "km too far outside allowed range [" +
                            Utility::str(m.MinHeight()/1000) + "km," +
                            Utility::str(m.MaxHeight()/1000) + "km]");
      if (verbose) {
        std::cerr << "Magnetic file: " << m.MagneticFile()      << "\n"
                  << "Name: "          << m.MagneticModelName() << "\n"
                  << "Description: "   << m.Description()       << "\n"
                  << "Date & Time: "   << m.DateTime()          << "\n"
                  << "Time range: ["
                  << m.MinTime() << ","
                  << m.MaxTime() << "]\n"
                  << "Height range: ["
                  << m.MinHeight()/1000 << "km,"
                  << m.MaxHeight()/1000 << "km]\n";
      }
      if ((timeset || circle) && (time < m.MinTime() || time > m.MaxTime()))
        std::cerr << "WARNING: Time " << time
                  << " outside allowed range ["
                  << m.MinTime() << "," << m.MaxTime() << "]\n";
      if (circle && (h < m.MinHeight() || h > m.MaxHeight()))
        std::cerr << "WARNING: Height " << h/1000
                  << "km outside allowed range ["
                  << m.MinHeight()/1000 << "km,"
                  << m.MaxHeight()/1000 << "km]\n";
      const MagneticCircle c(circle ? m.Circle(time, lat, h) :
                             MagneticCircle());
      std::string s, eol, stra, strb;
      std::istringstream str;
      while (std::getline(*input, s)) {
        try {
          eol = "\n";
          if (!cdelim.empty()) {
            std::string::size_type n = s.find(cdelim);
            if (n != std::string::npos) {
              eol = " " + s.substr(n) + "\n";
              s = s.substr(0, n);
            }
          }
          str.clear(); str.str(s);
          if (!(timeset || circle)) {
            if (!(str >> stra))
              throw GeographicErr("Incomplete input: " + s);
            time = Utility::fractionalyear<real>(stra);
            if (time < m.MinTime() - tguard || time > m.MaxTime() + tguard)
              throw GeographicErr("Time " + Utility::str(time) +
                                  " too far outside allowed range [" +
                                  Utility::str(m.MinTime()) + "," +
                                  Utility::str(m.MaxTime()) +
                                  "]");
            if (time < m.MinTime() || time > m.MaxTime())
              std::cerr << "WARNING: Time " << time
                        << " outside allowed range ["
                        << m.MinTime() << "," << m.MaxTime() << "]\n";
          }
          real lon;
          if (circle) {
            if (!(str >> strb))
              throw GeographicErr("Incomplete input: " + s);
            DMS::flag ind;
            lon = DMS::Decode(strb, ind);
            if (ind == DMS::LATITUDE)
              throw GeographicErr("Bad hemisphere letter on " + strb);
          } else {
            if (!(str >> stra >> strb))
              throw GeographicErr("Incomplete input: " + s);
            DMS::DecodeLatLon(stra, strb, lat, lon, longfirst);
            h = 0;              // h is optional
            if (str >> h) {
              if (h < m.MinHeight() - hguard || h > m.MaxHeight() + hguard)
                throw GeographicErr("Height " + Utility::str(h/1000) +
                                    "km too far outside allowed range [" +
                                    Utility::str(m.MinHeight()/1000) + "km," +
                                    Utility::str(m.MaxHeight()/1000) + "km]");
              if (h < m.MinHeight() || h > m.MaxHeight())
                std::cerr << "WARNING: Height " << h/1000
                          << "km outside allowed range ["
                          << m.MinHeight()/1000 << "km,"
                          << m.MaxHeight()/1000 << "km]\n";
            }
            else
              str.clear();
          }
          if (str >> stra)
            throw GeographicErr("Extra junk in input: " + s);
          real bx, by, bz, bxt, byt, bzt;
          if (circle)
            c(lon, bx, by, bz, bxt, byt, bzt);
          else
            m(time, lat, lon, h, bx, by, bz, bxt, byt, bzt);
          real H, F, D, I, Ht, Ft, Dt, It;
          MagneticModel::FieldComponents(bx, by, bz, bxt, byt, bzt,
                                         H, F, D, I, Ht, Ft, Dt, It);

          *output << DMS::Encode(D, prec + 1, DMS::NUMBER) << " "
                  << DMS::Encode(I, prec + 1, DMS::NUMBER) << " "
                  << Utility::str(H, prec) << " "
                  << Utility::str(by, prec) << " "
                  << Utility::str(bx, prec) << " "
                  << Utility::str(-bz, prec) << " "
                  << Utility::str(F, prec) << eol;
          if (rate)
            *output << DMS::Encode(Dt, prec + 1, DMS::NUMBER) << " "
                    << DMS::Encode(It, prec + 1, DMS::NUMBER) << " "
                    << Utility::str(Ht, prec) << " "
                    << Utility::str(byt, prec) << " "
                    << Utility::str(bxt, prec) << " "
                    << Utility::str(-bzt, prec) << " "
                    << Utility::str(Ft, prec) << eol;
        }
        catch (const std::exception& e) {
          *output << "ERROR: " << e.what() << "\n";
          retval = 1;
        }
      }
    }
    catch (const std::exception& e) {
      std::cerr << "Error reading " << model << ": " << e.what() << "\n";
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
    "    MagneticField [ -n name ] [ -d dir ] [ -N Nmax ] [ -M Mmax ] [ -t time\n"
    "    | -c time lat h ] [ -r ] [ -w ] [ -T tguard ] [ -H hguard ] [ -p prec ]\n"
    "    [ -v ] [ --comment-delimiter commentdelim ] [ --version | -h | --help ]\n"
    "    [ --input-file infile | --input-string instring ] [ --line-separator\n"
    "    linesep ] [ --output-file outfile ]\n"
    "\n"
    "For full documentation type:\n"
    "    MagneticField --help\n"
    "or visit:\n"
    "    https://geographiclib.sourceforge.io/C++/2.4/MagneticField.1.html\n";
  else
    (retval ? std::cerr : std::cout) << "Man page:\n"
    "NAME\n"
    "       MagneticField -- compute the earth's magnetic field\n"
    "\n"
    "SYNOPSIS\n"
    "       MagneticField [ -n name ] [ -d dir ] [ -N Nmax ] [ -M Mmax ] [ -t time\n"
    "       | -c time lat h ] [ -r ] [ -w ] [ -T tguard ] [ -H hguard ] [ -p prec ]\n"
    "       [ -v ] [ --comment-delimiter commentdelim ] [ --version | -h | --help ]\n"
    "       [ --input-file infile | --input-string instring ] [ --line-separator\n"
    "       linesep ] [ --output-file outfile ]\n"
    "\n"
    "DESCRIPTION\n"
    "       MagneticField reads in times and positions on standard input and prints\n"
    "       out the geomagnetic field on standard output and, optionally, its rate\n"
    "       of change.\n"
    "\n"
    "       The input line is of the form time lat lon h. time is a date of the\n"
    "       form 2012-11-25 (yyyy-mm-dd or yyyy-mm), a fractional year such as\n"
    "       2012.9, or the string \"now\".  lat and lon are the latitude and\n"
    "       longitude expressed as decimal degrees or degrees, minutes, and\n"
    "       seconds; for details on the allowed formats for latitude and longitude,\n"
    "       see the \"GEOGRAPHIC COORDINATES\" section of GeoConvert(1).  h is the\n"
    "       height above the ellipsoid in meters; this is optional and defaults to\n"
    "       zero.  Alternatively, time can be given on the command line as the\n"
    "       argument to the -t option, in which case it should not be included on\n"
    "       the input lines.  Finally, the magnetic field can be computed at\n"
    "       various points on a circle of latitude (constant time, lat, and h) via\n"
    "       the -c option; in this case only the longitude should be given on the\n"
    "       input lines.\n"
    "\n"
    "       The output consists of the following 7 items:\n"
    "\n"
    "         the declination (the direction of the horizontal component of\n"
    "           the magnetic field measured clockwise from north) in degrees,\n"
    "         the inclination (the direction of the magnetic field measured\n"
    "           down from the horizontal) in degrees,\n"
    "         the horizontal field in nanotesla (nT),\n"
    "         the north component of the field in nT,\n"
    "         the east component of the field in nT,\n"
    "         the vertical component of the field in nT (down is positive),\n"
    "         the total field in nT.\n"
    "\n"
    "       If the -r option is given, a second line is printed giving the rates of\n"
    "       change of these quantities in degrees/yr and nT/yr.\n"
    "\n"
    "       The WGS84 ellipsoid is used, a = 6378137 m, f = 1/298.257223563.\n"
    "\n"
    "OPTIONS\n"
    "       -n name\n"
    "           use magnetic field model name instead of the default \"wmm2020\".\n"
    "           See \"MODELS\".\n"
    "\n"
    "       -d dir\n"
    "           read magnetic models from dir instead of the default.  See\n"
    "           \"MODELS\".\n"
    "\n"
    "       -N Nmax\n"
    "           limit the degree of the model to Nmax.\n"
    "\n"
    "       -M Mmax\n"
    "           limit the order of the model to Mmax.\n"
    "\n"
    "       -t time\n"
    "           evaluate the field at time instead of reading the time from the\n"
    "           input lines.\n"
    "\n"
    "       -c time lat h\n"
    "           evaluate the field on a circle of latitude given by time, lat, h\n"
    "           instead of reading these quantities from the input lines.  In this\n"
    "           case, MagneticField can calculate the field considerably more\n"
    "           quickly.\n"
    "\n"
    "       -r  toggle whether to report the rates of change of the field.\n"
    "\n"
    "       -w  toggle the longitude first flag (it starts off); if the flag is on,\n"
    "           then on input and output, longitude precedes latitude (except that,\n"
    "           on input, this can be overridden by a hemisphere designator, N, S,\n"
    "           E, W).\n"
    "\n"
    "       -T tguard\n"
    "           signal an error if time lies tguard years (default 50 yr) beyond\n"
    "           the range for the model.\n"
    "\n"
    "       -H hguard\n"
    "           signal an error if h lies hguard meters (default 500000 m) beyond\n"
    "           the range for the model.\n"
    "\n"
    "       -p prec\n"
    "           set the output precision to prec (default 1).  Fields are printed\n"
    "           with precision with prec decimal places; angles use prec + 1\n"
    "           places.\n"
    "\n"
    "       -v  print information about the magnetic model on standard error before\n"
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
    "       -h  print usage, the default magnetic path and name, and exit.\n"
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
    "MODELS\n"
    "       MagneticField computes the geomagnetic field using one of the following\n"
    "       models\n"
    "\n"
    "           wmm2010, the World Magnetic Model 2010, which approximates the\n"
    "             main magnetic field for the period 2010-2015.  See\n"
    "             https://ngdc.noaa.gov/geomag/WMM/DoDWMM.shtml\n"
    "           wmm2015v2, the World Magnetic Model 2015, which approximates the\n"
    "             main magnetic field for the period 2015-2020.  See\n"
    "             https://ngdc.noaa.gov/geomag/WMM/DoDWMM.shtml\n"
    "           wmm2015, a deprecated version of wmm2015v2\n"
    "           wmm2020, the World Magnetic Model 2020, which approximates the\n"
    "             main magnetic field for the period 2020-2025.  See\n"
    "             https://ngdc.noaa.gov/geomag/WMM/DoDWMM.shtml\n"
    "           igrf11, the International Geomagnetic Reference Field (11th\n"
    "             generation), which approximates the main magnetic field for\n"
    "             the period 1900-2015.  See\n"
    "             https://ngdc.noaa.gov/IAGA/vmod/igrf.html\n"
    "           igrf12, the International Geomagnetic Reference Field (12th\n"
    "             generation), which approximates the main magnetic field for\n"
    "             the period 1900-2020.  See\n"
    "             https://ngdc.noaa.gov/IAGA/vmod/igrf.html\n"
    "           igrf13, the International Geomagnetic Reference Field (13th\n"
    "             generation), which approximates the main magnetic field for\n"
    "             the period 1900-2025.  See\n"
    "             https://ngdc.noaa.gov/IAGA/vmod/igrf.html\n"
    "           emm2010, the Enhanced Magnetic Model 2010, which approximates\n"
    "             the main and crustal magnetic fields for the period 2010-2015.\n"
    "             See https://ngdc.noaa.gov/geomag/EMM/index.html\n"
    "           emm2015, the Enhanced Magnetic Model 2015, which approximates\n"
    "             the main and crustal magnetic fields for the period 2000-2020.\n"
    "             See https://ngdc.noaa.gov/geomag/EMM/index.html\n"
    "           emm2017, the Enhanced Magnetic Model 2017, which approximates\n"
    "             the main and crustal magnetic fields for the period 2000-2022.\n"
    "             See https://ngdc.noaa.gov/geomag/EMM/index.html\n"
    "\n"
    "       These models approximate the magnetic field due to the earth's core and\n"
    "       (in the case of emm20xx) its crust.  They neglect magnetic fields due\n"
    "       to the ionosphere, the magnetosphere, nearby magnetized materials,\n"
    "       electrical machinery, etc.\n"
    "\n"
    "       By default, the \"wmm2020\" magnetic model is used.  This may changed by\n"
    "       setting the environment variable \"GEOGRAPHICLIB_MAGNETIC_NAME\" or with\n"
    "       the -n option.\n"
    "\n"
    "       The magnetic models will be loaded from a directory specified at\n"
    "       compile time.  This may changed by setting the environment variables\n"
    "       \"GEOGRAPHICLIB_MAGNETIC_PATH\" or \"GEOGRAPHICLIB_DATA\", or with the -d\n"
    "       option.  The -h option prints the default magnetic path and name.  Use\n"
    "       the -v option to ascertain the full path name of the data file.\n"
    "\n"
    "       Instructions for downloading and installing magnetic models are\n"
    "       available at\n"
    "       <https://geographiclib.sourceforge.io/C++/doc/magnetic.html#magneticinst>.\n"
    "\n"
    "ENVIRONMENT\n"
    "       GEOGRAPHICLIB_MAGNETIC_NAME\n"
    "           Override the compile-time default magnetic name of \"wmm2020\".  The\n"
    "           -h option reports the value of GEOGRAPHICLIB_MAGNETIC_NAME, if\n"
    "           defined, otherwise it reports the compile-time value.  If the -n\n"
    "           name option is used, then name takes precedence.\n"
    "\n"
    "       GEOGRAPHICLIB_MAGNETIC_PATH\n"
    "           Override the compile-time default magnetic path.  This is typically\n"
    "           \"/usr/local/share/GeographicLib/magnetic\" on Unix-like systems and\n"
    "           \"C:/ProgramData/GeographicLib/magnetic\" on Windows systems.  The -h\n"
    "           option reports the value of GEOGRAPHICLIB_MAGNETIC_PATH, if\n"
    "           defined, otherwise it reports the compile-time value.  If the -d\n"
    "           dir option is used, then dir takes precedence.\n"
    "\n"
    "       GEOGRAPHICLIB_DATA\n"
    "           Another way of overriding the compile-time default magnetic path.\n"
    "           If it is set (and if GEOGRAPHICLIB_MAGNETIC_PATH is not set), then\n"
    "           $GEOGRAPHICLIB_DATA/magnetic is used.\n"
    "\n"
    "ERRORS\n"
    "       An illegal line of input will print an error message to standard output\n"
    "       beginning with \"ERROR:\" and causes MagneticField to return an exit code\n"
    "       of 1.  However, an error does not cause MagneticField to terminate;\n"
    "       following lines will be converted.  If time or h are outside the\n"
    "       recommended ranges for the model (but inside the ranges increase by\n"
    "       tguard and hguard), a warning is printed on standard error and the\n"
    "       field (which may be inaccurate) is returned in the normal way.\n"
    "\n"
    "EXAMPLES\n"
    "       The magnetic field from WMM2020 in Timbuktu on 2020-12-25\n"
    "\n"
    "           echo 2020-12-25 16:46:33N 3:00:34W 300 | MagneticField -r\n"
    "           => -1.47 11.98 33994.9 33983.7 -871.7 7214.7 34752.1\n"
    "              0.13 -0.02 21.9 23.9 77.9 -8.4 19.7\n"
    "\n"
    "       The first two numbers returned are the declination and inclination of\n"
    "       the field.  The second line gives the annual change.\n"
    "\n"
    "SEE ALSO\n"
    "       GeoConvert(1), geographiclib-get-magnetic(8).\n"
    "\n"
    "AUTHOR\n"
    "       MagneticField was written by Charles Karney.\n"
    "\n"
    "HISTORY\n"
    "       MagneticField was added to GeographicLib,\n"
    "       <https://geographiclib.sourceforge.io>, in version 1.15.\n"
    ;
  return retval;
}
