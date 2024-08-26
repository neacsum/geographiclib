#include <utpp/utpp.h>

#include <GeographicLib/MGRS.hpp>
#include <GeographicLib/UTMUPS.hpp>
#include <GeographicLib/DMS.hpp>

using namespace std;

SUITE (GeoConvert)
{

using namespace GeographicLib;

TEST (test0)
{
  double lat = 33.3, lon = 44.4; // Baghdad
  int zone;
  bool northp;
  double x, y;
  UTMUPS::Forward (lat, lon, zone, northp, x, y);
  string mgrs;
  MGRS::Forward (zone, northp, x, y, lat, 2, mgrs);
  CHECK_EQUAL ("38SMB4484", mgrs);
}

TEST (test1)
{
  int zone, prec;
  bool northp;
  double x, y;
  double lat, lon;
  MGRS::Reverse ("38SMB", zone, northp, x, y, prec);
  UTMUPS::Reverse (zone, northp, x, y, lat, lon);
  auto lat_str = DMS::Encode (lat, 5);
  CHECK_EQUAL ("32d59'14.1\"", lat_str);
  auto lon_str = DMS::Encode (lon, 5);
  CHECK_EQUAL ("44d27'53.4\"", lon_str);
}

TEST (test8)
{
  int zone;
  bool northp;
  double x, y;
  UTMUPS::Forward (86, 0, zone, northp, x, y);
  CHECK_CLOSE (2000000.0, x, 1e-6);
  CHECK_CLOSE (1555731.570643, y, 1e-6);
}

TEST (test9_13)
{
  DMS::flag ind;
  CHECK_THROW (DMS::Decode ("5d70.0 10", ind), GeographicErr);
  CHECK_THROW (DMS::Decode ("5d60 10", ind), GeographicErr);

  CHECK_CLOSE (DMS::Decode ("5d59", ind), 5.98333, 1e-5);
  CHECK_EQUAL (DMS::Decode ("5d60.", ind), 6);
  CHECK_EQUAL (DMS::Decode ("5d60.0", ind), 6);
}




} //end suite GeoConvert