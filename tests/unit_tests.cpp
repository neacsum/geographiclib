#include <utpp/utpp.h>

#include <GeographicLib/Math.hpp>
#include <GeographicLib/Geodesic.hpp>
#include <GeographicLib/Intersect.hpp>
#include <GeographicLib/PolygonArea.hpp>
#include <GeographicLib/Utility.hpp>
#include <GeographicLib/DMS.hpp>
#include <limits>
#include <fstream>

using namespace GeographicLib;

TEST_MAIN (int /*argc*/, char** /*argv*/)
{
  auto ret = UnitTest::RunAllTests ();
  if (ret)
  {
    //Save results in case of failure
    std::ofstream xml ("tests_results.xml");
    UnitTest::ReporterXml rpt(xml);
    UnitTest::RunAllTests (rpt);
  }

  return ret;
}

typedef Math::real real;

SUITE (Geodesic)
{
  static const real testcases[][12] = {
    {35.60777, -139.44815, 111.098748429560326,
     -11.17491, -69.95921, 129.289270889708762,
     8935244.5604818305, 80.50729714281974, 6273170.2055303837,
     0.16606318447386067, 0.16479116945612937, 12841384694976.432},
    {55.52454, 106.05087, 22.020059880982801,
     77.03196, 197.18234, 109.112041110671519,
     4105086.1713924406, 36.892740690445894, 3828869.3344387607,
     0.80076349608092607, 0.80101006984201008, 61674961290615.615},
    {-21.97856, 142.59065, -32.44456876433189,
     41.84138, 98.56635, -41.84359951440466,
     8394328.894657671, 75.62930491011522, 6161154.5773110616,
     0.24816339233950381, 0.24930251203627892, -6637997720646.717},
    {-66.99028, 112.2363, 173.73491240878403,
     -12.70631, 285.90344, 2.512956620913668,
     11150344.2312080241, 100.278634181155759, 6289939.5670446687,
     -0.17199490274700385, -0.17722569526345708, -121287239862139.744},
    {-17.42761, 173.34268, -159.033557661192928,
     -15.84784, 5.93557, -20.787484651536988,
     16076603.1631180673, 144.640108810286253, 3732902.1583877189,
     -0.81273638700070476, -0.81299800519154474, 97825992354058.708},
    {32.84994, 48.28919, 150.492927788121982,
     -56.28556, 202.29132, 48.113449399816759,
     16727068.9438164461, 150.565799985466607, 3147838.1910180939,
     -0.87334918086923126, -0.86505036767110637, -72445258525585.010},
    {6.96833, 52.74123, 92.581585386317712,
     -7.39675, 206.17291, 90.721692165923907,
     17102477.2496958388, 154.147366239113561, 2772035.6169917581,
     -0.89991282520302447, -0.89986892177110739, -1311796973197.995},
    {-50.56724, -16.30485, -105.439679907590164,
     -33.56571, -94.97412, -47.348547835650331,
     6455670.5118668696, 58.083719495371259, 5409150.7979815838,
     0.53053508035997263, 0.52988722644436602, 41071447902810.047},
    {-58.93002, -8.90775, 140.965397902500679,
     -8.91104, 133.13503, 19.255429433416599,
     11756066.0219864627, 105.755691241406877, 6151101.2270708536,
     -0.26548622269867183, -0.27068483874510741, -86143460552774.735},
    {-68.82867, -74.28391, 93.774347763114881,
     -50.63005, -8.36685, 34.65564085411343,
     3956936.926063544, 35.572254987389284, 3708890.9544062657,
     0.81443963736383502, 0.81420859815358342, -41845309450093.787},
    {-10.62672, -32.0898, -86.426713286747751,
     5.883, -134.31681, -80.473780971034875,
     11470869.3864563009, 103.387395634504061, 6184411.6622659713,
     -0.23138683500430237, -0.23155097622286792, 4198803992123.548},
    {-21.76221, 166.90563, 29.319421206936428,
     48.72884, 213.97627, 43.508671946410168,
     9098627.3986554915, 81.963476716121964, 6299240.9166992283,
     0.13965943368590333, 0.14152969707656796, 10024709850277.476},
    {-19.79938, -174.47484, 71.167275780171533,
     -11.99349, -154.35109, 65.589099775199228,
     2319004.8601169389, 20.896611684802389, 2267960.8703918325,
     0.93427001867125849, 0.93424887135032789, -3935477535005.785},
    {-11.95887, -116.94513, 92.712619830452549,
     4.57352, 7.16501, 78.64960934409585,
     13834722.5801401374, 124.688684161089762, 5228093.177931598,
     -0.56879356755666463, -0.56918731952397221, -9919582785894.853},
    {-87.85331, 85.66836, -65.120313040242748,
     66.48646, 16.09921, -4.888658719272296,
     17286615.3147144645, 155.58592449699137, 2635887.4729110181,
     -0.90697975771398578, -0.91095608883042767, 42667211366919.534},
    {1.74708, 128.32011, -101.584843631173858,
     -11.16617, 11.87109, -86.325793296437476,
     12942901.1241347408, 116.650512484301857, 5682744.8413270572,
     -0.44857868222697644, -0.44824490340007729, 10763055294345.653},
    {-25.72959, -144.90758, -153.647468693117198,
     -57.70581, -269.17879, -48.343983158876487,
     9413446.7452453107, 84.664533838404295, 6356176.6898881281,
     0.09492245755254703, 0.09737058264766572, 74515122850712.444},
    {-41.22777, 122.32875, 14.285113402275739,
     -7.57291, 130.37946, 10.805303085187369,
     3812686.035106021, 34.34330804743883, 3588703.8812128856,
     0.82605222593217889, 0.82572158200920196, -2456961531057.857},
    {11.01307, 138.25278, 79.43682622782374,
     6.62726, 247.05981, 103.708090215522657,
     11911190.819018408, 107.341669954114577, 6070904.722786735,
     -0.29767608923657404, -0.29785143390252321, 17121631423099.696},
    {-29.47124, 95.14681, -163.779130441688382,
     -27.46601, -69.15955, -15.909335945554969,
     13487015.8381145492, 121.294026715742277, 5481428.9945736388,
     -0.51527225545373252, -0.51556587964721788, 104679964020340.318} };

  TEST (Geodesic_Inverse)
  {
    real lat1, lon1, azi1, lat2, lon2, azi2, s12, a12, m12, M12, M21, S12;
    real azi1a, azi2a, s12a, a12a, m12a, M12a, M21a, S12a;
    auto& g = Geodesic::WGS84 ();
    int ncases = sizeof (testcases) / sizeof (testcases[0]);
    for (int i = 0; i < ncases; ++i) {
      lat1 = testcases[i][0]; lon1 = testcases[i][1]; azi1 = testcases[i][2];
      lat2 = testcases[i][3]; lon2 = testcases[i][4]; azi2 = testcases[i][5];
      s12 = testcases[i][6]; a12 = testcases[i][7]; m12 = testcases[i][8];
      M12 = testcases[i][9]; M21 = testcases[i][10]; S12 = testcases[i][11];
      a12a = g.GenInverse (lat1, lon1, lat2, lon2, Geodesic::ALL,
                          s12a, azi1a, azi2a, m12a, M12a, M21a, S12a);
      CHECK_CLOSE_EX (azi1, azi1a, 1e-13, "test case %d", i);
      CHECK_CLOSE_EX (azi2, azi2a, 1e-13, "test case %d", i);
      CHECK_CLOSE_EX (s12, s12a, 1e-8, "test case %d", i);
      CHECK_CLOSE_EX (a12, a12a, 1e-13, "test case %d", i);
      CHECK_CLOSE_EX (m12, m12a, 1e-8, "test case %d", i);
      CHECK_CLOSE_EX (M12, M12a, 1e-15, "test case %d", i);
      CHECK_CLOSE_EX (M21, M21a, 1e-15, "test case %d", i);
      CHECK_CLOSE_EX (S12, S12a, 0.1, "test case %d", i);
    }
  }

  TEST (GeodesicExact_Inverse)
  {
    real lat1, lon1, azi1, lat2, lon2, azi2, s12, a12, m12, M12, M21, S12;
    real azi1a, azi2a, s12a, a12a, m12a, M12a, M21a, S12a;
    auto& g = GeodesicExact::WGS84 ();
    int ncases = sizeof (testcases) / sizeof (testcases[0]);
    for (int i = 0; i < ncases; ++i) {
      lat1 = testcases[i][0]; lon1 = testcases[i][1]; azi1 = testcases[i][2];
      lat2 = testcases[i][3]; lon2 = testcases[i][4]; azi2 = testcases[i][5];
      s12 = testcases[i][6]; a12 = testcases[i][7]; m12 = testcases[i][8];
      M12 = testcases[i][9]; M21 = testcases[i][10]; S12 = testcases[i][11];
      a12a = g.GenInverse (lat1, lon1, lat2, lon2, GeodesicExact::ALL,
                          s12a, azi1a, azi2a, m12a, M12a, M21a, S12a);
      // Allow 2x error with GeodesicExact calculations (for WGS84)
      CHECK_CLOSE_EX (azi1, azi1a, 2e-13, "test case %d", i);
      CHECK_CLOSE_EX (azi2, azi2a, 2e-13, "test case %d", i);
      CHECK_CLOSE_EX (s12, s12a, 2e-8, "test case %d", i);
      CHECK_CLOSE_EX (a12, a12a, 2e-13, "test case %d", i);
      CHECK_CLOSE_EX (m12, m12a, 2e-8, "test case %d", i);
      CHECK_CLOSE_EX (M12, M12a, 2e-15, "test case %d", i);
      CHECK_CLOSE_EX (M21, M21a, 2e-15, "test case %d", i);
      CHECK_CLOSE_EX (S12, S12a, 0.2, "test case %d", i);
    }
  }

  TEST (Geodesic_Direct)
  {
    real lat1, lon1, azi1, lat2, lon2, azi2, s12, a12, m12, M12, M21, S12;
    real lat2a, lon2a, azi2a, s12a, a12a, m12a, M12a, M21a, S12a;
    auto& g = Geodesic::WGS84 ();
    int ncases = sizeof (testcases) / sizeof (testcases[0]);
    for (int i = 0; i < ncases; ++i) {
      lat1 = testcases[i][0]; lon1 = testcases[i][1]; azi1 = testcases[i][2];
      lat2 = testcases[i][3]; lon2 = testcases[i][4]; azi2 = testcases[i][5];
      s12 = testcases[i][6]; a12 = testcases[i][7]; m12 = testcases[i][8];
      M12 = testcases[i][9]; M21 = testcases[i][10]; S12 = testcases[i][11];
      a12a = g.GenDirect (lat1, lon1, azi1, false, s12, Geodesic::ALL | Geodesic::LONG_UNROLL,
                lat2a, lon2a, azi2a, s12a, m12a, M12a, M21a, S12a);
      CHECK_CLOSE_EX (lat2, lat2a, 1e-13, "test case %d", i);
      CHECK_CLOSE_EX (lon2, lon2a, 1e-13, "test case %d", i);
      CHECK_CLOSE_EX (azi2, azi2a, 1e-13, "test case %d", i);
      CHECK_EQUAL_EX (s12, s12a, "test case %d", i);
      CHECK_CLOSE_EX (a12, a12a, 1e-13, "test case %d", i);
      CHECK_CLOSE_EX (m12, m12a, 1e-8, "test case %d", i);
      CHECK_CLOSE_EX (M12, M12a, 1e-15, "test case %d", i);
      CHECK_CLOSE_EX (M21, M21a, 1e-15, "test case %d", i);
      CHECK_CLOSE_EX (S12, S12a, 0.1, "test case %d", i);
    }
  }

  TEST (GeodesicExact_Direct)
  {
    real lat1, lon1, azi1, lat2, lon2, azi2, s12, a12, m12, M12, M21, S12;
    real lat2a, lon2a, azi2a, s12a, a12a, m12a, M12a, M21a, S12a;
    auto& g = GeodesicExact::WGS84 ();
    int ncases = sizeof (testcases) / sizeof (testcases[0]);
    for (int i = 0; i < ncases; ++i) {
      lat1 = testcases[i][0]; lon1 = testcases[i][1]; azi1 = testcases[i][2];
      lat2 = testcases[i][3]; lon2 = testcases[i][4]; azi2 = testcases[i][5];
      s12 = testcases[i][6]; a12 = testcases[i][7]; m12 = testcases[i][8];
      M12 = testcases[i][9]; M21 = testcases[i][10]; S12 = testcases[i][11];
      a12a = g.GenDirect (lat1, lon1, azi1, false, s12, GeodesicExact::ALL | GeodesicExact::LONG_UNROLL,
                lat2a, lon2a, azi2a, s12a, m12a, M12a, M21a, S12a);

      // Allow 2x error with GeodesicExact calculations (for WGS84)
      CHECK_CLOSE_EX (lat2, lat2a, 2e-13, "test case %d", i);
      CHECK_CLOSE_EX (lon2, lon2a, 2e-13, "test case %d", i);
      CHECK_CLOSE_EX (azi2, azi2a, 2e-13, "test case %d", i);
      CHECK_EQUAL_EX (s12, s12a, "test case %d", i);
      CHECK_CLOSE_EX (a12, a12a, 2e-13, "test case %d", i);
      CHECK_CLOSE_EX (m12, m12a, 2e-8, "test case %d", i);
      CHECK_CLOSE_EX (M12, M12a, 2e-15, "test case %d", i);
      CHECK_CLOSE_EX (M21, M21a, 2e-15, "test case %d", i);
      CHECK_CLOSE_EX (S12, S12a, 0.2, "test case %d", i);
    }
  }

  TEST (azi1)
  {
    // azimuth of geodesic line with points on equator determined by signs of
    // latitude
    // lat1 lat2 azi1/2
    real C[2][3] = {
      { +real (0), -real (0), 180 },
      { -real (0), +real (0),   0 }
    };
    const Geodesic& g = Geodesic::WGS84 ();
    const GeodesicExact& ge = GeodesicExact::WGS84 ();
    real azi1, azi2;
    for (int k = 0; k < 2; ++k) {
      g.Inverse (C[k][0], real (0), C[k][1], real (0), azi1, azi2);
      CHECK_EQUAL_EX (azi1, C[k][2], "test case %d", k);
      CHECK_EQUAL_EX (azi2, C[k][2], "test case %d", k);
      ge.Inverse (C[k][0], real (0), C[k][1], real (0), azi1, azi2);
      CHECK_EQUAL_EX (azi1, C[k][2], "test case %d", k);
      CHECK_EQUAL_EX (azi2, C[k][2], "test case %d", k);
    }
  }

  TEST (azi2)
  {
    // Does the nearly antipodal equatorial solution go north or south?
    // lat1 lat2 azi1 azi2
    real C[2][4] = {
      { +real (0), +real (0),  56, 124},
      { -real (0), -real (0), 124,  56}
    };
    const Geodesic& g = Geodesic::WGS84 ();
    const GeodesicExact& ge = GeodesicExact::WGS84 ();
    real azi1, azi2;
    for (int k = 0; k < 2; ++k) {
      g.Inverse (C[k][0], real (0), C[k][1], real (179.5), azi1, azi2);
      CHECK_CLOSE_EX (azi1, C[k][2], 1, "test case %d", k);
      CHECK_CLOSE_EX (azi2, C[k][3], 1, "test case %d", k);
      ge.Inverse (C[k][0], real (0), C[k][1], real (179.5), azi1, azi2);
      CHECK_CLOSE_EX (azi1, C[k][2], 1, "test case %d", k);
      CHECK_CLOSE_EX (azi2, C[k][3], 1, "test case %d", k);
    }
  }

  TEST (azi3)
  {
    // How does the exact antipodal equatorial path go N/S + E/W
    // lat1 lat2 lon2 azi1 azi2
    real C[4][5] = {
      { +real (0), +real (0), +180,   +real (0), +180},
      { -real (0), -real (0), +180, +180,   +real (0)},
      { +real (0), +real (0), -180,   -real (0), -180},
      { -real (0), -real (0), -180, -180,   -real (0)}
    };
    const Geodesic& g = Geodesic::WGS84 ();
    const GeodesicExact& ge = GeodesicExact::WGS84 ();
    real azi1, azi2;
    for (int k = 0; k < 4; ++k) {
      g.Inverse (C[k][0], real (0), C[k][1], C[k][2], azi1, azi2);
      CHECK_EQUAL_EX (azi1, C[k][3], "test case %d", k);
      CHECK_EQUAL_EX (azi2, C[k][4], "test case %d", k);
      ge.Inverse (C[k][0], real (0), C[k][1], C[k][2], azi1, azi2);
      CHECK_EQUAL_EX (azi1, C[k][3], "test case %d", k);
      CHECK_EQUAL_EX (azi2, C[k][4], "test case %d", k);
    }
  }

  TEST (azi4)
  {
    // Antipodal points on the equator with prolate ellipsoid
    // lon2 azi1/2
    real C[2][2] = {
      { +180, +90 },
      { -180, -90 }
    };
    const Geodesic g (real (6.4e6), -1 / real (300));
    const GeodesicExact ge (real (6.4e6), -1 / real (300));
    real azi1, azi2;
    for (int k = 0; k < 2; ++k) {
      g.Inverse (real (0), real (0), real (0), C[k][0], azi1, azi2);
      CHECK_EQUAL_EX (azi1, C[k][1], "test case %d", k);
      CHECK_EQUAL_EX (azi2, C[k][1], "test case %d", k);
      ge.Inverse (real (0), real (0), real (0), C[k][0], azi1, azi2);
      CHECK_EQUAL_EX (azi1, C[k][1], "test case %d", k);
      CHECK_EQUAL_EX (azi2, C[k][1], "test case %d", k);
    }
  }

  TEST (azi5)
  {
    // azimuths = +/-0 and +/-180 for the direct problem
    // azi1, lon2, azi2
    real C[4][3] = {
      { +real (0), +180, +180  },
      { -real (0), -180, -180  },
      { +180 , +180, +real (0) },
      { -180 , -180, -real (0) }
    };
    const Geodesic& g = Geodesic::WGS84 ();
    const GeodesicExact& ge = GeodesicExact::WGS84 ();
    real lon2, azi2;
    for (int k = 0; k < 4; ++k) {
      real t;
      g.GenDirect (real (0), real (0), C[k][0], false, real (15e6),
                  Geodesic::LONGITUDE | Geodesic::AZIMUTH |
                  Geodesic::LONG_UNROLL,
                  t, lon2, azi2,
                  t, t, t, t, t);
      CHECK_EQUAL_EX (lon2, C[k][1], "test case %d", k);
      CHECK_EQUAL_EX (azi2, C[k][2], "test case %d", k);
      ge.GenDirect (real (0), real (0), C[k][0], false, real (15e6),
                   Geodesic::LONGITUDE | Geodesic::AZIMUTH |
                   Geodesic::LONG_UNROLL,
                   t, lon2, azi2,
                   t, t, t, t, t);
      CHECK_EQUAL_EX (lon2, C[k][1], "test case %d", k);
      CHECK_EQUAL_EX (azi2, C[k][2], "test case %d", k);
    }
  }
}

void equatorialseg (Intersect& inter, 
  real lonx1, real lonx2, 
  real lony1, real lony2, 
  real xval, real yval, 
  int test)
{
  const real epsilon = 1e-6;
  int segmode = 0;
  auto p1 = inter.Segment (0, lonx1, 0, lonx2,
                          0, lony1, 0, lony2, segmode),
    p2 = inter.Segment (0, lony1, 0, lony2,
                       0, lonx1, 0, lonx2, segmode);
  CHECK_CLOSE_EX (p1.first, xval, epsilon, "test %d", test);
  CHECK_CLOSE_EX (p1.second, yval, epsilon, "test %d", test);
  CHECK_CLOSE_EX (p2.first, yval, epsilon, "test %d", test);
  CHECK_CLOSE_EX (p2.second, xval, epsilon, "test %d", test);
}

TEST (checkcoincident1)
{
  struct {
    real lonx1, lonx2, lony1, lony2, px, py;
  } test_data[] = {
    {0, 40, -20, -10, -5, 15},
    {0, 40, -10, 0, 0, 10},
    {0, 40, -8, 2, 1, 9},
    {0, 40, -2, 8, 4, 6},
    {0, 40, 0, 10, 5, 5},
    {0, 40, 2, 12, 7, 5},
    {0, 40, 15, 25, 20, 5},
    {0, 40, 30, 40, 35, 5},
    {0, 40, 32, 42, 36, 4},
    {0, 40, 38, 48, 39, 1},
    {0, 40, 40, 50, 40, 0},
    {0, 40, 50, 60, 45, -5},
    {40, 0, -20, -10, 40 + 5, 15},
    {40, 0, -10, -0, 40, 10},
    {40, 0, -8, 2, 40 - 1, 9},
    {40, 0, -2, 8, 40 - 4, 6},
    {40, 0, 0, 10, 40 - 5, 5},
    {40, 0, 2, 12, 40 - 7, 5},
    {40, 0, 15, 25, 40 - 20, 5},
    {40, 0, 30, 40, 40 - 35, 5},
    {40, 0, 32, 42, 40 - 36, 4},
    {40, 0, 38, 48, 40 - 39, 1},
    {40, 0, 40, 50, 40 - 40, 0},
    {40, 0, 50, 60, 40 - 45, -5}
  };
  int n = sizeof (test_data) / sizeof (test_data[0]);
  real a = real (180) / Math::pi ();
  for (int exact = 0; exact < 2; ++exact)
  {
    for (int fi = -1; fi <= 1; ++fi)
    {
      for (int i = 0; i < n; i++)
      {
        Geodesic geod (a, fi / real (10), exact != 0);
        Intersect inter (geod);
        equatorialseg (inter, test_data[i].lonx1, test_data[i].lonx2,
          test_data[i].lony1, test_data[i].lony2, test_data[i].px,
          test_data[i].py, i);
      }
    }
  }
}

SUITE (Polygons)
{
  TEST (Planimeter15)
  {
    // Coverage tests, includes Planimeter15 - Planimeter18 (combinations of
    // reverse and sign) + calls to testpoint, testedge, geod_polygonarea.
    const Geodesic& g = Geodesic::WGS84 ();
    PolygonArea polygon (g);
    real lat[] = { 2, 1, 3 }, lon[] = { 1, 2, 3 };
    real perim, area, s12, azi1, azi2;
    real r = 18454562325.45119,
      a0 = 510065621724088.5093;  // ellipsoid area
    polygon.AddPoint (lat[0], lon[0]);
    polygon.AddPoint (lat[1], lon[1]);
    polygon.TestPoint (lat[2], lon[2], false, true, perim, area);
    CHECK_CLOSE (area, r, 0.5);
    polygon.TestPoint (lat[2], lon[2], false, false, perim, area);
    CHECK_CLOSE (area, r, 0.5);
    polygon.TestPoint (lat[2], lon[2], true, true, perim, area);
    CHECK_CLOSE (area, -r, 0.5);
    polygon.TestPoint (lat[2], lon[2], true, false, perim, area);
    CHECK_CLOSE (area, a0 - r, 0.5);
    g.Inverse (lat[1], lon[1], lat[2], lon[2], s12, azi1, azi2);
    polygon.TestEdge (azi1, s12, false, true, perim, area);
    CHECK_CLOSE (area, r, 0.5);
    polygon.TestEdge (azi1, s12, false, false, perim, area);
    CHECK_CLOSE (area, r, 0.5);
    polygon.TestEdge (azi1, s12, true, true, perim, area);
    CHECK_CLOSE (area, -r, 0.5);
    polygon.TestEdge (azi1, s12, true, false, perim, area);
    CHECK_CLOSE (area, a0 - r, 0.5);
    polygon.AddPoint (lat[2], lon[2]);
    polygon.Compute (false, true, perim, area);
    CHECK_CLOSE (area, r, 0.5);
    polygon.Compute (false, false, perim, area);
    CHECK_CLOSE (area, r, 0.5);
    polygon.Compute (true, true, perim, area);
    CHECK_CLOSE (area, -r, 0.5);
    polygon.Compute (true, false, perim, area);
    CHECK_CLOSE (area, a0 - r, 0.5);
  }

  TEST (Planimeter1)
  {
    // Coverage tests, includes Planimeter19 - Planimeter20 (degenerate
    // polygons) + extra cases.
    const Geodesic& g = Geodesic::WGS84 ();
    PolygonArea polygon (g, false);
    PolygonArea polyline (g, true);
    real perim, area;
    polygon.Compute (false, true, perim, area);
    CHECK_EQUAL (area, 0);
    CHECK_EQUAL (perim, 0);
    polygon.TestPoint (1, 1, false, true, perim, area);
    CHECK_EQUAL (area, 0);
    CHECK_EQUAL (perim, 0);
    polygon.TestEdge (90, 1000, false, true, perim, area);
    CHECK_NAN (area);
    CHECK_NAN (perim);
    polygon.AddPoint (1, 1);
    polygon.Compute (false, true, perim, area);
    CHECK_EQUAL (area, 0);
    CHECK_EQUAL (perim, 0);
    polyline.Compute (false, true, perim, area);
    CHECK_EQUAL (perim, 0);
    polyline.TestPoint (1, 1, false, true, perim, area);
    CHECK_EQUAL (perim, 0);
    polyline.TestEdge (90, 1000, false, true, perim, area);
    CHECK_NAN (perim);
    polyline.AddPoint (1, 1);
    polyline.Compute (false, true, perim, area);
    CHECK_EQUAL (perim, 0);
    polyline.AddPoint (1, 1);
    polyline.TestEdge (90, 1000, false, true, perim, area);
    CHECK_CLOSE (perim, 1000, 1e-10);
    polyline.TestPoint (2, 2, false, true, perim, area);
    CHECK_CLOSE (perim, 156876.149, 0.5e-3);
  }

  TEST (Planimeter21)
  {
    // Some test to add code coverage: multiple circlings of pole (includes
    // Planimeter21 - Planimeter28) + invocations via testpoint and testedge.
    const Geodesic& g = Geodesic::WGS84 ();
    PolygonArea polygon (g);
    real perim, area, lat = 45,
      a = 39.2144607176828184218, s = 8420705.40957178156285,
      r = 39433884866571.4277,    // Area for one circuit
      a0 = 510065621724088.5093;  // Ellipsoid area

    polygon.AddPoint (lat, 60);
    polygon.AddPoint (lat, 180);
    polygon.AddPoint (lat, -60);
    polygon.AddPoint (lat, 60);
    polygon.AddPoint (lat, 180);
    polygon.AddPoint (lat, -60);
    for (int i = 3; i <= 4; ++i) {
      polygon.AddPoint (lat, 60);
      polygon.AddPoint (lat, 180);
      polygon.TestPoint (lat, -60, false, true, perim, area);
      CHECK_CLOSE (area, i * r, 0.5);
      polygon.TestPoint (lat, -60, false, false, perim, area);
      CHECK_CLOSE (area, i * r, 0.5);
      polygon.TestPoint (lat, -60, true, true, perim, area);
      CHECK_CLOSE (area, -i * r, 0.5);
      polygon.TestPoint (lat, -60, true, false, perim, area);
      CHECK_CLOSE (area, -i * r + a0, 0.5);
      polygon.TestEdge (a, s, false, true, perim, area);
      CHECK_CLOSE (area, i * r, 0.5);
      polygon.TestEdge (a, s, false, false, perim, area);
      CHECK_CLOSE (area, i * r, 0.5);
      polygon.TestEdge (a, s, true, true, perim, area);
      CHECK_CLOSE (area, -i * r, 0.5);
      polygon.TestEdge (a, s, true, false, perim, area);
      CHECK_CLOSE (area, -i * r + a0, 0.5);
      polygon.AddPoint (lat, -60);
      polygon.Compute (false, true, perim, area);
      CHECK_CLOSE (area, i * r, 0.5);
      polygon.Compute (false, false, perim, area);
      CHECK_CLOSE (area, i * r, 0.5);
      polygon.Compute (true, true, perim, area);
      CHECK_CLOSE (area, -i * r, 0.5);
      polygon.Compute (true, false, perim, area);
      CHECK_CLOSE (area, -i * r + a0, 0.5);
    }
  }

  TEST (Planimeter29)
  {
    // Check fix to transitdirect vs transit zero handling inconsistency
    const Geodesic& g = Geodesic::WGS84 ();
    PolygonArea polygon (g);
    real perim, area;
    polygon.AddPoint (0, 0);
    polygon.AddEdge (90, 1000);
    polygon.AddEdge (0, 1000);
    polygon.AddEdge (-90, 1000);
    polygon.Compute (false, true, perim, area);
    // The area should be 1e6.  Prior to the fix it was 1e6 - A/2, where
    // A = ellipsoid area.
    CHECK_CLOSE (area, 1000000.0, 0.01);
  }

  TEST (AddEdge)
  {
    // Implement check on AddEdge bug: AddPoint(0,0) +
    // AddEdge(90, 1000) + AddEdge(0, 1000) + AddEdge(-90, 0).  The area
    // should be 1e6.  Prior to the fix it was 1e6 - A/2, where A = ellipsoid
    // area.
    // add_test (NAME Planimeter29 COMMAND Planimeter ...)
    const Geodesic& g = Geodesic::WGS84 ();
    PolygonArea polygon (g);
    polygon.AddPoint (0, 0);
    polygon.AddEdge (90, 1000);
    polygon.AddEdge (0, 1000);
    polygon.AddEdge (-90, 1000);
    real perim, area;
    // The area should be 1e6.  Prior to the fix it was 1e6 - A/2, where
    // A = ellipsoid area.
    polygon.Compute (false, true, perim, area);
    CHECK_CLOSE (area, 1000000, 0.01);
  }
}
SUITE (AngleTest)
{
  const real inf = Math::infinity (),
    nan = Math::NaN (),
    eps = std::numeric_limits<real>::epsilon (),
    ovf = 1 / Math::sq (eps);

  TEST (AngRound)
  {
    CHECK_EQUAL (Math::AngRound (-eps / 32), -eps / 32);
    CHECK_EQUAL (Math::AngRound (-eps / 64), -0.0);
    CHECK_EQUAL (Math::AngRound (-real (0)), -0.0);
    CHECK_EQUAL (Math::AngRound (real (0)), +0.0);
    CHECK_EQUAL (Math::AngRound (eps / 64), +0.0);
    CHECK_EQUAL (Math::AngRound (eps / 32), +eps / 32);
    CHECK_EQUAL (Math::AngRound ((1 - 2 * eps) / 64), (1 - 2 * eps) / 64);
    CHECK_EQUAL (Math::AngRound ((1 - eps) / 64), real (1) / 64);
    CHECK_EQUAL (Math::AngRound ((1 - eps / 2) / 64), real (1) / 64);
    CHECK_EQUAL (Math::AngRound ((1 - eps / 4) / 64), real (1) / 64);
    CHECK_EQUAL (Math::AngRound (real (1) / 64), real (1) / 64);
    CHECK_EQUAL (Math::AngRound ((1 + eps / 2) / 64), real (1) / 64);
    CHECK_EQUAL (Math::AngRound ((1 + eps) / 64), real (1) / 64);
    CHECK_EQUAL (Math::AngRound ((1 + 2 * eps) / 64), (1 + 2 * eps) / 64);
    CHECK_EQUAL (Math::AngRound ((1 - eps) / 32), (1 - eps) / 32);
    CHECK_EQUAL (Math::AngRound ((1 - eps / 2) / 32), real (1) / 32);
    CHECK_EQUAL (Math::AngRound ((1 - eps / 4) / 32), real (1) / 32);
    CHECK_EQUAL (Math::AngRound (real (1) / 32), real (1) / 32);
    CHECK_EQUAL (Math::AngRound ((1 + eps / 2) / 32), real (1) / 32);
    CHECK_EQUAL (Math::AngRound ((1 + eps) / 32), (1 + eps) / 32);
    CHECK_EQUAL (Math::AngRound ((1 - eps) / 16), (1 - eps) / 16);
    CHECK_EQUAL (Math::AngRound ((1 - eps / 2) / 16), (1 - eps / 2) / 16);
    CHECK_EQUAL (Math::AngRound ((1 - eps / 4) / 16), real (1) / 16);
    CHECK_EQUAL (Math::AngRound (real (1) / 16), real (1) / 16);
    CHECK_EQUAL (Math::AngRound ((1 + eps / 4) / 16), real (1) / 16);
    CHECK_EQUAL (Math::AngRound ((1 + eps / 2) / 16), real (1) / 16);
    CHECK_EQUAL (Math::AngRound ((1 + eps) / 16), (1 + eps) / 16);
    CHECK_EQUAL (Math::AngRound ((1 - eps) / 8), (1 - eps) / 8);
    CHECK_EQUAL (Math::AngRound ((1 - eps / 2) / 8), (1 - eps / 2) / 8);
    CHECK_EQUAL (Math::AngRound ((1 - eps / 4) / 8), real (1) / 8);
    CHECK_EQUAL (Math::AngRound ((1 + eps / 2) / 8), real (1) / 8);
    CHECK_EQUAL (Math::AngRound ((1 + eps) / 8), (1 + eps) / 8);
    CHECK_EQUAL (Math::AngRound (1 - eps), 1 - eps);
    CHECK_EQUAL (Math::AngRound (1 - eps / 2), 1 - eps / 2);
    CHECK_EQUAL (Math::AngRound (1 - eps / 4), 1);
    CHECK_EQUAL (Math::AngRound (real (1)), 1);
    CHECK_EQUAL (Math::AngRound (1 + eps / 4), 1);
    CHECK_EQUAL (Math::AngRound (1 + eps / 2), 1);
    CHECK_EQUAL (Math::AngRound (1 + eps), 1 + eps);
    CHECK_EQUAL (Math::AngRound (real (90) - 64 * eps), 90 - 64 * eps);
    CHECK_EQUAL (Math::AngRound (real (90) - 32 * eps), 90);
    CHECK_EQUAL (Math::AngRound (real (90)), 90);
  }

  TEST (sind)
  {
    CHECK_NAN (Math::sind (-inf));
    CHECK_EQUAL (Math::sind (-real (720)), -0.0);
    CHECK_EQUAL (Math::sind (-real (540)), -0.0);
    CHECK_EQUAL (Math::sind (-real (360)), -0.0);
    CHECK_EQUAL (Math::sind (-real (180)), -0.0);
    CHECK_EQUAL (Math::sind (-real (0)), -0.0);
    CHECK_EQUAL (Math::sind (+real (0)), +0.0);
    CHECK_EQUAL (Math::sind (+real (180)), +0.0);
    CHECK_EQUAL (Math::sind (+real (360)), +0.0);
    CHECK_EQUAL (Math::sind (+real (540)), +0.0);
    CHECK_EQUAL (Math::sind (+real (720)), +0.0);
    CHECK_NAN (Math::sind (+inf));
  }

  TEST (cosd)
  {
    CHECK_NAN (Math::cosd (-inf));
    CHECK_EQUAL (Math::cosd (-real (810)), +0.0);
    CHECK_EQUAL (Math::cosd (-real (630)), +0.0);
    CHECK_EQUAL (Math::cosd (-real (450)), +0.0);
    CHECK_EQUAL (Math::cosd (-real (270)), +0.0);
    CHECK_EQUAL (Math::cosd (-real (90)), +0.0);
    CHECK_EQUAL (Math::cosd (+real (90)), +0.0);
    CHECK_EQUAL (Math::cosd (+real (270)), +0.0);
    CHECK_EQUAL (Math::cosd (+real (450)), +0.0);
    CHECK_EQUAL (Math::cosd (+real (630)), +0.0);
    CHECK_EQUAL (Math::cosd (+real (810)), +0.0);
    CHECK_NAN (Math::cosd (+inf));
  }

  TEST (sincosd)
  {
    real sx, cx;

    Math::sincosd (-inf, sx, cx);
    CHECK_NAN (sx);
    CHECK_NAN (cx);

    Math::sincosd (+inf, sx, cx);
    CHECK_NAN (sx);
    CHECK_NAN (cx);

    Math::sincosd (nan, sx, cx);
    CHECK_NAN (sx);
    CHECK_NAN (cx);

#define checksincosd(x, s, c) do {                  \
    Math::sincosd(x, sx, cx);                       \
    CHECK_EQUAL (s, sx);                            \
    CHECK_EQUAL (c, cx);                            \
  } while (false)

    checksincosd (-real (810), -1.0, +0.0);
    checksincosd (-real (720), -0.0, +1.0);
    checksincosd (-real (630), +1.0, +0.0);
    checksincosd (-real (540), -0.0, -1.0);
    checksincosd (-real (450), -1.0, +0.0);
    checksincosd (-real (360), -0.0, +1.0);
    checksincosd (-real (270), +1.0, +0.0);
    checksincosd (-real (180), -0.0, -1.0);
    checksincosd (-real (90), -1.0, +0.0);
    checksincosd (-real (0), -0.0, +1.0);
    checksincosd (+real (0), +0.0, +1.0);
    checksincosd (+real (90), +1.0, +0.0);
    checksincosd (+real (180), +0.0, -1.0);
    checksincosd (+real (270), -1.0, +0.0);
    checksincosd (+real (360), +0.0, +1.0);
    checksincosd (+real (450), +1.0, +0.0);
    checksincosd (+real (540), +0.0, -1.0);
    checksincosd (+real (630), -1.0, +0.0);
    checksincosd (+real (720), +0.0, +1.0);
    checksincosd (+real (810), +1.0, +0.0);
  }

  TEST (tand)
  {
    CHECK_NAN (Math::tand (-inf));
    CHECK_EQUAL (Math::tand (-real (810)), -ovf);
    CHECK_EQUAL (Math::tand (-real (720)), -0.0);
    CHECK_EQUAL (Math::tand (-real (630)), +ovf);
    CHECK_EQUAL (Math::tand (-real (540)), +0.0);
    CHECK_EQUAL (Math::tand (-real (450)), -ovf);
    CHECK_EQUAL (Math::tand (-real (360)), -0.0);
    CHECK_EQUAL (Math::tand (-real (270)), +ovf);
    CHECK_EQUAL (Math::tand (-real (180)), +0.0);
    CHECK_EQUAL (Math::tand (-real (90)), -ovf);
    CHECK_EQUAL (Math::tand (-real (0)), -0.0);
    CHECK_EQUAL (Math::tand (+real (0)), +0.0);
    CHECK_EQUAL (Math::tand (+real (90)), +ovf);
    CHECK_EQUAL (Math::tand (+real (180)), -0.0);
    CHECK_EQUAL (Math::tand (+real (270)), -ovf);
    CHECK_EQUAL (Math::tand (+real (360)), +0.0);
    CHECK_EQUAL (Math::tand (+real (450)), +ovf);
    CHECK_EQUAL (Math::tand (+real (540)), -0.0);
    CHECK_EQUAL (Math::tand (+real (630)), -ovf);
    CHECK_EQUAL (Math::tand (+real (720)), +0.0);
    CHECK_EQUAL (Math::tand (+real (810)), +ovf);
    CHECK_NAN (Math::tand (+inf));
  }

  TEST (atan2d)
  {
    CHECK_EQUAL (Math::atan2d (+real (0), -real (0)), +180);
    CHECK_EQUAL (Math::atan2d (-real (0), -real (0)), -180);
    CHECK_EQUAL (Math::atan2d (+real (0), +real (0)), +0.0);
    CHECK_EQUAL (Math::atan2d (-real (0), +real (0)), -0.0);
    CHECK_EQUAL (Math::atan2d (+real (0), -real (1)), +180);
    CHECK_EQUAL (Math::atan2d (-real (0), -real (1)), -180);
    CHECK_EQUAL (Math::atan2d (+real (0), +real (1)), +0.0);
    CHECK_EQUAL (Math::atan2d (-real (0), +real (1)), -0.0);
    CHECK_EQUAL (Math::atan2d (-real (1), +real (0)), -90);
    CHECK_EQUAL (Math::atan2d (-real (1), -real (0)), -90);
    CHECK_EQUAL (Math::atan2d (+real (1), +real (0)), +90);
    CHECK_EQUAL (Math::atan2d (+real (1), -real (0)), +90);
    CHECK_EQUAL (Math::atan2d (+real (1), -inf), +180);
    CHECK_EQUAL (Math::atan2d (-real (1), -inf), -180);
    CHECK_EQUAL (Math::atan2d (+real (1), +inf), +0.0);
    CHECK_EQUAL (Math::atan2d (-real (1), +inf), -0.0);
    CHECK_EQUAL (Math::atan2d (+inf, +real (1)), +90);
    CHECK_EQUAL (Math::atan2d (+inf, -real (1)), +90);
    CHECK_EQUAL (Math::atan2d (-inf, +real (1)), -90);
    CHECK_EQUAL (Math::atan2d (-inf, -real (1)), -90);
    CHECK_EQUAL (Math::atan2d (+inf, -inf), +135);
    CHECK_EQUAL (Math::atan2d (-inf, -inf), -135);
    CHECK_EQUAL (Math::atan2d (+inf, +inf), +45);
    CHECK_EQUAL (Math::atan2d (-inf, +inf), -45);
    CHECK_NAN (Math::atan2d (nan, +real (1)));
    CHECK_NAN (Math::atan2d (+real (1), nan));

    real x = 7e-16;
    CHECK_EQUAL (Math::atan2d (x, -real (1)), 180 - Math::atan2d (x, real (1)));
  }

  TEST (sum)
  {
    real x;
    CHECK_EQUAL (Math::sum (+real (9), -real (9), x), +0.0);
    CHECK_EQUAL (Math::sum (-real (9), +real (9), x), +0.0);
    CHECK_EQUAL (Math::sum (-real (0), +real (0), x), +0.0);
    CHECK_EQUAL (Math::sum (+real (0), -real (0), x), +0.0);
    CHECK_EQUAL (Math::sum (-real (0), -real (0), x), -0.0);
    CHECK_EQUAL (Math::sum (+real (0), +real (0), x), +0.0);
  }

  TEST (AngNormalize)
  {
    CHECK_EQUAL (Math::AngNormalize (-real (900)), -180);
    CHECK_EQUAL (Math::AngNormalize (-real (720)), -0.0);
    CHECK_EQUAL (Math::AngNormalize (-real (540)), -180);
    CHECK_EQUAL (Math::AngNormalize (-real (360)), -0.0);
    CHECK_EQUAL (Math::AngNormalize (-real (180)), -180);
    CHECK_EQUAL (Math::AngNormalize (-real (0)), -0.0);
    CHECK_EQUAL (Math::AngNormalize (+real (0)), +0.0);
    CHECK_EQUAL (Math::AngNormalize (real (180)), +180);
    CHECK_EQUAL (Math::AngNormalize (real (360)), +0.0);
    CHECK_EQUAL (Math::AngNormalize (real (540)), +180);
    CHECK_EQUAL (Math::AngNormalize (real (720)), +0.0);
    CHECK_EQUAL (Math::AngNormalize (real (900)), +180);
  }

  TEST (AngDiff)
  {
    real x;
    CHECK_EQUAL (Math::AngDiff (+real (0), +real (0), x), +0.0);
    CHECK_EQUAL (Math::AngDiff (+real (0), -real (0), x), -0.0);
    CHECK_EQUAL (Math::AngDiff (-real (0), +real (0), x), +0.0);
    CHECK_EQUAL (Math::AngDiff (-real (0), -real (0), x), +0.0);
    CHECK_EQUAL (Math::AngDiff (+real (5), +real (365), x), +0.0);
    CHECK_EQUAL (Math::AngDiff (+real (365), +real (5), x), -0.0);
    CHECK_EQUAL (Math::AngDiff (+real (5), +real (185), x), +180.0);
    CHECK_EQUAL (Math::AngDiff (+real (185), +real (5), x), -180.0);
    CHECK_EQUAL (Math::AngDiff (+eps, +real (180), x), +180.0);
    CHECK_EQUAL (Math::AngDiff (-eps, +real (180), x), -180.0);
    CHECK_EQUAL (Math::AngDiff (+eps, -real (180), x), +180.0);
    CHECK_EQUAL (Math::AngDiff (-eps, -real (180), x), -180.0);

    x = 138 + 128 * eps;
    real y = -164;
    CHECK_EQUAL (Math::AngDiff (x, y), 58 - 128 * eps);
  }

  TEST (val_str)
  {
    CHECK_EQUAL (Utility::val<real> ("+0"), +0.0);
    CHECK_EQUAL (Utility::val<real> ("-0"), -0.0);
    CHECK_NAN (Utility::val<real> ("nan"));
    CHECK_EQUAL (Utility::val<real> ("+inf"), +inf);
    CHECK_EQUAL (Utility::val<real> ("inf"), +inf);
    CHECK_EQUAL (Utility::val<real> ("-inf"), -inf);

    CHECK_EQUAL (Utility::str (nan, 0), "nan");
    CHECK_EQUAL (Utility::str (-inf, 0), "-inf");
    CHECK_EQUAL (Utility::str (+inf, 0), "inf");
    CHECK_EQUAL (Utility::str (-real (3.5), 0), "-4");

    CHECK_EQUAL (Utility::str (nan, 0), "nan");
    CHECK_EQUAL (Utility::str (-inf, 0), "-inf");
    CHECK_EQUAL (Utility::str (+inf, 0), "inf");
    CHECK_EQUAL (Utility::str (-real (3.5), 0), "-4");
    CHECK_EQUAL (Utility::str (-real (2.5), 0), "-2");
    CHECK_EQUAL (Utility::str (-real (1.5), 0), "-2");
    CHECK_EQUAL (Utility::str (-real (0.5), 0), "-0");
    CHECK_EQUAL (Utility::str (-real (0), 0), "-0");
    CHECK_EQUAL (Utility::str (+real (0), 0), "0");
    CHECK_EQUAL (Utility::str (+real (0.5), 0), "0");
    CHECK_EQUAL (Utility::str (+real (1.5), 0), "2");
    CHECK_EQUAL (Utility::str (+real (2.5), 0), "2");
    CHECK_EQUAL (Utility::str (+real (3.5), 0), "4");
    CHECK_EQUAL (Utility::str (-real (1.75), 1), "-1.8");
    CHECK_EQUAL (Utility::str (-real (1.25), 1), "-1.2");
    CHECK_EQUAL (Utility::str (-real (0.75), 1), "-0.8");
    CHECK_EQUAL (Utility::str (-real (0.25), 1), "-0.2");
    CHECK_EQUAL (Utility::str (-real (0), 1), "-0.0");
    CHECK_EQUAL (Utility::str (+real (0), 1), "0.0");
    CHECK_EQUAL (Utility::str (+real (0.25), 1), "0.2");
    CHECK_EQUAL (Utility::str (+real (0.75), 1), "0.8");
    CHECK_EQUAL (Utility::str (+real (1.25), 1), "1.2");
    CHECK_EQUAL (Utility::str (+real (1.75), 1), "1.8");
  }

  TEST (DMS_Encode)
  {
    CHECK_EQUAL (DMS::Encode (nan, DMS::DEGREE, 0), "nan");
    CHECK_EQUAL (DMS::Encode (-inf, DMS::DEGREE, 0), "-inf");
    CHECK_EQUAL (DMS::Encode (+inf, DMS::DEGREE, 0), "inf");
    CHECK_EQUAL (DMS::Encode (-real (3.5), DMS::DEGREE, 0), "-4");
    CHECK_EQUAL (DMS::Encode (-real (2.5), DMS::DEGREE, 0), "-2");
    CHECK_EQUAL (DMS::Encode (-real (1.5), DMS::DEGREE, 0), "-2");
    CHECK_EQUAL (DMS::Encode (-real (0.5), DMS::DEGREE, 0), "-0");
    CHECK_EQUAL (DMS::Encode (-real (0), DMS::DEGREE, 0), "-0");
    CHECK_EQUAL (DMS::Encode (+real (0), DMS::DEGREE, 0), "0");
    CHECK_EQUAL (DMS::Encode (+real (0.5), DMS::DEGREE, 0), "0");
    CHECK_EQUAL (DMS::Encode (+real (1.5), DMS::DEGREE, 0), "2");
    CHECK_EQUAL (DMS::Encode (+real (2.5), DMS::DEGREE, 0), "2");
    CHECK_EQUAL (DMS::Encode (+real (3.5), DMS::DEGREE, 0), "4");
    CHECK_EQUAL (DMS::Encode (-real (1.75), DMS::DEGREE, 1), "-1.8");
    CHECK_EQUAL (DMS::Encode (-real (1.25), DMS::DEGREE, 1), "-1.2");
    CHECK_EQUAL (DMS::Encode (-real (0.75), DMS::DEGREE, 1), "-0.8");
    CHECK_EQUAL (DMS::Encode (-real (0.25), DMS::DEGREE, 1), "-0.2");
    CHECK_EQUAL (DMS::Encode (-real (0), DMS::DEGREE, 1), "-0.0");
    CHECK_EQUAL (DMS::Encode (+real (0), DMS::DEGREE, 1), "0.0");
    CHECK_EQUAL (DMS::Encode (+real (0.25), DMS::DEGREE, 1), "0.2");
    CHECK_EQUAL (DMS::Encode (+real (0.75), DMS::DEGREE, 1), "0.8");
    CHECK_EQUAL (DMS::Encode (+real (1.25), DMS::DEGREE, 1), "1.2");
    CHECK_EQUAL (DMS::Encode (+real (1.75), DMS::DEGREE, 1), "1.8");
    CHECK_EQUAL (DMS::Encode (real (1e20), DMS::DEGREE, 0), "100000000000000000000");
    CHECK_EQUAL (DMS::Encode (real (1e21), DMS::DEGREE, 0), "1000000000000000000000");

    real x = -(1 + 2 / real (60) + real (2.99) / 3600);

    CHECK_EQUAL (DMS::Encode (x, DMS::DEGREE, 0, DMS::NONE), "-1");
    CHECK_EQUAL (DMS::Encode (x, DMS::DEGREE, 0, DMS::LATITUDE), "01S");
    CHECK_EQUAL (DMS::Encode (x, DMS::DEGREE, 0, DMS::LONGITUDE), "001W");
    CHECK_EQUAL (DMS::Encode (-x, DMS::DEGREE, 0, DMS::AZIMUTH), "001");
    CHECK_EQUAL (DMS::Encode (x, DMS::DEGREE, 1, DMS::NONE), "-1.0");
    CHECK_EQUAL (DMS::Encode (x, DMS::DEGREE, 1, DMS::LATITUDE), "01.0S");
    CHECK_EQUAL (DMS::Encode (x, DMS::DEGREE, 1, DMS::LONGITUDE), "001.0W");
    CHECK_EQUAL (DMS::Encode (-x, DMS::DEGREE, 1, DMS::AZIMUTH), "001.0");
    CHECK_EQUAL (DMS::Encode (x, DMS::MINUTE, 0, DMS::NONE), "-1d02'");
    CHECK_EQUAL (DMS::Encode (x, DMS::MINUTE, 0, DMS::LATITUDE), "01d02'S");
    CHECK_EQUAL (DMS::Encode (x, DMS::MINUTE, 0, DMS::LONGITUDE), "001d02'W");
    CHECK_EQUAL (DMS::Encode (-x, DMS::MINUTE, 0, DMS::AZIMUTH), "001d02'");
    CHECK_EQUAL (DMS::Encode (x, DMS::MINUTE, 1, DMS::NONE), "-1d02.0'");
    CHECK_EQUAL (DMS::Encode (x, DMS::MINUTE, 1, DMS::LATITUDE), "01d02.0'S");
    CHECK_EQUAL (DMS::Encode (x, DMS::MINUTE, 1, DMS::LONGITUDE), "001d02.0'W");
    CHECK_EQUAL (DMS::Encode (-x, DMS::MINUTE, 1, DMS::AZIMUTH), "001d02.0'");
    CHECK_EQUAL (DMS::Encode (x, DMS::SECOND, 0, DMS::NONE), "-1d02'03\"");
    CHECK_EQUAL (DMS::Encode (x, DMS::SECOND, 0, DMS::LATITUDE), "01d02'03\"S");
    CHECK_EQUAL (DMS::Encode (x, DMS::SECOND, 0, DMS::LONGITUDE), "001d02'03\"W");
    CHECK_EQUAL (DMS::Encode (-x, DMS::SECOND, 0, DMS::AZIMUTH), "001d02'03\"");
    CHECK_EQUAL (DMS::Encode (x, DMS::SECOND, 1, DMS::NONE), "-1d02'03.0\"");
    CHECK_EQUAL (DMS::Encode (x, DMS::SECOND, 1, DMS::LATITUDE), "01d02'03.0\"S");
    CHECK_EQUAL (DMS::Encode (x, DMS::SECOND, 1, DMS::LONGITUDE), "001d02'03.0\"W");
    CHECK_EQUAL (DMS::Encode (-x, DMS::SECOND, 1, DMS::AZIMUTH), "001d02'03.0\"");
  }

  TEST (Decode)
  {
    DMS::flag ind;
    CHECK_EQUAL (DMS::Decode (" +0 ", ind), +0.0);
    CHECK_EQUAL (DMS::Decode ("-0  ", ind), -0.0);
    CHECK_NAN (DMS::Decode (" nan", ind));
    CHECK_EQUAL (DMS::Decode ("+inf", ind), +inf);
    CHECK_EQUAL (DMS::Decode (" inf", ind), +inf);
    CHECK_EQUAL (DMS::Decode ("-inf", ind), -inf);
    CHECK_EQUAL (DMS::Decode (" +0N", ind), +0.0);
    CHECK_EQUAL (DMS::Decode ("-0N ", ind), -0.0);
    CHECK_EQUAL (DMS::Decode ("+0S ", ind), -0.0);
    CHECK_EQUAL (DMS::Decode (" -0S", ind), +0.0);
  }
}

