/**
 * \file Geoid.hpp
 * \brief Header for GeographicLib::Geoid class
 *
 * Copyright (c) Charles Karney (2009-2022) <karney@alum.mit.edu> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_GEOID_HPP)
#define GEOGRAPHICLIB_GEOID_HPP 1

#include <vector>
#include <fstream>
#include <GeographicLib/Constants.hpp>

#if defined(_MSC_VER)
// Squelch warnings about dll vs vector and constant conditional expressions
#  pragma warning (push)
#  pragma warning (disable: 4251 4127)
#endif

#if !defined(GEOGRAPHICLIB_GEOID_PGM_PIXEL_WIDTH)
/**
 * The size of the pixel data in the pgm data files for the geoids.  2 is the
 * standard size corresponding to a maxval 2<sup>16</sup>&minus;1.  Setting it
 * to 4 uses a maxval of 2<sup>32</sup>&minus;1 and changes the extension for
 * the data files from .pgm to .pgm4.  Note that the format of these pgm4 files
 * is a non-standard extension of the pgm format.
 **********************************************************************/
#  define GEOGRAPHICLIB_GEOID_PGM_PIXEL_WIDTH 2
#endif

namespace GeographicLib {

  /**
   * \brief Looking up the height of the geoid above the ellipsoid
   *
   * This class evaluates the height of one of the standard geoids, EGM84,
   * EGM96, or EGM2008 by bilinear or cubic interpolation into a rectangular
   * grid of data.  These geoid models are documented in
   * - EGM84:
   *   https://earth-info.nga.mil/index.php?dir=wgs84&action=wgs84#tab_egm84
   * - EGM96:
   *   https://earth-info.nga.mil/index.php?dir=wgs84&action=wgs84#tab_egm96
   * - EGM2008:
   *   https://earth-info.nga.mil/index.php?dir=wgs84&action=wgs84#tab_egm2008
   *
   * The geoids are defined in terms of spherical harmonics.  However in order
   * to provide a quick and flexible method of evaluating the geoid heights,
   * this class evaluates the height by interpolation into a grid of
   * precomputed values.
   *
   * The height of the geoid above the ellipsoid, \e N, is sometimes called the
   * geoid undulation.  It can be used to convert a height above the ellipsoid,
   * \e h, to the corresponding height above the geoid (the orthometric height,
   * roughly the height above mean sea level), \e H, using the relations
   *
   * &nbsp;&nbsp;&nbsp;\e h = \e N + \e H;
   * &nbsp;&nbsp;\e H = &minus;\e N + \e h.
   *
   * See \ref geoid for details of how to install the data sets, the data
   * format, estimates of the interpolation errors, and how to use caching.
   *
   * This class is typically \e not thread safe in that a single instantiation
   * cannot be safely used by multiple threads because of the way the object
   * reads the data set and because it maintains a single-cell cache.  If
   * multiple threads need to calculate geoid heights they should all construct
   * thread-local instantiations.  Alternatively, set the optional \e
   * threadsafe parameter to true in the constructor.  This causes the
   * constructor to read all the data into memory and to turn off the
   * single-cell caching which results in a Geoid object which \e is thread
   * safe.
   *
   * Example of use:
   * \include examples/Geoid.cpp
   *
   * <a href="GeoidEval.1.html">GeoidEval</a> is a command-line utility
   * providing access to the functionality of Geoid.
   **********************************************************************/

  class GEOGRAPHICLIB_EXPORT Geoid {
  private:
#if GEOGRAPHICLIB_GEOID_PGM_PIXEL_WIDTH != 4
    typedef unsigned short pixel_t;
    static const unsigned pixel_size_ = 2;
    static const unsigned pixel_max_ = 0xffffu;
#else
    typedef unsigned pixel_t;
    static const unsigned pixel_size_ = 4;
    static const unsigned pixel_max_ = 0xffffffffu;
#endif
    static const unsigned stencilsize_ = 12;
    static const unsigned nterms_ = ((3 + 1) * (3 + 2))/2; // for a cubic fit
    static const int c0_;
    static const int c0n_;
    static const int c0s_;
    static const int c3_[stencilsize_ * nterms_];
    static const int c3n_[stencilsize_ * nterms_];
    static const int c3s_[stencilsize_ * nterms_];

    std::string name_, dir_, filename_;
    const bool cubic_;
    const real a_, e2_, degree_, eps_;
    mutable std::ifstream file_;
    real rlonres_, rlatres_;
    std::string description_, datetime_;
    real offset_, scale_, maxerror_, rmserror_;
    int width_, height_;
    unsigned long long datastart_, swidth_;
    bool threadsafe_;
    // Area cache
    mutable std::vector< std::vector<pixel_t> > data_;
    mutable bool cache_;
    // NE corner and extent of cache
    mutable int xoffset_, yoffset_, xsize_, ysize_;
    // Cell cache
    mutable int ix_, iy_;
    mutable real v00_, v01_, v10_, v11_;
    mutable real t_[nterms_];
    void filepos(int ix, int iy) const {
      file_.seekg(std::streamoff
                  (datastart_ +
                   pixel_size_ * (unsigned(iy)*swidth_ + unsigned(ix))));
    }
    real rawval(int ix, int iy) const {
      if (ix < 0)
        ix += width_;
      else if (ix >= width_)
        ix -= width_;
      if (cache_ && iy >= yoffset_ && iy < yoffset_ + ysize_ &&
          ((ix >= xoffset_ && ix < xoffset_ + xsize_) ||
           (ix + width_ >= xoffset_ && ix + width_ < xoffset_ + xsize_))) {
        return real(data_[iy - yoffset_]
                    [ix >= xoffset_ ? ix - xoffset_ : ix + width_ - xoffset_]);
      } else {
        if (iy < 0 || iy >= height_) {
          iy = iy < 0 ? -iy : 2 * (height_ - 1) - iy;
          ix += (ix < width_/2 ? 1 : -1) * width_/2;
        }
        try {
          filepos(ix, iy);
          // initial values to suppress warnings in case get fails
          char a = 0, b = 0;
          file_.get(a);
          file_.get(b);
          unsigned r = ((unsigned char)(a) << 8) | (unsigned char)(b);
          // for C++17 use if constexpr
          if (pixel_size_ == 4) {
            file_.get(a);
            file_.get(b);
            r = (r << 16) | ((unsigned char)(a) << 8) | (unsigned char)(b);
          }
          return real(r);
        }
        catch (const std::exception& e) {
          // throw GeographicErr("Error reading " + filename_ + ": "
          //                      + e.what());
          // triggers complaints about the "binary '+'" under Visual Studio.
          // So use '+=' instead.
          std::string err("Error reading ");
          err += filename_;
          err += ": ";
          err += e.what();
          throw GeographicErr(err);
        }
      }
    }
    real height(real lat, real lon) const;
    Geoid(const Geoid&) = delete;            // copy constructor not allowed
    Geoid& operator=(const Geoid&) = delete; // copy assignment not allowed
  public:

    /**
     * Flags indicating conversions between heights above the geoid and heights
     * above the ellipsoid.
     **********************************************************************/
    enum convertflag {
      /**
       * The multiplier for converting from heights above the ellipsoid to
       * heights above the geoid.
       **********************************************************************/
      ELLIPSOIDTOGEOID = -1,
      /**
       * No conversion.
       **********************************************************************/
      NONE = 0,
      /**
       * The multiplier for converting from heights above the geoid to heights
       * above the ellipsoid.
       **********************************************************************/
      GEOIDTOELLIPSOID = 1,
    };

    /** \name Setting up the geoid
     **********************************************************************/
    ///@{
    /**
     * Construct a geoid.
     *
     * @param[in] name the name of the geoid.
     * @param[in] path (optional) directory for data file.
     * @param[in] cubic (optional) interpolation method; false means bilinear,
     *   true (the default) means cubic.
     * @param[in] threadsafe (optional), if true, construct a thread safe
     *   object.  The default is false
     * @exception GeographicErr if the data file cannot be found, is
     *   unreadable, or is corrupt.
     * @exception GeographicErr if \e threadsafe is true but the memory
     *   necessary for caching the data can't be allocated.
     *
     * The data file is formed by appending ".pgm" to the name.  If \e path is
     * specified (and is non-empty), then the file is loaded from directory, \e
     * path.  Otherwise the path is given by DefaultGeoidPath().  If the \e
     * threadsafe parameter is true, the data set is read into memory, the data
     * file is closed, and single-cell caching is turned off; this results in a
     * Geoid object which \e is thread safe.
     **********************************************************************/
    explicit Geoid(const std::string& name, const std::string& path = "",
                   bool cubic = true, bool threadsafe = false);

    /**
     * Set up a cache.
     *
     * @param[in] south latitude (degrees) of the south edge of the cached
     *   area.
     * @param[in] west longitude (degrees) of the west edge of the cached area.
     * @param[in] north latitude (degrees) of the north edge of the cached
     *   area.
     * @param[in] east longitude (degrees) of the east edge of the cached area.
     * @exception GeographicErr if the memory necessary for caching the data
     *   can't be allocated (in this case, you will have no cache and can try
     *   again with a smaller area).
     * @exception GeographicErr if there's a problem reading the data.
     * @exception GeographicErr if this is called on a threadsafe Geoid.
     *
     * Cache the data for the specified "rectangular" area bounded by the
     * parallels \e south and \e north and the meridians \e west and \e east.
     * \e east is always interpreted as being east of \e west, if necessary by
     * adding 360&deg; to its value.  \e south and \e north should be in
     * the range [&minus;90&deg;, 90&deg;].
     **********************************************************************/
    void CacheArea(real south, real west, real north, real east) const;

    /**
     * Cache all the data.
     *
     * @exception GeographicErr if the memory necessary for caching the data
     *   can't be allocated (in this case, you will have no cache and can try
     *   again with a smaller area).
     * @exception GeographicErr if there's a problem reading the data.
     * @exception GeographicErr if this is called on a threadsafe Geoid.
     *
     * On most computers, this is fast for data sets with grid resolution of 5'
     * or coarser.  For a 1' grid, the required RAM is 450MB; a 2.5' grid needs
     * 72MB; and a 5' grid needs 18MB.
     **********************************************************************/
    void CacheAll() const { CacheArea(real(-90), real(0),
                                      real(90), real(360)); }

    /**
     * Clear the cache.  This never throws an error.  (This does nothing with a
     * thread safe Geoid.)
     **********************************************************************/
    void CacheClear() const;

    ///@}

    /** \name Compute geoid heights
     **********************************************************************/
    ///@{
    /**
     * Compute the geoid height at a point
     *
     * @param[in] lat latitude of the point (degrees).
     * @param[in] lon longitude of the point (degrees).
     * @exception GeographicErr if there's a problem reading the data; this
     *   never happens if (\e lat, \e lon) is within a successfully cached
     *   area.
     * @return the height of the geoid above the ellipsoid (meters).
     *
     * The latitude should be in [&minus;90&deg;, 90&deg;].
     **********************************************************************/
    real operator()(real lat, real lon) const {
      return height(lat, lon);
    }

    /**
     * Convert a height above the geoid to a height above the ellipsoid and
     * vice versa.
     *
     * @param[in] lat latitude of the point (degrees).
     * @param[in] lon longitude of the point (degrees).
     * @param[in] h height of the point (degrees).
     * @param[in] d a Geoid::convertflag specifying the direction of the
     *   conversion; Geoid::GEOIDTOELLIPSOID means convert a height above the
     *   geoid to a height above the ellipsoid; Geoid::ELLIPSOIDTOGEOID means
     *   convert a height above the ellipsoid to a height above the geoid.
     * @exception GeographicErr if there's a problem reading the data; this
     *   never happens if (\e lat, \e lon) is within a successfully cached
     *   area.
     * @return converted height (meters).
     **********************************************************************/
    real ConvertHeight(real lat, real lon, real h,
                             convertflag d) const {
      return h + real(d) * height(lat, lon);
    }

    ///@}

    /** \name Inspector functions
     **********************************************************************/
    ///@{
    /**
     * @return geoid description, if available, in the data file; if
     *   absent, return "NONE".
     **********************************************************************/
    const std::string& Description() const { return description_; }

    /**
     * @return date of the data file; if absent, return "UNKNOWN".
     **********************************************************************/
    const std::string& DateTime() const { return datetime_; }

    /**
     * @return full file name used to load the geoid data.
     **********************************************************************/
    const std::string& GeoidFile() const { return filename_; }

    /**
     * @return "name" used to load the geoid data (from the first argument of
     *   the constructor).
     **********************************************************************/
    const std::string& GeoidName() const { return name_; }

    /**
     * @return directory used to load the geoid data.
     **********************************************************************/
    const std::string& GeoidDirectory() const { return dir_; }

    /**
     * @return interpolation method ("cubic" or "bilinear").
     **********************************************************************/
    const std::string Interpolation() const
    { return std::string(cubic_ ? "cubic" : "bilinear"); }

    /**
     * @return estimate of the maximum interpolation and quantization error
     *   (meters).
     *
     * This relies on the value being stored in the data file.  If the value is
     * absent, return &minus;1.
     **********************************************************************/
    real MaxError() const { return maxerror_; }

    /**
     * @return estimate of the RMS interpolation and quantization error
     *   (meters).
     *
     * This relies on the value being stored in the data file.  If the value is
     * absent, return &minus;1.
     **********************************************************************/
    real RMSError() const { return rmserror_; }

    /**
     * @return offset (meters).
     *
     * This in used in converting from the pixel values in the data file to
     * geoid heights.
     **********************************************************************/
    real Offset() const { return offset_; }

    /**
     * @return scale (meters).
     *
     * This in used in converting from the pixel values in the data file to
     * geoid heights.
     **********************************************************************/
    real Scale() const { return scale_; }

    /**
     * @return true if the object is constructed to be thread safe.
     **********************************************************************/
    bool ThreadSafe() const { return threadsafe_; }

    /**
     * @return true if a data cache is active.
     **********************************************************************/
    bool Cache() const { return cache_; }

    /**
     * @return west edge of the cached area; the cache includes this edge.
     **********************************************************************/
    real CacheWest() const {
      return cache_ ? ((xoffset_ + (xsize_ == width_ ? 0 : cubic_)
                        + width_/2) % width_ - width_/2) / rlonres_ :
        0;
    }

    /**
     * @return east edge of the cached area; the cache excludes this edge.
     **********************************************************************/
    real CacheEast() const {
      return  cache_ ?
        CacheWest() +
        (xsize_ - (xsize_ == width_ ? 0 : 1 + 2 * cubic_)) / rlonres_ :
        0;
    }

    /**
     * @return north edge of the cached area; the cache includes this edge.
     **********************************************************************/
    real CacheNorth() const {
      return cache_ ? real(90) - (yoffset_ + cubic_) / rlatres_ : 0;
    }

    /**
     * @return south edge of the cached area; the cache excludes this edge
     *   unless it's the south pole.
     **********************************************************************/
    real CacheSouth() const {
      return cache_ ?
        real(90) - ( yoffset_ + ysize_ - 1 - cubic_) / rlatres_ :
        0;
    }

    /**
     * @return \e a the equatorial radius of the WGS84 ellipsoid (meters).
     *
     * (The WGS84 value is returned because the supported geoid models are all
     * based on this ellipsoid.)
     **********************************************************************/
    real EquatorialRadius() const
    { return Constants::WGS84_a(); }

    /**
     * @return \e f the flattening of the WGS84 ellipsoid.
     *
     * (The WGS84 value is returned because the supported geoid models are all
     * based on this ellipsoid.)
     **********************************************************************/
    real Flattening() const { return Constants::WGS84_f(); }
    ///@}

    /**
     * @return the default path for geoid data files.
     *
     * This is the value of the environment variable GEOGRAPHICLIB_GEOID_PATH,
     * if set; otherwise, it is $GEOGRAPHICLIB_DATA/geoids if the environment
     * variable GEOGRAPHICLIB_DATA is set; otherwise, it is a compile-time
     * default (/usr/local/share/GeographicLib/geoids on non-Windows systems
     * and C:/ProgramData/GeographicLib/geoids on Windows systems).
     **********************************************************************/
    static std::string DefaultGeoidPath();

    /**
     * @return the default name for the geoid.
     *
     * This is the value of the environment variable GEOGRAPHICLIB_GEOID_NAME,
     * if set; otherwise, it is "egm96-5".  The Geoid class does not use this
     * function; it is just provided as a convenience for a calling program
     * when constructing a Geoid object.
     **********************************************************************/
    static std::string DefaultGeoidName();

  };

} // namespace GeographicLib

#if defined(_MSC_VER)
#  pragma warning (pop)
#endif

#endif  // GEOGRAPHICLIB_GEOID_HPP
