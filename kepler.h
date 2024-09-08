
// ---------------------------------------------------------------------------
//  Copyright (c) 2009-2024, Andre Caceres Carrilho. All rights reserved. 
//
//   kepler.h --Main library header
// ---------------------------------------------------------------------------

#ifndef KEPLER_h
#define KEPLER_h

#if (__cplusplus < 201703L)
  #error C++17 or later is required 
  /// C++17 or later is needed in order to have proper alignment of 
  /// SIMD data inside std::vector<> containers (otherwise segfault)
#endif

#if defined(__AVX2__)
  #include <immintrin.h>
  #define SIMD_x86 1 /// Assuming we have both AVX2 and FMA3 
#endif 

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cstring>

#include <cmath>
#include <cstdio>
#include <cstdint>
#include <cstdlib>

#ifdef DEBUG
  void error_(const char *from, const char *msg);
  #define error(x) error_(__PRETTY_FUNCTION__, x) // & exit 1  
  void warn_(const char *from, const char *msg);
  #define warn(x) warn_(__PRETTY_FUNCTION__, x)
#endif 

//////////////////////////////////////////////////////////////////////
// Analogous to hypot() but for 3D values
//  sqrt(x*x + y*y + z*z)
template<typename real>
real pythag(real x, real y, real z)
{
  real w;
  x= fabs(x);
  y= fabs(y);
  z= fabs(z);
  if(x > y){ w=x; x=y; y=w; }
  if(x > z){ w=x; x=z; z=w; }
  if(y > z){ w=y; y=z; z=w; }
  w= real(1) / z;
  x*= w;
  y*= w;
  return sqrt(x*x + y*y + real(1)) * z;
}


//////////////////////////////////////////////////////////////////////
//  Unix Time point

class alignas(16) Time{
public:
  int64_t t_sec;
  double t_frac;
public:
  Time();
  Time(const Time& t);
  template<typename FT> Time(FT t);
  Time& from_cal(const int *cal); // y m d h m s
  Time& from_rnx(const char *rnx); // RINEX epoch string
  Time& from_gps(int gpsweek, double gpstow);
  int gps_week() const;
  double gps_tow() const;
  double to_double() const;
  bool operator>(const Time& t) const; // for sorting
  bool operator<(const Time& t) const;
  bool operator>=(const Time& t) const;
  bool operator<=(const Time& t) const;
  template<typename FT> bool operator>(FT t) const;
  template<typename FT> bool operator<(FT t) const;
  template<typename FT> bool operator>=(FT t) const;
  template<typename FT> bool operator<=(FT t) const;
  Time& operator=(const Time& t);
  Time& operator+=(const Time& t);
  Time& operator-=(const Time& t);
  Time operator+(const Time& t) const;
  Time operator-(const Time& t) const;
  template<typename FT> Time& operator=(FT t);
  template<typename FT> Time& operator+=(FT t);
  template<typename FT> Time& operator-=(FT t);
  template<typename FT> Time operator+(FT t) const;
  template<typename FT> Time operator-(FT t) const;
  friend std::ostream& operator<<(std::ostream& os, const Time& t);
  static int civ2day(int y, int m, int d);
  static void day2civ(int z, int& y, int& m, int& d);
  static int64_t cal2unx(const int *cal);
  static void unx2cal(int64_t secs, int *cal);
private:
  void normalize();
};

#include "timesys.h" // implementation for templated Time methods 


//////////////////////////////////////////////////////////////////////
//  Vector 3D

class alignas(32) Vec3{
public:
  union{
    // last element is ignored throughout the code
    // it is only defined to match the size of AVX2 register
    double v[4]; 
#ifdef SIMD_x86  
    __m256d r;   
#endif 
  };
public:
  Vec3(double X=0.0,double Y=0.0,double Z=0.0);
  Vec3(const Vec3& b);
  Vec3& operator=(const Vec3& b);
  Vec3& operator+=(const Vec3& b);
  Vec3& operator-=(const Vec3& b);
  Vec3& operator*=(const Vec3& b); // element-wise product
  Vec3& operator*=(double s);
  Vec3& operator/=(double s);
  Vec3 operator+(const Vec3& b) const;
  Vec3 operator-(const Vec3& b) const;
  Vec3 operator*(const Vec3& b) const;
  Vec3 operator*(double s) const;
  Vec3 operator/(double s) const;
  Vec3& unit();
  double norm() const;
        double* data()      { return &v[0]; }
  const double* data() const{ return &v[0]; }
  double x() const{ return v[0]; }
  double y() const{ return v[1]; }
  double z() const{ return v[2]; }
  double& x() { return v[0]; }
  double& y() { return v[1]; }
  double& z() { return v[2]; }
  friend std::ostream& operator<<(std::ostream& os, const Vec3& b);
};

double dist(const Vec3& a, const Vec3& b);
double dot(const Vec3& a, const Vec3& b);
Vec3 cross(const Vec3& a, const Vec3& b);
void cross(Vec3& c, const Vec3& a, const Vec3& b);


#define MATRIX_EPS 1e-11


//////////////////////////////////////////////////////////////////////
//  Matrix (row-major layout)

class Mat{
private:
  int m,n;
  std::vector<double> mat;
public:
  Mat(int m_=3, int n_=3);
  Mat(const Mat& b);
  Mat& zero();
  Mat& eye(); 
  Mat t() const;
  Mat inv() const;
  //Mat invspd() const;
  //double det() const;
  double tr() const;
  int rows() const{ return m; }
  int cols() const{ return n; }
  //Mat row(int i) const;
  //Mat col(int j) const;
        double* data()      { return &mat[0]; }
  const double* data() const{ return &mat[0]; }
        double& operator()(int i, int j);
  const double& operator()(int i, int j) const;
  Mat operator*(const Mat& b) const;
  Mat operator+(const Mat& b) const;
  Mat operator-(const Mat& b) const;
  Mat operator*(double s) const;
  Mat operator/(double s) const;
  Mat operator-() const;
  Mat& operator=(const Mat& b);
  Mat& operator+=(const Mat& b);
  Mat& operator-=(const Mat& b);
  Mat& operator*=(double s);
  Mat& operator/=(double s);
  bool issquare() const{ return m==n; }
  bool iszero(double tol=MATRIX_EPS) const;
  bool issymmetric(double tol=MATRIX_EPS) const;
  bool isidentity(double tol=MATRIX_EPS) const;
  bool compare(const Mat& b, double tol=MATRIX_EPS) const;
  bool operator==(const Mat& b) const{ return compare(b); }
  friend std::ostream& operator<<(std::ostream& os, const Mat& b);
};


//////////////////////////////////////////////////////////////////////
//  Spheroid 

class Spheroid{
private:
  double a,b,f,e,e2,rr;
  double alpha[8],beta[8];  
public:
  enum eSPHEROID{
    WGS66,
    WGS72,
    WGS84,
    GRS67,
    GRS80,
    PZ90 ,
    SAD69
  };
public:
  Spheroid(eSPHEROID ellps=WGS84);
  void set(eSPHEROID ellps);
  void set(double a_, double f_);  
  double geodesic(double lat1deg, double lon1deg, 
                  double lat2deg, double lon2deg) const;
  void ecf2geo(const double *xyz, double *geo) const;
  void geo2ecf(const double *geo, double *xyz) const;
  void utm2geo(const double *utm, double *geo, int zone, char h) const;
  void geo2utm(const double *geo, double *utm, int *zone, char *h) const;
  static double utmscale(double lon_deg);
private:
  static void utmzone(double lon_deg, int *zone, double *cm_deg); // or public?
  double cnflat(double phi) const;
};



//////////////////////////////////////////////////////////////////////
//  Broadcast ephemeris

struct Nav{
public:
  char prn[4];
//int sat;            // satellite number 
  int iode,iodc;      // IODE,IODC 
  int sva;            // SV accuracy (URA index) 
  int svh;            // SV health (0:ok) 
  int week;           // GPS/QZS: gps week, GAL: galileo week 
  int code;           // GPS/QZS: code on L2, GAL/CMP: data sources 
  int flag;           // GPS/QZS: L2 P data flag, CMP: nav type 
  Time toe,toc,ttr;   // Toe,Toc,T_trans 
                      // SV orbit parameters 
  double A,e,i0,OMG0,omg,M0,deln,OMGd,idot;
  double crc,crs,cuc,cus,cic,cis;
  double toes;        // Toe (s) in week 
  double fit;         // fit interval (h) 
  double f0,f1,f2;    // SV clock parameters (af0,af1,af2) 
  double tgd[4];      // group delay parameters 
                      // GPS/QZS:tgd[0]=TGD 
                      // GAL    :tgd[0]=BGD E5a/E1,tgd[1]=BGD E5b/E1 
                      // CMP    :tgd[0]=BGD1,tgd[1]=BGD2 
  double Adot,ndot;   // Adot,ndot for CNAV 
public:
  Nav();
  Nav(const char *rnx);
  void rnx2nav(const char *rnx);
  void nav2ecf(const Time& t, double *xyz, double *clock_bias) const;
  Vec3 nav2ecf(const Time& t, double *clock_bias) const;
  std::string nav2rnx() const;
};

#endif 
