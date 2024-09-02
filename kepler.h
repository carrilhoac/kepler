
// ---------------------------------------------------------------------------
//  Copyright (c) 2009-2024, Andre Caceres Carrilho. All rights reserved. 
//
//   kepler.h --Main library header
// ---------------------------------------------------------------------------

#ifndef KEPLER_h
#define KEPLER_h

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <cstdint>
#include <cstdlib>

void error_(const char *from, const char *msg);
void warn_(const char *from, const char *msg);

#define error(x) error_(__PRETTY_FUNCTION__, x)
#define warn(x) warn_(__PRETTY_FUNCTION__, x)

// ---------------------------------------------------------------------------
//  Time point
// ---------------------------------------------------------------------------
struct epoch_t{
  int64_t t_sec;
  double t_frac;
};

epoch_t rnx2unx(const char *rnx);
epoch_t gps2unx(int gps_week, double gps_tow);
double unx2gps(const epoch_t& utc_ts, int *gps_week);
double timesub(const epoch_t& a, const epoch_t& b);

// ---------------------------------------------------------------------------
//  Spheroid 
// ---------------------------------------------------------------------------
enum ELLPS{
  ELLPS_WGS84,
  ELLPS_GRS80
};

class Spheroid{
private:
  double a,b,f,e,e2,pr;
  double alpha[8],beta[8];
  
public:
  Spheroid(ELLPS ellps=ELLPS_WGS84);
  void set(ELLPS ellps);
  void set(double a_, double f_);
  
  double geodesic(double lat1deg, double lon1deg, 
                  double lat2deg, double lon2deg) const;
  
  void ecf2geo(const double *xyz, double *geo) const;
  void geo2ecf(const double *geo, double *xyz) const;
  
  void utm2geo(const double *utm, double *geo, int zone, char h) const;
  void geo2utm(const double *geo, double *utm, int *zone, char *h) const;
  
  static double utmscale(double lon_deg);
  
private:
  static void utmzone(double lon_deg, int *zone, double *cm_deg);
  double cnflat(double phi) const;
};


#define MATRIX_EPS 1e-11

// ---------------------------------------------------------------------------
//  Matrix 
// ---------------------------------------------------------------------------
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
  friend std::ostream& operator << (std::ostream& os, const Mat& b);
};
  
#endif 
