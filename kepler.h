
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

//////////////////////////////////////////////////////////////////////
//  Time point (needs refactoring)

struct Epoch{
  int64_t t_sec;
  double t_frac;
};

Epoch rnx2unx(const char *rnx);
Epoch gps2unx(int gps_week, double gps_tow);
double unx2gps(const Epoch& utc_ts, int *gps_week);
double timesub(const Epoch& a, const Epoch& b);

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
  static void utmzone(double lon_deg, int *zone, double *cm_deg);
  double cnflat(double phi) const;
};


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
  Epoch toe,toc,ttr;// Toe,Toc,T_trans 
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
  Nav(const char *rnx);
  void rnx2nav(const char *rnx);
  void nav2ecf(const Epoch& t, double *xyz, double *clk_bias) const;
};

#endif 
