
#include "test.h"

#define D2R      0.017453292519943295 // pi/180
#define R2D      57.29577951308232087 // 180/pi

/*
G12 20096650.251 -6284389.016 16054576.329
G18 16612563.445 -20451984.604 1114099.112
G29 21230934.099 -5654129.427 -14829027.098
G30 -2939357.310 -18880281.743 -18814949.538
G21 7902410.931 -17362539.372 -17739506.489
G15 26409819.803 3695186.959 735104.310
G22 5895222.365 -22160371.672 13479675.997
G25 19173618.669 -17646407.492 4923460.234
G31 -6334182.339 -25491495.029 3097282.676
G16 -10299229.490 -10599464.909 -22065484.500
*/

static void test_azel_ex(double *rec, double *sat, double *ref_azel)
{
  double r,los[3],geo[3],azel[2];
  
  Spheroid wgs84(Spheroid::WGS84);
  wgs84.ecf2geo(rec,geo);
  
  geo[0]*=D2R;
  geo[1]*=D2R;
  
  r=geomdist(sat, rec, los);
  satazel(geo, los, azel);
  
  azel[0]*=R2D;
  azel[1]*=R2D;
  
  std::cout << std::fixed;
  std::cout << r << " " << azel[0] << " " << azel[1] << std::endl;
}

static void test_azel()
{
  double rec[3]={ 3687624.367, -4620818.683, -2386880.382}; // PPTE
 // double g12[3]={20096650.251, -6284389.016, 16054576.329};
  double g30[3]={-2939357.310, -18880281.743, -18814949.538};
  
  test_azel_ex(rec, g30, 0);
  //fail("incorrect ECEF coordinates");
}

void test_atmosphere()
{
  test_azel();
}
