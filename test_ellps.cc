
#include "test.h"

static int test_geodesic()
{
  double d;
  Spheroid wgs84(ELLPS_WGS84);
  
  // distance between PPTE and BRAZ (RBMC)
  d=wgs84.geodesic(
    -22.119904740399434, -51.408534025148890, // PPTE
    -15.947475342062091, -47.877868983050924);// BRAZ
    
  if(fabs(777678.312-d)>0.0005)
    fail("incorrect geodesic distance");
  
  return 0;
}

static int test_geo2ecf()
{
  // PPTE
  double ref[3]={3687624.3674,-4620818.6827,-2386880.3805};
  double geo[3]={-22.119904740399434, -51.408534025148890, 431.049};
  double xyz[3]={0};
  
  Spheroid grs80(ELLPS_GRS80);
  grs80.geo2ecf(geo, xyz);
  
  xyz[0]-=ref[0];
  xyz[1]-=ref[1];
  xyz[2]-=ref[2];
  
  if(fabs(xyz[0])>0.00025
   ||fabs(xyz[1])>0.00025
   ||fabs(xyz[2])>0.00025)
    fail("incorrect ECEF coordinates");
    
  return 0;
}

static int test_ecf2geo()
{
  double ref[3]={-22.119904740399434,-51.408534025148890,431.049};
  double geo[3]={0};
  double xyz[3]={3687624.3672, -4620818.6825, -2386880.3804};
  
  Spheroid grs80(ELLPS_GRS80);
  grs80.ecf2geo(xyz, geo);
  
  geo[0]-=ref[0];
  geo[1]-=ref[1];
  geo[2]-=ref[2];
  
  if(fabs(geo[0])>1e-9
   ||fabs(geo[1])>1e-9
   ||fabs(geo[2])>1e-4)
    fail("incorrect Geodetic coordinates");
    
  return 0;
}

static int test_geo2utm()
{
  double ref[2]={457866.057,7553844.609};
  double geo[3]={-22.119904740399434, -51.408534025148890, 431.049};
  double utm[2]={0};
  int zone;
  char h;
  
  Spheroid grs80(ELLPS_GRS80);
  grs80.geo2utm(geo, utm,&zone,&h);
  
  if(zone!=22||h!='S')
    fail("incorrect UTM zone and/or hemisphere");
  
  utm[0]-=ref[0];
  utm[1]-=ref[1];
  
  if(fabs(utm[0])>1e-3
   ||fabs(utm[1])>1e-3)
    fail("incorrect UTM coordinates");
    
  return 0;
}

static int test_utm2geo()
{
  double ref[2]={-22.119904740399434,-51.408534025148890};
  double geo[2]={0};
  double utm[2]={457866.057,7553844.609};
  
  Spheroid grs80(ELLPS_GRS80);
  grs80.utm2geo(utm, geo, 22,'S');
  
  //printf("%.9f %.9f\n",geo[0],geo[1]);
  
  geo[0]-=ref[0];
  geo[1]-=ref[1];
  
  if(fabs(geo[0])>1e-4
   ||fabs(geo[1])>1e-4)
    fail("incorrect Geodetic coordinates");
    
  //printf("%.9f %.9f\n",geo[0],geo[1]);
  return 0;
}

void test_ellps()
{ 
  test_geodesic();
  test_geo2ecf();
  test_ecf2geo();
  test_geo2utm();
  test_utm2geo();
  //std::cout<<"[ELLPS] all tests run successfully"<<std::endl;
}
