
#include "test.h"

void test_from_rnx(const char *epoch, int gps_week, 
                   double gps_tow, double unix_timestamp)
{
  Time t;
  t.from_rnx(epoch);
  
  if(t.gps_week()!=gps_week)
    fail("different GPS Week");
  if(fabs(t.gps_tow()-gps_tow)>0.001)
    fail("different GPS Time of Week");
  if(fabs(t.to_double()-unix_timestamp)>0.001)
    fail("different Unix timestamp");
}

void test_from_rnx_unix(const char *epoch, double unix_timestamp)
{
  Time t;
  t.from_rnx(epoch);
  
  if(fabs(t.to_double()-unix_timestamp)>0.001)
    fail("different Unix timestamp");
}

void test_time()
{  
  { // from_rnx()
    Time rnx; // no leap seconds considered
    
    // TODO: read from text file 
    test_from_rnx_unix("1969 12 31 23 59 59",-1.0); // 1s before Unix epoch
    test_from_rnx_unix("1970 01 01 00 00 00", 0.0); // Unix epoch
    test_from_rnx_unix("1980 01 06 00 00 00", 315964800.0);// GPS epoch
    
    // this first one is the DooM release date
    test_from_rnx("1993 12 10 00 00 00",  726, 432000.0,  755481600.0);
    test_from_rnx("1999 01 01 00 00 00",  990, 432000.0,  915148800.0);
    test_from_rnx("1999 12 31 23 59 59", 1042, 518399.0,  946684799.0);
    test_from_rnx("2000 01 01 00 00 00", 1042, 518400.0,  946684800.0);
    test_from_rnx("2024 07 15 22 00 00", 2323, 165600.0, 1721080800.0);
    test_from_rnx("2030 08 08 14 19 57", 2639, 397197.0, 1912429197.0);
  }
}
