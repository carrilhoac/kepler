
#include "test.h"

void test_time()
{  
  { // from_rnx()
    Time rnx; // no leap seconds considered
    
    rnx.from_rnx("1980 01 06 00 00 00"); // GPS epoch
    if(fabs(rnx.to_double()-315964800.0)>0.001)
      fail("different timestamp (315964800 s)");
    
    rnx.from_rnx("1970 01 01 00 00 00"); // Unix epoch
    if(fabs(rnx.to_double())>0.001)
      fail("different timestamp (0 s)");
      
    rnx.from_rnx("1999 01 01 00 00 00"); 
    if(rnx.gps_week()!=990)
      fail("different GPS Week (990)");
    if(fabs(rnx.gps_tow()-432000.0)>0.001)
      fail("different GPS Time of Week (432000 s)");
    
    rnx.from_rnx("2000 01 01 00 00 00"); 
    if(rnx.gps_week()!=1042)
      fail("different GPS Week (1042)");
    if(fabs(rnx.gps_tow()-518400.0)>0.001)
      fail("different GPS Time of Week (518400 s)");
      
    rnx.from_rnx("1993 12 10 00 00 00"); // DooM release date
    if(rnx.gps_week()!=726)
      fail("different GPS Week (726)");
    if(fabs(rnx.gps_tow()-432000.0)>0.001)
      fail("different GPS Time of Week (432000 s)");
    if(fabs(rnx.to_double()-755481600.0)>0.001)
      fail("different timestamp (755481600 s)");
    
    rnx.from_rnx("2024 07 15 22 00 00");
    if(rnx.gps_week()!=2323)
      fail("different GPS Week (2323)");
    if(fabs(rnx.gps_tow()-165600.0)>0.001)
      fail("different GPS Time of Week (165600 s)");
    if(fabs(rnx.to_double()-1721080800.0)>0.001)
      fail("different timestamp (1721080800 s)");
    
    rnx.from_rnx("2030 08 08 14 19 57"); 
    if(rnx.gps_week()!=2639)
      fail("different GPS Week (2639)");
    if(fabs(rnx.gps_tow()-397197.0)>0.001)
      fail("different GPS Time of Week (397197 s)");
    if(fabs(rnx.to_double()-1912429197.0)>0.001)
      fail("different timestamp (1912429197 s)");
  
    //std::cout << std::fixed;
    //std::cout << rnx.to_double() << std::endl;
  }
}
