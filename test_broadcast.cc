
#include "test.h"
#include <cstring>

void test_broadcast()
{
  {
  char rinex[]=   // mixed exponents
"G 1 2024 07 15 22 00 00 2.417615614831D-04-7.162270776462d-12 0.000000000000E+00\n"
"     8.100000000000e+01 7.109375000000D+01 5.544516665660d-09-2.734409949984E+00\n"
"     3.712251782417e-06 1.337930315640D-02 9.505078196526d-06 5.153784959793E+03\n"
"     1.656000000000e+05 8.754432201385D-08-9.331363423370d-01-7.450580596924E-08\n"
"     9.540004770577e-01 1.908750000000D+02 1.030326911217d+00-7.669248026393E-09\n"
"    -1.371485699313e-10 1.000000000000D+00 2.323000000000d+03 0.000000000000E+00\n"
"     2.000000000000e+00 6.300000000000D+01-1.955777406693d-08 8.100000000000E+01\n"
"     1.584180000000e+05 4.000000000000D+00";

    Nav eph(rinex);
    
    if(strncmp(eph.prn,"G01",3))
      fail("wrong PRN");
    if(fabs(eph.toc.to_double()-1721080800.0)>0.001)
      fail("wrong Time of Clock (as Unix Timestamp)");
    if(eph.week!=2323)
      fail("wrong GPS Week");
    if(eph.week!=eph.toc.gps_week())
      fail("wrong GPS Week parsed");
    if(fabs(eph.toc.gps_tow()-165600.0)>0.001)
      fail("wrong GPS Time of Week for Time of Clock");
    
    //std::cout << std::fixed;
  }
}
