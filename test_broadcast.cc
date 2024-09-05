
#include "test.h"
#include <cstring>

void test_broadcast()
{
  {
  // a valid (but nasty) GPS nav record
  //
  // missing minutes and seconds (both zero)
  // missing leading 0 on sat PRN
  // missing leading 0 on some values (starts at the dot)
  // missing some elements (which are zero, btw)
  // mixed exponents (e,E,d,D), mixed EOL (Win32:'\r\n'; Linux:'\n'; and Mac:'\r')
  char rinex[]= 
"G 1 2024  7 15 22       2.417615614831D-04-7.162270776462d-12\r"
"      .810000000000e+02 7.109375000000D+01 5.544516665660d-09-2.734409949984E+00\n"
"     3.712251782417e-06  .133793031564D-01 9.505078196526d-06 5.153784959793E+03\r\n"
"      .165600000000e+06 8.754432201385D-08-9.331363423370d-01-7.450580596924E-08 \t \n"
"     9.540004770577e-01  .190875000000D+03 1.030326911217d+00-7.669248026393E-09 \r"
"    -1.371485699313e-10 1.000000000000D+00  .232300000000d+04\n"
"     2.000000000000e+00 6.300000000000D+01-1.955777406693d-08 8.100000000000E+01\r\n"
"      .158418000000e+06 4.000000000000D+00";

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
  }
}
