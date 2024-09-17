
#include "test.h"
#include <cstring>

void test_ephemeris()
{
  {
    // a valid (but nasty) GPS nav record
    //
    // missing minutes and seconds (both zero)
    // missing leading 0 on sat PRN
    // missing leading 0 on some values (starts at the dot)
    // missing some elements (which are zero, btw)
    // mixed exponents (e,E,d,D), mixed EOL (Win32:'\r\n'; Unix:'\n'; and old MacOS:'\r')
    const char rinex[]= // "G1","G2"..."G9" would be valid? give wawrning perhaps?
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
        
    const std::string ref= // this is missing some values yet
"G01 2024 07 15 22 00 00 2.417615614831E-04-7.162270776462E-12 0.000000000000E+00\n"
"     8.100000000000E+01 7.109375000000E+01 5.544516665660E-09-2.734409949984E+00\n"
"     3.712251782417E-06 1.337930315640E-02 9.505078196526E-06 5.153784959793E+03\n"
"     1.656000000000E+05 8.754432201385E-08-9.331363423370E-01-7.450580596924E-08\n"
"     9.540004770577E-01 1.908750000000E+02 1.030326911217E+00-7.669248026393E-09\n"
"    -1.371485699313E-10 1.000000000000E+00 2.323000000000E+03 0.000000000000E+00\n"
"     0.000000000000E+00 6.300000000000E+01-1.955777406693E-08 8.100000000000E+01\n"
"     0.000000000000E+00 4.000000000000E+00 0.000000000000E+00 0.000000000000E+00\n";

    if(ref.compare(eph.nav2rnx()))
      fail("incorrect RINEX Nav output");
    
    #if 0
    Time t;
    double sat_xyz[3];
    double clock_bias;
    
    t.from_rnx("2024 07 15 22 00 00");
    eph.nav2ecf(t, sat_xyz, &clock_bias);
    
    sat_xyz[0]/=1e3; // from meters to kilometers
    sat_xyz[1]/=1e3;
    sat_xyz[2]/=1e3;
    clock_bias *= 1e6; // convert from seconds to microseconds
    
    std::cout<<std::fixed;
    std::cout<<sat_xyz[0]<<'\t'<<sat_xyz[1]<<'\t'<<sat_xyz[2]<<'\t';
    std::cout<<clock_bias<<std::endl;
   #endif 
   
//*  2024  7 15 22  0  0.00000000
//PG01 -10071.415123 -12260.750456 -21707.551434    241.257286

// From SP3  
// P: position
// G01 satellite PRN
// -10071.415123 -12260.750456 -21707.551434: ECEF coordinates (km)
// 241.257286 clock microseconds
  }
  
  
#if 1
  {
    /*
    const char rinex[]= 
"G12 2012 07 15 12 00 00 8.004158735280D-05 2.046363078990D-12 0.000000000000D+00\n"
"     5.800000000000D+01-1.050937500000D+02 3.911234347180D-09 1.104899710370D+00\n"
"    -5.422160029410D-06 4.072798416020D-03 7.687136530880D-06 5.153766078950D+03\n"
"     4.320000000000D+04 0.000000000000D+00 2.108448110950D+00 9.685754776000D-08\n"
"     9.803929996810D-01 2.392812500000D+02 1.127125627250D-01-7.783181343600D-09\n"
"    -2.389385241770D-10 1.000000000000D+00 1.697000000000D+03 0.000000000000D+00\n"
"     2.000000000000D+00 0.000000000000D+00-1.210719347000D-08 5.800000000000D+01\n"
"     3.957000000000D+04 4.000000000000D+00";
*/

    const char rinex[]= 
"G12 2012 07 15 10 00 00 8.002668619160D-05 2.046363078990D-12 0.000000000000D+00\n"
"     2.500000000000D+01-1.110000000000D+02 3.935878230840D-09 5.493558895060D-02\n"
"    -5.858018994330D-06 4.072350449860D-03 7.569789886470D-06 5.153764709470D+03\n"
"     3.600000000000D+04 7.264316082000D-08 2.108504061720D+00 2.421438694000D-08\n"
"     9.803942387720D-01 2.425625000000D+02 1.125559178460D-01-7.852827101770D-09\n"
"    -1.825076021740D-10 1.000000000000D+00 1.697000000000D+03 0.000000000000D+00\n"
"     2.000000000000D+00 0.000000000000D+00-1.210719347000D-08 2.500000000000D+01\n"
"     3.399000000000D+04 4.000000000000D+00";
    Nav eph(rinex);
    
    Time t;
    double sat_xyz[3];
    double clock_bias;
    
    t.from_rnx("2012 07 15 11 13 45");
    eph.nav2ecf(t-0.057, sat_xyz, &clock_bias);

    clock_bias *= 1e6; // convert from seconds to microseconds
    
    std::cout<<std::fixed;
    std::cout<<sat_xyz[0]<<'\t'<<sat_xyz[1]<<'\t'<<sat_xyz[2]<<'\t';
    std::cout<<clock_bias<<std::endl;
  }
#endif 

}
