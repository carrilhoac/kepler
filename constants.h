
#ifndef KEPLER_constants_h
#define KEPLER_constants_h

#define CLIGHT   299792458.0      // speed of light (m/s)
#define PI       3.1415926535897932384626
#define PI2      1.5707963267948966192313
#define D2R      0.017453292519943295 // pi/180
#define R2D      57.29577951308232087 // 180/pi
#define EPS      1e-12            // machine epsilon IEEE 64-bit
#define MU_GPS   3.9860050E14     // gravitational constant         ref [1] 
#define MU_GLO   3.9860044E14     // gravitational constant         ref [2] 
#define MU_GAL   3.986004418E14   // earth gravitational constant   ref [7] 
#define MU_CMP   3.986004418E14   // earth gravitational constant   ref [9] 
#define OMGE_GLO 7.292115E-5      // earth angular velocity (rad/s) ref [2] 
#define OMGE_GAL 7.2921151467E-5  // earth angular velocity (rad/s) ref [7] 
#define OMGE_CMP 7.292115E-5      // earth angular velocity (rad/s) ref [9] 
#define OMGE_GPS 7.2921151467E-5  // earth angular velocity (rad/s) (IS-GPS)

#endif 
