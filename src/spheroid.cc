
// ---------------------------------------------------------------------------
//  Copyright (C) 2009-2024, All rights reserved. Andre Caceres Carrilho
//
//   spheroid.cc --Spheroid class, Geodetic computations & UTM projection
//   The term 'Spheroid' is what Open Geospatial Consortium (OGC) uses.
// ---------------------------------------------------------------------------

#include "kepler.h"
#include "constants.h"

#define POW2(x) (x*x)

#ifdef DEBUG
static bool check_lat(double lat_deg)
{
  return lat_deg >= -90.0 && lat_deg <= 90.0;
}
static bool check_lon(double lon_deg)
{
  return lon_deg >= -180.0 && lon_deg <= 180.0;
}
#endif 

Spheroid::Spheroid(eSPHEROID ellps)
{
  set(ellps);
}

void Spheroid::set(eSPHEROID ellps)
{
  switch(ellps)
  {
  case WGS66:
    set(6378145.0,1.0/298.25); break;
  case WGS72:
    set(6378135.0,1.0/298.26); break;
  case WGS84:
    set(6378137.0,1.0/298.257223563); break;
  case GRS67:
    set(6378160.0,1.0/298.247167427); break;
  case GRS80:
    set(6378137.0,1.0/298.257222100882711); break;
  case PZ90:
    set(6378136.0,1.0/298.257839303); break;
  case SAD69:
    set(6378160.0,1.0/298.25); break;
  default: // to WGS84
    set(6378137.0,1.0/298.257223563); break;
  }
}
  
// ---------------------------------------------------------------------------
//  f=0 constructs a Sphere
// ---------------------------------------------------------------------------
void Spheroid::set(double a_, double f_)
{
  double n1,n2,n3,n4,n5,n6,n7,n8;
 
#ifdef DEBUG
  if(a_<=0.0)
    error("Invalid equatorial radius");
#endif

  a=a_; // equatorial radius (m)
  f=f_; // flattening
  e2=f*(2.0-f); // eccentricity^2
  e=sqrt(e2); // eccentricity
  b=a*sqrt(1.0-e2); // polar semi axis (m)
  
  // Third flattening (up to the 8th power). Equation (9)
  // This is only needed to compute alpha, beta and the rectifying radius.
  n1=f/(2.0-f);
  n2=POW2(n1);
  n3=n1*n2;
  n4=POW2(n2);
  n5=n2*n3;
  n6=POW2(n3);
  n7=n3*n4;
  n8=POW2(n4);
  
  // Equation (41) (rectifying radius)
  rr=a/(1.0+n1)*(1.0+0.25*n2+0.015625*n4+0.00390625*n6+0.00006103515625*n8);
  
  // Coefficients for Karney-Kruger UTM equations: 
  // Deakin, R.E.; Hunter, M.N.; Karney, C.F.F. A fresh look at the UTM projection.
  
  // Equation (62) Krueger eq, but extended to order n^8
  alpha[0]=(n1*0.5-n2*2/3+n3*5/16+n4*41/180-n5*127/288+
            n6*(7891/37800)+n7*(72161/387072)-n8*(18975107/50803200));
  alpha[1]=(n2*13/48-n3*3/5+n4*557/1440+n5*281/630-
            n6*(1983433/1935360)+n7*(13769/28800)
            +n8*(148003883/174182400));
  alpha[2]=(n3*61/240-n4*103/140+n5*(15061/26880)+n6*(167603/181440)
            -n7*(67102379/29030400)+n8*(79682431/79833600));
  alpha[3]=(n4*(49561/161280)-n5*179/168+n6*(6601661/7257600)+
            n7*(97445/49896)-n8*(40176129013/7664025600));
  alpha[4]=(n5*(34729/80640)-n6*(3418889/1995840) +
            n7*(14644087/9123840)+n8*(2605413599/622702080));
  alpha[5]=(n6*(212378941/319334400)-n7*(30705481/10378368)+
            n8*(175214326799/58118860800));
  alpha[6]=n7*(1522256789/1383782400)-n8*(16759934899/3113510400);
  alpha[7]=n8*(1424729850961/743921418240);
  
  // Equation (64) Krueger eq, but extended to order n^8
  beta[0] =(-n1*0.5+n2*(2/3)-n3*(37/96)+n4*(1/360)+n5*(81/512)
             -n6*(96199/604800)+n7*(5406467/38707200)
             -n8*(7944359/67737600));
  beta[1]=(-n2*1/48-n3*1/15+n4*437/1140-n5*46/105+
             n6*(1118711/3870720)-n7*(51841/1209600)-
             n8*(24749483/348364800));
  beta[2]=(-n3*17/480+n4*37/840+n5*209/4480-
             n6*(5569/90720)-n7*(9261899/58060800)+
             n8*(6457463/17740800));
  beta[3]=(-n4*(4397/161280)+n5*11/504+
             n6*(830251/7257600)-n7*(466511/2494800)-
             n8*(324154477/7664025600));
  beta[4]=(-n5*(4583/161280)+n6*(108847/3991680)+
             n7*(8005831/63866880)-n8*(22894433/124540416));
  beta[5]=(-n6*(20648693/638668800)+n7*(16363163/518918400)+
             n8*(2204645983/12915302400));
  beta[6]=(-n7*(219941297/5535129600)+n8*(497323811/12454041600));
  beta[7]=(-n8*(191773887257/3719607091200));
}

// ---------------------------------------------------------------------------
//  Inverse Geodesic problem using Thadeus Vincenty formula (1975-76).
//  Might fail for antipodal points (or take thousands of iterations). 
//  
//  Args:
//      - (lat1, lon1): Coordinates of first point (decimal degrees)
//      - (lat2, lon2): Coordinates of second point (decimal degrees)
//      
//  Returns:
//      Distance between the two coordinates (in meters). Error < 0.5 mm 
// ---------------------------------------------------------------------------
double Spheroid::geodesic(
  double lat1deg, double lon1deg,
  double lat2deg, double lon2deg) const
{
#ifdef DEBUG
  if(!check_lat(lat1deg)||!check_lat(lat2deg))
    warn("latitude out of range [-90; 90]");
  if(!check_lon(lon1deg)||!check_lon(lon2deg))
    warn("longitude out of range [-180; 180]");
#endif

  int i;
  double dlmb,lmb,lmb0;
  double u1,u2,h0,h1;
  double k0,k1,k2,ss,sc,sg;
  double va,vb,vc,w,gg,ds;
  double alphas,alphac2;
  double s2c,s2c2,clmb,slmb;
  double cu1,cu2,su1,su2;
  
  // reduced latitudes (latitude on the auxiliary sphere)
  u1=atan((1.0-f)*tan(lat1deg*D2R));
  u2=atan((1.0-f)*tan(lat2deg*D2R));
  
  // iterative point longitude difference on the auxiliary sphere 
  dlmb=(lon2deg-lon1deg)*D2R; // delta lambda 
  lmb0=dlmb; // initial approximation 
  lmb=dlmb; // current approximation 
  
  cu1=cos(u1);
  su1=sin(u1);
  cu2=cos(u2);
  su2=sin(u2);
  
  // usually takes less than 25 iterations 
  for(i=0; i<1000; i++){
    clmb=cos(lmb);
    slmb=sin(lmb);
    
    h0=cu2*slmb;
    h1=cu1*su2-su1*cu2*clmb;
    ss=hypot(h0,h1);
    
    // angular separation between points 
    sc=su1*su2+cu1*cu2*clmb;
    sg=atan2(ss,sc);
    
    // forward azimuth of the geodesic at the 
    // equator if it were extended that far 
    alphas=cu1*cu2*slmb/ss;
    alphac2=1.0-POW2(alphas);
    s2c=sc-(2.0*su1*su2/alphac2);
    s2c2=POW2(s2c);
    vc=f*0.0625*alphac2*(4.0+f*(4.0-3.0*alphac2));
    
    // new approximation for lambda 
    lmb0=lmb;
    w=sg+vc*ss*(s2c+vc*sc*(2.0*s2c2-1.0));
    lmb=dlmb+(1.0-vc)*f*alphas*w;
    
    // convergence test 
    if(fabs(lmb0-lmb)<EPS)
      break;
  }

#ifdef DEBUG
  if(i>25)
    warn("more than 25 iterations");
#endif

  k0=a/b;
  k0=alphac2*(POW2(k0)-1.0);
  k2=sqrt(1.0+k0);
  k1=(k2-1.0)/(k2+1.0);

  // simpler formulation given in Vincenty 1976
  va=(1.0+POW2(k1)*0.25)/(1.0-k1);
  vb=k1*(1.0-0.375*POW2(k1));
  gg=vb/6.0*s2c*(4.0*ss*ss-3.0)*(4.0*s2c2-3.0);
  ds=vb*ss*(s2c+vb*0.25*((sc*(2.0*s2c2-1.0))-gg));
  
  // distance in meters 
  return b*va*(sg-ds);
}

// ---------------------------------------------------------------------------
//  Conversion between Geodetic and Geocentric Cartesian coordinates,
//  which is also known as Earth Centered, Earth Fixed (ECEF).
// 
//  Args:
//      (lat,lon,h) Coordinates of the point in decimal degrees.
//      The height must be given in meters.
// 
//  Returns:
//      (X,Y,Z) - ECEF coordinates (in meters).
// ---------------------------------------------------------------------------
void Spheroid::geo2ecf(const double *geo, double *xyz) const
{
#ifdef DEBUG
  if(!geo||!xyz)
    error("null pointers");
  if(!check_lat(geo[0]))
    warn("latitude out of range [-90; 90]");
  if(!check_lon(geo[1]))
    warn("longitude out of range [-180; 180]");
#endif
  double phi=geo[0]*D2R;
  double lmb=geo[1]*D2R;
  double cf=cos(phi);
  double sf=sin(phi);
  //radius of curvature in prime vertical
  double n=a/sqrt(1.0-e2*POW2(sf));
  double nh=n+geo[2];
  //geocentric (ecef)
  xyz[0]=nh*cf*cos(lmb);
  xyz[1]=nh*cf*sin(lmb);
  xyz[2]=(n*(1.0-e2)+geo[2])*sf;
}

// ---------------------------------------------------------------------------
//  Conversion from ECEF  to Geodetic coordinates.
//
//   Args:
//       (X,Y,Z) - ECEF coordinates (in meters).
//   
//   Returns:
//       (lat,lon,h) Coordinates of the point in decimal degrees.
//       The height is given in meters.
//   
//   Note:
//       Adapted from RTKlib.
// ---------------------------------------------------------------------------
void Spheroid::ecf2geo(const double *xyz, double *geo) const
{
#ifdef DEBUG
  if(!xyz||!geo)
    error("null pointers");
  if(pythag(xyz[0],xyz[1],xyz[2])<EPS)
    warn("invalid ECEF coordinates");
#endif
  int i;
  double rho,n;
  double zi,zk,sf;
  double phi,lmb;
  n=a;zi=xyz[2];zk=0;
  rho=hypot(xyz[0],xyz[1]);
  for(i=0; i<25; i++){
    zk=zi;
    sf=zi/hypot(rho,zi);
    n=a/sqrt(1.0-e2*POW2(sf));
    zi=xyz[2]+n*e2*sf;
    if(fabs(zi-zk)<EPS)
      break;
  }
#ifdef DEBUG
  if(i>10)
    warn("more than 10 iterations");
#endif
  lmb=atan2(xyz[1],xyz[0]);
  if(rho>EPS)
    phi=atan(zi/rho);
  else 
    phi=xyz[2]<0.0?-PI2:PI2;
  geo[0]=phi*R2D;
  geo[1]=lmb*R2D;
  geo[2]=hypot(rho,zi)-n;
}

// ---------------------------------------------------------------------------
//  Forward UTM projection
//    
//  Args:
//      (lat,lon) Coordinates of the point in decimal degrees.
//  
//  Returns:
//      (E,N) Coordinates in UTM projection in meters; 
//      UTM Zone; and
//      Hemisphere ('N' or 'S').
// ---------------------------------------------------------------------------
void Spheroid::geo2utm(const double *geo, double *utm, int *zone, char *h) const
{
#ifdef DEBUG
  if(!geo||!utm)
    error("null pointers");
  if(!check_lat(geo[0]))
    warn("latitude out of range [-90; 90]");
  if(!check_lon(geo[1]))
    warn("longitude out of range [-180; 180]");
#endif

  int i;
  double lmb0;
  double x,y,z,w,u,v;
  utmzone(geo[1],zone,&lmb0);// zone is only written if not null
  w=(geo[1]-lmb0)*D2R;// omega
  z=cnflat(geo[0]*D2R);// conformal latitude
  u=atan(z/cos(w));
  v=asinh(sin(w)/hypot(z,cos(w)));
  for(i=0,x=0,y=0,w=2.; i<8; i++,w+=2.0){
    x+=alpha[i]*cos(w*u)*sinh(w*v);
    y+=alpha[i]*sin(w*u)*cosh(w*v);
  }
  x=0.9996*rr*(v+x);
  y=0.9996*rr*(u+y);
  utm[0]=x+5E+05;
  utm[1]=geo[0]>0.0?y:y+1E+07;
  if(h)
    *h=geo[0]>0.0?'N':'S';
}

// ---------------------------------------------------------------------------
//  Backward UTM projection
//    
//  Args:
//      - (E,N) Coordinates in UTM projection in meters; 
//      - UTM Zone;
//      - Hemisphere ('N','n' or 'S','s').
//  
//  Returns:
//      - (lat,lon) Coordinates of the point in decimal degrees.
// ---------------------------------------------------------------------------
void Spheroid::utm2geo(const double *utm, double *geo, int zone, char h) const
{
#ifdef DEBUG
  if(!utm||!geo)
    error("null pointers");
  if(utm[0]<0.0||utm[1]<0.0)
    warn("negative UTM coordinates");
  if(zone<1||zone>60)
    warn("UTM zone out of range [1;60]");
  if(h!='N'&&h!='n'&&h!='S'&&h!='s')
    warn("invalid hemisphere label [N,S]");
#endif

  int i;
  double x,y,w,u,v;
  double ta,ti,re,rt,sg;
  double rs,tf,dr,nt;
  double lmb0;
  x=utm[0]-5E+05;
  y=(h=='N')||(h=='n')?utm[1]:utm[1]-1E+07;
  x*=1.0/(rr*0.9996);
  y*=1.0/(rr*0.9996);
  for(i=0,u=0,v=0,w=2.0; i<8; i++,w+=2.0){
    v+=beta[i]*cos(w*y)*sinh(w*x);
    u+=beta[i]*sin(w*y)*cosh(w*x);
  }
  v+=x;
  u+=y;
  w=atan(sinh(v)/cos(u));
  ta=sin(u)/hypot(sinh(v),cos(u));
  ti=ta;
  re=1.0-e2;
  for(i=0;i<25;i++){
    rt=sqrt(1.0+POW2(ta));
    sg=sinh(e*atanh(e*ta/rt));
    rs=sqrt(1.0+POW2(sg));
    tf=ta*rs-sg*rt-ti;
    dr=(rs*rt-sg*ta)*(re*rt)/(1.0+re*POW2(ta));
    nt=ta-tf/dr;
    dr=fabs(ta-nt);
    ta=nt;
    if(dr<EPS)
      break;
  }
#ifdef DEBUG
  if(i>10)
    warn("more than 10 iterations");
#endif
  lmb0=(zone-1)*6-177;
  geo[0]=R2D*atan(ta);
  geo[1]=R2D*w+lmb0;
}

// ---------------------------------------------------------------------------
//  Conformal latitude
//  Equations 88 & 89, used on the Forward projection (geo2utm)
// ---------------------------------------------------------------------------
double Spheroid::cnflat(double phi) const
{
  double x,y,z,w,r,s;
  x=tan(phi); y=x*x;
  z=e*x/sqrt(1.0+y);
  w=sinh(e*atanh(z));
  r=x*sqrt(1.0+w*w);
  s=w*sqrt(1.0+y);
  return r-s;
}

// ---------------------------------------------------------------------------
//  Computes UTM Zone and Central Meridian longitude.
//
//  Args:
//      - lon: Longitude of point in decimal degrees.
// 
//  Returns (if ptr is not null):
//      - UTM Zone (int) (1-60)
//      - Central Meridian longitude (in decimal degrees).
// ---------------------------------------------------------------------------
void Spheroid::utmzone(double lon_deg, int *zone, double *cm_deg)
{
#ifdef DEBUG
  if(!check_lon(lon_deg))
    warn("longitude out of range [-180; 180]");
#endif

  int z,w;
  z=1+ (int)((lon_deg+180)/6);
  w=(z-1)*6-177;
  if(zone) 
    *zone=z;
  if(cm_deg) // central meridian (deg)
    *cm_deg=(double)w;
}

// ---------------------------------------------------------------------------
//  UTM Scale Factor 
//  abs(err) < 5.923E-07 vs complete model
// ---------------------------------------------------------------------------
double Spheroid::utmscale(double lon_deg)
{
  double cm,w;
  utmzone(lon_deg,0,&cm);
  w=(cm-lon_deg)*D2R;// omega
  return 0.9996+0.5035348161*w*w; 
}
