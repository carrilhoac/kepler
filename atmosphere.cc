
#include "kepler.h"
#include "constants.h"

//////////////////////////////////////////////////////////////////////
// Utility

static void mat_mulv3(double *x, const double *A, const double *b)
{
  x[0]=A[0]*b[0]+A[1]*b[1]+A[2]*b[2];
  x[1]=A[3]*b[0]+A[4]*b[1]+A[5]*b[2];
  x[2]=A[6]*b[0]+A[7]*b[1]+A[8]*b[2];
}

static void xyz2enu(const double *geo, double *E)
{
  double sf,cf,sl,cl;

  sf=sin(geo[0]); cf=cos(geo[0]);
  sl=sin(geo[1]); cl=cos(geo[1]);
  E[0]=-sl;    E[1]=cl;     E[2]=0.0; // row-major
  E[3]=-sf*cl; E[4]=-sf*sl; E[5]=cf;
  E[6]=cf*cl;  E[7]=cf*sl;  E[8]=sf;
}

static void enu2xyz(const double *geo, double *E)
{
  double sf,cf,sl,cl;

  sf=sin(geo[0]); cf=cos(geo[0]);
  sl=sin(geo[1]); cl=cos(geo[1]);
  E[0]=-sl; E[1]=-sf*cl; E[2]=cf*cl; 
  E[3]=cl;  E[4]=-sf*sl; E[5]=cf*sl;
  E[6]=0.0; E[7]=cf;     E[8]=sf;
}

void ecf2enu(const double *geo, const double *xyz, double *enu) 
{
  double E[9];
  xyz2enu(geo,E);
  mat_mulv3(enu,E,xyz);
}

void enu2ecf(const double *geo, const double *enu, double *xyz) 
{
  double E[9];
  enu2xyz(geo,E);
  mat_mulv3(xyz,E,enu);
}

//los: line of sight
double geomdist(const double *sat, const double *rec, double *los)
{
  double r,s;

  los[0]=sat[0]-rec[0];
  los[1]=sat[1]-rec[1];
  los[2]=sat[2]-rec[2];

  r=pythag(los[0],los[1],los[2]);
  
#if DEBUG
  if(fabs(r)<MATRIX_EPS)
    warn("division by zero");
#endif 
  s=1.0/r;
  los[0]*=s;
  los[1]*=s;
  los[2]*=s;
  
  s=sat[0]*rec[1]-sat[1]*rec[0];
  return r+OMGE_GPS*s/CLIGHT;
}

double satazel(const double *geo, const double *los, double *azel)
{
  double az,el,enu[3]={0};

  ecf2enu(geo,los,enu);
  az=atan2(enu[0],enu[1]);
  el=asin(enu[2]);
  if(az<0.0)
    az+=2*PI;

  if(azel) {
    azel[0]=az;
    azel[1]=el;
  }
  return el;
}

// Saastamoinen
double tropmod(const double *geo, const double *azel, double humi)
{
  const double temp0=15.0; // temperature at sea level (Celsius)
  double hgt,pres,temp,e,z,trph,trpw;

  // standard atmosphere
  hgt=geo[2]<0.0?0.0:geo[2];

  pres=1013.25*pow(1.0-2.2557E-5*hgt,5.2568);
  temp=temp0-6.5E-3*hgt+273.16;
  e=6.108*humi*exp((17.15*temp-4684.0)/(temp-38.45));

  // saastamoinen model
  z=PI/2.0-azel[1];
  trph=0.0022768*pres/(1.0-0.00266*cos(2.0*geo[0])-0.00028*hgt/1E3)/cos(z);
  trpw=0.002277*(1255.0/temp+0.05)*e/cos(z);
  return trph+trpw;
}

// 2004/1/1
static const double ion_default[8] = 
{ 
  0.1118E-07,-0.7451E-08,-0.5961E-07, 0.1192E-06,
  0.1167E+06,-0.2294E+06,-0.1311E+06, 0.1049E+07
};
  
#if 0
// Move this to a class
// Klobuchar
double ionmod(const Time& t, const double *ion, const double *geo, 
  const double *azel)
{
  double tt,tow,f,psi,phi,lam,amp,per,x;
  int week;

  //if (geo[2]<-1E3||azel[1]<=0) return 0.0;
  //if (norm(ion,8)<=0.0) ion=ion_default;
  if(!ion)
    ion=ion_default;

  // earth centered angle (semi-circle)
  psi=0.0137/(azel[1]/PI+0.11)-0.022;

  // subionospheric latitude/longitude (semi-circle)
  phi=geo[0]/PI+psi*cos(azel[0]);
  if (phi> 0.416) 
    phi= 0.416;
  else if (phi<-0.416) 
    phi=-0.416;
  lam=geo[1]/PI+psi*sin(azel[0])/cos(phi*PI);

  // geomagnetic latitude (semi-circle)
  phi+=0.064*cos((lam-1.617)*PI);

  // local time (s)
  ///tow=unx2gps(t,&week);
  tt=43200.0*lam+tow;//time2gpst(t,&week);
  tt-=floor(tt/86400.0)*86400.0; // 0<=tt<86400

  // slant factor
  f=1.0+16.0*pow(0.53-azel[1]/PI,3.0);

  // ionospheric delay
  amp=ion[0]+phi*(ion[1]+phi*(ion[2]+phi*ion[3]));
  per=ion[4]+phi*(ion[5]+phi*(ion[6]+phi*ion[7]));
  amp=amp<    0.0?    0.0:amp;
  per=per<72000.0?72000.0:per;
  x=2.0*PI*(tt-50400.0)/per;

  return CLIGHT*f*(fabs(x)<1.57?5E-9+amp*(1.0+x*x*(-0.5+x*x/24.0)):5E-9);
}
#endif 

