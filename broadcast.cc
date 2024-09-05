
#include "kepler.h"
#include <cstring>

#define MU_GPS   3.9860050E14     // gravitational constant         ref [1] 
#define MU_GLO   3.9860044E14     // gravitational constant         ref [2] 
#define MU_GAL   3.986004418E14   // earth gravitational constant   ref [7] 
#define MU_CMP   3.986004418E14   // earth gravitational constant   ref [9] 
#define OMGE_GLO 7.292115E-5      // earth angular velocity (rad/s) ref [2] 
#define OMGE_GAL 7.2921151467E-5  // earth angular velocity (rad/s) ref [7] 
#define OMGE_CMP 7.292115E-5      // earth angular velocity (rad/s) ref [9] 
#define OMGE_GPS 7.2921151467E-5  // earth angular velocity (rad/s) (IS-GPS)

#define CLIGHT      299792458.0   // speed of light (m/s)

#define POW2(x) (x*x)

// ura values (ref [3] 20.3.3.3.1.1)
static const double ura_eph[16] =
{
  2.4   ,3.4    , 4.85   , 6.85 ,
  9.65  ,13.65  , 24.0   , 48.0 ,
  96.0  ,192.0  , 384.0  , 768.0,
  1536.0,3072.0 , 6144.0 , 0.0
};

// ura value (m) to ura index 
static int uraindex(double value)
{
  for(int i=0;i<15;i++) 
    if(ura_eph[i]>=value) 
      return i;
  return 0;
}

// FORTRAN to C double exponent 
// We need to replace the exponent character from 'D' to 'E'  
// due to legacy FORTRAN compatibility on the RINEX format 
// specification (for writing double precision floating point) 
static void fixflt(char *s)
{
  for(;*s;s++)
    *s=(*s=='D'||*s=='d')?'E':*s;
}

// adjust time considering week handover 
static void adjweek(Time& t, const Time& t0)
{
  double tt=(t-t0).to_double();
  if(tt<-302400.0)
    t.t_sec+=604800;
  if(tt> 302400.0)
    t.t_sec-=604800;
}

Nav::Nav()
{}

Nav::Nav(const char *rnx)
{
  rnx2nav(rnx);
}

#ifdef DEBUG
void check_prn(const char *prn)
{
  int satnum;
  
  satnum=strtol(&prn[1],0,10);
  switch(prn[0]){
  case 'G': // GPS
    if(satnum<1||satnum>32)
      error("invalid GPS satellite PRN");
  break; 
  //case 'E': // Galileo
  //break;
  //case 'C':
  //break;
  //case 'I':
  //break;
  default:
    error("invalid satellite system");
  break;
  }
}
#endif

void Nav::rnx2nav(const char *rnx)
{
  int i,j,k,n,w;
  const char *lines[8];
  char buf[20]={0};
  double par[32]={0};

  // parsing lines
  lines[0]=rnx;
  for(i=1,j=1;i<8;j++){ 
#ifdef DEBUG
    if(!rnx[j]) // we expect eight lines in the nav record
      error("incomplete GPS Nav record");
#endif 
    if(rnx[j]=='\n' // handles both Win32 & Linux
    ||(rnx[j]=='\r'&&rnx[j+1]!='\n')){// for MacOS
      lines[i++]=rnx+(++j); // new line
    }
  }

  // some RINEX files might omit the leading zero on the
  // satellite PRN, we standardize it here for later checks
  prn[0]=lines[0][0]; 
  prn[1]=lines[0][1]==' '?'0':lines[0][1]; 
  prn[2]=lines[0][2];
  prn[3]=0;
  
#ifdef DEBUG
  check_prn(prn);
#endif

  // reading parameters
  for(k=0,i=0;i<8;i++){
    n=strlen(lines[i]);    
    for(j=0;j<4;j++){
      w=4+19*j;      
      if(w+19>n){
      // some entries might not have the entire line populated
      // thus, we need to check the line length before reading 
        k++; // skip parameter
        continue; 
      }      
      strncpy(buf,&lines[i][w],19);      
      if(!i&&!j){ // first parameter is the time of clock (toc)
        toc.from_rnx(buf);
      } else {
        fixflt(buf);
        par[k++]=strtod(buf,0);       
      }      
    }
  }

  f0=par[0]; 
  f1=par[1]; 
  f2=par[2];

  A=POW2(par[10]);e=par[ 8]; i0  =par[15]; OMG0=par[13]; 
  omg =par[17]; M0 =par[ 6]; deln=par[ 5]; OMGd=par[18]; 
  idot=par[19]; crc=par[16]; crs =par[ 4]; cuc =par[ 7];
  cus =par[ 9]; cic=par[12]; cis =par[14];

  iode=(int)par[ 3];
  iodc=(int)par[26];
  toes=     par[11];
  week=(int)par[21];
  
  toe.from_gps(week, par[11]);
  ttr.from_gps(week, par[27]);
  adjweek(toe,toc);
  adjweek(ttr,toc);

  code=(int)par[20];
  svh =(int)par[24];
  sva =uraindex(par[23]);
  flag=(int)par[22];
  tgd[0]=   par[25];
  fit =     par[28];
}

void Nav::nav2ecf(const Time& t, double *xyz, double *clock_bias) const
{
  int itr;
  double tk,M,E,Ek,sinE,cosE,Q,sinQ,cosQ;
  double u,r,i,x,y,dts,cosi,sin2u,cos2u;
  double mu,omge;

  mu=MU_GPS;
  omge=OMGE_GPS;

  tk=(t-toe).to_double();
  M=M0+(sqrt(mu/(A*A*A))+deln)*tk;

  for(E=M,Ek=0, itr=0;itr<30 && fabs(E-Ek)<1e-14;itr++){
    Ek=E; E-=(E-e*sin(E)-M)/(1.0-e*cos(E));
  }

  sinE=sin(E);
  cosE=cos(E);
  u=atan2(sqrt(1.0-e*e)*sinE,cosE-e)+omg;
  r=A*(1.0-e*cosE);
  i=i0+idot*tk;
  sin2u=sin(2.0*u); 
  cos2u=cos(2.0*u);
  u+=cus*sin2u+cuc*cos2u;
  r+=crs*sin2u+crc*cos2u;
  i+=cis*sin2u+cic*cos2u;
  x=r*cos(u); 
  y=r*sin(u); 
  cosi=cos(i);

  // GPS
  Q=OMG0+(OMGd-omge)*tk-omge*toes;
  sinQ=sin(Q); 
  cosQ=cos(Q);
  xyz[0]=x*cosQ-y*cosi*sinQ;
  xyz[1]=x*sinQ+y*cosi*cosQ;
  xyz[2]=y*sin(i);

  // relativity correction
  tk=(t-toc).to_double();
  dts=f0+f1*tk+f2*tk*tk;
  dts-=2.0*sqrt(mu*A)*e*sinE/POW2(CLIGHT);

  if(clock_bias)
    *clock_bias=dts;
}
