
#include "kepler.h"

#define SEC_MIN   60
#define SEC_HOUR  3600
#define SEC_DAY   86400
#define SEC_WEEK  604800

static const int gps0[6] = {1980,1,6,0,0,0}; // GPS ref epoch

// Adapted from Howard Hinnant paper 
static int civ2day(int y, int m, int d)
{
  int era,yoe,doy,doe;

  y-=m<=2;
  era=(y>=0?y:y-399)/400;
  yoe=y-era*400;
  doy=(153*(m>2?m-3:m+9)+2)/5+d-1;
  doe=yoe*365+yoe/4-yoe/100+doy;

  // 146097 days in 400 Gregorian years (era) 
  // 719468 days from 0000-03-01 to 1970-01-01 (Unix ref epoch) 
  return era*146097+doe-719468;
}

static int64_t cal2unx(const int *cal)
{
  int64_t days,secs;

  days=civ2day(cal[0],cal[1],cal[2]);
  secs=cal[3]*SEC_HOUR+cal[4]*SEC_MIN+cal[5];
  return days*SEC_DAY+secs;
}

// offsets into RINEX (version 3 and later) epoch string with 
// the follwing format: YYYY MM DD HH mm SS.sssssss
static const int r_offs[6] = {0,5,8,11,14,17};

// This function assumes standard locale: setlocale(LC_ALL,"C"); 
extern Epoch rnx2unx(const char *rnx)
{  
  int i,cal[6];
  double sec;
  Epoch r;

  for(i=0;i<6;i++)
    cal[i]=strtol(rnx+r_offs[i],0,10);
  sec=strtod(rnx+r_offs[5],0);
  
  r.t_sec=cal2unx(cal);
  r.t_frac=sec-floor(sec);
  return r;
}

extern double unx2gps(const Epoch& utc_ts, int *gps_week)
{
  int w;
  double s;
  int64_t t0, dt;

  t0=cal2unx(gps0);
  dt=utc_ts.t_sec-t0;
  w=dt/SEC_WEEK;
  s=dt-w*SEC_WEEK+utc_ts.t_frac;

  if(gps_week)
    *gps_week=w;
  return s;
}

extern Epoch gps2unx(int gps_week, double gps_tow)
{
  int64_t t0;
  Epoch r;

  t0=cal2unx(gps0);
  r.t_sec=t0+gps_week*SEC_WEEK+(int64_t)(gps_tow);
  r.t_frac=gps_tow-floor(gps_tow);
  return r;
}

extern double timesub(const Epoch& a, const Epoch& b)
{
  return (double)(a.t_sec-b.t_sec)+(a.t_frac-b.t_frac);
}
