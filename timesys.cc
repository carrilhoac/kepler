
#include "kepler.h"

#define SEC_MIN   60
#define SEC_HOUR  3600
#define SEC_DAY   86400
#define SEC_WEEK  604800

static const int gps0[6] = {1980,1,6,0,0,0}; // GPS ref epoch

Time::Time()
{}

Time::Time(const Time& t)
  : t_sec(t.t_sec)
  , t_frac(t.t_frac)
{}

Time& Time::from_cal(const int *cal)
{
  t_sec=Time::cal2unx(cal);
  t_frac=0.0;
  return *this;
}

// offsets into RINEX (version 3 and later) epoch string with 
// the follwing format: YYYY MM DD HH mm SS.sssssss
static const int r_offs[6] = {0,5,8,11,14,17};

// This function assumes standard locale: setlocale(LC_ALL,"C"); 
Time& Time::from_rnx(const char *rnx)
{
  int cal[6];
  double sec;
  
  for(int i=0;i<6;i++)
    cal[i]=strtol(rnx+r_offs[i],0,10);
  sec=strtod(rnx+r_offs[5],0);
  
  t_sec=Time::cal2unx(cal);
  t_frac=sec-floor(sec);
  return *this;
}

Time& Time::from_gps(int gpsweek, double gpstow)
{
  int64_t t0;

  t0=cal2unx(gps0);
  t_sec=t0+gpsweek*SEC_WEEK+(int64_t)(gpstow);
  t_frac=gpstow-floor(gpstow);
  return *this;
}

int Time::gps_week() const
{
  int64_t t0,dt;

  t0=Time::cal2unx(gps0);
  dt=t_sec-t0;
  return dt/SEC_WEEK;
}

double Time::gps_tow() const
{
  double s;
  int64_t t0,dt,w;

  t0=Time::cal2unx(gps0);
  dt=t_sec-t0;
  w=dt/SEC_WEEK;
  s=dt-w*SEC_WEEK+t_frac;
  return s;
}

// might have loss of precision if t_sec >= 1e8
// but it is fine otherwise (gps tow, for instance)
double Time::to_double() const
{
  return double(t_sec)+t_frac;
}

Time& Time::operator=(const Time& t)
{
  t_sec=t.t_sec;
  t_frac=t.t_frac;  
  return *this;
}

Time& Time::operator+=(const Time& t)
{
  t_sec+=t.t_sec;
  t_frac+=t.t_frac;
  normalize();
  return *this;
}

Time& Time::operator-=(const Time& t)
{
  t_sec-=t.t_sec;
  t_frac-=t.t_frac;
  normalize();
  return *this;
}

Time Time::operator+(const Time& t) const
{
  Time tt(*this);
  tt+=t;
  return tt;
}

Time Time::operator-(const Time& t) const
{
  Time tt(*this);
  tt-=t;
  return tt;
}

// Adapted from Howard Hinnant paper 
int Time::civ2day(int y, int m, int d)
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

int64_t Time::cal2unx(const int *cal)
{
  int64_t days,secs;

  days=Time::civ2day(cal[0],cal[1],cal[2]);
  secs=cal[3]*SEC_HOUR+cal[4]*SEC_MIN+cal[5];
  return days*SEC_DAY+secs;
}

std::ostream& operator<<(std::ostream& os, const Time& t)
{
  os<<t.to_double();
  return os;
}

void Time::normalize()
{
}
