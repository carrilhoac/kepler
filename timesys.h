
#ifndef KEPLER_TIMESYS_h
#define KEPLER_TIMESYS_h 1

template<typename FT> 
Time::Time(FT t)
{
  (*this)=t;
}

template<typename FT> 
Time& Time::operator=(FT t)
{
  double td=double(t);
  double tf=floor(td); // or trunc?
  t_sec=int64_t(tf);
  t_frac=td-tf;
  return *this;
}

template<typename FT> 
Time& Time::operator+=(FT t)
{
  Time tt(t);
  return (*this)+=tt;
}

template<typename FT> 
Time& Time::operator-=(FT t)
{
  Time tt(t);
  return (*this)-=tt;
}

template<typename FT> 
Time Time::operator+(FT t) const
{
  Time tt(*this);
  tt+=t;
  return tt;
}

template<typename FT> 
Time Time::operator-(FT t) const
{
  Time tt(*this);
  tt-=t;
  return tt;
}

#endif 
