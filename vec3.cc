
#include "kepler.h"

Vec3::Vec3(double X, double Y, double Z)
#ifdef SIMD_x86
{ 
  r = _mm256_setr_pd(X, Y, Z, 0.0); 
}
#else 
{
  v[0]=X;
  v[1]=Y;
  v[2]=Z;
}
#endif 

Vec3::Vec3(const Vec3 &a)
#ifdef SIMD_x86
{ 
  r = a.r; 
}
#else 
{
  v[0]=a.v[0];
  v[1]=a.v[1];
  v[2]=a.v[2];
}
#endif 

Vec3& Vec3::operator=(const Vec3& a) 
{
#ifdef SIMD_x86
  r = a.r;
#else 
  v[0]=a.v[0];
  v[1]=a.v[1];
  v[2]=a.v[2];
#endif
  return *this;
}

Vec3& Vec3::operator+=(const Vec3& a) 
{
#ifdef SIMD_x86
  r = _mm256_add_pd(r, a.r);
#else
  v[0]+=a.v[0];
  v[1]+=a.v[1];
  v[2]+=a.v[2];
#endif
  return *this;
}

Vec3& Vec3::operator-=(const Vec3& a) 
{
#ifdef SIMD_x86
  r = _mm256_sub_pd(r, a.r);
#else
  v[0]-=a.v[0];
  v[1]-=a.v[1];
  v[2]-=a.v[2];
#endif
  return *this;
}

///  Element-wise product 
///  a Vector with elements to be multiplied
Vec3& Vec3::operator*=(const Vec3& a) 
{
#ifdef SIMD_x86
  r = _mm256_mul_pd(r, a.r);
#else
  v[0]*=a.v[0];
  v[1]*=a.v[1];
  v[2]*=a.v[2];
#endif
  return *this;
}

Vec3& Vec3::operator*=(double s) 
{
#ifdef SIMD_x86
  r = _mm256_mul_pd(r, _mm256_set1_pd(s));
#else 
  v[0]*=s;
  v[1]*=s;
  v[2]*=s;
#endif
  return *this;
}
Vec3& Vec3::operator/=(double s) 
{
#if DEBUG
  if(fabs(s)<1e-11)
    warn("division by zero");
#endif 
  s=1.0/s;
#ifdef SIMD_x86
  r = _mm256_mul_pd(r, _mm256_set1_pd(s));
#else 
  v[0]*=s;
  v[1]*=s;
  v[2]*=s;
#endif
  return *this;
}

Vec3 Vec3::operator+(const Vec3& a) const
{
  Vec3 r(*this);
  r+=a;
  return r;
}
Vec3 Vec3::operator-(const Vec3& a) const
{
  Vec3 r(*this);
  r-=a;
  return r;
}
Vec3 Vec3::operator*(const Vec3& a) const
{
  Vec3 r(*this);
  r*=a;
  return r;
}
Vec3 Vec3::operator*(double s) const
{
  Vec3 r(*this);
  r*=s;
  return r;
}
Vec3 Vec3::operator/(double s) const
{
  Vec3 r(*this);
  r/=s;
  return r;
}

///  Computes the 3D norm of the Vector
///  \f$ sqrt(x*x + y*y + z*z) \f$
double Vec3::norm() const 
{
#if 1
  return pythag(v[0],v[1],v[2]);
#else 
  return sqrt(
    v[0]*v[0]+ 
    v[1]*v[1]+ 
    v[2]*v[2]);
#endif 
}

///  Sets 'this' Vector length to one, i.e. 
/// \f$ sqrt(x*x + y*y + z*z) = 1 \f$
Vec3& Vec3::unit()
{
  double n=norm();
#if DEBUG
  if(fabs(n)<0.0001)
    warn("division by zero");
#endif 
  (*this)*=1.0/n;
  return *this;
}

///  Euclidian distance between \p a and \p b 
///  a  3D point
///  b  3D point
///  \f$ sqrt( dx^2 + dy^2 + dz^2 ) \f$
///  only (X,Y,Z) are used in the formula!
double dist(const Vec3 &a, const Vec3 &b)
{
  Vec3 v(a);
  v -= b;
  return sqrt(dot(v,v));
}

///  Dot Product
///  ax*bx + ay*by + az*bz
///  only (X,Y,Z) are used in the formula!
double dot(const Vec3 &a, const Vec3 &b)
{
  return 
    a.v[0]*b.v[0]+ 
    a.v[1]*b.v[1]+ 
    a.v[2]*b.v[2];
}

///  Cross Product
///  \f$ c = a * b \f$
Vec3 cross(const Vec3 &a, const Vec3 &b)
{
  Vec3 c;
  cross(c,a,b);
  return c;
}

///  Cross Product
///  \f$ c = a * b \f$
void cross(Vec3 &c, const Vec3 &a, const Vec3 &b)
{
#ifdef SIMD_x86
  c.r = _mm256_permute4x64_pd(
    _mm256_sub_pd(
			_mm256_mul_pd(a.r, _mm256_permute4x64_pd(b.r, _MM_SHUFFLE(3, 0, 2, 1))),
			_mm256_mul_pd(b.r, _mm256_permute4x64_pd(a.r, _MM_SHUFFLE(3, 0, 2, 1)))
    ), _MM_SHUFFLE(3, 0, 2, 1)
  );
#else 
  c.v[0]=a.v[1]*b.v[2]-a.v[2]*b.v[1];
  c.v[1]=a.v[2]*b.v[0]-a.v[0]*b.v[2];
  c.v[2]=a.v[0]*b.v[1]-a.v[1]*b.v[0];
#endif 
}

std::ostream& operator<<(std::ostream& os, const Vec3& b)
{
  os<<std::setw(14)<<b.v[0];
  os<<std::setw(14)<<b.v[1];
  os<<std::setw(14)<<b.v[2];
  return os;
}
  