
// ---------------------------------------------------------------------------
//  Copyright (C) 2009-2024, All rights reserved. Andre Caceres Carrilho
//
//   matrix.cc --Mat class, Linear Algebra and Matrix computations
// ---------------------------------------------------------------------------

#include "kepler.h"

void error_(const char *from, const char *msg)
{
  std::cerr<<"[ERROR] "<<from<<": "<<msg<<std::endl;
  exit(1);
}

void warn_(const char *from, const char *msg)
{
  std::cerr<<"[WARNING] "<<from<<": "<<msg<<std::endl;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

#if 0
// Cholesky factorization 
static void chol(double *A, int n)
{
  int i,j,k;
  double x, r;
  for(j=0;j<n;j++){
    x=A[j*n+j];
    for(k=0;k<j;k++)
      x-=A[j*n+k]*A[j*n+k];
#ifdef DEBUG
    if(x<0.0){
      warn("matrix is not S.P.D.");
      return;
    }
#endif 
    x=sqrt(x);
    A[j*n+j]=x;
    r=1.0/x;
    for(i=j+1;i<n;i++){
      x=A[i*n+j];
      for(k=0;k<j;k++)
        x-=A[i*n+k]*A[j*n+k];
      A[i*n+j]=x*r;
    }
  }
}

// Cholesky back-substitution
static void chbksb(const double *L, int n, double *x, int j)
{
  int i,k;

  // Solve L*y = b 
  for(k=j;k<n;++k){
    for(i=j;i<k;++i)
      x[k]-=x[i]*L[k*n+i];
    x[k]/=L[k*n+k];
  }

  // Solve L'*X = Y 
  for(k=n-1;k>=j;--k){
    for(i=k+1;i<n;++i)
      x[k]-=x[i]*L[i*n+k];
    x[k]/=L[k*n+k];
  }
}
#endif 

//LU decomposition
static double ludcmp(double *A, int m, int n, int *piv)
{
  int i,j,k,p,q;
  double s,t,d,w;
  w=1.0;//det
  q=m<n?m:n;
  for(i=0;i<m;i++) 
    piv[i]=i;//initialize array for partial pivoting
  for(i=0;i<q;i++){
    p=i;//find pivot
    s=fabs(A[i*n+i]);
    for(j=i+1;j<m;j++){
      t=fabs(A[j*n+i]);
      if(t>s){
        p=j;
        s=t;
      }
    }
    if(p!=i){ //swap if necessary
      w=-w;//swap det sign
      k=piv[i];//swap pivot indexes
      piv[i]=piv[p];
      piv[p]=k;      
      for(j=0;j<n;j++){//swap rows
        t=A[i*n+j];
        A[i*n+j]=A[p*n+j];
        A[p*n+j]=t;
      }
    }
    if(i<m){//Gaussian elimination
      d=1.0/A[i*n+i];
      for(j=i+1;j<m;j++)
        A[j*n+i]*=d;
    }
    if(i<(q-1)){
      for(j=i+1;j<m;j++){
        d=A[j*n+i];
        for(k=i+1;k<n;k++)
          A[j*n+k]-=d*A[i*n+k];
      }
    }
  }
  for(i=0;i<n;i++)
    w*=A[i*n+i];//computing the determinant
  return w;
}

//LU back-substitution (inverse matrix only)
static void lubksb(double *LU, int n, double *Bt, int *piv)
{
  int i,j,k;
  double b;
  for(k=0;k<n;k++){
    for(i=0;i<n;i++){
      b=piv[i]==k?1.0:0.0;// this is genius, btw
      for(j=0;j<i;j++)
        b-=LU[i*n+j]*Bt[k*n+j];
      Bt[k*n+i]=b;
    }
    for(i=n-1;i>=0;i--){
      b=Bt[k*n+i];
      for(j=i+1;j<n;j++)
        b-=LU[i*n+j]*Bt[k*n+j];
      Bt[k*n+i]=b/LU[i*n+i];
    }
  }
}


// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

#define EPS 1e-12

Mat::Mat(int m_, int n_)
  : m(m_), n(n_), mat(m_*n_)
{
#ifdef DEBUG
  if(m<=0||n<=0)
    error("invalid matrix size");
  if(m>64&&n>64)
    warn("large matrix allocated");
#endif
}

Mat::Mat(const Mat& b)
  : m(b.m), n(b.n), mat(b.mat)
{}

Mat& Mat::zero()
{
  double *A=&mat[0];
  for(int k=m*n,i=0; i<k; i++)
    A[i]=0.0;
  return *this;
}

Mat& Mat::eye()
{
#ifdef DEBUG
  if(m!=n)
    warn("matrix is not square");
#endif
  double *A=&mat[0];
  for(int i=0; i<m; i++)
    for(int j=0; j<n; j++)
      A[i*n +j]= i==j?1.0:0.0;
  return *this;
}

Mat Mat::t() const
{
  Mat At(n,m);
  const double *A=&mat[0];
  double *T=&At.mat[0];
  for(int i=0; i<m; i++)
    for(int j=0; j<n; j++)
      T[j*m+i]=A[i*n+j];
  return At;
}

double Mat::tr() const
{
#ifdef DEBUG
  if(m!=n)
    warn("matrix is not square");
#endif
  const double *A=&mat[0];
  double trace=A[0];
  for(int i=1; i<n; i++)
    trace+=A[i*n +i];
  return trace;
}

static void mat_inv2x2(const double *A, double *B)
{
  double det, inv_det;
  det= A[0]*A[3]-A[1]*A[2];
#ifdef DEBUG
  if(fabs(det)<EPS)
    warn("matrix is singular or near singular");
#endif
  inv_det = 1.0 / det;
  B[0]= inv_det*A[3];
  B[1]=-inv_det*A[1];
  B[2]=-inv_det*A[2];
  B[3]= inv_det*A[0];
}

static void mat_inv3x3(const double *A, double *B)
{
  double det, inv_det;  
  B[0]= A[4]*A[8] - A[5]*A[7];
  B[1]= A[2]*A[7] - A[1]*A[8];
  B[2]= A[1]*A[5] - A[2]*A[4];  
  B[3]= A[5]*A[6] - A[3]*A[8];
  B[4]= A[0]*A[8] - A[2]*A[6];
  B[5]= A[2]*A[3] - A[0]*A[5];  
  B[6]= A[3]*A[7] - A[4]*A[6];
  B[7]= A[1]*A[6] - A[0]*A[7];
  B[8]= A[0]*A[4] - A[1]*A[3];
  det= A[0]*B[0] + A[1]*B[3] +A[2]*B[6];
  
#ifdef DEBUG
  if(fabs(det)<EPS)
    warn("matrix is singular or near singular");
#endif

  inv_det = 1.0 / det;  
  for(int i=0; i<9; i++)
    B[i]*=inv_det;
}

static void mat_inv4x4(const double *A, double *B)
{
  double det, inv_det;
  double det2_01_01=A[0]*A[ 5]-A[1]*A[ 4];
  double det2_01_02=A[0]*A[ 6]-A[2]*A[ 4];
  double det2_01_03=A[0]*A[ 7]-A[3]*A[ 4];
  double det2_01_12=A[1]*A[ 6]-A[2]*A[ 5];
  double det2_01_13=A[1]*A[ 7]-A[3]*A[ 5];
  double det2_01_23=A[2]*A[ 7]-A[3]*A[ 6];
  double det2_03_01=A[0]*A[13]-A[1]*A[12];
  double det2_03_02=A[0]*A[14]-A[2]*A[12];
  double det2_03_03=A[0]*A[15]-A[3]*A[12];
  double det2_03_12=A[1]*A[14]-A[2]*A[13];
  double det2_03_13=A[1]*A[15]-A[3]*A[13];
  double det2_03_23=A[2]*A[15]-A[3]*A[14];
  double det2_13_01=A[4]*A[13]-A[5]*A[12];
  double det2_13_02=A[4]*A[14]-A[6]*A[12];
  double det2_13_03=A[4]*A[15]-A[7]*A[12];
  double det2_13_12=A[5]*A[14]-A[6]*A[13];
  double det2_13_13=A[5]*A[15]-A[7]*A[13];
  double det2_13_23=A[6]*A[15]-A[7]*A[14];
  B[ 0]=-(A[ 9]*det2_13_23-A[10]*det2_13_13+A[11]*det2_13_12);
  B[ 1]= (A[ 9]*det2_03_23-A[10]*det2_03_13+A[11]*det2_03_12);
  B[ 2]= (A[13]*det2_01_23-A[14]*det2_01_13+A[15]*det2_01_12);
  B[ 3]=-(A[ 9]*det2_01_23-A[10]*det2_01_13+A[11]*det2_01_12);
  B[ 4]= (A[ 8]*det2_13_23-A[10]*det2_13_03+A[11]*det2_13_02);
  B[ 5]=-(A[ 8]*det2_03_23-A[10]*det2_03_03+A[11]*det2_03_02);
  B[ 6]=-(A[12]*det2_01_23-A[14]*det2_01_03+A[15]*det2_01_02);
  B[ 7]= (A[ 8]*det2_01_23-A[10]*det2_01_03+A[11]*det2_01_02);
  B[ 8]=-(A[ 8]*det2_13_13-A[ 9]*det2_13_03+A[11]*det2_13_01);
  B[ 9]= (A[ 8]*det2_03_13-A[ 9]*det2_03_03+A[11]*det2_03_01);
  B[10]= (A[12]*det2_01_13-A[13]*det2_01_03+A[15]*det2_01_01);
  B[11]=-(A[ 8]*det2_01_13-A[ 9]*det2_01_03+A[11]*det2_01_01);
  B[12]= (A[ 8]*det2_13_12-A[ 9]*det2_13_02+A[10]*det2_13_01);
  B[13]=-(A[ 8]*det2_03_12-A[ 9]*det2_03_02+A[10]*det2_03_01);
  B[14]=-(A[12]*det2_01_12-A[13]*det2_01_02+A[14]*det2_01_01);
  B[15]= (A[ 8]*det2_01_12-A[ 9]*det2_01_02+A[10]*det2_01_01);
  det= B[7]*A[13] +B[3]*A[12] +B[11]*A[14] +B[15]*A[15];
  
#ifdef DEBUG
  if(fabs(det)<EPS)
    warn("matrix is singular or near singular");
#endif

  inv_det = 1.0 / det;  
  for(int i=0; i<16; i++)
    B[i]*=inv_det;
}

#define SWAP(a,b) t=a; a=b; b=t;
Mat Mat::inv() const
{
#ifdef DEBUG
  if(m!=n)
    error("matrix is not square");
#endif
  Mat Q(n,n);
  const double *A=&mat[0];
  double *B=&Q.mat[0];
  switch(n)
  {
  case 1:
    B[0]=1.0 / A[0]; break;
  case 2:
    mat_inv2x2(A,B); break;
  case 3:
    mat_inv3x3(A,B); break;
  case 4:
    mat_inv4x4(A,B); break;
  default:
    Mat lu(*this);
    Q.eye();
    double det,t;
    double *LU=&lu.mat[0];
    std::vector<int> pivot(n);
    int *piv=&pivot[0];
    det=ludcmp(LU,n,n,piv);
#ifdef DEBUG
    if(fabs(det)<EPS)
      warn("matrix is singular or near singular");
#endif
    lubksb(LU,n,B,piv);
    for(int i=0; i<n; i++){
      for(int j=0; j<i; j++){
        SWAP(B[i*n+j], B[j*n+i])
      }
    }
    break;
  }
  return Q;
}
#undef SWAP

Mat Mat::operator*(const Mat& b) const
{
#ifdef DEBUG
  if(n!=b.m)
    error("incompatible matrix sizes");
#endif
  int q=b.n;
  Mat Q(m,q);
  
  const double *A=&mat[0];
  const double *B=&b.mat[0];
  double *C=&Q.mat[0];
  
  for(int k=m*q,i=0; i<k; i++)
    C[i]=0.0;
  
  for(int i=0; i<m; i++)
    for(int k=0; k<n; k++)
      for(int j=0; j<q; j++)
        C[i*q+j] += A[i*n+k] * B[k*q+j];
      
  return Q;
}
Mat Mat::operator+(const Mat& b) const
{
  Mat Q(*this);
  Q+=b;
  return Q;
}
Mat Mat::operator-(const Mat& b) const
{
  Mat Q(*this);
  Q-=b;
  return Q;
}
Mat Mat::operator*(double s) const
{
  Mat Q(*this);
  Q*=s;
  return Q;
}
Mat Mat::operator/(double s) const
{
  Mat Q(*this);
  Q/=s;
  return Q;
}

Mat Mat::operator-() const
{
  Mat Q(m,n);
  const double *A=&mat[0];
  double *B=&Q.mat[0];
  for(int k=m*n,i=0; i<k; i++)
    B[i]=-A[i];
  return Q;
}
  
double& Mat::operator()(int i, int j)
{
#ifdef DEBUG
  if(i<0||i>=m||j<0||j>=n)
    error("index out of range");
#endif
  return mat[i*n + j];
}

const double& Mat::operator()(int i, int j) const
{
#ifdef DEBUG
  if(i<0||i>=m||j<0||j>=n)
    error("index out of range");
#endif
  return mat[i*n + j];
}
  
Mat& Mat::operator=(const Mat& b)
{
  int n1=m*n;
  int n2=b.m*b.n;
  double *A=&mat[0];
  const double *B=&b.mat[0];
  if(n1!=n2)
    mat.resize(n2);
  for(int i=0; i<n2; i++)
    A[i]=B[i];
  m=b.m;
  n=b.n;
  return *this;
}

Mat& Mat::operator+=(const Mat& b)
{
#ifdef DEBUG
  if(m!=b.m||n!=b.n)
    error("incompatible matrix sizes");
#endif
  int k=m*n;
  double *A=&mat[0];
  const double *B=&b.mat[0];
  for(int i=0; i<k; i++)
    A[i]+=B[i];
  return *this;
}

Mat& Mat::operator-=(const Mat& b)
{
#ifdef DEBUG
  if(m!=b.m||n!=b.n)
    error("incompatible matrix sizes");
#endif
  int k=m*n;
  double *A=&mat[0];
  const double *B=&b.mat[0];
  for(int i=0; i<k; i++)
    A[i]-=B[i];
  return *this;
}

Mat& Mat::operator*=(double s)
{
  int k=m*n;
  double *A=&mat[0];
  for(int i=0; i<k; i++)
    A[i]*=s;
  return *this;
}

Mat& Mat::operator/=(double s)
{
  int k=m*n;
  double *A=&mat[0];
#ifdef DEBUG
  if(fabs(s)<EPS)
    warn("division by zero");
#endif
  s=1.0/s;
  for(int i=0; i<k; i++)
    A[i]*=s;
  return *this;
}

bool Mat::iszero(double tol) const
{
  const double *A=&mat[0];
  for(int k=m*n,i=0; i<k; i++){
    if(fabs(A[i])>tol){
      return false;
    }
  }
  return true;
}

bool Mat::issymmetric(double tol) const
{
  if(m!=n)
    return false;
  const double *A=&mat[0];
  for(int i=0; i<n; i++){
    for(int j=0; j<i; j++){
      if(fabs(A[i*n+j]-A[j*n+i])>tol){
        return false;
      }
    }
  }
  return true;
}

bool Mat::isidentity(double tol) const
{
  if(m!=n)
    return false;
  const double *A=&mat[0];
  for(int i=0; i<m; i++){
    for(int j=0; j<n; j++){
      if(fabs(A[i*n+j]-(i==j?1.0:0.0))>tol){
        return false;
      }
    }
  }
  return true;
}
  
bool Mat::compare(const Mat& b, double tol) const
{
  if(m!=b.m||n!=b.n)
    return false;
  const double *A=&mat[0];
  const double *B=&b.mat[0];
  for(int k=m*n,i=0; i<k; i++)
    if(fabs(A[i]-B[i])>tol)
      return false;
  return true;
}
  
std::ostream& operator << (std::ostream& os, const Mat& b)
{
  for(int i=0; i<b.m; i++){
    for(int j=0; j<b.n; j++){
      os<<std::setw(14)<<b.mat[i*b.n+j];
    }
    os<<std::endl;
  }
  return os;
}
