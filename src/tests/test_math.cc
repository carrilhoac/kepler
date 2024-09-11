
#include "test.h"

static int test_constr()
{
  {
    Mat a(4,2);
    
    if(a.rows()!=4||a.cols()!=2)
      fail("incorrect matrix size");
  }
  
  {
    Mat a(4,4);
    a.eye();
    
    Mat b(a);
    
    if(!b.isidentity())
      fail("should be identity matrix");
  }
  
  {
    Mat a(3,3);
    a.eye();
    const double *A=a.data();
    const double B[9]=
    {
      1.0,0.0,0.0, 
      0.0,1.0,0.0, 
      0.0,0.0,1.0
    };
    for(int i=0; i<9; i++)
      if(A[i]!=B[i])
        fail("eye() didnt result in identity matrix");
  }
  
  {
    Mat a, b;
    
    if(!a.iszero())
      fail("constructor didnt init elements as zero");
    
    b.eye();
    if(!b.isidentity())
      fail("eye() didnt result in identity matrix");
    
    b.zero();
    if(!a.iszero())
      fail("zero() didnt clear matrix elements");
    
    a.eye();
    if(a==b)
      fail("matrices should be different");
  }
  
  return 0;
}

static int test_ops()
{
  {
    Mat a(4,4), b(4,4);
    
    a.eye();
    b.eye();
    
    a *= 2;
    if(a.isidentity())
      fail("should be different from identity");
    
    a -= b;
    
    if(!a.isidentity())
      fail("not identity matrix");
    
    if(!a.issymmetric())
      fail("not symmetric matrix");
    
    if(a.tr() != 4.0)
      fail("tr() should be 4.0");
  }
  
  {
    Mat a(2,3);
    Mat b(3,2);
    
    a(0,0)=1; a(0,1)=2; a(0,2)=3;
    a(1,0)=4; a(1,1)=5; a(1,2)=6;
    
    b(0,0)=1; b(0,1)=4;
    b(1,0)=2; b(1,1)=5;
    b(2,0)=3; b(2,1)=6;
    
    if(a.t()!=b)
      fail("transposition failed");
  }
  
  {
    Mat a(4,4);
    Mat b(4,4);
    
    a.eye();
    b.eye();
    a = -a;
    
    a += b;
    
    if(!a.iszero())
      fail("matrix should be zero");
  }
  
  return 0;
}

static int test_mult()
{
  { // test 1
    Mat a(2,3);
    Mat b(3,2);
    Mat c(2,2);
    Mat d(3,3);
    
    a(0,0)=1; a(0,1)=2; a(0,2)=3;
    a(1,0)=4; a(1,1)=5; a(1,2)=6;
    
    b(0,0)= 7; b(0,1)= 8;
    b(1,0)= 9; b(1,1)=10;
    b(2,0)=11; b(2,1)=12;
    
    c(0,0)= 58; c(0,1)= 64;
    c(1,0)=139; c(1,1)=154;
    
    d(0,0)=39; d(0,1)=54; d(0,2)=69;
    d(1,0)=49; d(1,1)=68; d(1,2)=87;
    d(2,0)=59; d(2,1)=82; d(2,2)=105;
    
    if(c!=(a*b))
      fail("matrices should be equal");
    if(d!=(b*a))
      fail("matrices should be equal");
  }
  
  { // test 2
    Mat a(1,3);
    Mat b(3,1);
    Mat c(1,1);
    Mat d(3,3);
    
    a(0,0)=1; a(0,1)=2; a(0,2)=3;
    
    b(0,0)= 4;
    b(1,0)= 5;
    b(2,0)= 6;
    
    c(0,0)= 32;
    
    d(0,0)= 4; d(0,1)= 8; d(0,2)=12;
    d(1,0)= 5; d(1,1)=10; d(1,2)=15;
    d(2,0)= 6; d(2,1)=12; d(2,2)=18;
    
    if(c!=(a*b))
      fail("matrices should be equal");
    if(d!=(b*a))
      fail("matrices should be equal");
  }
  
  return 0;
}

static int test_inv()
{
  { // 2x2
    Mat a(2,2);
    Mat b(2,2);
    Mat c(2,2);
    
    c.eye();
    
    a(0,0)=4; a(0,1)=7;
    a(1,0)=2; a(1,1)=6;
    
    b(0,0)= 6; b(0,1)=-7;
    b(1,0)=-2; b(1,1)= 4;
    b /= 10;
    
    if(b!=a.inv())
      fail("matrices should be equal");
    
    if(c!=(a*a.inv()))
      fail("product between A and inv(A) should be identity");
  }
  { // 2x2 linear system
    Mat a(2,2);
    Mat b(2,1);
    Mat c(2,1);
    
    a(0,0)=3.0; a(0,1)=3.2;
    a(1,0)=3.5; a(1,1)=3.6;
    
    b(0,0)= 118.4; 
    b(1,0)= 135.2; 
    
    c(0,0)=16; 
    c(1,0)=22; 
    
    if(c!=(a.inv()*b))
      fail("solution is different from reference");   
  }
  
  { // 3x3 linear system
    Mat a(3,3);
    Mat b(3,1);
    Mat c(3,1);
    Mat d(3,3);
    
    a(0,0)=2; a(0,1)=1; a(0,2)=1;
    a(1,0)=1; a(1,1)=2; a(1,2)=1;
    a(2,0)=1; a(2,1)=3; a(2,2)=3;
    
    b(0,0)=13;
    b(1,0)=11;
    b(2,0)=19;
    
    c(0,0)=4;
    c(1,0)=2;
    c(2,0)=3;
    
    d(0,0)= 3; d(0,1)= 0; d(0,2)=-1;
    d(1,0)=-2; d(1,1)= 5; d(1,2)=-1;
    d(2,0)= 1; d(2,1)=-5; d(2,2)= 3;
    d /= 5;
    
    if(c!=(a.inv()*b))
      fail("solution is different from reference");   
    
    if(d!=a.inv())
      fail("matrices should be equal");  
    
    if(!(a.inv()*d.inv()).isidentity())
      fail("matrix should be identity");  
  }
  
  { // 4x4 
    Mat a(4,4);
    Mat b(4,4);
    
    a(0,0)=2; a(0,1)=1; a(0,2)=3; a(0,3)= 1;
    a(1,0)=4; a(1,1)=3; a(1,2)=1; a(1,3)= 8;
    a(2,0)=6; a(2,1)=2; a(2,2)=7; a(2,3)= 1;
    a(3,0)=2; a(3,1)=1; a(3,2)=1; a(3,3)=-1;
    
    b(0,0)=-61; b(0,1)= 4; b(0,2)= 26; b(0,3)= -3;
    b(1,0)= 88; b(1,1)= 0; b(1,2)=-44; b(1,3)= 44;
    b(2,0)= 28; b(2,1)=-4; b(2,2)= -4; b(2,3)= -8;
    b(3,0)= -6; b(3,1)= 4; b(3,2)=  4; b(3,3)=-14;
    b /= 44;
    
    if(b!=a.inv())
      fail("matrices should be equal");
    
    if(!(a.inv()*b.inv()).isidentity())
      fail("matrix should be identity");
  }
  
  { // 5x5 (uses LU)
    Mat a(5,5);
    Mat b(5,5);
    
    a(0,0)= 0; a(0,1)= 6; a(0,2)= -2; a(0,3)=-1; a(0,4)= 5;
    a(1,0)= 0; a(1,1)= 0; a(1,2)=  0; a(1,3)=-9; a(1,4)=-7;
    a(2,0)= 0; a(2,1)=15; a(2,2)= 35; a(2,3)= 0; a(2,4)= 0;
    a(3,0)= 0; a(3,1)=-1; a(3,2)=-11; a(3,3)=-2; a(3,4)= 1;
    a(4,0)=-2; a(4,1)=-2; a(4,2)=  3; a(4,3)= 0; a(4,4)=-2;
    
    b(0,0)=  305; b(0,1)= 335; b(0,2)= -398; b(0,3)=-1660; b(0,4)=-1240;
    b(1,0)=-1610; b(1,1)=-630; b(1,2)= 1052; b(1,3)= 3640; b(1,4)=    0;
    b(2,0)=  690; b(2,1)= 270; b(2,2)= -380; b(2,3)=-1560; b(2,4)=    0;
    b(3,0)=-1820; b(3,1)=-820; b(3,2)=  952; b(3,3)= 3360; b(3,4)=    0;
    b(4,0)= 2340; b(4,1)= 700; b(4,2)=-1224; b(4,3)=-4320; b(4,4)=    0;
    b /= 2480;
    
    if(b!=a.inv())
      fail("matrices should be equal");
    
    if(!(a.inv()*b.inv()).isidentity())
      fail("matrix should be identity");
    
    if(a.tr()!=31.0)
      fail("tr() should be 31.0");
    
    //if(a.det()!=2480.0)
    //  fail("matrix determinant should be 2480.0");
  }
  
  return 0;
}

void test_math()
{
  test_constr();
  test_ops();
  test_mult();
  test_inv();
}
