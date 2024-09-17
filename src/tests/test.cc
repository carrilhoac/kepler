
#include "test.h"

void fail_(const char *from, const char *msg)
{
  std::cerr<<"[FAILED] "<<from<<": "<<msg<<std::endl;
  exit(1);
}

int main(int argc, char **argv)
{ 
  std::cout<<"[CORE] ";
  test_core();
  std::cout<<"all tests run successfully"<<std::endl;
  
  std::cout<<"[MATH] ";
  test_math();
  std::cout<<"all tests run successfully"<<std::endl;
  
  std::cout<<"[TIME] ";
  test_time();
  std::cout<<"all tests run successfully"<<std::endl;
  
  std::cout<<"[SPHEROID] ";
  test_spheroid();
  std::cout<<"all tests run successfully"<<std::endl;
  
  std::cout<<"[EPHEMERIS] ";
  test_ephemeris();
  std::cout<<"all tests run successfully"<<std::endl;
  
  std::cout<<"[AMOSPHERE] ";
  test_atmosphere();
  std::cout<<"all tests run successfully"<<std::endl;
  
  return 0;
}
