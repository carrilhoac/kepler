
#include "test.h"

void fail_(const char *from, const char *msg)
{
  std::cerr<<"[FAILED] "<<from<<": "<<msg<<std::endl;
  exit(1);
}

int main(int argc, char **argv)
{ 
  std::cout<<"[SPHEROID] ";
  test_ellps();
  std::cout<<"all tests run successfully"<<std::endl;
  
  std::cout<<"[MATRIX] ";
  test_linalg();
  std::cout<<"all tests run successfully"<<std::endl;
  
  std::cout<<"[TIME] ";
  test_time();
  std::cout<<"all tests run successfully"<<std::endl;
  
  std::cout<<"[BROADCAST] ";
  test_broadcast();
  std::cout<<"all tests run successfully"<<std::endl;
  
  return 0;
}
