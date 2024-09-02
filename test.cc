
#include "test.h"

void fail_(const char *from, const char *msg)
{
  std::cerr<<"[FAILED] "<<from<<": "<<msg<<std::endl;
  exit(1);
}

int main(int argc, char **argv)
{ 
  test_ellps();
  std::cout<<"[ELLPS] all tests run successfully"<<std::endl;
  
  test_linalg();
  std::cout<<"[LINALG] all tests run successfully"<<std::endl;
  
  test_time();
  std::cout<<"[TIME] all tests run successfully"<<std::endl;
  
  return 0;
}
