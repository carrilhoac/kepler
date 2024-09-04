
#include "test.h"

void test_time()
{  
  Time rnx;
  rnx.from_rnx("2024 07 15 22 00 00");
  
  std::cout << std::fixed;
  std::cout << rnx.to_double()+10800.0 << std::endl;
  std::cout << 1721091600.0 << std::endl;
}
