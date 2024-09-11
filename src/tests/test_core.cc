
#include "test.h"

void test_core()
{  
  {
    Mem f;
    f.load("./data/dummy.txt");
    
    if(f.size()!=44)
      fail("different file size");
  }
}
