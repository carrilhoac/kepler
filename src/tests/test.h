
#ifndef KEPLER_test_h
#define KEPLER_test_h

#include "../kepler.h"

void fail_(const char *from, const char *msg);
#define fail(x) fail_(__PRETTY_FUNCTION__, x)

void test_ellps();
void test_linalg();
void test_time();
void test_broadcast();

#endif 