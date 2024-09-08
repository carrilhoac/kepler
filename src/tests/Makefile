
CC=g++
CFLAGS= -Wall -O3 -mavx2 -pedantic -std=c++20 -DDEBUG

OBJS_TEST = test_math.o test_time.o test_spheroid.o test_ephemeris.o

all: tests
	./tests
	
# ---------------------------------------------------------------------------
# TESTS
# ---------------------------------------------------------------------------
test_spheroid.o: test_spheroid.cc 
	${CC} ${CFLAGS} -c test_spheroid.cc
	
test_math.o: test_math.cc 
	${CC} ${CFLAGS} -c test_math.cc
	
test_time.o: test_time.cc 
	${CC} ${CFLAGS} -c test_time.cc
	
test_ephemeris.o: test_ephemeris.cc
	${CC} ${CFLAGS} -c test_ephemeris.cc
	
tests: ../kepler.h ../libkepler.a test.h test.cc ${OBJS_TEST}
	${CC} ${CFLAGS} -o tests test.cc ${OBJS_TEST} ../libkepler.a
	
# ---------------------------------------------------------------------------
# CLEAN
# ---------------------------------------------------------------------------
clean:
	rm -f *.o tests
