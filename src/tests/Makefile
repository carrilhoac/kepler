
CC=g++
CFLAGS= -Wall -O3 -mavx2 -pedantic -std=c++20 -DDEBUG

OBJS_TEST = test_core.o test_math.o test_time.o test_spheroid.o test_ephemeris.o test_atmosphere.o

all: test

# ---------------------------------------------------------------------------
# TESTS
# ---------------------------------------------------------------------------

test_core.o: test_core.cc 
	${CC} ${CFLAGS} -c test_core.cc
	
test_math.o: test_math.cc 
	${CC} ${CFLAGS} -c test_math.cc
	
test_time.o: test_time.cc 
	${CC} ${CFLAGS} -c test_time.cc
	
test_spheroid.o: test_spheroid.cc 
	${CC} ${CFLAGS} -c test_spheroid.cc
	
test_ephemeris.o: test_ephemeris.cc
	${CC} ${CFLAGS} -c test_ephemeris.cc
	
test_atmosphere.o: test_atmosphere.cc
	${CC} ${CFLAGS} -c test_atmosphere.cc
	
test_all: ../kepler.h ../libkepler.a test.h test.cc ${OBJS_TEST}
	${CC} ${CFLAGS} -o test_all test.cc ${OBJS_TEST} ../libkepler.a
	
test: test_all
	./test_all

# ---------------------------------------------------------------------------
# CLEAN
# ---------------------------------------------------------------------------
clean:
	rm -f *.o test_all
