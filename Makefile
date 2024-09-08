
CC=g++
CFLAGS= -Wall -O3 -mavx2 -pedantic -std=c++20 -DDEBUG

OBJS_LIB = spheroid.o math.o time.o ephemeris.o atmosphere.o
OBJS_TEST= test_spheroid.o test_math.o test_time.o test_ephemeris.o

all: libkepler.a

# ---------------------------------------------------------------------------
#  STATIC LIB
# ---------------------------------------------------------------------------
libkepler.a: kepler.h constants.h ${OBJS_LIB}
	ar rcs libkepler.a ${OBJS_LIB}

# ---------------------------------------------------------------------------
# OBJ files
# ---------------------------------------------------------------------------
spheroid.o: kepler.h spheroid.cc
	${CC} ${CFLAGS} -c spheroid.cc
math.o: kepler.h math.cc
	${CC} ${CFLAGS} -c math.cc
time.o: kepler.h time.cc
	${CC} ${CFLAGS} -c time.cc
ephemeris.o: kepler.h ephemeris.cc 
	${CC} ${CFLAGS} -c ephemeris.cc
atmosphere.o: kepler.h atmosphere.cc 
	${CC} ${CFLAGS} -c atmosphere.cc

# ---------------------------------------------------------------------------
# TESTS
# ---------------------------------------------------------------------------
test_spheroid.o: test_spheroid.cc libkepler.a
	${CC} ${CFLAGS} -c test_spheroid.cc
test_math.o: test_math.cc libkepler.a
	${CC} ${CFLAGS} -c test_math.cc
test_time.o: test_time.cc libkepler.a
	${CC} ${CFLAGS} -c test_time.cc
test_ephemeris.o: test_ephemeris.cc libkepler.a 
	${CC} ${CFLAGS} -c test_ephemeris.cc
tests: kepler.h libkepler.a test.h test.cc ${OBJS_TEST}
	${CC} ${CFLAGS} -o tests test.cc ${OBJS_TEST} ./libkepler.a
test: tests
	./tests

# ---------------------------------------------------------------------------
# CLEAN
# ---------------------------------------------------------------------------
clean:
	rm -f *.o libkepler.a tests
