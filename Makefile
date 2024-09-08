
CC=g++
CFLAGS= -Wall -O3 -mavx2 -pedantic -std=c++20 -DDEBUG

OBJS_LIB = spheroid.o vec3.o matrix.o timesys.o broadcast.o atmosphere.o
OBJS_TEST= test_spheroid.o test_matrix.o test_timesys.o test_broadcast.o

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
vec3.o: kepler.h vec3.cc
	${CC} ${CFLAGS} -c vec3.cc
matrix.o: kepler.h matrix.cc
	${CC} ${CFLAGS} -c matrix.cc
timesys.o: kepler.h timesys.cc
	${CC} ${CFLAGS} -c timesys.cc
broadcast.o: kepler.h broadcast.cc 
	${CC} ${CFLAGS} -c broadcast.cc
atmosphere.o: kepler.h atmosphere.cc 
	${CC} ${CFLAGS} -c atmosphere.cc

# ---------------------------------------------------------------------------
# TESTS
# ---------------------------------------------------------------------------
test_spheroid.o: test_spheroid.cc libkepler.a
	${CC} ${CFLAGS} -c test_spheroid.cc
test_matrix.o: test_matrix.cc libkepler.a
	${CC} ${CFLAGS} -c test_matrix.cc
test_timesys.o: test_timesys.cc libkepler.a
	${CC} ${CFLAGS} -c test_timesys.cc
test_broadcast.o: test_broadcast.cc libkepler.a 
	${CC} ${CFLAGS} -c test_broadcast.cc
tests: kepler.h libkepler.a test.h test.cc ${OBJS_TEST}
	${CC} ${CFLAGS} -o tests test.cc ${OBJS_TEST} ./libkepler.a
test: tests
	./tests

# ---------------------------------------------------------------------------
# CLEAN
# ---------------------------------------------------------------------------
clean:
	rm -f *.o libkepler.a tests
