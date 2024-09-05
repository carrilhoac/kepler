
CC=g++
CFLAGS= -Wall -O3 -pedantic -std=c++20 -DDEBUG

OBJS_LIB = spheroid.o matrix.o timesys.o broadcast.o
OBJS_TEST= test_spheroid.o test_matrix.o test_timesys.o

all: libkepler.a

# ---------------------------------------------------------------------------
#  STATIC LIB
# ---------------------------------------------------------------------------
libkepler.a: kepler.h ${OBJS_LIB}
	ar rcs libkepler.a ${OBJS_LIB}

# ---------------------------------------------------------------------------
# OBJ files
# ---------------------------------------------------------------------------
spheroid.o: kepler.h spheroid.cc
	${CC} ${CFLAGS} -c spheroid.cc
matrix.o: kepler.h matrix.cc
	${CC} ${CFLAGS} -c matrix.cc
timesys.o: kepler.h timesys.cc
	${CC} ${CFLAGS} -c timesys.cc
broadcast.o: kepler.h broadcast.cc 
	${CC} ${CFLAGS} -c broadcast.cc

# ---------------------------------------------------------------------------
# TESTS
# ---------------------------------------------------------------------------
test_spheroid.o: test_spheroid.cc libkepler.a
	${CC} ${CFLAGS} -c test_spheroid.cc
test_matrix.o: test_matrix.cc libkepler.a
	${CC} ${CFLAGS} -c test_matrix.cc
test_timesys.o: test_timesys.cc libkepler.a
	${CC} ${CFLAGS} -c test_timesys.cc
tests: kepler.h libkepler.a test.h test.cc ${OBJS_TEST}
	${CC} ${CFLAGS} -o tests test.cc ${OBJS_TEST} ./libkepler.a
test: tests
	./tests

# ---------------------------------------------------------------------------
# CLEAN
# ---------------------------------------------------------------------------
clean:
	rm -f *.o libkepler.a tests
