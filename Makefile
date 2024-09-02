
CC=g++
CFLAGS= -Wall -O3 -pedantic -std=c++20 -DDEBUG

OBJS_LIB = ellps.o linalg.o time.o
OBJS_TEST= test_ellps.o test_linalg.o test_time.o

all: libkepler.a

# ---------------------------------------------------------------------------
#  STATIC LIB
# ---------------------------------------------------------------------------
libkepler.a: kepler.h ${OBJS_LIB}
	ar rcs libkepler.a ${OBJS_LIB}

# ---------------------------------------------------------------------------
# OBJ files
# ---------------------------------------------------------------------------
ellps.o: kepler.h ellps.cc
	${CC} ${CFLAGS} -c ellps.cc
linalg.o: kepler.h linalg.cc
	${CC} ${CFLAGS} -c linalg.cc
time.o: kepler.h time.cc
	${CC} ${CFLAGS} -c time.cc


# ---------------------------------------------------------------------------
# TESTS
# ---------------------------------------------------------------------------
test_ellps.o: test_ellps.cc libkepler.a
	${CC} ${CFLAGS} -c test_ellps.cc
test_linalg.o: test_linalg.cc libkepler.a
	${CC} ${CFLAGS} -c test_linalg.cc
test_time.o: test_time.cc libkepler.a
	${CC} ${CFLAGS} -c test_time.cc
tests: kepler.h libkepler.a test.h test.cc ${OBJS_TEST}
	${CC} ${CFLAGS} -o tests test.cc ${OBJS_TEST} ./libkepler.a
test: tests
	./tests

# ---------------------------------------------------------------------------
# CLEAN
# ---------------------------------------------------------------------------
clean:
	rm -f *.o libkepler.a tests
