
CC=g++
CFLAGS= -Wall -O3 -mavx2 -pedantic -std=c++20 -DDEBUG

OBJS_LIB = core.o spheroid.o math.o time.o ephemeris.o atmosphere.o

all: libkepler.a

# ---------------------------------------------------------------------------
#  STATIC LIB
# ---------------------------------------------------------------------------
libkepler.a: kepler.h constants.h ${OBJS_LIB}
	ar rcs libkepler.a ${OBJS_LIB}

# ---------------------------------------------------------------------------
# OBJ files
# ---------------------------------------------------------------------------
core.o: kepler.h core.cc
	${CC} ${CFLAGS} -c core.cc
	
math.o: kepler.h math.cc
	${CC} ${CFLAGS} -c math.cc
	
time.o: kepler.h time.cc
	${CC} ${CFLAGS} -c time.cc
	
spheroid.o: kepler.h spheroid.cc
	${CC} ${CFLAGS} -c spheroid.cc
	
ephemeris.o: kepler.h ephemeris.cc 
	${CC} ${CFLAGS} -c ephemeris.cc
	
atmosphere.o: kepler.h atmosphere.cc 
	${CC} ${CFLAGS} -c atmosphere.cc


test: libkepler.a
	cd tests && $(MAKE) test

# ---------------------------------------------------------------------------
# CLEAN
# ---------------------------------------------------------------------------
clean:
	rm -f *.o libkepler.a && cd tests && $(MAKE) clean
