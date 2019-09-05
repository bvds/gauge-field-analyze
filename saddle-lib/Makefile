#
#   Stand-alone program, maybe called from Mathematica
#
MINRESDIR = ../../linear-algebra/minres-qlp
TRLANDIR = ../../linear-algebra/nutrlan
#
#  To use librsb for the sparse matrices, add this flag:
#                                           -DUSE_LIBRSB
CFLAGS = -fPIC -Wall -g -Wextra -pedantic -std=gnu11 -I$(TRLANDIR)/INC

MINRESOBJS = $(MINRESDIR)/minresqlpModule.o $(MINRESDIR)/minresqlpDataModule.o  $(MINRESDIR)/minresqlpBlasModule.o
SHIFTOBJS = shifts.o dynamic.o eigensystem.o cutoff.o linear.o

shifts: $(SHIFTOBJS)
	$(CC) $(CFLAGS) -o shifts $(SHIFTOBJS) $(MINRESOBJS) -L$(TRLANDIR) -L/usr/local/lib -lnutrlan -lgfortran -llapack -lblas -lrsb -lm -lcjson

%.o: %.c shifts.h
	$(CC) -c -o $@ $< $(CFLAGS)

clean:
	rm shifts *.o
