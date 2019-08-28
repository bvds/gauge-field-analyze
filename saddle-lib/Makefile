#
#   Stand-alone program, maybe called from Mathematica
#
MINRESDIR = ../../linear-algebra/minresqlp_f90_ver27
TRLANDIR = ../../linear-algebra/nutrlan-0.2
CFLAGS = -fPIC -Wall -g -Wextra -pedantic -std=gnu11

MINRESOBJS = $(MINRESDIR)/minresqlpModule.o $(MINRESDIR)/minresqlpDataModule.o  $(MINRESDIR)/minresqlpBlasModule.o

shifts: shifts.o dynamic.o
	$(CC) $(CFLAGS) -o shifts shifts.o dynamic.o $(MINRESOBJS) -I$(TRLANDIR)/INC -L$(TRLANDIR) -L/usr/local/lib -lnutrlan -lgfortran -lm -lcjson

%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS)

clean:
	rm shifts *.o
