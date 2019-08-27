#
#   Stand-alone program, maybe called from Mathematica
#
MINRESLIB = ../../linear-algebra/minresqlp_f90_ver27
TRLANLIB = ../../linear-algebra/nutrlan-0.2
CFLAGS = -fPIC -Wall -g -Wextra -pedantic -std=gnu11

shifts: shifts.o dynamic.o
	$(CC) $(CFLAGS) -o shifts shifts.o dynamic.o $(MINRESLIB)/minresqlpModule.o  $(MINRESLIB)/minresqlpDataModule.o  $(MINRESLIB)/minresqlpBlasModule.o -L/usr/local/lib -lgfortran -lm -lcjson

%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS)

clean:
	rm shifts *.o
