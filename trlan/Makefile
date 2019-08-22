MINSTALLDIR = /usr/local/Wolfram/Mathematica/11.2
INCMMA = $(MINSTALLDIR)/SystemFiles/IncludeFiles/C

CFLAGS = -fPIC

lib:
	$(CC) -shared $(CFLAGS) -L ../../trlan-201009/ -ltrlan -lgfortran -I $(INCMMA) -o libtrlan-wrap.so trlan-wrap.c

