# Install jJSON from github
#
#  Had trouble with runtime finding library in /usr/local/lib
#  Ran "sudo ldconfig" to update cache.
#
CFLAGS = -fPIC -Wall

shifts: shifts.o
	$(CC) $(CFLAGS)  -o shifts shifts.o -L/usr/local/lib -lm -lcjson

%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS)

clean:
	rm shifts *.o
