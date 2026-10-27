CC= gcc
CFLAGS=  -Wall -lm

default: md.4.c
	$(CC) -o md.x md.4.c $(CFLAGS) -O3

debug: md.4.c
	$(CC) -o md.x md.4.c $(CFLAGS) -pg -g

clean:
	rm -f *.x
