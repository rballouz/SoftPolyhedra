#QHULLDIR = ./qhull/qhull
QHULLDIR = /usr/local/Cellar/qhull/2012.1/include/libqhull/
CFLAGS = -O2 -g
#CFLAGS = -g
LDLIBS = -lm
CC = gcc

gjkdemo:	gjkdemo.o gjk.o
	$(CC) $(CFLAGS) gjkdemo.o gjk.o -o gjkdemo -lm 

gjkqhull:	gjkqhull.o gjk.o sac.o
	$(CC) $(CFLAGS) gjkqhull.o gjk.o sac.o -o gjkqhull -L$(QHULLDIR) -lqhull -lm 

gjkqhull.o:	gjkdemo.c
	$(CC) -c $(CFLAGS) -DQHULL gjkdemo.c -o gjkqhull.o

sac.o:	sac.c gjk.h
	$(CC) -c $(CFLAGS) -I$(QHULLDIR) sac.c -o sac.o

clean:
	rm -f *.o core

.c.o:	gjk.h
	$(CC) -c $(CFLAGS) $*.c
