SHELL := /bin/bash

CC=gcc
CCFLAGS=-g -I/home/rballouz/SoftPolyhedra/Dynamics
#OMP=-assume nocc-omp -openmp -openmp_report1

PROGS = 
OBJS = simple.o

all: $(PROGS) MAIN
   
MAIN: simple.c
	$(CC) $(CCFLAGS) -o simple simple.c -lm
#	$(CC) $(CCFLAGS) $(OBJS) -o simple -lm

#Utility Targets
.PHONY: all clean 

clean:
	rm -rf *.o *.temp simple
