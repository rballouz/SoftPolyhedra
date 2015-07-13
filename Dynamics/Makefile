SHELL := /bin/bash

CC=gcc
CCFLAGS=-g -I/home/rballouz/SoftPolyhedra/Dynamics
#OMP=-assume nocc-omp -openmp -openmp_report1

PROGS = 
OBJS = simple.o

all: $(PROGS) MAIN
   
MAIN: simple.c
	$(CC) $(CCFLAGS) -o simple simple.c -lm

poly: simple_poly.c
	$(CC) $(CCFLAGS) -o simple_poly simple_poly.c -lm

#Utility Targets
.PHONY: all clean 

clean:
	rm -rf *.o *.temp simple simply_poly
