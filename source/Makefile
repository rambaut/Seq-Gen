CLINKER = cc 
CLIBS = -lm

RM = /bin/rm

LIBS = $(CLIBS)

CFLAGS = -O3
DEBUG_FLAGS = -g -O0

OBJ = seq-gen.o global.o evolve.o model.o nucmodels.o aamodels.o eigen.o treefile.o progress.o gamma.o twister.o

.SUFFIXES: .o .c .h

.c.o:
	$(CLINKER) $(CFLAGS) -c $*.c 

dummy:	default

clean:
	$(RM) -f *~ *.o core

default: all

all: seq-gen

seq-gen: $(OBJ)
	$(CLINKER) $(CFLAGS) -o seq-gen $(OBJ) $(LIBS)

debug: CFLAGS = $(DEBUG_FLAGS)
debug: seq-gen

%.o: %.c
	$(CLINKER) $(CFLAGS) -c $<

%.o: %.c %.h
	$(CLINKER) $(CFLAGS) -c $<
