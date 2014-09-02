# Makefile for depSolver v2.0
#
# Copyright 2006, 2008 Carlos Rosales Fernandez and IHPC (A*STAR).             *
#                                                                              *
# This file is part of depSolver.                                              *
#                                                                              *
# depSolver is free software: you can redistribute it and/or modify it under   *
# the terms of the GNU GPL version 3 or (at your option) any later version.    *
#                                                                              *
# depSolver is distributed in the hope that it will be useful, but WITHOUT ANY *
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS    *
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for details.    *
#                                                                              *
# You should have received a copy of the GNU General Public License along with *
# depSolver. If not, see <http://gnu.org/licenses/gpl.html>.                   *

SHELL = /bin/sh

# Compiler
CC = gcc

# Flags
# In Debian and Debian variants like the Ubuntu family, please add the flag
# -fno-stack-protector because libssp0-dev has a bug and will post an error
# message on return from main. This error is bogus. The code DOES NOT smash the stack.
LFLAGS =  -lm -msse2 -mfpmath=sse -funroll-loops -fprefetch-loop-arrays -O3
CFLAGS =  -msse2 -mfpmath=sse -funroll-loops -fprefetch-loop-arrays -O3

# Linking Step
depSolver: assembly_el.o depolarization.o depSolver.o errorHandler.o field.o  \
           force.o forceIntegrals.o gaussBksb.o globalToLocal.o gmres_el.o    \
           kernelIntegrals.o memHandler.o normals.o postProcess.o potential.o \
           preProcess.o shapeFunctions.o
	$(CC) $(LFLAGS) assembly_el.o depolarization.o depSolver.o errorHandler.o   \
  field.o force.o forceIntegrals.o gaussBksb.o globalToLocal.o gmres_el.o     \
  kernelIntegrals.o memHandler.o normals.o postProcess.o potential.o          \
  preProcess.o shapeFunctions.o -o depSolver


depSolver-sft: assembly_el.o depolarization.o depSolver-shift.o errorHandler.o\
               field.o force.o forceIntegrals.o gaussBksb.o globalToLocal.o   \
               gmres_el.o kernelIntegrals.o memHandler.o normals.o            \
               postProcess.o potential.o preProcess.o shapeFunctions.o
	$(CC) $(LFLAGS) assembly_el.o depolarization.o depSolver-shift.o            \
  errorHandler.o field.o force.o forceIntegrals.o gaussBksb.o globalToLocal.o \
  gmres_el.o kernelIntegrals.o memHandler.o normals.o postProcess.o           \
  potential.o preProcess.o shapeFunctions.o -o depSolver-sft


# Compiling Step
assembly_el.o: assembly_el.c constants.h kernelIntegrals.h memHandler.h
	$(CC) -c $(CFLAGS) assembly_el.c -o assembly_el.o

depolarization.o: depolarization.c errorHandler.h globalToLocal.h             \
                  shapeFunctions.h
	$(CC) -c $(CFLAGS) depolarization.c -o depolarization.o

depSolver.o: depSolver.c assembly_el.h errorHandler.h gaussBksb.h gaussData.h \
             gmres_el.h memHandler.h postProcess.h preProcess.h
	$(CC) -c $(CFLAGS) depSolver.c -o depSolver.o

depSolver-shift.o: depSolver-shift.c assembly_el.h errorHandler.h gaussBksb.h \
                   gaussData.h gmres_el.h memHandler.h postProcess.h preProcess.h
	$(CC) -c $(CFLAGS) depSolver-shift.c -o depSolver-shift.o

errorHandler.o: errorHandler.c
	$(CC) -c $(CFLAGS) errorHandler.c -o errorHandler.o

field.o: field.c constants.h kernelIntegrals.h normals.h
	$(CC) -c $(CFLAGS) field.c -o field.o

force.o: force.c constants.h depolarization.h field.h forceIntegrals.h        \
         globalToLocal.h shapeFunctions.h memHandler.h
	$(CC) -c $(CFLAGS) force.c -o force.o

forceIntegrals.o: forceIntegrals.c constants.h globalToLocal.h                \
                  shapeFunctions.h normals.h
	$(CC) -c $(CFLAGS) forceIntegrals.c -o forceIntegrals.o

gaussBksb.o: gaussBksb.c
	$(CC) -c $(CFLAGS) gaussBksb.c -o gaussBksb.o

globalToLocal.o: globalToLocal.c
	$(CC) -c $(CFLAGS) globalToLocal.c -o globalToLocal.o

gmres_el.o: gmres_el.c errorHandler.h memHandler.h
	$(CC) -c $(CFLAGS) gmres_el.c -o gmres_el.o

kernelIntegrals.o: kernelIntegrals.c constants.h globalToLocal.h              \
                   shapeFunctions.h
	$(CC) -c $(CFLAGS) kernelIntegrals.c -o kernelIntegrals.o

memHandler.o: memHandler.c errorHandler.h
	$(CC) -c $(CFLAGS) memHandler.c -o memHandler.o

normals.o: normals.c 
	$(CC) -c $(CFLAGS) normals.c -o normals.o

postProcess.o: postProcess.c constants.h memHandler.h field.h force.h         \
               potential.h
	$(CC) -c $(CFLAGS) postProcess.c -o postProcess.o

potential.o: potential.c constants.h kernelIntegrals.h normals.h
	$(CC) -c $(CFLAGS) potential.c -o potential.o

preProcess.o: preProcess.c
	$(CC) -c $(CFLAGS) preProcess.c -o preProcess.o

shapeFunctions.o: shapeFunctions.c 
	$(CC) -c $(CFLAGS) shapeFunctions.c -o shapeFunctions.o

# Build options
all: depSolver depSolver-sft

build: clean depSolver depSolver-sft install

install:
	mv ./depSolver ../bin/depSolver
	mv ./depSolver-sft ../bin/depSolver-sft

clean:
	rm -f ./*.o

