# Makefile for depSolver
#
# COPYRIGHT & LICENSE INFORMATION
# 
# Copyright 2006 Carlos Rosales Fernandez and The Institute of High
# Performance Computing (A*STAR)
#
# This file is part of depSolver.
#
# depSolver is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# depSolver is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with depSolver; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

SHELL = /bin/sh

# Compiler
CC = gcc

# Flags
LFLAGS = -lm -msse2 -mfpmath=sse -funroll-loops -fprefetch-loop-arrays -O3
CFLAGS = -msse2 -mfpmath=sse -funroll-loops -fprefetch-loop-arrays -O3

# Linking Step
depSolver: comFilter.o depSolver.o depolarization.o dotProd.o doubleMatrix.o\
     doublePointer.o doubleVector.o electricFormA_tria3.o                     \
     electricFormFullA_tria3.o elemType.o errorHandler.o field_tria3.o        \
     forceEllipsoid_tria3.o forceMST_tria3.o forceMultipole_tria3.o           \
     freeDoubleMatrix.o freeDoublePointer.o freeUintMatrix.o gaussBksb.o      \
     getLocalNormal_tria3.o getNormal_tria3.o initRes_el.o                    \
     intDE_tria3.o intF_tria3.o intG_tria3.o intH_tria3.o intSingularG_tria3.o\
     intSingularH_tria3.o iterGMRES_el.o L2Norm.o matVectProd_el.o            \
     postProcess_tria3.o potential_tria3.o shape_line2.o shape_tria3.o        \
     solverGMRES_el.o uintMatrix.o uintVector.o X2L_line2.o X2L_tria3.o       \
     X2L_tria6.o electricFormA_tria6.o electricFormFullA_tria6.o field_tria6.o\
     potential_tria6.o forceEllipsoid_tria6.o forceMST_tria6.o                \
     forceMultipole_tria6.o getLocalNormal_tria6.o getNormal_tria6.o          \
     intDE_tria6.o intF_tria6.o intG_tria6.o intH_tria6.o intSingularG_tria6.o\
     intSingularH_tria6.o postProcess_tria6.o postProcessCol_tria3.o          \
     postProcessCol_tria6.o shape_tria6.o
	$(CC) $(LFLAGS) comFilter.o depSolver.o depolarization.o dotProd.o      \
     doubleMatrix.o doublePointer.o doubleVector.o electricFormA_tria3.o      \
     electricFormFullA_tria3.o elemType.o errorHandler.o field_tria3.o        \
     forceEllipsoid_tria3.o forceMST_tria3.o forceMultipole_tria3.o           \
     freeDoubleMatrix.o freeDoublePointer.o freeUintMatrix.o gaussBksb.o      \
     getLocalNormal_tria3.o getNormal_tria3.o initRes_el.o intDE_tria3.o      \
     intF_tria3.o intG_tria3.o intH_tria3.o intSingularG_tria3.o              \
     intSingularH_tria3.o iterGMRES_el.o L2Norm.o matVectProd_el.o            \
     postProcess_tria3.o potential_tria3.o shape_line2.o shape_tria3.o        \
     solverGMRES_el.o uintMatrix.o uintVector.o X2L_line2.o X2L_tria3.o       \
     X2L_tria6.o electricFormA_tria6.o electricFormFullA_tria6.o field_tria6.o\
     potential_tria6.o forceEllipsoid_tria6.o forceMST_tria6.o                \
     forceMultipole_tria6.o getLocalNormal_tria6.o getNormal_tria6.o          \
     intDE_tria6.o intF_tria6.o intG_tria6.o intH_tria6.o intSingularG_tria6.o\
     intSingularH_tria6.o postProcess_tria6.o postProcessCol_tria3.o          \
     postProcessCol_tria6.o shape_tria6.o -o depSolver
	mv depSolver ../bin/depSolver

depSolver-sft: comFilter.o depSolver-shift.o depolarization.o dotProd.o     \
     doubleMatrix.o doublePointer.o doubleVector.o electricFormA_tria3.o      \
     electricFormFullA_tria3.o elemType.o errorHandler.o field_tria3.o        \
     forceEllipsoid_tria3.o forceMST_tria3.o forceMultipole_tria3.o           \
     freeDoubleMatrix.o freeDoublePointer.o freeUintMatrix.o gaussBksb.o      \
     getLocalNormal_tria3.o getNormal_tria3.o initRes_el.o                    \
     intDE_tria3.o intF_tria3.o intG_tria3.o intH_tria3.o intSingularG_tria3.o\
     intSingularH_tria3.o iterGMRES_el.o L2Norm.o matVectProd_el.o            \
     postProcess_tria3.o potential_tria3.o shape_line2.o shape_tria3.o        \
     solverGMRES_el.o uintMatrix.o uintVector.o X2L_line2.o X2L_tria3.o       \
     X2L_tria6.o electricFormA_tria6.o electricFormFullA_tria6.o field_tria6.o\
     potential_tria6.o forceEllipsoid_tria6.o forceMST_tria6.o                \
     forceMultipole_tria6.o getLocalNormal_tria6.o getNormal_tria6.o          \
     intDE_tria6.o intF_tria6.o intG_tria6.o intH_tria6.o intSingularG_tria6.o\
     intSingularH_tria6.o postProcess_tria6.o postProcessCol_tria3.o          \
     postProcess_tria6.o postProcessCol_tria3.o postProcessCol_tria6.o        \
     shape_tria6.o
	$(CC) $(LFLAGS) comFilter.o depSolver-shift.o depolarization.o dotProd.o\
     doubleMatrix.o doublePointer.o doubleVector.o electricFormA_tria3.o      \
     electricFormFullA_tria3.o elemType.o errorHandler.o field_tria3.o        \
     forceEllipsoid_tria3.o forceMST_tria3.o forceMultipole_tria3.o           \
     freeDoubleMatrix.o freeDoublePointer.o freeUintMatrix.o gaussBksb.o      \
     getLocalNormal_tria3.o getNormal_tria3.o initRes_el.o intDE_tria3.o      \
     intF_tria3.o intG_tria3.o intH_tria3.o intSingularG_tria3.o              \
     intSingularH_tria3.o iterGMRES_el.o L2Norm.o matVectProd_el.o            \
     postProcess_tria3.o potential_tria3.o shape_line2.o shape_tria3.o        \
     solverGMRES_el.o uintMatrix.o uintVector.o X2L_line2.o X2L_tria3.o       \
     X2L_tria6.o electricFormA_tria6.o electricFormFullA_tria6.o field_tria6.o\
     potential_tria6.o forceEllipsoid_tria6.o forceMST_tria6.o                \
     forceMultipole_tria6.o  getLocalNormal_tria6.o getNormal_tria6.o         \
     intDE_tria6.o intF_tria6.o intG_tria6.o intH_tria6.o intSingularG_tria6.o\
     intSingularH_tria6.o postProcess_tria6.o postProcessCol_tria3.o          \
     postProcessCol_tria6.o shape_tria6.o -o depSolver-sft
	mv depSolver-sft ../bin/depSolver-sft

# Compiling Step
comFilter.o: comFilter.c
	$(CC) -c $(CFLAGS) comFilter.c -o comFilter.o

depSolver.o: depSolver.c depSolver.h gaussData.h
	$(CC) -c $(CFLAGS) depSolver.c -o depSolver.o

depSolver-shift.o: depSolver-shift.c depSolver.h gaussData.h
	$(CC) -c $(CFLAGS) depSolver-shift.c -o depSolver-shift.o

depolarization.o: depolarization.c depolarization.h gaussData.h
	$(CC) -c $(CFLAGS) depolarization.c -o depolarization.o

dotProd.o: dotProd.c
	$(CC) -c $(CFLAGS) dotProd.c -o dotProd.o

doubleMatrix.o: doubleMatrix.c 
	$(CC) -c $(CFLAGS) doubleMatrix.c -o doubleMatrix.o

doublePointer.o: doublePointer.c
	$(CC) -c $(CFLAGS) doublePointer.c -o doublePointer.o

doubleVector.o: doubleVector.c
	$(CC) -c $(CFLAGS) doubleVector.c -o doubleVector.o

electricFormA_tria3.o: electricFormA_tria3.c electricFormA_tria3.h constants.h
	$(CC) -c $(CFLAGS) electricFormA_tria3.c -o electricFormA_tria3.o

electricFormA_tria6.o: electricFormA_tria6.c electricFormA_tria6.h constants.h
	$(CC) -c $(CFLAGS) electricFormA_tria6.c -o electricFormA_tria6.o

electricFormFullA_tria3.o: electricFormFullA_tria3.c electricFormA_tria3.h constants.h
	$(CC) -c $(CFLAGS) electricFormFullA_tria3.c -o electricFormFullA_tria3.o

electricFormFullA_tria6.o: electricFormFullA_tria6.c electricFormA_tria6.h constants.h
	$(CC) -c $(CFLAGS) electricFormFullA_tria6.c -o electricFormFullA_tria6.o

elemType.o: elemType.c
	$(CC) -c $(CFLAGS) elemType.c -o elemType.o

errorHandler.o: errorHandler.c
	$(CC) -c $(CFLAGS) errorHandler.c -o errorHandler.o

field_tria3.o: field_tria3.c field_tria3.h constants.h
	$(CC) -c $(CFLAGS) field_tria3.c -o field_tria3.o

field_tria6.o: field_tria6.c field_tria6.h constants.h
	$(CC) -c $(CFLAGS) field_tria6.c -o field_tria6.o

forceEllipsoid_tria3.o: forceEllipsoid_tria3.c forceEllipsoid_tria3.h constants.h
	$(CC) -c $(CFLAGS) forceEllipsoid_tria3.c -o forceEllipsoid_tria3.o

forceEllipsoid_tria6.o: forceEllipsoid_tria6.c forceEllipsoid_tria6.h constants.h
	$(CC) -c $(CFLAGS) forceEllipsoid_tria6.c -o forceEllipsoid_tria6.o

forceMST_tria3.o: forceMST_tria3.c forceMST_tria3.h
	$(CC) -c $(CFLAGS) forceMST_tria3.c -o forceMST_tria3.o

forceMST_tria6.o: forceMST_tria6.c forceMST_tria6.h
	$(CC) -c $(CFLAGS) forceMST_tria6.c -o forceMST_tria6.o

forceMultipole_tria3.o: forceMultipole_tria3.c forceMultipole_tria3.h constants.h
	$(CC) -c $(CFLAGS) forceMultipole_tria3.c -o forceMultipole_tria3.o

forceMultipole_tria6.o: forceMultipole_tria6.c forceMultipole_tria6.h constants.h
	$(CC) -c $(CFLAGS) forceMultipole_tria6.c -o forceMultipole_tria6.o

freeDoubleMatrix.o: freeDoubleMatrix.c
	$(CC) -c $(CFLAGS) freeDoubleMatrix.c -o freeDoubleMatrix.o

freeDoublePointer.o: freeDoublePointer.c
	$(CC) -c $(CFLAGS) freeDoublePointer.c -o freeDoublePointer.o

freeUintMatrix.o: freeUintMatrix.c
	$(CC) -c $(CFLAGS) freeUintMatrix.c -o freeUintMatrix.o

gaussBksb.o: gaussBksb.c
	$(CC) -c $(CFLAGS) gaussBksb.c -o gaussBksb.o

getLocalNormal_tria3.o: getLocalNormal_tria3.c
	$(CC) -c $(CFLAGS) getLocalNormal_tria3.c -o getLocalNormal_tria3.o

getLocalNormal_tria6.o: getLocalNormal_tria6.c
	$(CC) -c $(CFLAGS) getLocalNormal_tria6.c -o getLocalNormal_tria6.o

getNormal_tria3.o: getNormal_tria3.c
	$(CC) -c $(CFLAGS) getNormal_tria3.c -o getNormal_tria3.o

getNormal_tria6.o: getNormal_tria6.c
	$(CC) -c $(CFLAGS) getNormal_tria6.c -o getNormal_tria6.o

initRes_el.o: initRes_el.c initRes_el.h
	$(CC) -c $(CFLAGS) initRes_el.c -o initRes_el.o

intDE_tria3.o: intDE_tria3.c integral_tria3.h gaussData.h
	$(CC) -c $(CFLAGS) intDE_tria3.c -o intDE_tria3.o

intDE_tria6.o: intDE_tria6.c integral_tria6.h gaussData.h
	$(CC) -c $(CFLAGS) intDE_tria6.c -o intDE_tria6.o

intF_tria3.o: intF_tria3.c integral_tria3.h gaussData.h
	$(CC) -c $(CFLAGS) intF_tria3.c -o intF_tria3.o

intF_tria6.o: intF_tria6.c integral_tria6.h gaussData.h
	$(CC) -c $(CFLAGS) intF_tria6.c -o intF_tria6.o

intG_tria3.o: intG_tria3.c integral_tria3.h gaussData.h
	$(CC) -c $(CFLAGS) intG_tria3.c -o intG_tria3.o

intG_tria6.o: intG_tria6.c integral_tria6.h gaussData.h
	$(CC) -c $(CFLAGS) intG_tria6.c -o intG_tria6.o

intH_tria3.o: intH_tria3.c integral_tria3.h gaussData.h
	$(CC) -c $(CFLAGS) intH_tria3.c -o intH_tria3.o

intH_tria6.o: intH_tria6.c integral_tria6.h gaussData.h
	$(CC) -c $(CFLAGS) intH_tria6.c -o intH_tria6.o

intSingularG_tria3.o: intSingularG_tria3.c integral_tria3.h gaussData.h
	$(CC) -c $(CFLAGS) intSingularG_tria3.c -o intSingularG_tria3.o

intSingularG_tria6.o: intSingularG_tria6.c integral_tria6.h gaussData.h
	$(CC) -c $(CFLAGS) intSingularG_tria6.c -o intSingularG_tria6.o

intSingularH_tria3.o: intSingularH_tria3.c integral_tria3.h gaussData.h
	$(CC) -c $(CFLAGS) intSingularH_tria3.c -o intSingularH_tria3.o

intSingularH_tria6.o: intSingularH_tria6.c integral_tria6.h
	$(CC) -c $(CFLAGS) intSingularH_tria6.c -o intSingularH_tria6.o

iterGMRES_el.o: iterGMRES_el.c iterGMRES_el.h
	$(CC) -c $(CFLAGS) iterGMRES_el.c -o iterGMRES_el.o

L2Norm.o: L2Norm.c
	$(CC) -c $(CFLAGS) L2Norm.c -o L2Norm.o

matVectProd_el.o: matVectProd_el.c
	$(CC) -c $(CFLAGS) matVectProd_el.c -o matVectProd_el.o

postProcessCol_tria3.o: postProcessCol_tria3.c postProcessCol_tria3.h
	$(CC) -c $(CFLAGS) postProcessCol_tria3.c -o postProcessCol_tria3.o

postProcessCol_tria6.o: postProcessCol_tria6.c postProcessCol_tria6.h
	$(CC) -c $(CFLAGS) postProcessCol_tria6.c -o postProcessCol_tria6.o

postProcess_tria3.o: postProcess_tria3.c postProcess_tria3.h
	$(CC) -c $(CFLAGS) postProcess_tria3.c -o postProcess_tria3.o

postProcess_tria6.o: postProcess_tria6.c postProcess_tria6.h
	$(CC) -c $(CFLAGS) postProcess_tria6.c -o postProcess_tria6.o

potential_tria3.o: potential_tria3.c potential_tria3.h constants.h
	$(CC) -c $(CFLAGS) potential_tria3.c -o potential_tria3.o

potential_tria6.o: potential_tria6.c potential_tria6.h
	$(CC) -c $(CFLAGS) potential_tria6.c -o potential_tria6.o

shape_line2.o: shape_line2.c
	$(CC) -c $(CFLAGS) shape_line2.c -o shape_line2.o

shape_tria3.o: shape_tria3.c
	$(CC) -c $(CFLAGS) shape_tria3.c -o shape_tria3.o

shape_tria6.o: shape_tria6.c
	$(CC) -c $(CFLAGS) shape_tria6.c -o shape_tria6.o

solverGMRES_el.o: solverGMRES_el.c solverGMRES_el.h
	$(CC) -c $(CFLAGS) solverGMRES_el.c -o solverGMRES_el.o

uintMatrix.o: uintMatrix.c
	$(CC) -c $(CFLAGS) uintMatrix.c -o uintMatrix.o

uintVector.o: uintVector.c
	$(CC) -c $(CFLAGS) uintVector.c -o uintVector.o

X2L_line2.o: X2L_line2.c
	$(CC) -c $(CFLAGS) X2L_line2.c -o X2L_line2.o

X2L_tria3.o: X2L_tria3.c
	$(CC) -c $(CFLAGS) X2L_tria3.c -o X2L_tria3.o

X2L_tria6.o: X2L_tria6.c
	$(CC) -c $(CFLAGS) X2L_tria6.c -o X2L_tria6.o

# Eliminate all object files
clean:
	rm -f *.o

