/******************************************************************************
* File      : depSolver.h                                                     *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
******************************************************************************/

/******************************************************************************
* COPYRIGHT & LICENSE INFORMATION                                             *
*                                                                             *
* Copyright 2006 Carlos Rosales Fernandez and The Institute of High           *
* Performance Computing (A*STAR)                                              *
*                                                                             *
* This file is part of depSolver.                                             *
*                                                                             *
* depSolver is free software; you can redistribute it and/or modify           *
* it under the terms of the GNU General Public License as published by        *
* the Free Software Foundation; either version 2 of the License, or           *
* (at your option) any later version.                                         *
*                                                                             *
* depSolver is distributed in the hope that it will be useful,                *
* but WITHOUT ANY WARRANTY; without even the implied warranty of              *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               *
* GNU General Public License for more details.                                *
*                                                                             *
* You should have received a copy of the GNU General Public License           *
* along with depSolver; if not, write to the Free Software                    *
* Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA  *
******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>

/* Function prototypes */
void errorHandler(char errorText[]);

void comFilter(char *FileName);

void freeDoubleMatrix(double **M, unsigned int ROWS);

void freeUintMatrix(unsigned int **M, unsigned int ROWS);

void solverGMRES_el(unsigned int preCond, unsigned int nInit, 
                    double **Amatrix, double *Rhs);

unsigned int *uintVector(unsigned int LENGTH, unsigned int INIT);

unsigned int **uintMatrix(unsigned int ROWS, unsigned int COLS, 
                            unsigned int INIT);
                            
int elemType();

int electricFormA_tria3(double **mNodes, unsigned int **mElems, 
                        unsigned int *vBCType, unsigned int **vInterfaces, 
                        double **vMatParam, double *vProbParam, double **mA, 
                        double **mBC, double *vB);
                        
int electricFormA_tria6(double **mNodes, unsigned int **mElems, 
                        unsigned int *vBCType, unsigned int **vInterfaces, 
                        double **vMatParam, double *vProbParam, double **mA, 
                        double **mBC, double *vB);
                            
int electricFormFullA_tria3(double **mNodes, unsigned int **mElems, 
                            unsigned int *vBCType, 
                            unsigned int **vInterfaces, double **vMatParam, 
                            double *vProbParam, double **mA, double **mBC, 
                            double *vB);  
                            
int electricFormFullA_tria6(double **mNodes, unsigned int **mElems, 
                            unsigned int *vBCType, 
                            unsigned int **vInterfaces, double **vMatParam, 
                            double *vProbParam, double **mA, double **mBC, 
                            double *vB);             
                      
int postProcessCol_tria3(unsigned int ANALYSIS, double *axis, double **XF, 
                         double **Xinner, double **mNodes, 
                         unsigned int **mElems, double **vMatParam, 
                         double *vProbParam, unsigned int *vBCType, 
                         double **xCols, double *vB);
                         
int postProcessCol_tria6(unsigned int ANALYSIS, double *axis, double **XF, 
                         double **Xinner, double **mNodes, 
                         unsigned int **mElems, double **vMatParam, 
                         double *vProbParam, unsigned int *vBCType, 
                         double **xCols, double *vB);

int postProcess_tria3(unsigned int ANALYSIS, double *axis, double **XF, 
                      double **Xinner, double **mNodes, 
                      unsigned int **mElems, double **vMatParam, 
                      double *vProbParam, unsigned int *vBCType, double *vB);
                      
int postProcess_tria6(unsigned int ANALYSIS, double *axis, double **XF, 
                      double **Xinner, double **mNodes, 
                      unsigned int **mElems, double **vMatParam, 
                      double *vProbParam, unsigned int *vBCType, double *vB);     
                      
int gaussBksb(unsigned int N, double **A, double *B);

double *doubleVector(unsigned int LENGTH, unsigned int INIT);

double **doubleMatrix(unsigned int ROWS, unsigned int COLS, 
                       unsigned int INIT);


/* Global variables */
FILE *file_log;
char cElemType[6];
unsigned int nDim, nElems, nElemType, nFPoints, nInternalPoints, nInterfaces, 
               nMats, nNodes, nNodesInElem, nCols, nColType;
