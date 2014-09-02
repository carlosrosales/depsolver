/******************************************************************************
* File      : postProcess_tria6.h                                             *
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
#include <math.h>
#include "constants.h"

/* Function prototypes */
void freeDoubleMatrix(double **M, unsigned int ROWS);

double **doubleMatrix(unsigned int ROWS, unsigned int COLS, unsigned int INIT);

int field_tria6(double *Xin, double **mNodes, unsigned int **mElems, double *vB, 
                double *mE);

int forceEllipsoid_tria6(unsigned int ANALYSIS, unsigned int nOrder, 
                         double *axis, double *Xeval, double **mNodes, 
                         unsigned int **mElems, double *vB, double **vMatParam, 
                         double *vProbParam, double *F);
                         
int forceMST_tria6(double **mNodes, unsigned int **mElems, 
                   unsigned int *vBCType, double *vB, double Eps, double **Fce, 
                   double **Xce, double *Fcm, double *Xcm);
                   
int forceMultipole_tria6(unsigned int nOrder, double R, double *Xeval, 
                         double **mNodes, unsigned int **mElems, double *vB, 
                         double **vMatParam, double *vProbParam, double *F);

int potential_tria6(double *Xin, double **mNodes, unsigned int **mElems, 
                    double *vB, double *mPot);

/* Global variables */
extern FILE *file_log;
extern unsigned int nElems, nInternalPoints, nFPoints;
