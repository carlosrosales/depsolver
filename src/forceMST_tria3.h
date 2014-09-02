/******************************************************************************
* File      : forceMST_tria3.h                                                *
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Function prototypes */
double **doubleMatrix(unsigned int ROWS, unsigned int COLS, unsigned int INIT);

int field_tria3(double *Xin, double **mNodes, unsigned int **mElems, double *vB, 
                double *mE);
                
void freeDoubleMatrix(double **M, unsigned int ROWS);

int intF_tria3(double X[3][3], double E[3][3], double *Int);

int shape_tria3(double *L, double *N);

double X2L_tria3(double X[3][3], double *XL, double *L, double *N);

/* Global variables */
extern unsigned int nElems, nNodes;
