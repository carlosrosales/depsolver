/******************************************************************************
* File      : iterGMRES_el.h                                                  *
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

#define GMRES_RESTART   0   /* 0 = not activated; 1 = activated */
#if GMRES_RESTART
    #define MAX_OUTER_LOOP  10
    #define MAX_ITERS       20      
#else
    #define MAX_OUTER_LOOP  1
    #define MAX_ITERS       600     
#endif


/* Function prototypes */
void errorHandler(char errorText[]);

void freeDoubleMatrix(double **M, unsigned int ROWS);

void freeDoublePointer(double **P, unsigned int ROWS);

void initRes_el(unsigned int size, double **A, double *x, double *B, double *r);   
                                       
void matVectProd_el(double **A, double *x, double *y);

double dotProd(unsigned int size, double *x, double *y);

double L2Norm(unsigned int size, double *x);

double *doubleVector(unsigned int LENGTH, unsigned int INIT);

double **doubleMatrix(unsigned int ROWS, unsigned int COLS, unsigned int INIT);
                       
double **doublePointer(unsigned int LENGTH, unsigned int INIT);

/* Global variables */
extern FILE *fg;
