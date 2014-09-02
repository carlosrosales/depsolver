/******************************************************************************
* File      : matVectProd_el.c                                                *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* Computes the matrix-vector product y = A*x.This "_el" version has been      *
* modified to use a [nNodes]x[2xnNodes] matrix.                               *
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

extern unsigned int nNodes;

void matVectProd_el(double **A, double *x, double *y)
{
	int i, j;
   
  	for(i = 0; i < nNodes; i++){
		y[i] = 0.0;
        y[i + nNodes] = 0.0;
        for(j = 0; j < nNodes; j++){
            y[i] += A[i][j]*x[j] + A[i][j + nNodes]*x[j + nNodes];
            y[i + nNodes] += A[i][j]*x[j + nNodes] - A[i][j + nNodes]*x[j];
        }
  	}
}
