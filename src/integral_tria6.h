/******************************************************************************
* File      : integral_tria6.h                                                *
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

#define NGAUSS  8
#define TNGAUSS 7
#define TSNGAUSS    7
#define NSUBDIVISIONS   8

#define shift(a,b,c){double d = (a); (a) = (b); (b) = (c); (c) = d;}
#define swap(a,b){double c = (a); (a) = (b); (b) = c;}

/* Function prototypes */
int shape_tria6(double *L, double *N);
int getLocalNormal_tria6(double *L, double X[6][3], double *normal);
double X2L_tria6(double X[6][3], double *XL, double *L, double *N);

/* Global variables */
extern const double TGauss[TNGAUSS][3], TSGauss[TSNGAUSS][3];
extern const double Gauss[NGAUSS][2], GaussJacobi[NGAUSS][2];
