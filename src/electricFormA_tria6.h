/******************************************************************************
* File      : electricFormA_tria6.h                                           *
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

/* Function prototypes */
double *doubleVector(unsigned int LENGTH, unsigned int INIT);
int getNormal_tria6(unsigned int nNodeID, double **mNodes, 
                    unsigned int **mElems, double *normal);
int intG_tria6(double X[6][3], double *Xeval, double *Int);
int intH_tria6(const unsigned int FLAG, double X[6][3], double *Xeval, 
               double *Int, double *normal);
int intSingularG_tria6(unsigned int SinNode, double X[6][3], double *Xeval, 
                       double *Int);
int intSingularH_tria6(const unsigned int FLAG, unsigned int SinNode, 
                       double X[6][3], double *Xeval, double *Int, 
                       double *normal);

/* Global variables */
extern unsigned int nElems, nInterfaces, nNodes;
