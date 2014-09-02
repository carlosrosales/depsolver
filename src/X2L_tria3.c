/******************************************************************************
* File      : X2L_tria3.c                                                     *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* Takes a vector N with the shape functions, a vector L with the local        *
* variables, and a matrix X with the global variables of the element and      *
* returns a vector XL that is the transformed L in global space, and the      *
* Jacobian of the transformation.                                             *
* Works for linear interpolation in triangular elements (3-noded triangles).  *
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

double X2L_tria3(double X[3][3], double *XL, double *L, double *N)
{
    double dx1, dx2, dy1, dy2, dz1, dz2, j1, j2, j3;

    /* Global coordinates */
    XL[0] = N[0]*X[0][0] + N[1]*X[1][0] + N[2]*X[2][0];
    XL[1] = N[0]*X[0][1] + N[1]*X[1][1] + N[2]*X[2][1];
    XL[2] = N[0]*X[0][2] + N[1]*X[1][2] + N[2]*X[2][2];
    
    /* Derivatives of (x,y,z) in global space with respect to L1 and L2 */
    dx1 = X[0][0] - X[2][0];
    dx2 = X[1][0] - X[2][0];
    dy1 = X[0][1] - X[2][1];
    dy2 = X[1][1] - X[2][1];
    dz1 = X[0][2] - X[2][2];
    dz2 = X[1][2] - X[2][2];

    /* Jacobian */
    j1 = dy1*dz2 - dy2*dz1;
    j2 = dx2*dz1 - dx1*dz2;
    j3 = dx1*dy2 - dx2*dy1;

    return sqrt(j1*j1 + j2*j2 + j3*j3);
}
