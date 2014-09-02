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
* Works for quadratic interpolation in triangular elements (6-noded triangles)*
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

double X2L_tria6(double X[6][3], double *XL, double *L, double *N)
{
    double ax, ay, az, dx1, dx2, dy1, dy2, dz1, dz2, j1, j2, j3, L1, L2;

    /* Global coordinates */
    XL[0] = N[0]*X[0][0] + N[1]*X[1][0] + N[2]*X[2][0] + N[3]*X[3][0] + N[4]*X[4][0] + N[5]*X[5][0];
    XL[1] = N[0]*X[0][1] + N[1]*X[1][1] + N[2]*X[2][1] + N[3]*X[3][1] + N[4]*X[4][1] + N[5]*X[5][1];
    XL[2] = N[0]*X[0][2] + N[1]*X[1][2] + N[2]*X[2][2] + N[3]*X[3][2] + N[4]*X[4][2] + N[5]*X[5][2];

    /* Auxiliar quantities */
    L1 = L[0];
    L2 = L[1];
    ax = X[1][0] - X[3][0] + X[4][0] - X[5][0];
    ay = X[1][1] - X[3][1] + X[4][1] - X[5][1];
    az = X[1][2] - X[3][2] + X[4][2] - X[5][2];

    /* Derivatives of (x,y,z) in global space with respect to L1 and L2 */
    dx1 = 4.0*(L1*(X[0][0] + X[4][0] - 2.0*X[5][0]) + L2*ax + X[5][0]) - X[0][0] - 3.0*X[4][0];
    dx2 = 4.0*(L1*ax + L2*(X[2][0] + X[4][0] - 2.0*X[3][0]) + X[3][0]) - X[2][0] - 3.0*X[4][0];
    dy1 = 4.0*(L1*(X[0][1] + X[4][1] - 2.0*X[5][1]) + L2*ay + X[5][1]) - X[0][1] - 3.0*X[4][1];
    dy2 = 4.0*(L1*ay + L2*(X[2][1] + X[4][1] - 2.0*X[3][1]) + X[3][1]) - X[2][1] - 3.0*X[4][1];
    dz1 = 4.0*(L1*(X[0][2] + X[4][2] - 2.0*X[5][2]) + L2*az + X[5][2]) - X[0][2] - 3.0*X[4][2];
    dz2 = 4.0*(L1*az + L2*(X[2][2] + X[4][2] - 2.0*X[3][2]) + X[3][2]) - X[2][2] - 3.0*X[4][2];

    /* Jacobian */
    j1 = dy1*dz2 - dy2*dz1;
    j2 = dx2*dz1 - dx1*dz2;
    j3 = dx1*dy2 - dx2*dy1;

    return sqrt(j1*j1 + j2*j2 + j3*j3);
}
