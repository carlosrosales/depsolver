/******************************************************************************
* File      : x2L_line2.c                                                     *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* Takes a vector N with the shape functions at the local variables and a      *
* matrix X with the global variables of the element and returns a vector XL   * 
* XL that is the transformed L in global space, and the Jacobian.             *
* Works for linear interpolation in line elements (2-noded lines).            *
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

double X2L_line2(double X[2][2], double *XL, double *N)
{
    double dx, dy;

    /* Global coordinates */
    XL[0] = N[0]*X[0][0] + N[1]*X[1][0];
    XL[1] = N[0]*X[0][1] + N[1]*X[1][1];

    /* Jacobian */
    dx = 0.5*(X[1][0] - X[0][0]);
    dy = 0.5*(X[1][1] - X[0][1]); 
    return sqrt(dx*dx+dy*dy);
}
