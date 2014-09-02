/******************************************************************************
* File      : shape_tria6.c                                                   *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* Returns a vector N with the shape functions       3                         *
* for quadratic interpolation in a isoparametric    |\                        *
* triangular element (6-noded triangle).           4| \2                      *
* Node order as in figure.                          |__\                      *
*                                                  5  6  1                    *
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

int shape_tria6(double *L, double *N)
{
    double L1, L2, L3;

    /* Auxiliar quantities */
    L1 = L[0];
    L2 = L[1];
    L3 = 1.0 - L1 - L2;

    /* Shape function definition for node order 1-2-3-4-5-6 */
    N[0] = L1*(2.0*L1 - 1.0);
    N[1] = 4.0*L1*L2;
    N[2] = L2*(2.0*L2 - 1.0);
    N[3] = 4.0*L2*L3;
    N[4] = L3*(2.0*L3 - 1.0);
    N[5] = 4.0*L1*L3;

    return 0;
}

