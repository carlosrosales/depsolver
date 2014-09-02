/******************************************************************************
* File      : intG_tria3.c                                                    *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* Computes the integral of the Green's function G for one element and returns *
* the integrals in vector Int[] for all shape functions.                      *
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

#include "integral_tria3.h"

int intG_tria3(double X[3][3], double *Xeval, double *Int)
{
    unsigned int i;
    double A, dx, dy, dz, J, W;
    double L[2], N[3], XL[3];

    /* Initialize */
    Int[0] = Int[1] = Int[2] = 0.0; 

    /* Integrate using Gauss-Legendre quadrature */
    for(i = 0; i < TNGAUSS; i++){
        L[0] = TGauss[i][0];
        L[1] = TGauss[i][1];
        W = TGauss[i][2];

        /* Shape functions N and Jacobian J for these (L1,L2) values */
        shape_tria3(L,N);       
        J = X2L_tria3(X,XL,L,N);

        /* Auxiliar quantities */
        dx = Xeval[0] - XL[0];
        dy = Xeval[1] - XL[1];
        dz = Xeval[2] - XL[2];
        A = W*J/sqrt(dx*dx + dy*dy + dz*dz);

        /* Add contribution for each node in this element */
        Int[0] += A*N[0];
        Int[1] += A*N[1];
        Int[2] += A*N[2];  
    }

    return 0;
}
