/******************************************************************************
* File      : intH_tria6.c                                                    *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* Computes the integral of the Green's function normal derivative H for one   *
* element and returns the integrals in vector Int[] for all shape functions.  *
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

#include "integral_tria6.h"

int intH_tria6(const unsigned int FLAG, double X[6][3], double *Xeval, 
               double *Int, double *normal)
{
    const unsigned int NODES_IN_ELEM = 6;
    unsigned int i, j;
    double A, B, C, dx, dy, dz, J, nxdx, rsq, W;
    double L[2], N[6], XL[3];

    /* Initialize */
    if(!FLAG) for(i = 0; i < NODES_IN_ELEM; i++) Int[i] = 0.0;
    else for(i = 0; i < 18; i++) Int[i] = 0.0;

    /* Integrate using Gauss-Legendre quadrature */
    for(i = 0; i < TNGAUSS; i++){
        L[0] = TGauss[i][0];
        L[1] = TGauss[i][1];
        W = TGauss[i][2];

        /* Shape functions N and Jacobian J for these (L1,L2) values */
        shape_tria6(L,N);       
        J = X2L_tria6(X,XL,L,N);

        /* Auxiliar quantities */
        dx = Xeval[0] - XL[0];
        dy = Xeval[1] - XL[1];
        dz = Xeval[2] - XL[2];
        rsq = dx*dx + dy*dy + dz*dz;
        A = J*W/pow(rsq,1.5);

        /* Add contribution for each node in this element */
        if(!FLAG){
            nxdx = dx*normal[0] + dy*normal[1] + dz*normal[2];
            A = A*nxdx;
            for(j = 0; j < NODES_IN_ELEM; j++) Int[j] -= A*N[j];
        }
        else{
            C = A*dz;
            B = A*dy;
            A = A*dx;
            for(j = 0; j < NODES_IN_ELEM; j++){
                Int[j*3] += A*N[j];
                Int[j*3+1] += B*N[j];
                Int[j*3+2] += C*N[j];
            }          
        }
    }

    return 0;
}
