/******************************************************************************
* File      : intSingularG_tria3.c                                            *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* Computes the integral of the Green's function G when a collocation point is *
* one of the nodes of the integration element and returns the values for all  *
* shape functions in Int[]. Uses a regularization transformation that changes *
* the triangle into a degenerated square and Gauss-Jacobi quadrature.         *
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

int intSingularG_tria3(unsigned int SinNode, double X[3][3], double *Xeval, 
                       double *Int)
{
    unsigned int i, j;
    double A, dx, dy, dz, J, U1, U2, W1, W2;
    double L[2], N[3], XL[3];
 
    /* Initialize */
    Int[0] = Int[1] = Int[2] = 0.0;

    /* Integrate using Gauss-Legendre and Gauss-Jacobi quadratures */
    for(i = 0; i < NGAUSS; i++){
        U1 = Gauss[i][0];
        W1 = Gauss[i][1];
    
        for(j = 0; j < NGAUSS; j++){
            U2 = GaussJacobi[j][0];
            W2 = GaussJacobi[j][1];

            /* Linear transformation to get the singularity at node 1 */
            if(SinNode == 1){
                L[0] = 0.5*(1.0 - U2);
                L[1] = 0.25*(1.0 + U1)*(1.0 + U2);
            }
            else if(SinNode == 2){
                L[0] = 0.25*(1.0 - U1)*(1.0 + U2);
                L[1] = 0.5*(1.0 - U2);
            }
            else{
                L[0] = 0.25*(1.0 + U1)*(1.0 + U2);
                L[1] = 0.25*(1.0 - U1)*(1.0 + U2);
            }

            /* Shape functions N and Jacobian J for these (L1,L2) values */
            shape_tria3(L,N);
            J = X2L_tria3(X,XL,L,N);
            
            /* Contribution from point (i,j) */
            dx = Xeval[0] - XL[0];
            dy = Xeval[1] - XL[1];
            dz = Xeval[2] - XL[2];
            A = W1*W2*J*0.125/sqrt(dx*dx + dy*dy + dz*dz);              
            Int[0] += A*N[0];
            Int[1] += A*N[1];
            Int[2] += A*N[2];
        }
    }

    return 0;
}
