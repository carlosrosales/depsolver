/******************************************************************************
* File      : intSingularG_tria6.c                                            *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* Computes the integral of the Green's function G when a collocation point is *
* one of the nodes of the integration element and returns the values for all  *
* shape functions in Int[]. Uses a regularization transformation that changes *
* the triangle into a degenerated square and Gauss-Jacobi quadrature.         *
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

int intSingularG_tria6(unsigned int SinNode, double X[6][3], double *Xeval, 
                       double *Int)
{
    const unsigned int NODES_IN_ELEM = 6;
    unsigned int i, j, k, T;
    double A, B, dx, dy, dz, J, U1, U2, W1, W2;
    double L1[2], L2[2], N1[6], N2[6], XL[3];
 
    /* Initialize */
    for(i = 0; i < NODES_IN_ELEM; i++) Int[i] = 0.0;

    /* Choose number of subelements */
    if(SinNode%2 == 1){
        T = 1;
        B = 0.125;
    }
    else{
        T = 2;
        B = 0.0625;
    }

    /* Integrate using Gauss-Legendre and Gauss-Jacobi quadratures */
    for(i = 0; i < NGAUSS; i++){
        U1 = Gauss[i][0];
        W1 = Gauss[i][1];
    
        for(j = 0; j < NGAUSS; j++){
            U2 = GaussJacobi[j][0];
            W2 = GaussJacobi[j][1];

            /* Linear transformation to get the singularity at node 1 */
            if(SinNode == 1){
                L1[0] = 0.5*(1.0 - U2);
                L1[1] = 0.25*(1.0 + U1)*(1.0 + U2);
            }
            else if(SinNode == 2){
                L1[0] = 0.25*(1.0 - U2);
                L1[1] = 0.5 + 0.25*U1*(1.0 + U2);
                L2[0] = 0.5 - 0.25*U1*(1.0 + U2);
                L2[1] = 0.25*(1.0 - U2);
            }
            else if(SinNode == 3){
                L1[0] = 0.25*(1.0 - U1)*(1.0 + U2);
                L1[1] = 0.5*(1.0 - U2);
            }
            else if(SinNode == 4){
                L1[0] = 0.25*(1.0 + U1)*(1.0 + U2);
                L1[1] = 0.5 - 0.25*U1*(1.0 + U2);
                L2[0] = 0.25*(1.0 - U1)*(1.0 + U2);
                L2[1] = 0.25*(1.0 - U2);
            }
            else if(SinNode == 5){
                L1[0] = 0.25*(1.0 + U1)*(1.0 + U2);
                L1[1] = 0.25*(1.0 - U1)*(1.0 + U2);
            }
            else{
                L1[0] = 0.5 + 0.25*U1*(1.0 + U2);
                L1[1] = 0.25*(1.0 - U1)*(1.0 + U2);
                L2[0] = 0.25*(1 - U2);
                L2[1] = 0.25*(1.0 + U1)*(1.0 + U2);
            }

            /* Contribution from triangle 1 */
            shape_tria6(L1,N1);
            J = X2L_tria6(X,XL,L1,N1);
            dx = Xeval[0] - XL[0];
            dy = Xeval[1] - XL[1];
            dz = Xeval[2] - XL[2];
            A = W1*W2*J*B/sqrt(dx*dx + dy*dy + dz*dz);
            for(k = 0; k < NODES_IN_ELEM; k++) Int[k] += A*N1[k];

            /* Contribution from triangle 2 (only if necessary) */
            if(T == 2){
                shape_tria6(L2,N2);
                J = X2L_tria6(X,XL,L2,N2);
                dx = Xeval[0] - XL[0];
                dy = Xeval[1] - XL[1];
                dz = Xeval[2] - XL[2];
                A = W1*W2*J*B/sqrt(dx*dx + dy*dy + dz*dz);  
                for(k = 0; k < NODES_IN_ELEM; k++) Int[k] += A*N2[k];
            }
        }
    }

    return 0;
}
