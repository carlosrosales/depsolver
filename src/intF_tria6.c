/******************************************************************************
* File      : intF_tria6.c                                                    *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* Computes the integral of the Maxwell Stress Tensor for one element and      *
* returns the integrals in vector Int[] for the three spatial dimentions.     *
* Works for quadratic interpolation in triangular elements (6-noded triangles)*
* Int[0] = Fx    |    Int[1] = Fy    |    Int[2] = Fz                         *
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

int intF_tria6(double X[6][3], double E[6][3], double *Int)
{
    unsigned int i;
    double A, Ex, Ey, Ez, J, W, Esq, NE;
    double L[2], N[6], normal[3], XL[3];

    /* Initialize */
    Int[0] = Int[1] = Int[2] = 0.0;  

    /* Integrate using Gauss-Legendre quadrature */
    for(i = 0; i < TNGAUSS; i++){
        L[0] = TGauss[i][0];
        L[1] = TGauss[i][1];
        W = TGauss[i][2];
    
        /* Shape functions N and Jacobian J for these (L1,L2) values */
        shape_tria6(L,N);       
        J = X2L_tria6(X,XL,L,N);
        getLocalNormal_tria6(L,X,normal);

        /* Electric field value at this integration point */
        Ex = E[0][0]*N[0] + E[1][0]*N[1] + E[2][0]*N[2] + E[3][0]*N[3] + E[4][0]*N[4] + E[5][0]*N[5];
        Ey = E[0][1]*N[0] + E[1][1]*N[1] + E[2][1]*N[2] + E[3][1]*N[3] + E[4][1]*N[4] + E[5][1]*N[5];
        Ez = E[0][2]*N[0] + E[1][2]*N[1] + E[2][2]*N[2] + E[3][2]*N[3] + E[4][2]*N[4] + E[5][2]*N[5];
    
        /* Auxiliar quantities */
        NE = Ex*normal[0] + Ey*normal[1] + Ez*normal[2];
        Esq = Ex*Ex + Ey*Ey + Ez*Ez;
        A = J*W;

        /* Add contribution from this integration point */
        Int[0] += A*(Ex*NE - 0.5*Esq*normal[0]);
        Int[1] += A*(Ey*NE - 0.5*Esq*normal[1]);
        Int[2] += A*(Ez*NE - 0.5*Esq*normal[2]);
    }

    return 0;
}
