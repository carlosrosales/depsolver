/******************************************************************************
* File      : depolarization.c                                                *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* Takes the semiaxis lengths of an ellipsoid in array axis, and calculates    *
* the depolarization factors numerically using a high quality Gauss-Legendre  *
* quadrature, returning them in array Ld.                                     *
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

#include "depolarization.h"

int depolarization(double *axis, double *Ld)
{
    unsigned int  i;
    double  A, a2, b2, c2, J, L, R, W, x;
    double  N[2], XL[2];
    double  X[2][2];

    /* Initialize */
    Ld[0] = Ld[1] = Ld[2] = 0.0;
    X[0][0] = X[0][1] = X[1][1] = 0.0;
    A = 0.5*axis[0]*axis[1]*axis[2];
    a2 = axis[0]*axis[0];
    b2 = axis[1]*axis[1];
    c2 = axis[2]*axis[2];

    /* Integrate for Lx using Gaussian quadrature */
    X[1][0] = 100*a2;
    for(i = 0; i < 64; i ++){
        L = GaussH[i][0];
        W = GaussH[i][1];
    
        shape_line2(L,N);
        J = X2L_line2(X,XL,N);
        x = XL[0];
        R = J*W/sqrt((x+a2)*(x+b2)*(x+c2));
    
        Ld[0] += R/(x+a2);
    }

    /* Integrate for Ly using Gaussian quadrature */
    X[1][0] = 100*b2;
    for(i = 0; i < 64; i ++){
        L = GaussH[i][0];
        W = GaussH[i][1];
    
        shape_line2(L,N);
        J = X2L_line2(X,XL,N);
        x = XL[0];
        R = J*W/sqrt((x+a2)*(x+b2)*(x+c2));
    
        Ld[1] += R/(x+b2);
    }

    /* Integrate for Lz using Gaussian quadrature */
    X[1][0] = 100*c2;
    for(i = 0; i < 64; i ++){
        L = GaussH[i][0];
        W = GaussH[i][1];
    
        shape_line2(L,N);
        J = X2L_line2(X,XL,N);
        x = XL[0];
        R = J*W/sqrt((x+a2)*(x+b2)*(x+c2));
    
        Ld[2] += R/(x+c2);
    }


    Ld[0] = A*Ld[0];
    Ld[1] = A*Ld[1];
    Ld[2] = A*Ld[2];

    /* Safety checks */
    for(i = 0; i < 3; i++){
        if(Ld[i] < 0.0) errorHandler("Error: depolarization() - factors must be positive!");
    }
    if(fabs(1.0 - Ld[0] - Ld[1] - Ld[2]) > 0.05)
        printf("Warning: depolarization() - error is larger than 0.05");

    return 0;
}
