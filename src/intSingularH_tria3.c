/******************************************************************************
* File      : intSingularH_tria3.c                                            *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* Computes the integral of the Green's function normal derivative H for one   *
* element when the collocation point is one of the nodes of the integration   *
* element, and returns the values for all shape functions in vector Int[].    *
* Uses element subdivision and standard Gauss-Legendre quadrature integration.*
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

int intSingularH_tria3(const unsigned int FLAG, unsigned int SinNode, 
                       double X[3][3], double *Xeval, double *Int, 
                       double *normal)
{
    const unsigned int NODES_IN_ELEM = 3;
    unsigned int i, j, k, m;
    double A, B, C, dx, dy, dz, J, JL, nxdx, rsq, W, M1, M2, M3;
    double L[2], N[3], XL[3];
    double xs[3][3][2];
 
    /* Initialize */
    if(!FLAG) Int[0] = Int[1] = Int[2] = 0.0;
    else for(i = 0; i < 9; i++) Int[i] = 0.0;
    for(i = 0; i < 3; i++)
        for(j = 0; j < 3; j++)
            for(k = 0; k <2; k++) xs[i][j][k] = 0.0;

    /* Reorder so that singularity is at node 3 */
    if(SinNode == 1){
        for(i = 0; i < 3; i++) shift(X[2][i],X[0][i],X[1][i]);
    }
    else if(SinNode == 2){
        for(i = 0; i < 3; i++) shift(X[2][i],X[1][i],X[0][i]);
    }

    /* Setup 3 initial triangles */
    xs[0][0][0] = 1.0;  xs[0][0][1] = 0.0;
    xs[0][1][0] = 0.5;  xs[0][1][1] = 0.5;
    xs[0][2][0] = 0.5;  xs[0][2][1] = 0.0;

    xs[1][0][0] = 0.5;  xs[1][0][1] = 0.0;
    xs[1][1][0] = 0.5;  xs[1][1][1] = 0.5;
    xs[1][2][0] = 0.0;  xs[1][2][1] = 0.5;
    
    xs[2][0][0] = 0.5;  xs[2][0][1] = 0.5;
    xs[2][1][0] = 0.0;  xs[2][1][1] = 1.0;
    xs[2][2][0] = 0.0;  xs[2][2][1] = 0.5;
    
    /* Start subdivision loop */
    J = 1.0;
    for(m = 0; m < NSUBDIVISIONS; m++){
        J = J*0.25;
    
        /* Integrate over the 3 triangles furthest away from the singularity */
        for(i = 0; i < 3; i++){
            
            /* Integrate over each triangle using Gauss-Legendre quadrature */
            for(j = 0; j < TSNGAUSS; j++){
                M1 = TSGauss[j][0];
                M2 = TSGauss[j][1];
                M3 = 1.0 - M1 - M2;
                W = TSGauss[j][2];
                
                /**************************************************** 
                * Get L and N as a function of (M1,M2). This is the *
                * transformation between the unit triangle (M1,M2)  *
                * and the smaller subtriangle (L1,L2)               *
                ****************************************************/
                L[0] = xs[i][0][0]*M1 + xs[i][1][0]*M2 + xs[i][2][0]*M3;
                L[1] = xs[i][0][1]*M1 + xs[i][1][1]*M2 + xs[i][2][1]*M3;
                shape_tria3(L,N);
                JL = X2L_tria3(X,XL,L,N);

                /* Auxiliar quantities */
                dx = Xeval[0] - XL[0];
                dy = Xeval[1] - XL[1];
                dz = Xeval[2] - XL[2];
                rsq = dx*dx + dy*dy + dz*dz;
                A = J*JL*W/pow(rsq,1.5);

                if(!FLAG){
                    nxdx = dx*normal[0]+dy*normal[1]+dz*normal[2];
                    A = A*nxdx;
                    Int[0] -= A*N[0];
                    Int[1] -= A*N[1];
                    Int[2] -= A*N[2];
                }
                else{
                    C = A*dz;
                    B = A*dy;
                    A = A*dx;
                    for(k = 0; k < NODES_IN_ELEM; k++){
                        Int[k*3] += A*N[k];
                        Int[k*3+1] += B*N[k];
                        Int[k*3+2] += C*N[k];
                    }
                }
            }
        }
    
        /* Set triangles for new subdivision */
        for(i = 0; i < 3; i++) {
            for(j = 0; j < 3; j++){
                xs[i][j][0] = xs[i][j][0]*0.5;
                xs[i][j][1] = xs[i][j][1]*0.5;
            }
        }
    }

    /* Unscramble node contributions to return them in their right places */    
    if(FLAG == 0){
        if(SinNode == 1){
            shift(Int[1],Int[0],Int[2]);
        }
        else if(SinNode == 2){
            shift(Int[0],Int[1],Int[2]);
        }
    }
    else{
        if(SinNode == 1){
            for(i = 0; i < 3; i++) shift(Int[3+i],Int[i],Int[6+i]);
        }
        else if(SinNode == 2){
            for(i = 0; i < 3; i++) shift(Int[i],Int[3+i],Int[6+i]);
        }
    }

    return 0;
}



















