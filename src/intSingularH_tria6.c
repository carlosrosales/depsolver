/******************************************************************************
* File      : intSingularH_tria6.c                                            *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* Computes the integral of the Green's function normal derivative H for one   * 
* element when the collocation point is one of the nodes of the integration   *
* element, and returns the values for all shape functions in vector Int[].    *
* Uses element subdivision and standard Gauss-Legendre quadrature integration.*
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

int intSingularH_tria6(const unsigned int FLAG, unsigned int SinNode, 
                       double X[6][3], double *Xeval, double *Int, 
                       double *normal)
{
    const unsigned int NODES_IN_ELEM = 6;
    unsigned int i, j, k, m;
    double A, A1, A2, B, B1, B2, C, C1, C2, dx, dx1, dx2, dy, dy1, dy2, dz, 
            dz1, dz2, J, JL, JL1, JL2, dr, dr1, dr2, tmp, W, M1, M2, M3, U1, U2, U3;
    double L[2], L1[2], L2[2], N[6], N1[6], N2[6], XL[3], XL1[3], XL2[3];
    double xs[3][3][2];
 
    /* Initialize */
    if(FLAG == 0) for(i = 0; i < NODES_IN_ELEM; i++) Int[i] = 0.0;
    else for(i = 0; i < 18; i++) Int[i] = 0.0;

    /* Reorder so that singularity is at node 5 */
    if(SinNode == 1){
        for(i = 0; i < 3; i++){
            shift(X[4][i],X[0][i],X[2][i]);
            shift(X[5][i],X[1][i],X[3][i]);         
        }
    }
    else if(SinNode == 3){
        for(i = 0; i < 3; i++){
            shift(X[4][i],X[2][i],X[0][i]);
            shift(X[3][i],X[1][i],X[5][i]);
        }
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
                if(SinNode%2 == 1){
                    L[0] = xs[i][0][0]*M1 + xs[i][1][0]*M2 + xs[i][2][0]*M3;
                    L[1] = xs[i][0][1]*M1 + xs[i][1][1]*M2 + xs[i][2][1]*M3;
                    shape_tria6(L,N);
                    JL = X2L_tria6(X,XL,L,N);
                    
                    /* Auxiliar quantities */
                    dx = Xeval[0] - XL[0];
                    dy = Xeval[1] - XL[1];
                    dz = Xeval[2] - XL[2];
                    dr = sqrt(dx*dx + dy*dy + dz*dz);
                    A = J*JL*W/(dr*dr*dr);
                    
                    if(FLAG == 0){
                        A = A*(dx*normal[0] + dy*normal[1] + dz*normal[2]);
                        for(k = 0; k < NODES_IN_ELEM; k++) Int[k] -= A*N[k];
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
                else{
                    if(SinNode == 2){
                        U1 = 0.5*M3; U2 = M1 + 0.5*M3; U3 = 1.0 - U1 - U2;
                        L1[0] = xs[i][0][0]*U1 + xs[i][1][0]*U2 + xs[i][2][0]*U3;
                        L1[1] = xs[i][0][1]*U1 + xs[i][1][1]*U2 + xs[i][2][1]*U3;
                        U1 = M1 + 0.5*M3; U2 = M3; U3 = 1.0 - U1 - U2;
                        L2[0] = xs[i][0][0]*U1 + xs[i][1][0]*U2 + xs[i][2][0]*U3;
                        L2[1] = xs[i][0][1]*U1 + xs[i][1][1]*U2 + xs[i][2][1]*U3;
                    }
                    else if(SinNode == 4){
                        U1 = M1; U2 = 0.5*M3; U3 = 1.0 - U1 - U2;
                        L1[0] = xs[i][0][0]*U1 + xs[i][1][0]*U2 + xs[i][2][0]*U3;
                        L1[1] = xs[i][0][1]*U1 + xs[i][1][1]*U2 + xs[i][2][1]*U3;
                        U1 = M1; U2 = M2 + 0.5*M3; U3 = 1.0 - U1 - U2;
                        L2[0] = xs[i][0][0]*U1 + xs[i][1][0]*U2 + xs[i][2][0]*U3;
                        L2[1] = xs[i][0][1]*U1 + xs[i][1][1]*U2 + xs[i][2][1]*U3;
                    }
                    else{
                        U1 = 0.5*M3; U2 = M1; U3 = 1.0 - U1 - U2;
                        L1[0] = xs[i][0][0]*U1 + xs[i][1][0]*U2 + xs[i][2][0]*U3;
                        L1[1] = xs[i][0][1]*U1 + xs[i][1][1]*U2 + xs[i][2][1]*U3;
                        U1 = M1 + 0.5*M3; U2 = M2; U3 = 1.0 - U1 - U2;
                        L2[0] = xs[i][0][0]*U1 + xs[i][1][0]*U2 + xs[i][2][0]*U3;
                        L2[1] = xs[i][0][1]*U1 + xs[i][1][1]*U2 + xs[i][2][1]*U3;
                    }
                    shape_tria6(L1,N1);
                    shape_tria6(L2,N2);
                    JL1 = X2L_tria6(X,XL1,L1,N1);
                    JL2 = X2L_tria6(X,XL2,L2,N2);   

                    /* Auxiliar quantities */
                    dx1 = Xeval[0] - XL1[0];
                    dy1 = Xeval[1] - XL1[1];
                    dz1 = Xeval[2] - XL1[2];
                    dx2 = Xeval[0] - XL2[0];
                    dy2 = Xeval[1] - XL2[1];
                    dz2 = Xeval[2] - XL2[2];
                    dr1 = sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1);
                    dr2 = sqrt(dx2*dx2 + dy2*dy2 + dz2*dz2);
                    A1 = 0.5*J*JL1*W/(dr1*dr1*dr1);
                    A2 = 0.5*J*JL2*W/(dr2*dr2*dr2);
                    
                    if(FLAG == 0){
                        A1 = A1*(dx1*normal[0] + dy1*normal[1] + dz1*normal[2]);
                        A2 = A2*(dx2*normal[0] + dy2*normal[1] + dz2*normal[2]);

                        for(k = 0; k < NODES_IN_ELEM; k++) Int[k] -= A1*N1[k] + A2*N2[k];
                    }
                    else{
                        C1 = A1*dz1;
                        B1 = A1*dy1;
                        A1 = A1*dx1;
                        C2 = A2*dz2;
                        B2 = A2*dy2;
                        A2 = A2*dx2;
                        for(k = 0; k < NODES_IN_ELEM; k++){
                            Int[k*3] += A1*N1[k] + A2*N2[k];
                            Int[k*3+1] += B1*N1[k] + B2*N2[k];
                            Int[k*3+2] += C1*N1[k] + B2*N2[k];
                        }
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
            shift(Int[2],Int[0],Int[4]);
            shift(Int[3],Int[1],Int[5]);            
        }
        else if(SinNode == 3){
            shift(Int[0],Int[2],Int[4]);
            shift(Int[5],Int[1],Int[3]);
        }
    }
    else{
        if(SinNode == 1){
            for(i = 0; i < 3; i++){
                shift(Int[6+i],Int[i],Int[12+i]);
                shift(Int[9+i],Int[3+i],Int[15+i]);         
            }
        }
        else if(SinNode == 3){
            for(i = 0; i < 3; i++){
                shift(Int[i],Int[6+i],Int[12+i]);
                shift(Int[15+i],Int[3+i],Int[9+i]);
            }
        }
    }

    return 0;
}

