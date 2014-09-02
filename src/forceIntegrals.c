/******************************************************************************
* File     : forceIntegrals.c                                                 *
* Author   : Carlos Rosales Fernandez (carlos.rosales.fernandez(at)gmail.com) *
* Revision : 2.0 (2008-06-16)                                                 *
******************************************************************************/
/**
 * @brief Integrals necessary in the calculation of the DEP force
 *
 * @file
 * Definitions of the four functions necessary for the integrals in the force
 * calculations of depSolver: intDE_tria3, intDE_tria6, intF_tria3, intF_tria6.
 */

/*******************************************************************************
* Copyright 2006, 2008 Carlos Rosales Fernandez and IHPC (A*STAR).             *
*                                                                              *
* This file is part of depSolver.                                              *
*                                                                              *
* depSolver is free software: you can redistribute it and/or modify it under   *
* the terms of the GNU GPL version 3 or (at your option) any later version.    *
*                                                                              *
* depSolver is distributed in the hope that it will be useful, but WITHOUT ANY *
* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS    *
* FOR A PARTICULAR PURPOSE. See the GNU General Public License for details.    *
*                                                                              *
* You should have received a copy of the GNU General Public License along with *
* depSolver. If not, see <http://gnu.org/licenses/gpl.html>.                   *
*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constants.h"
#include "forceIntegrals.h"
#include "globalToLocal.h"
#include "shapeFunctions.h"
#include "normals.h"

extern const double TGauss[TNGAUSS][3];

/**
 * Computes the integral of the 2nd, 3rd and 4th order derivatives of the 3D
 * Green's function for one element and returns the integrals in matrix
 * \a E for the 3 dimensions and each node. These are some examples:
 *
 * \a Efield[0][1][0][0] : Ex @ node 1 \n
 * \a Efield[0][0][1][0] : Ey @ node 1 \n
 * \a Efield[0][2][0][0] : dEx/dx @ node 1 \n
 * \a Efield[0][1][0][1] : dEx/dz = dEz/dx @ node 1 \n
 * \a Efield[2][0][1][3] : d3Ez/dydz2 = d3Ey/dz3 @ node 3 \n
 *
 * Where Ex = dG/dx, Ey = dG/dy and Ez = dG/dz and G is Green's function.
 *
 * Works for linear interpolation in triangular elements (3-noded triangles).
 *
 * @param nOrder : [ Input ]  Multipolar approximation order (implemented 1, 2, 3 only)
 * @param X      : [ Input ]  Coordinates of the nodes defining the element
 * @param Xeval  : [ Input ]  Coordinates of the evaluation point
 * @param Efield : [ Output ] Integrals of the differentials of the Green's function
 */
int intDE_tria3(unsigned int nOrder,
                double X[3][3],
                double *Xeval,
                double Efield[][7][7][7])
{
    const unsigned int NODES_IN_ELEM = 3;
    unsigned int i, j, k, m;
    double A, B, C, D, E, F, dx, dx2, dx4, dy, dy2, dy4, dz, dz2, dz4, J, rsq, W;
    double DE[80], L[2], N[3], XL[3];

    /* Initialize */
    for(i = 0; i < 80; i++) DE[i] = 0.0;
    for(i = 0; i < NODES_IN_ELEM; i++)
        for(j = 0; j < 7; j++)
            for(k = 0; k < 7; k++)
                for(m = 0; m < 7; m++) Efield[i][j][k][m] = 0.0;

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
        dx2 = dx*dx;
        dy2 = dy*dy;
        dz2 = dz*dz;
        rsq = dx2 + dy2 + dz2;
        A = J*W/pow(rsq,1.5);
        B = 3.0*A/rsq;

        /***********************************************/
        /*** GET CONTRIBUTION FROM EACH ELEMENT NODE ***/
        /***********************************************/
    
        /* First order derivatives of the electric field */
        DE[0] = A - B*dx2;  /* dEx/dx */
        DE[1] = -B*dx*dy;   /* dEx/dy */
        DE[2] = -B*dx*dz;   /* dEx/dz */
    
        DE[3] = A - B*dy2;  /* dEy/dy */
        DE[4] = -B*dy*dz;   /* dEy/dz */

        DE[5] = A - B*dz2;  /* dEz/dz */
    
        /* Second order derivatives of the electric field */
        if(nOrder > 1){
            C = 5.0*B/rsq;

            DE[6] = (C*dx2 - 3.0*B)*dx;     /* d2Ex/dx2  */
            DE[7] = (C*dx2 - B)*dy;         /* d2Ex/dxdy */
            DE[8] = (C*dx2 - B)*dz;         /* d2Ex/dxdz */
            DE[9] = (C*dy2 - B)*dx;         /* d2Ex/dy2  */
            DE[10] = C*dx*dy*dz;            /* d2Ex/dydz */
            DE[11] = (C*dz2 - B)*dx;        /* d2Ex/dzdz */
        
            DE[12] = (C*dy2 - 3.0*B)*dy;    /* d2Ey/dy2  */
            DE[13] = (C*dy2 - B)*dz;        /* d2Ey/dydz */
            DE[14] = (C*dz2 - B)*dy;        /* d2Ey/dz2  */
        
            DE[15] = (C*dz2 - 3.0*B)*dz;    /* d2Ez/dz2  */
        
            /* Third order derivatives of the electric field */
            if(nOrder > 2){
                D = 7.0*C/rsq;
            
                DE[16] = (-D*dx2 + 6.0*C)*dx2 - 3.0*B;      /* d3Ex/dx3    */
                DE[17] = (-D*dx2 + 3.0*C)*dx*dy;            /* d3Ex/dx2dy  */
                DE[18] = (-D*dx2 + 3.0*C)*dx*dz;            /* d3Ex/dx2dz  */
                DE[19] = -D*dx2*dy2 + C*(dx2 + dy2) - B;    /* d3Ex/dxdy2  */
                DE[20] = (-D*dx2 + C)*dy*dz;                /* d3Ex/dxdydz */
                DE[21] = -D*dx2*dz2 + C*(dx2 + dz2) - B;    /* d3Ex/dxdz2  */
                DE[22] = (-D*dy2 + 3.0*C)*dx*dy;            /* d3Ex/dy3    */
                DE[23] = (-D*dy2 + C)*dx*dz;                /* d3Ex/dy2dz  */
                DE[24] = (-D*dz2 + C)*dx*dy;                /* d3Ex/dydz2  */
                DE[25] = (-D*dz2 + 3.0*C)*dx*dz;            /* d3Ex/dz3    */
        
                DE[26] = (-D*dy2 + 6.0*C)*dy2 - 3.0*B;      /* d3Ey/dy3    */
                DE[27] = (-D*dy2 + 3.0*C)*dy*dz;            /* d3Ey/dy2dz  */
                DE[28] = -D*dy2*dz2 + C*(dy2 + dz2) - B;    /* d3Ey/dydz2  */
                DE[29] = (-D*dz2 + 3.0*C)*dy*dz;            /* d3Ey/dz3    */
        
                DE[30] = (-D*dz2 + 6.0*C)*dz2 - 3.0*B;      /* d3Ez/dz3    */
            
                /* Fourth order derivatives of the electric field */
                if(nOrder > 3){
                    E=9.0*D/rsq;

                    dx4 = dx2*dx2;
                    dy4 = dy2*dy2;
                    dz4 = dz2*dz2;
    
                    DE[31] = (E*dx4 - 10.0*D*dx2 + 15.0*C)*dx;              /* d4Ex/dx4     */
                    DE[32] = (E*dx4 - 6.0*D*dx2 + 3.0*C)*dy;                /* d4Ex/dx3dy   */
                    DE[33] = (E*dx4 - 6.0*D*dx2 + 3.0*C)*dz;                /* d4Ex/dx3dz   */
                    DE[34] = (E*dx2*dy2-D*(dx2 + 3.0*dy2) + 3.0*C)*dx;      /* d4Ex/dx2dy2  */
                    DE[35] = (E*dx2 - 3.0*D)*dx*dy*dz;                      /* d4Ex/dx2dydz */
                    DE[36] = (E*dx2*dz2 - D*(dx2 + 3.0*dz2) + 3.0*C)*dx;    /* d4Ex/dx2dz2  */
                    DE[37] = (E*dx2*dy2 - D*(dy2 + 3.0*dx2) + 3.0*C)*dy;    /* d4Ex/dxdy3   */
                    DE[38] = (E*dx2*dy2 - D*(dx2 + dy2) + C)*dz;            /* d4Ex/dxdy2dz */
                    DE[39] = (E*dx2*dz2 - D*(dx2 + dz2) + C)*dy;            /* d4Ex/dxdydz2 */
                    DE[40] = (E*dx2*dz2 - D*(dz2 + 3.0*dx2) + 3.0*C)*dz;    /* d4Ex/dxdz3   */
                    DE[41] = (E*dy4 - 6.0*D*dy2 + 3.0*C)*dx;                /* d4Ex/dy4     */
                    DE[42] = (E*dy2 - 3.0*D)*dx*dy*dz;                      /* d4Ex/dy3dz   */
                    DE[43] = (E*dy2*dz2 - D*(dy2 + dz2) + C)*dx;            /* d4Ex/dy2dz2  */
                    DE[44] = (E*dz2 - 3.0*D)*dx*dy*dz;                      /* d4Ex/dydz3   */
                    DE[45] = (E*dz4 - 6.0*D*dz2 + 3.0*C)*dx;                /* d4E/dz4      */
                
                    DE[46] = (E*dy4 - 10.0*D*dy2 + 15.0*C)*dy;              /* d4Ey/dy4     */
                    DE[47] = (E*dy4 - 6.0*D*dy2 + 3.0*C)*dz;                /* d4Ey/dy3dz   */
                    DE[48] = (E*dy2*dz2 - D*(dy2 + 3.0*dz2) + 3.0*C)*dy;    /* d4Ey/dy2dz2  */
                    DE[49] = (E*dy2*dz2 - D*(dz2 + 3.0*dy2) + 3.0*C)*dz;    /* d4Ey/dydz3   */
                    DE[50] = (E*dz4 - 6.0*D*dz2 + 3.0*C)*dy;                /* d4Ey/dz4     */
                
                    DE[51] = (E*dz4 - 10.0*D*dz2 + 15.0*C)*dz;              /* d4Ez/dz4     */
                
                    /* Fifth order derivatives of the electric field */
                    if(nOrder > 4){
                        F = 11.0*E/rsq;
                
                        DE[52] = -F*dx4*dx2 + 15.0*E*dx4 - 45.0*D*dx2 + 15.0*C;                         /* d5Ex/dx5      */
                        DE[53] = (-F*dx4 + 10.0*E*dx2 - 15.0*D)*dx*dy;                                  /* d5Ex/dx4dy    */
                        DE[54] = (-F*dx4 + 10.0*E*dx2 - 15.0*D)*dx*dz;                                  /* d5Ex/dx4dz    */
                        DE[55] = -F*dx4*dy2 + E*(dx4 + 6.0*dx2*dy2) - D*(3.0*dy2 + 6.0*dx2) + 3.0*C;    /* d5Ex/dx3dy2   */
                        DE[56] = (-F*dx4 + 6.0*E*dx2 - 3.0*D)*dy*dz;                                    /* d5Ex/dx3dydz  */
                        DE[57] = -F*dx4*dz2 + E*(dx4 + 6.0*dx2*dz2) - D*(3.0*dz2 + 6.0*dx2) + 3.0*C;    /* d5Ex/dx3dz2   */
                        DE[58] = (-F*dx2*dy2 + 3.0*E*(dx2 + dy2) - 9.0*D)*dx*dy;                        /* d5Ex/dx2dy3   */
                        DE[59] = (-F*dx2*dy2 + E*(dx2 + 3.0*dy2) - 3.0*D)*dx*dz;                        /* d5Ex/dx2dy2dz */
                        DE[60] = (-F*dx2*dz2 + E*(dx2 + 3.0*dz2) - 3.0*D)*dx*dy;                        /* d5Ex/dx2dydz2 */
                        DE[61] = (-F*dx2*dz2 + 3.0*E*(dx2 + dz2) - 9.0*D)*dx*dz;                        /* d5Ex/dx2dz3   */
                        DE[62] = -F*dx2*dy4 + E*(dy4 + 6.0*dy2*dx2) - D*(3.0*dx2 + 6.0*dy2) + 3.0*C;    /* d5Ex/dxdy4    */
                        DE[63] = (-F*dx2*dy2 + E*(dy2 + 3.0*dx2) - 3.0*D)*dy*dz;                        /* d5Ex/dxdy3dz  */
                        DE[64] = -F*dx2*dy2*dz2 + E*(dx2*dy2 + dx2*dz2 + dy2*dz2) - D*rsq + C;          /* d5Ex/dxdy2dz2 */
                        DE[65] = (-F*dx2*dz2 + E*(dz2 + 3.0*dx2) - 3.0*D)*dy*dz;                        /* d5Ex/dxdydz3  */
                        DE[66] = -F*dx2*dz4 + E*(dz4 + 6.0*dx2*dz2) - D*(3.0*dx2 + 6.0*dz2) + 3.0*C;    /* d5Ex/dxdz4    */
                        DE[67] = (-F*dy4 + 10.0*E*dy2 - 15.0*D)*dx*dy;                                  /* d5Ex/dy5      */
                        DE[68] = (-F*dy4 + 6.0*E*dy2 - 3.0*D)*dx*dz;                                    /* d5Ex/dy4dz    */
                        DE[69] = (-F*dy2*dz2 + E*(dy2 + 3.0*dz2) - 3.0*D)*dx*dy;                        /* d5Ex/dy3dz2   */
                        DE[70] = (-F*dy2*dz2 + E*(dz2 + 3.0*dy2) - 3.0*D)*dx*dz;                        /* d5Ex/dy2dz3   */
                        DE[71] = (-F*dz4 + 6.0*E*dz2 - 3.0*D)*dx*dy;                                    /* d5Ex/dydz4    */
                        DE[72] = (-F*dz4 + 10.0*E*dz2 - 15.0*D)*dx*dz;                                  /* d5Ex/dz5      */
    
                        DE[73] = -F*dy4*dy2 + 15.0*E*dy4 - 45.0*D*dy2 + 15.0*C;                         /* d5Ey/dy5      */
                        DE[74] = (-F*dy4 + 10.0*E*dy2 - 15.0*D)*dy*dz;                                  /* d5Ey/dy4dz    */
                        DE[75] = -F*dy4*dz2 + E*(dy4 + 6.0*dz2*dy2) - D*(3.0*dz2 + 6.0*dy2) + 3.0*C;    /* d5Ey/dy3dz2   */
                        DE[76] = (-F*dy2*dz2 + 3.0*E*(dy2 + dz2) - 9.0*D)*dy*dz;                        /* d5Ey/dy2dz3   */
                        DE[77] = -F*dy2*dz4 + E*(dz4 + 6.0*dy2*dz2) - D*(3.0*dy2 + 6.0*dz2) + 3.0*C;    /* d5Ey/dydz4    */
                        DE[78] = (-F*dz4 + 10.0*E*dz2 - 15.0*D)*dy*dz;                                  /* d5Ey/dz5      */
                    
                        DE[79] = -F*dz4*dz2 + 15.0*E*dz4 - 45.0*D*dz2 + 15.0*C;                         /* d5Ez/dz5      */
                    }
                }
            }
        }
    
        /* Add contributions from all nodes */
        for(j = 0; j < NODES_IN_ELEM; j++){
            Efield[j][2][0][0] += DE[0]*N[j];
            Efield[j][1][1][0] += DE[1]*N[j];
            Efield[j][1][0][1] += DE[2]*N[j];
            Efield[j][0][2][0] += DE[3]*N[j];
            Efield[j][0][1][1] += DE[4]*N[j];
            Efield[j][0][0][2] += DE[5]*N[j];
            
            if(nOrder > 1){
                Efield[j][3][0][0] += DE[6]*N[j];
                Efield[j][2][1][0] += DE[7]*N[j];
                Efield[j][2][0][1] += DE[8]*N[j];
                Efield[j][1][2][0] += DE[9]*N[j];
                Efield[j][1][1][1] += DE[10]*N[j];
                Efield[j][1][0][2] += DE[11]*N[j];
            
                Efield[j][0][3][0] += DE[12]*N[j];
                Efield[j][0][2][1] += DE[13]*N[j];
                Efield[j][0][1][2] += DE[14]*N[j];
            
                Efield[j][0][0][3] += DE[15]*N[j];  

                if(nOrder > 2){
                    Efield[j][4][0][0] += DE[16]*N[j];
                    Efield[j][3][1][0] += DE[17]*N[j];
                    Efield[j][3][0][1] += DE[18]*N[j];
                    Efield[j][2][2][0] += DE[19]*N[j];
                    Efield[j][2][1][1] += DE[20]*N[j];
                    Efield[j][2][0][2] += DE[21]*N[j];
                    Efield[j][1][3][0] += DE[22]*N[j];
                    Efield[j][1][2][1] += DE[23]*N[j];
                    Efield[j][1][1][2] += DE[24]*N[j];
                    Efield[j][1][0][3] += DE[25]*N[j];
            
                    Efield[j][0][4][0] += DE[26]*N[j];
                    Efield[j][0][3][1] += DE[27]*N[j];
                    Efield[j][0][2][2] += DE[28]*N[j];
                    Efield[j][0][1][3] += DE[29]*N[j];
            
                    Efield[j][0][0][4] += DE[30]*N[j];
                
                    if(nOrder > 3){
                        Efield[j][5][0][0] += DE[31]*N[j];
                        Efield[j][4][1][0] += DE[32]*N[j];
                        Efield[j][4][0][1] += DE[33]*N[j];
                        Efield[j][3][2][0] += DE[34]*N[j];
                        Efield[j][3][1][1] += DE[35]*N[j];
                        Efield[j][3][0][2] += DE[36]*N[j];
                        Efield[j][2][3][0] += DE[37]*N[j];
                        Efield[j][2][2][1] += DE[38]*N[j];
                        Efield[j][2][1][2] += DE[39]*N[j];
                        Efield[j][2][0][3] += DE[40]*N[j];
                        Efield[j][1][4][0] += DE[41]*N[j];
                        Efield[j][1][3][1] += DE[42]*N[j];
                        Efield[j][1][2][2] += DE[43]*N[j];
                        Efield[j][1][1][3] += DE[44]*N[j];
                        Efield[j][1][0][4] += DE[45]*N[j];
        
                        Efield[j][0][5][0] += DE[46]*N[j];
                        Efield[j][0][4][1] += DE[47]*N[j];
                        Efield[j][0][3][2] += DE[48]*N[j];
                        Efield[j][0][2][3] += DE[49]*N[j];
                        Efield[j][0][1][4] += DE[50]*N[j];
                    
                        Efield[j][0][0][5] += DE[51]*N[j];
                    
                        if(nOrder > 4){
                            Efield[j][6][0][0] += DE[52]*N[j];
                            Efield[j][5][1][0] += DE[53]*N[j];
                            Efield[j][5][0][1] += DE[54]*N[j];
                            Efield[j][4][2][0] += DE[55]*N[j];
                            Efield[j][4][1][1] += DE[56]*N[j];
                            Efield[j][4][0][2] += DE[57]*N[j];
                            Efield[j][3][3][0] += DE[58]*N[j];
                            Efield[j][3][2][1] += DE[59]*N[j];
                            Efield[j][3][1][2] += DE[60]*N[j];
                            Efield[j][3][0][3] += DE[61]*N[j];
                            Efield[j][2][4][0] += DE[62]*N[j];
                            Efield[j][2][3][1] += DE[63]*N[j];
                            Efield[j][2][2][2] += DE[64]*N[j];
                            Efield[j][2][1][3] += DE[65]*N[j];
                            Efield[j][2][0][4] += DE[66]*N[j];
                            Efield[j][1][5][0] += DE[67]*N[j];
                            Efield[j][1][4][1] += DE[68]*N[j];
                            Efield[j][1][3][2] += DE[69]*N[j];
                            Efield[j][1][2][3] += DE[70]*N[j];
                            Efield[j][1][1][4] += DE[71]*N[j];
                            Efield[j][1][0][5] += DE[72]*N[j];

                            Efield[j][0][6][0] += DE[73]*N[j];
                            Efield[j][0][5][1] += DE[74]*N[j];
                            Efield[j][0][4][2] += DE[75]*N[j];
                            Efield[j][0][3][3] += DE[76]*N[j];
                            Efield[j][0][2][4] += DE[77]*N[j];
                            Efield[j][0][1][5] += DE[78]*N[j];
                        
                            Efield[j][0][0][6] += DE[79]*N[j];
                        }
                    }
                }
            }
        }          
    }

    return 0;
}


/**
 * Computes the integral of the 2nd, 3rd and 4th order derivatives of the 3D
 * Green's function for one element and returns the integrals in matrix
 * \a E for the 3 dimensions and each node. These are some examples:
 *
 * \a Efield[0][1][0][0] : Ex @ node 1 \n
 * \a Efield[0][0][1][0] : Ey @ node 1 \n
 * \a Efield[0][2][0][0] : dEx/dx @ node 1 \n
 * \a Efield[0][1][0][1] : dEx/dz = dEz/dx @ node 1 \n
 * \a Efield[2][0][1][3] : d3Ez/dydz2 = d3Ey/dz3 @ node 3 \n
 *
 * Where Ex = dG/dx, Ey = dG/dy and Ez = dG/dz and G is Green's function.
 *
 * Works for quadratic interpolation in triangular elements (6-noded triangles).
 *
 * @param nOrder : [ Input ]  Multipolar approximation order (implemented 1, 2, 3 only)
 * @param X      : [ Input ]  Coordinates of the nodes defining the element
 * @param Xeval  : [ Input ]  Coordinates of the evaluation point
 * @param Efield : [ Output ] Integrals of the differentials of the Green's function
 */
int intDE_tria6(unsigned int nOrder,
                double X[6][3],
                double *Xeval,
                double Efield[][7][7][7])
{
    const unsigned int NODES_IN_ELEM = 6;
    unsigned int i, j, k, m;
    double A, B, C, D, E, F, dx, dx2, dx4, dy, dy2, dy4, dz, dz2, dz4, J, rsq, W;
    double DE[80], L[2], N[6], XL[3];

    /* Initialize */
    for(i = 0; i < 80; i++) DE[i] = 0.0;
    for(i = 0; i < NODES_IN_ELEM; i++)
        for(j = 0; j < 7; j++)
            for(k = 0; k < 7; k++)
                for(m = 0; m < 7; m++) Efield[i][j][k][m] = 0.0;

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
        dx2 = dx*dx;
        dy2 = dy*dy;
        dz2 = dz*dz;
        rsq = dx2 + dy2 + dz2;
        A = J*W/pow(rsq,1.5);
        B = 3.0*A/rsq;

        /***********************************************/
        /*** GET CONTRIBUTION FROM EACH ELEMENT NODE ***/
        /***********************************************/
    
        /* First order derivatives of the electric field */
        DE[0] = A - B*dx2;  /* dEx/dx */
        DE[1] = -B*dx*dy;   /* dEx/dy */
        DE[2] = -B*dx*dz;   /* dEx/dz */
    
        DE[3] = A - B*dy2;  /* dEy/dy */
        DE[4] = -B*dy*dz;   /* dEy/dz */

        DE[5] = A - B*dz2;  /* dEz/dz */
    
        /* Second order derivatives of the electric field */
        if(nOrder > 1){
            C = 5.0*B/rsq;

            DE[6] = (C*dx2 - 3.0*B)*dx;     /* d2Ex/dx2  */
            DE[7] = (C*dx2 - B)*dy;         /* d2Ex/dxdy */
            DE[8] = (C*dx2 - B)*dz;         /* d2Ex/dxdz */
            DE[9] = (C*dy2 - B)*dx;         /* d2Ex/dy2  */
            DE[10] = C*dx*dy*dz;            /* d2Ex/dydz */
            DE[11] = (C*dz2 - B)*dx;        /* d2Ex/dzdz */
        
            DE[12] = (C*dy2 - 3.0*B)*dy;    /* d2Ey/dy2  */
            DE[13] = (C*dy2 - B)*dz;        /* d2Ey/dydz */
            DE[14] = (C*dz2 - B)*dy;        /* d2Ey/dz2  */
        
            DE[15] = (C*dz2 - 3.0*B)*dz;    /* d2Ez/dz2  */
        
            /* Third order derivatives of the electric field */
            if(nOrder > 2){
                D = 7.0*C/rsq; 

                DE[16] = (-D*dx2 + 6.0*C)*dx2 - 3.0*B;      /* d3Ex/dx3    */
                DE[17] = (-D*dx2 + 3.0*C)*dx*dy;            /* d3Ex/dx2dy  */
                DE[18] = (-D*dx2 + 3.0*C)*dx*dz;            /* d3Ex/dx2dz  */
                DE[19] = -D*dx2*dy2 + C*(dx2 + dy2) - B;    /* d3Ex/dxdy2  */
                DE[20] = (-D*dx2 + C)*dy*dz;                /* d3Ex/dxdydz */
                DE[21] = -D*dx2*dz2 + C*(dx2 + dz2) - B;    /* d3Ex/dxdz2  */
                DE[22] = (-D*dy2 + 3.0*C)*dx*dy;            /* d3Ex/dy3    */
                DE[23] = (-D*dy2 + C)*dx*dz;                /* d3Ex/dy2dz  */
                DE[24] = (-D*dz2 + C)*dx*dy;                /* d3Ex/dydz2  */
                DE[25] = (-D*dz2 + 3.0*C)*dx*dz;            /* d3Ex/dz3    */
        
                DE[26] = (-D*dy2 + 6.0*C)*dy2 - 3.0*B;      /* d3Ey/dy3    */
                DE[27] = (-D*dy2 + 3.0*C)*dy*dz;            /* d3Ey/dy2dz  */
                DE[28] = -D*dy2*dz2 + C*(dy2 + dz2) - B;    /* d3Ey/dydz2  */
                DE[29] = (-D*dz2 + 3.0*C)*dy*dz;            /* d3Ey/dz3    */
        
                DE[30] = (-D*dz2 + 6.0*C)*dz2 - 3.0*B;      /* d3Ez/dz3    */
            
                /* Fourth order derivatives of the electric field */
                if(nOrder > 3){
                    E = 9.0*D/rsq;
                    dx4 = dx2*dx2;
                    dy4 = dy2*dy2;
                    dz4 = dz2*dz2;
                
                    DE[31] = (E*dx4 - 10.0*D*dx2 + 15.0*C)*dx;              /* d4Ex/dx4     */
                    DE[32] = (E*dx4 - 6.0*D*dx2 + 3.0*C)*dy;                /* d4Ex/dx3dy   */
                    DE[33] = (E*dx4 - 6.0*D*dx2 + 3.0*C)*dz;                /* d4Ex/dx3dz   */
                    DE[34] = (E*dx2*dy2 - D*(dx2 + 3.0*dy2) + 3.0*C)*dx;    /* d4Ex/dx2dy2  */
                    DE[35] = (E*dx2 - 3.0*D)*dx*dy*dz;                      /* d4Ex/dx2dydz */
                    DE[36] = (E*dx2*dz2 - D*(dx2 + 3.0*dz2) + 3.0*C)*dx;    /* d4Ex/dx2dz2  */
                    DE[37] = (E*dx2*dy2 - D*(dy2 + 3.0*dx2) + 3.0*C)*dy;    /* d4Ex/dxdy3   */
                    DE[38] = (E*dx2*dy2 - D*(dx2 + dy2) + C)*dz;            /* d4Ex/dxdy2dz */
                    DE[39] = (E*dx2*dz2 - D*(dx2 + dz2) + C)*dy;            /* d4Ex/dxdydz2 */
                    DE[40] = (E*dx2*dz2 - D*(dz2 + 3.0*dx2) + 3.0*C)*dz;    /* d4Ex/dxdz3   */
                    DE[41] = (E*dy4 - 6.0*D*dy2 + 3.0*C)*dx;                /* d4Ex/dy4     */
                    DE[42] = (E*dy2 - 3.0*D)*dx*dy*dz;                      /* d4Ex/dy3dz   */
                    DE[43] = (E*dy2*dz2 - D*(dy2 + dz2) + C)*dx;            /* d4Ex/dy2dz2  */
                    DE[44] = (E*dz2 - 3.0*D)*dx*dy*dz;                      /* d4Ex/dydz3   */
                    DE[45] = (E*dz4 - 6.0*D*dz2 + 3.0*C)*dx;                /* d4E/dz4      */
                
                    DE[46] = (E*dy4 - 10.0*D*dy2 + 15.0*C)*dy;              /* d4Ey/dy4     */
                    DE[47] = (E*dy4 - 6.0*D*dy2 + 3.0*C)*dz;                /* d4Ey/dy3dz   */
                    DE[48] = (E*dy2*dz2 - D*(dy2 + 3.0*dz2) + 3.0*C)*dy;    /* d4Ey/dy2dz2  */
                    DE[49] = (E*dy2*dz2 - D*(dz2 + 3.0*dy2) + 3.0*C)*dz;    /* d4Ey/dydz3   */
                    DE[50] = (E*dz4 - 6.0*D*dz2 + 3.0*C)*dy;                /* d4Ey/dz4     */
                
                    DE[51] = (E*dz4 - 10.0*D*dz2 + 15.0*C)*dz;              /* d4Ez/dz4     */
                
                    /* Fifth order derivatives of the electric field */
                    if(nOrder > 4){
                        F = 11.0*E/rsq;
                    
                        DE[52] = -F*dx4*dx2 + 15.0*E*dx4 - 45.0*D*dx2 + 15.0*C;                         /* d5Ex/dx5      */
                        DE[53] = (-F*dx4 + 10.0*E*dx2 - 15.0*D)*dx*dy;                                  /* d5Ex/dx4dy    */
                        DE[54] = (-F*dx4 + 10.0*E*dx2 - 15.0*D)*dx*dz;                                  /* d5Ex/dx4dz    */
                        DE[55] = -F*dx4*dy2 + E*(dx4 + 6.0*dx2*dy2) - D*(3.0*dy2 + 6.0*dx2) + 3.0*C;    /* d5Ex/dx3dy2   */
                        DE[56] = (-F*dx4 + 6.0*E*dx2 - 3.0*D)*dy*dz;                                    /* d5Ex/dx3dydz  */
                        DE[57] = -F*dx4*dz2 + E*(dx4 + 6.0*dx2*dz2) - D*(3.0*dz2 + 6.0*dx2) + 3.0*C;    /* d5Ex/dx3dz2   */
                        DE[58] = (-F*dx2*dy2 + 3.0*E*(dx2 + dy2) - 9.0*D)*dx*dy;                        /* d5Ex/dx2dy3   */
                        DE[59] = (-F*dx2*dy2 + E*(dx2 + 3.0*dy2) - 3.0*D)*dx*dz;                        /* d5Ex/dx2dy2dz */
                        DE[60] = (-F*dx2*dz2 + E*(dx2 + 3.0*dz2) - 3.0*D)*dx*dy;                        /* d5Ex/dx2dydz2 */
                        DE[61] = (-F*dx2*dz2 + 3.0*E*(dx2 + dz2) - 9.0*D)*dx*dz;                        /* d5Ex/dx2dz3   */
                        DE[62] = -F*dx2*dy4 + E*(dy4 + 6.0*dy2*dx2) - D*(3.0*dx2 + 6.0*dy2) + 3.0*C;    /* d5Ex/dxdy4    */
                        DE[63] = (-F*dx2*dy2 + E*(dy2 + 3.0*dx2) - 3.0*D)*dy*dz;                        /* d5Ex/dxdy3dz  */
                        DE[64] = -F*dx2*dy2*dz2 + E*(dx2*dy2 + dx2*dz2 + dy2*dz2) - D*rsq + C;          /* d5Ex/dxdy2dz2 */
                        DE[65] = (-F*dx2*dz2 + E*(dz2 + 3.0*dx2) - 3.0*D)*dy*dz;                        /* d5Ex/dxdydz3  */
                        DE[66] = -F*dx2*dz4 + E*(dz4 + 6.0*dx2*dz2) - D*(3.0*dx2 + 6.0*dz2) + 3.0*C;    /* d5Ex/dxdz4    */
                        DE[67] = (-F*dy4 + 10.0*E*dy2 - 15.0*D)*dx*dy;                                  /* d5Ex/dy5      */
                        DE[68] = (-F*dy4 + 6.0*E*dy2 - 3.0*D)*dx*dz;                                    /* d5Ex/dy4dz    */
                        DE[69] = (-F*dy2*dz2 + E*(dy2 + 3.0*dz2) - 3.0*D)*dx*dy;                        /* d5Ex/dy3dz2   */
                        DE[70] = (-F*dy2*dz2 + E*(dz2 + 3.0*dy2) - 3.0*D)*dx*dz;                        /* d5Ex/dy2dz3   */
                        DE[71] = (-F*dz4 + 6.0*E*dz2 - 3.0*D)*dx*dy;                                    /* d5Ex/dydz4    */
                        DE[72] = (-F*dz4 + 10.0*E*dz2 - 15.0*D)*dx*dz;                                  /* d5Ex/dz5      */
                    
                        DE[73] = -F*dy4*dy2 + 15.0*E*dy4 - 45.0*D*dy2 + 15.0*C;                         /* d5Ey/dy5      */
                        DE[74] = (-F*dy4 + 10.0*E*dy2 - 15.0*D)*dy*dz;                                  /* d5Ey/dy4dz    */
                        DE[75] = -F*dy4*dz2 + E*(dy4 + 6.0*dz2*dy2) - D*(3.0*dz2 + 6.0*dy2) + 3.0*C;    /* d5Ey/dy3dz2   */
                        DE[76] = (-F*dy2*dz2 + 3.0*E*(dy2 + dz2) - 9.0*D)*dy*dz;                        /* d5Ey/dy2dz3   */
                        DE[77] = -F*dy2*dz4 + E*(dz4 + 6.0*dy2*dz2) - D*(3.0*dy2 + 6.0*dz2) + 3.0*C;    /* d5Ey/dydz4    */
                        DE[78] = (-F*dz4 + 10.0*E*dz2 - 15.0*D)*dy*dz;                                  /* d5Ey/dz5      */
                    
                        DE[79] = -F*dz4*dz2 + 15.0*E*dz4 - 45.0*D*dz2 + 15.0*C;                         /* d5Ez/dz5      */
                    }
                }
            }
        }
    
        /* Add contributions from all nodes */
        for(j = 0; j < NODES_IN_ELEM; j++){
            Efield[j][2][0][0] += DE[0]*N[j];
            Efield[j][1][1][0] += DE[1]*N[j];
            Efield[j][1][0][1] += DE[2]*N[j];
            Efield[j][0][2][0] += DE[3]*N[j];
            Efield[j][0][1][1] += DE[4]*N[j];
            Efield[j][0][0][2] += DE[5]*N[j];
            
            if(nOrder > 1){
                Efield[j][3][0][0] += DE[6]*N[j];
                Efield[j][2][1][0] += DE[7]*N[j];
                Efield[j][2][0][1] += DE[8]*N[j];
                Efield[j][1][2][0] += DE[9]*N[j];
                Efield[j][1][1][1] += DE[10]*N[j];
                Efield[j][1][0][2] += DE[11]*N[j];
            
                Efield[j][0][3][0] += DE[12]*N[j];
                Efield[j][0][2][1] += DE[13]*N[j];
                Efield[j][0][1][2] += DE[14]*N[j];
            
                Efield[j][0][0][3] += DE[15]*N[j];  

                if(nOrder > 2){
                    Efield[j][4][0][0] += DE[16]*N[j];
                    Efield[j][3][1][0] += DE[17]*N[j];
                    Efield[j][3][0][1] += DE[18]*N[j];
                    Efield[j][2][2][0] += DE[19]*N[j];
                    Efield[j][2][1][1] += DE[20]*N[j];
                    Efield[j][2][0][2] += DE[21]*N[j];
                    Efield[j][1][3][0] += DE[22]*N[j];
                    Efield[j][1][2][1] += DE[23]*N[j];
                    Efield[j][1][1][2] += DE[24]*N[j];
                    Efield[j][1][0][3] += DE[25]*N[j];
            
                    Efield[j][0][4][0] += DE[26]*N[j];
                    Efield[j][0][3][1] += DE[27]*N[j];
                    Efield[j][0][2][2] += DE[28]*N[j];
                    Efield[j][0][1][3] += DE[29]*N[j];
            
                    Efield[j][0][0][4] += DE[30]*N[j];
                
                    if(nOrder > 3){
                        Efield[j][5][0][0] += DE[31]*N[j];
                        Efield[j][4][1][0] += DE[32]*N[j];
                        Efield[j][4][0][1] += DE[33]*N[j];
                        Efield[j][3][2][0] += DE[34]*N[j];
                        Efield[j][3][1][1] += DE[35]*N[j];
                        Efield[j][3][0][2] += DE[36]*N[j];
                        Efield[j][2][3][0] += DE[37]*N[j];
                        Efield[j][2][2][1] += DE[38]*N[j];
                        Efield[j][2][1][2] += DE[39]*N[j];
                        Efield[j][2][0][3] += DE[40]*N[j];
                        Efield[j][1][4][0] += DE[41]*N[j];
                        Efield[j][1][3][1] += DE[42]*N[j];
                        Efield[j][1][2][2] += DE[43]*N[j];
                        Efield[j][1][1][3] += DE[44]*N[j];
                        Efield[j][1][0][4] += DE[45]*N[j];
                    
                        Efield[j][0][5][0] += DE[46]*N[j];
                        Efield[j][0][4][1] += DE[47]*N[j];
                        Efield[j][0][3][2] += DE[48]*N[j];
                        Efield[j][0][2][3] += DE[49]*N[j];
                        Efield[j][0][1][4] += DE[50]*N[j];
                    
                        Efield[j][0][0][5] += DE[51]*N[j];
                        
                        if(nOrder > 4){
                            Efield[j][6][0][0] += DE[52]*N[j];
                            Efield[j][5][1][0] += DE[53]*N[j];
                            Efield[j][5][0][1] += DE[54]*N[j];
                            Efield[j][4][2][0] += DE[55]*N[j];
                            Efield[j][4][1][1] += DE[56]*N[j];
                            Efield[j][4][0][2] += DE[57]*N[j];
                            Efield[j][3][3][0] += DE[58]*N[j];
                            Efield[j][3][2][1] += DE[59]*N[j];
                            Efield[j][3][1][2] += DE[60]*N[j];
                            Efield[j][3][0][3] += DE[61]*N[j];
                            Efield[j][2][4][0] += DE[62]*N[j];
                            Efield[j][2][3][1] += DE[63]*N[j];
                            Efield[j][2][2][2] += DE[64]*N[j];
                            Efield[j][2][1][3] += DE[65]*N[j];
                            Efield[j][2][0][4] += DE[66]*N[j];
                            Efield[j][1][5][0] += DE[67]*N[j];
                            Efield[j][1][4][1] += DE[68]*N[j];
                            Efield[j][1][3][2] += DE[69]*N[j];
                            Efield[j][1][2][3] += DE[70]*N[j];
                            Efield[j][1][1][4] += DE[71]*N[j];
                            Efield[j][1][0][5] += DE[72]*N[j];
                        
                            Efield[j][0][6][0] += DE[73]*N[j];
                            Efield[j][0][5][1] += DE[74]*N[j];
                            Efield[j][0][4][2] += DE[75]*N[j];
                            Efield[j][0][3][3] += DE[76]*N[j];
                            Efield[j][0][2][4] += DE[77]*N[j];
                            Efield[j][0][1][5] += DE[78]*N[j];
                        
                            Efield[j][0][0][6] += DE[79]*N[j];
                        }
                    }
                }
            }
        }          
    }

    return 0;
}


/**
 * Computes the integral of the Maxwell Stress Tensor for one element and
 * returns the integrals in vector \a Int for the three spatial dimentions.
 *
 * Works for linear interpolation in triangular elements (3-noded triangles).
 *
 * @param X   : [ Input ]  Coordinates of the nodes defining the element
 * @param E   : [ Input ]  Electric field at the nodes defining the element
 * @param Int : [ Output ] Integral of the Maxwell Stress tensor for the element
 */
int intF_tria3(double X[3][3],
               double E[3][3],
               double *Int)
{
    unsigned int i;
    double A, Ex, Ey, Ez, J, W, Esq, NE;
    double L[2], N[3], normal[3], XL[3];

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
        getLocalNormal_tria3(X,normal);

        /* Electric field value at this integration point */
        Ex = E[0][0]*N[0] + E[1][0]*N[1] + E[2][0]*N[2];
        Ey = E[0][1]*N[0] + E[1][1]*N[1] + E[2][1]*N[2];
        Ez = E[0][2]*N[0] + E[1][2]*N[1] + E[2][2]*N[2];
    
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


/**
 * Computes the integral of the Maxwell Stress Tensor for one element and
 * returns the integrals in vector \a Int for the three spatial dimentions.
 *
 * Works for quadratic interpolation in triangular elements (6-noded triangles)
 *
 * @param X   : [ Input ]  Coordinates of the nodes defining the element
 * @param E   : [ Input ]  Electric field at the nodes defining the element
 * @param Int : [ Output ] Integral of the Maxwell Stress tensor for the element
 */
int intF_tria6(double X[6][3],
               double E[6][3],
               double *Int)
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


