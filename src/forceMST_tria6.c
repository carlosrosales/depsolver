/******************************************************************************
* File      : forceMST_tria6.c                                                *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* Calculates the DEP force at dielectric interfaces using the Maxwell Stress  *
* Tensor method. Returns the total force and the application point in arrays  *
* F[] and Xcm[].                                                              *
* Works for quadratic interpolation in triangular elements (6-noded triangles)*
* F[0] = Fx    |    F[1] = Fy    |    F[2] = Fz                               *
* Xcm[0] = x   |    Xcm[1] = y   |    Xcm[2] = z                              *
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

#include "forceMST_tria6.h"

int forceMST_tria6(double **mNodes, unsigned int **mElems, 
                   unsigned int *vBCType, double *vB, double Eps, double **Fce,
                   double **Xce, double *Fcm, double *Xcm)
{
    const unsigned int ELEMS = nElems, NODES = nNodes, NODES_IN_ELEM = 6;
    unsigned int currentNode, i, j, N;
    double A, J;
    double L[2], M[6], mF[3], mE[6], Xin[3];
    double X[6][3], E[6][3];
    double **Efield;

    /* Memory allocation for real part of electric field */
    Efield = doubleMatrix(NODES,3,1);

    /* Initialize */
    A = 0.5*Eps;
    L[0] = 0.25;
    L[1] = 0.25;
    N = 0;
    Xcm[0] = Xcm[1] = Xcm[2] = 0.0;
    Fcm[0] = Fcm[1] = Fcm[2] = 0.0;

    /* Electric field values in the surface nodes */
    for(i = 0; i < NODES; i++){
        if(vBCType[i] == 6){
            Xin[0] = mNodes[i][0];
            Xin[1] = mNodes[i][1];
            Xin[2] = mNodes[i][2];
            
            Xcm[0] += Xin[0];
            Xcm[1] += Xin[1];
            Xcm[2] += Xin[2];
            N++;

            field_tria6(Xin,mNodes,mElems,vB,mE);

            /* Save Re(E) for later use */
            Efield[i][0] = mE[0];
            Efield[i][1] = mE[2];
            Efield[i][2] = mE[4];
        }
    }
    Xcm[0] = Xcm[0]/N;
    Xcm[1] = Xcm[1]/N;
    Xcm[2] = Xcm[2]/N;

    /* Integrate the Maxwell's Stress Tensor over the surface of the particle */
    for(i = 0; i < ELEMS; i++){
        if(vBCType[mElems[i][0]-1] == 6){
            for(j = 0; j < NODES_IN_ELEM; j++){
                currentNode = mElems[i][j] - 1;
                X[j][0] = mNodes[currentNode][0];
                X[j][1] = mNodes[currentNode][1];
                X[j][2] = mNodes[currentNode][2];
                E[j][0] = Efield[currentNode][0];
                E[j][1] = Efield[currentNode][1];
                E[j][2] = Efield[currentNode][2];
            }
        
            /* Contribution from this element */
            intF_tria6(X,E,mF);
        
            /* Element centre */
            shape_tria6(L,M);
            J = X2L_tria6(X,Xin,L,M);

            Fce[i][0] = A*mF[0];
            Fce[i][1] = A*mF[1];
            Fce[i][2] = A*mF[2];
            Xce[i][0] = Xin[0];
            Xce[i][1] = Xin[1];
            Xce[i][2] = Xin[2];
        
            /* Add contribution to the total force at the centre of mass */
            Fcm[0] += Fce[i][0];
            Fcm[1] += Fce[i][1];
            Fcm[2] += Fce[i][2];
        }
    }

    freeDoubleMatrix(Efield,NODES);

    return 0;
}
