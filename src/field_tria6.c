/******************************************************************************
* File      : field_tria6.c                                                   *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* Calculates the Electric Field at the required point Xin[] and returns its   *
* real and imaginary parts in mE[].                                           *
* Works for quadratic interpolation in triangular elements (6-noded triangles)*
* mE[0] -> Re(Ex)    |    mE[2] -> Re(Ey)    |    mE[4] -> Re(Ez)             *
* mE[1] -> Im(Ex)    |    mE[3] -> Im(Ey)    |    mE[5] -> Im(Ez)             *
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

#include "constants.h"
#include "field_tria6.h"

int field_tria6(double *Xin, double **mNodes, unsigned int **mElems, double *vB,
                double *mE)
{
    const unsigned int ELEMS = nElems, DIM = 3, FLAG = 1, NODES_IN_ELEM = 6;
    unsigned int i, j, k, currentNode, SinNode, test;
    double A, B, dx, dy, dz;
    double normal[3], SubE[18];
    double X[6][3];

    /* Initialize */
    B = 0.5/eps0;
    A = B/pi2;
    for(i = 0; i < 6; i++) mE[i] = 0.0;

    for(i = 0; i < ELEMS; i++){
        
        /* Get coordinates of nodes in this element */
        for(j = 0; j < NODES_IN_ELEM; j++){
            currentNode = mElems[i][j] - 1;
            X[j][0] = mNodes[currentNode][0];
            X[j][1] = mNodes[currentNode][1];
            X[j][2] = mNodes[currentNode][2];
        }

        /* Test for singular case */
        test = 0;
        for(j = 0; j < NODES_IN_ELEM; j++){
            dx = Xin[0] - X[j][0];
            dy = Xin[1] - X[j][1];
            dz = Xin[2] - X[j][2];
            if(dx == 0 && dy == 0.0 && dz == 0.0 && test == 0){
                test = 1;
                SinNode = j+1;
                getNormal_tria6(mElems[i][j],mNodes,mElems,normal);
            }
        }

        /* Integrate to get contributions from all j nodes in element i */
        if(test == 1) intSingularH_tria6(FLAG,SinNode,X,Xin,SubE,normal);
        else intH_tria6(FLAG,X,Xin,SubE,normal);

        /* Add contribution to field from element i */
        for(j = 0; j < NODES_IN_ELEM; j++){
            currentNode = mElems[i][j] - 1;
            for(k = 0; k < DIM; k++){
                mE[k*2] += A*SubE[j*DIM+k]*vB[currentNode];             /* Re(E) */
                mE[k*2+1] += A*SubE[j*DIM+k]*vB[currentNode + nNodes];  /* Im(E) */
            }
            if(test == 1){
                for(k = 0; k < DIM; k++){
                    mE[k*2] += B*vB[currentNode]*normal[k];
                    mE[k*2+1] += B*vB[currentNode + nNodes]*normal[k];
                }
                test = 2;
            }
        }
    }

    return 0;
}
