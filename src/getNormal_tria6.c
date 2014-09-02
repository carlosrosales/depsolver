/******************************************************************************
* File      : getNormal_tria6.c                                               *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* Returns the unit normal[] vector at the point nNodeID. The average normal   *
* is used if the node is shared by several elements.                          *
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Global variables */
extern unsigned int nElems;

int getNormal_tria6(unsigned int nNodeID, double **mNodes, 
                    unsigned int **mElems, double *normal)
{
    const unsigned int NODES_IN_ELEM = 6;
    unsigned int i, j, nIsInElem;
    int nCurrentElem, nNodePosition;
    double ax, ay, az, dx1, dx2, dy1, dy2, dz1, dz2, L1, L2, normalmod;
    double tmp[3];
    double X[6][3];

    /* Initialize */
    tmp[0] = tmp[1] = tmp[2] = 0.0;
    nIsInElem = 0;
    nCurrentElem = -1;
    nNodePosition = -1;

    /* Loop through all elements and average their normals at nNodeID */
    do{
        
        /* Find element containing the evaluation point nNodeID */
        for(i = nCurrentElem+1; i < nElems; i++){ 
            for(j = 0; j < NODES_IN_ELEM; j++){
                if(mElems[i][j] == nNodeID) {
                    nCurrentElem = i;
                    nNodePosition = j;
                    nIsInElem = 1;
                    break;
                    }
            }
            if(nIsInElem) break;
        }

        if(nIsInElem){
            
            /* (L1,L2) for evaluation point in this element */
            if(nNodePosition == 0){
                L1 = 1.0;
                L2 = 0.0;
            }
            else if(nNodePosition == 1){
                L1 = 0.5;
                L2 = 0.5;
            }
            else if(nNodePosition == 2){
                L1 = 0.0;
                L2 = 1.0;
            }
            else if(nNodePosition == 3){
                L1 = 0.0;
                L2 = 0.5;
            }
            else if(nNodePosition == 4){
                L1 = 0.0;
                L2 = 0.0;
            }
            else{
                L1 = 0.5;
                L2 = 0.0;
            }

            /* Get coordinates of nodes of this element */
            for(j = 0; j < NODES_IN_ELEM; j++){
                X[j][0] = mNodes[mElems[nCurrentElem][j]-1][0];
                X[j][1] = mNodes[mElems[nCurrentElem][j]-1][1];
                X[j][2] = mNodes[mElems[nCurrentElem][j]-1][2];
            }

            /* Auxiliar quantities */
            ax = X[1][0] - X[3][0] + X[4][0] - X[5][0];
            ay = X[1][1] - X[3][1] + X[4][1] - X[5][1];
            az = X[1][2] - X[3][2] + X[4][2] - X[5][2];

            /* Derivatives of (x,y,z) with respect to L1 and L2 */
            dx1 = 4.0*(L1*(X[0][0] + X[4][0] - 2.0*X[5][0]) + L2*ax + X[5][0]) - X[0][0] - 3.0*X[4][0];
            dx2 = 4.0*(L1*ax + L2*(X[2][0] + X[4][0] - 2.0*X[3][0]) + X[3][0]) - X[2][0] - 3.0*X[4][0];
            dy1 = 4.0*(L1*(X[0][1] + X[4][1] - 2.0*X[5][1]) + L2*ay + X[5][1]) - X[0][1] - 3.0*X[4][1];
            dy2 = 4.0*(L1*ay + L2*(X[2][1] + X[4][1] - 2.0*X[3][1]) + X[3][1]) - X[2][1] - 3.0*X[4][1];
            dz1 = 4.0*(L1*(X[0][2] + X[4][2] - 2.0*X[5][2]) + L2*az + X[5][2]) - X[0][2] - 3.0*X[4][2];
            dz2 = 4.0*(L1*az + L2*(X[2][2] + X[4][2] - 2.0*X[3][2]) + X[3][2]) - X[2][2] - 3.0*X[4][2];

            /* Normal components and modulus (inverted) */
            normal[0] = dy1*dz2 - dy2*dz1;
            normal[1] = dx2*dz1 - dx1*dz2;
            normal[2] = dx1*dy2 - dx2*dy1;
            normalmod = 1.0/sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);

            /* Gunit normal */
            normal[0] = normal[0]*normalmod;
            normal[1] = normal[1]*normalmod;
            normal[2] = normal[2]*normalmod;

            /* Add contribution to average normal */
            tmp[0] += normal[0];
            tmp[1] += normal[1];
            tmp[2] += normal[2];

            /* Reset flag */
            nIsInElem = 0;
        }
    }while(i < nElems);
    
    /* Modulus of the average normal (inverted) */
    normalmod = 1.0/sqrt(tmp[0]*tmp[0] + tmp[1]*tmp[1] + tmp[2]*tmp[2]);

    /* Get unit average normal */
    normal[0] = tmp[0]*normalmod;
    normal[1] = tmp[1]*normalmod;
    normal[2] = tmp[2]*normalmod;

    return 0;
}

