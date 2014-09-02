/******************************************************************************
* File      : potential_tria3.c                                               *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* Calculates the electric potential at the required point Xin and returns its *
* real and imaginary parts in array mPot.                                     *
* Works for quadratic interpolation in triangular elements (6-noded triangles)*
*                                                                             *
* mPot[0] = Re(V)    |    mPot[1] = Im(V)                                     *
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
#include "potential_tria6.h"

int potential_tria6(double *Xin, double **mNodes, unsigned int **mElems, 
                    double *vB, double *mPot)
{
    const unsigned int ELEMS = nElems, NODES_IN_ELEM = 6;
    unsigned int currentNode, i, j, test, SinNode; 
    double A, dx, dy, dz;
    double normal[3], SubP[6];
    double X[6][3];

    /* Initialize */
    A = 1.0/(pi4*eps0);
    mPot[0] = mPot[1] = 0.0;
 
    for(i = 0; i < ELEMS; i++){
        
        /* Coordinates of nodes of this element */
        for(j = 0; j < NODES_IN_ELEM; j++){
            currentNode = mElems[i][j] - 1;
            X[j][0] = mNodes[currentNode][0];
            X[j][1] = mNodes[currentNode][1];
            X[j][2] = mNodes[currentNode][2];
        }

        /* Test for singular case */
        test = 0;
        for(j = 0; j < NODES_IN_ELEM; j++){
            dx = X[j][0] - Xin[0];
            dy = X[j][1] - Xin[1];
            dz = X[j][2] - Xin[2];
            if(dx == 0 && dy == 0.0 && dz == 0.0 && test == 0){
                test = 1;
                SinNode = j+1;
                getNormal_tria6(mElems[i][j],mNodes,mElems,normal);
            }
        }

        /* Contribution from all j nodes in element i */
        if(test == 1) intSingularG_tria6(SinNode,X,Xin,SubP);
        else intG_tria6(X,Xin,SubP);

        /* Add contribution to potential from element i */
        for(j = 0; j < NODES_IN_ELEM; j++){
            currentNode = mElems[i][j] - 1;
            mPot[0] += A*SubP[j]*vB[currentNode];
            mPot[1] += A*SubP[j]*vB[currentNode + nNodes];
        }
    }

    return 0;
}
