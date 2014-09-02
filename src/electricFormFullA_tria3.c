/******************************************************************************
* File      : electricFormFullA_tria3.c                                       *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* Forms the coefficient matrix mA by assembling submatrices of the element,   *
* SubA, for the IBEM AC electrostatic case.                                   *
* Works for linear interpolation in triangular elements (3-noded triangles).  **
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

#include "electricFormA_tria3.h"
#include "constants.h"

int electricFormFullA_tria3(double **mNodes, unsigned int **mElems, 
                            unsigned int *vBCType, unsigned int **vInterfaces, 
                            double **vMatParam, double *vProbParam, 
                            double **mA, double **mBC, double *vB)
{
    const unsigned int FLAG = 0, NODES_IN_ELEM = 3;
    unsigned int i, interfaceID, j, k, mat1, mat2, currentNode, SinNode;  
    double C, D, EpsDif, EpsTot, F, G, SigDif, SigTot, W;
    double Xeval[3], SubA[3], normal[3];
    double X[3][3];
    double *A, *B;

    /* Auxiliar quantities for eqs at conductors */
    C = 1.0/(pi4*eps0);

    /* Auxiliar quantities for eqs at dielectric interfaces */
    if(nInterfaces){
        W = vProbParam[0]*pi2*eps0;
        A = doubleVector(nInterfaces,0);
        B = doubleVector(nInterfaces,0);
 
        for(i = 0; i < nInterfaces; i++){
            mat1 = vInterfaces[i][0] - 1;
            mat2 = vInterfaces[i][1] - 1;
            SigTot = vMatParam[mat1][0] + vMatParam[mat2][0];
            SigDif = vMatParam[mat1][0] - vMatParam[mat2][0];
            EpsTot = W*(vMatParam[mat1][1] + vMatParam[mat2][1]);
            EpsDif = W*(vMatParam[mat1][1] - vMatParam[mat2][1]);
            D = C/(SigTot*SigTot + EpsTot*EpsTot);
            A[i] = D*(SigTot*SigDif + EpsTot*EpsDif);
            B[i] = D*(SigDif*EpsTot - SigTot*EpsDif);
        }
    }

    /* Find contributions to each node and assemble coefficient matrix */
    for(i = 0; i < nNodes; i++){
        Xeval[0] = mNodes[i][0];
        Xeval[1] = mNodes[i][1];
        Xeval[2] = mNodes[i][2];

        /* Operations necessary only for dielectric interfaces */
        if(vBCType[i] == 0 || vBCType[i] == 6){
            
            /* Get interface and normal at evaluation point */
            interfaceID = (unsigned int)mBC[i][1]-1;
            getNormal_tria3(i+1,mNodes,mElems,normal);

            /* Add diagonal terms due to jump across interface */
            mA[i][i] -= 0.5/eps0;
            mA[i + nNodes][i + nNodes] -= 0.5/eps0;
        }

        /* Loop through elements doing integrals for each pair of nodes */
        for(j = 0; j < nElems; j++){
               
            /* Test for singular case */
            SinNode = 0;
            for(k = 0; k < NODES_IN_ELEM; k++){
                currentNode = mElems[j][k] - 1;
                X[k][0] = mNodes[currentNode][0];
                X[k][1] = mNodes[currentNode][1];
                X[k][2] = mNodes[currentNode][2];
                if(currentNode == i) SinNode = k + 1;
            }

            /* Potential given */
            if(vBCType[i] == 1){
                if(!SinNode) intG_tria3(X,Xeval,SubA);
                else intSingularG_tria3(SinNode,X,Xeval,SubA);
            
                /* Assign contributions to the right places in the matrix */
                for(k = 0; k < NODES_IN_ELEM; k++){
                    F = C*SubA[k];
                    currentNode = mElems[j][k] - 1;
                    mA[i][currentNode] += F;
                    mA[i + nNodes][currentNode + nNodes] += F;
                }
            }
        
            /* Dielectric interface */
            else{
                if(!SinNode) intH_tria3(FLAG,X,Xeval,SubA,normal);
                else intSingularH_tria3(FLAG,SinNode,X,Xeval,SubA,normal);
            
                /* Assign contributions to the right places in the matrix */
                for(k = 0; k < NODES_IN_ELEM; k++){
                    F = A[interfaceID]*SubA[k];
                    G = B[interfaceID]*SubA[k];
                    currentNode = mElems[j][k] - 1;
                    mA[i][currentNode] += F;
                    mA[i][currentNode + nNodes] += G;
                    mA[i + nNodes][currentNode] -= G;
                    mA[i + nNodes][currentNode + nNodes] += F;
                }
            }
        }

        /* Assemble right hand side vector */
        vB[i] = mBC[i][0];
        vB[i + nNodes] = mBC[i + nNodes][0];
    }

    /* Free only when necessary */
    if(nInterfaces){
        free(A);
        free(B);
    }

    return 0;
}
