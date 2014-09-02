/******************************************************************************
* File      : forceMultipole_tria3.c                                          *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* Calculates the DEP force at dielectric interfaces using a multipolar        *
* approximation of order nOrder for a sphere of radius R. The force is        *
* returned in array F[].                                                      *
* Works for linear interpolation in triangular elements (3-noded triangles).  *
* F[0] = Fx    |    F[1] = Fy    |    F[2] = Fz                               *
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
#include "forceMultipole_tria3.h"

int forceMultipole_tria3(unsigned int nOrder, double R, double *Xeval, 
                         double **mNodes, unsigned int **mElems, double *vB, 
                         double **vMatParam, double *vProbParam, double *F)
{
    const unsigned int ELEMS = nElems, NODES_IN_ELEM = 3;
    unsigned int currentNode, i, j, p, q, r;
    double B, CM1, CM2, CM3, CM4, CM5, Eps, EpsDif, EpsPlus, SigDif, SigPlus, w;
    double C[6], mE[6];
    double X[3][3];
    double E[7][7][7];
    double intE[3][7][7][7];

    /* Initialize */
    w = vProbParam[0]*pi2*eps0;
    Eps = vMatParam[0][1]*eps0;
    B = 1.0/(pi4*eps0);
    F[0] = F[1] = F[2] = 0.0;
    for(p = 0; p < 7; p++)
        for(q = 0; q < 7 ; q++)
            for(r = 0; r < 7; r++) E[p][q][r] = 0.0;        

    /* Electric field value at the required position */
    field_tria3(Xeval,mNodes,mElems,vB,mE);

    /* Store Re(E) for later use */
    E[1][0][0] = mE[0];
    E[0][1][0] = mE[2];
    E[0][0][1] = mE[4];

    /* Electric field derivatives at the required position */
    for(i = 0; i < ELEMS; i++){
        for(j = 0; j < NODES_IN_ELEM; j++){
            currentNode = mElems[i][j] - 1;
            X[j][0] = mNodes[currentNode][0];
            X[j][1] = mNodes[currentNode][1];
            X[j][2] = mNodes[currentNode][2];
        }
        
        intDE_tria3(nOrder,X,Xeval,intE);
        for(j = 0; j < NODES_IN_ELEM; j++)
            for(p = 0; p < 7; p++)
                for(q = 0; q < 7 ; q++)
                    for(r = 0; r < 7; r++)
                        E[p][q][r] += B*intE[j][p][q][r]*vB[mElems[i][j] - 1];
    }

    /********************************************************/
    /***  CALCULATE FORCE USING MULTIPOLAR APPROXIMATION  ***/
    /********************************************************/

    /* Clausius-Mossotti factors for different aproximation orders */
    EpsDif = w*(vMatParam[1][1] - vMatParam[0][1]);
    EpsPlus = w*(vMatParam[1][1] + 2.0*vMatParam[0][1]);
    SigDif = vMatParam[1][0] - vMatParam[0][0];
    SigPlus = vMatParam[1][0] + 2.0*vMatParam[0][0];
    CM1 = (EpsDif*EpsPlus + SigDif*SigPlus)/(EpsPlus*EpsPlus + SigPlus*SigPlus);
    C[1] = pi2*Eps*R*R*R;
    if(nOrder > 1){
        EpsPlus = w*(2.0*vMatParam[1][1] + 3.0*vMatParam[0][1]);
        SigPlus = 2.0*vMatParam[1][0] + 3.0*vMatParam[0][0];
        CM2 = (EpsDif*EpsPlus + SigDif*SigPlus)/(EpsPlus*EpsPlus + SigPlus*SigPlus);
        C[2] = C[1]*R*R/3.0;
        
        if(nOrder > 2){
            EpsPlus = w*(3.0*vMatParam[1][1] + 4.0*vMatParam[0][1]);
            SigPlus = 3.0*vMatParam[1][0] + 4.0*vMatParam[0][0];
            CM3 = (EpsDif*EpsPlus + SigDif*SigPlus)/(EpsPlus*EpsPlus + SigPlus*SigPlus);
            C[3] = C[2]*R*R/10.0;
            
            if(nOrder > 3){
                EpsPlus = w*(4.0*vMatParam[1][1] + 5.0*vMatParam[0][1]);
                SigPlus = 4.0*vMatParam[1][0] + 5.0*vMatParam[0][0];
                CM4 = (EpsDif*EpsPlus + SigDif*SigPlus)/(EpsPlus*EpsPlus + SigPlus*SigPlus);
                C[4] = C[3]*R*R/21.0;
                
                if(nOrder > 4){
                    EpsPlus = w*(5.0*vMatParam[1][1] + 6.0*vMatParam[0][1]);
                    SigPlus = 5.0*vMatParam[1][0] + 6.0*vMatParam[0][0];
                    CM5 = (EpsDif*EpsPlus + SigDif*SigPlus)/(EpsPlus*EpsPlus + SigPlus*SigPlus);
                    C[5] = C[4]*R*R*CM5/36.0;
                }
                C[4] = C[4]*CM4;
            }
            C[3] = C[3]*CM3;
        }
        C[2] = C[2]*CM2;
    }
    C[1] = C[1]*CM1;

    /* Contribution to the force from dipolar term (n = 1) */
    F[0] = C[1]*(E[1][0][0]*E[2][0][0] + E[0][1][0]*E[1][1][0] + E[0][0][1]*E[1][0][1]); /* Fx (n=1) */
    F[1] = C[1]*(E[1][0][0]*E[1][1][0] + E[0][1][0]*E[0][2][0] + E[0][0][1]*E[0][1][1]); /* Fy (n=1) */
    F[2] = C[1]*(E[1][0][0]*E[1][0][1] + E[0][1][0]*E[0][1][1] + E[0][0][1]*E[0][0][2]); /* Fz (n=1) */

    /* Contribution to the force from quadrupolar term (n = 2) */    
    if(nOrder > 1){
        F[0] += C[2]*(E[2][0][0]*E[3][0][0] + E[0][2][0]*E[1][2][0] + E[0][0][2]*E[1][0][2] + 2.0*(E[1][1][0]*E[2][1][0] + E[1][0][1]*E[2][0][1] + E[0][1][1]*E[1][1][1])); /* Fx (n=2) */
        F[1] += C[2]*(E[2][0][0]*E[2][1][0] + E[0][2][0]*E[0][3][0] + E[0][0][2]*E[0][1][2] + 2.0*(E[1][1][0]*E[1][2][0] + E[1][0][1]*E[1][1][1] + E[0][1][1]*E[0][2][1])); /* Fy (n=2) */
        F[2] += C[2]*(E[2][0][0]*E[2][0][1] + E[0][2][0]*E[0][2][1] + E[0][0][2]*E[0][0][3] + 2.0*(E[1][1][0]*E[1][1][1] + E[1][0][1]*E[1][0][2] + E[0][1][1]*E[0][1][2])); /* Fz (n=2) */

        /* Contribution to the force from octupolar term (n = 3) */
        if(nOrder > 2){
            F[0] += C[3]*(E[3][0][0]*E[4][0][0] + E[0][3][0]*E[1][3][0] + E[0][0][3]*E[1][0][3] + 3.0*(E[2][1][0]*E[3][1][0] + E[2][0][1]*E[3][0][1] + E[1][2][0]*E[2][2][0] + E[1][0][2]*E[2][0][2] + E[0][2][1]*E[1][2][1] + E[0][1][2]*E[1][1][2]) + 6.0*E[1][1][1]*E[2][1][1]); /* Fx (n=3) */
            F[1] += C[3]*(E[3][0][0]*E[3][1][0] + E[0][3][0]*E[0][4][0] + E[0][0][3]*E[0][1][3] + 3.0*(E[2][1][0]*E[2][2][0] + E[2][0][1]*E[2][1][1] + E[1][2][0]*E[1][3][0] + E[1][0][2]*E[1][1][2] + E[0][2][1]*E[0][3][1] + E[0][1][2]*E[0][2][2]) + 6.0*E[1][1][1]*E[1][2][1]); /* Fy (n=3) */
            F[2] += C[3]*(E[3][0][0]*E[3][0][1] + E[0][3][0]*E[0][3][1] + E[0][0][3]*E[0][0][4] + 3.0*(E[2][1][0]*E[2][1][1] + E[2][0][1]*E[2][0][2] + E[1][2][0]*E[1][2][1] + E[1][0][2]*E[1][0][3] + E[0][2][1]*E[0][2][2] + E[0][1][2]*E[0][1][3]) + 6.0*E[1][1][1]*E[1][1][2]); /* Fz (n=3) */
        }
    }

    return 0;
}
