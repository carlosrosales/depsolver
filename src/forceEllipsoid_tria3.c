/******************************************************************************
* File      : forceEllipsoid_tria3.c                                          *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* Calculates the DEP force at dielectric interfaces using a Dipolar           *
* Approximation for an ellipsoid of semiaxis given by input array axis[].     *
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

#include "constants.h"
#include "forceEllipsoid_tria3.h"

int forceEllipsoid_tria3(unsigned int ANALYSIS, unsigned int nOrder, 
                         double *axis, double *Xeval, double **mNodes, 
                         unsigned int **mElems, double *vB, double **vMatParam,
                         double *vProbParam, double *F)
{
    const unsigned int ELEMS = nElems, NODES_IN_ELEM = 3;
    unsigned int currentNode, i, j, k, p, q, r;
    double A, B, C, Eps, EpsX, EpsY, EpsZ, SigX, SigY, SigZ, w;
    double K[3], EpsDif[3], EpsPlus[3], SigDif[3],SigPlus[3], mE[6], Ld[3];
    double X[3][3];
    double E[7][7][7];
    double intE[3][7][7][7];

    /* Initialize */
    w = vProbParam[0]*pi2*eps0;
    Eps = vMatParam[0][1]*eps0;
    A = pi2*Eps*axis[0]*axis[1]*axis[2]/3.0;
    B = 1.0/(pi4*eps0);
    C = axis[0]*axis[2]/(2.0*axis[2] + axis[0]);
    F[0] = F[1] = F[2] = 0.0;
    for(i = 0; i < 7; i++)
        for(j = 0; j < 7 ; j++)
            for(k = 0; k < 7; k++) E[i][j][k] = 0.0;

    /* Electric field value at the required position */
    field_tria3(Xeval,mNodes,mElems,vB,mE);

    /* Store only real part of the electric field */
    E[1][0][0] = mE[0];
    E[0][1][0] = mE[2];
    E[0][0][1] = mE[4];

    /* Electric field gradient at the required position */
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

    /***************************************************/
    /*** CALCULATE FORCE USING DIPOLAR APPROXIMATION ***/
    /***************************************************/

    /* Obtain depolarization factors */
    depolarization(axis,Ld);

    /* Generalized Clausius-Mosotti factors for homogeneous ellipsoid */
    if(ANALYSIS == 10){
        EpsDif[0] = w*(vMatParam[1][1] - vMatParam[0][1]);
        EpsDif[1] = w*(vMatParam[1][1] - vMatParam[0][1]);
        EpsDif[2] = w*(vMatParam[1][1] - vMatParam[0][1]);
        SigDif[0] = vMatParam[1][0] - vMatParam[0][0];
        SigDif[1] = vMatParam[1][0] - vMatParam[0][0];
        SigDif[2] = vMatParam[1][0] - vMatParam[0][0];
        EpsPlus[0] = w*(vMatParam[0][1] + (vMatParam[1][1] - vMatParam[0][1])*Ld[0]);
        EpsPlus[1] = w*(vMatParam[0][1] + (vMatParam[1][1] - vMatParam[0][1])*Ld[1]);
        EpsPlus[2] = w*(vMatParam[0][1] + (vMatParam[1][1] - vMatParam[0][1])*Ld[2]);
        SigPlus[0] = (vMatParam[0][0] + (vMatParam[1][0] - vMatParam[0][0])*Ld[0]);
        SigPlus[1] = (vMatParam[0][0] + (vMatParam[1][0] - vMatParam[0][0])*Ld[1]);
        SigPlus[2] = (vMatParam[0][0] + (vMatParam[1][0] - vMatParam[0][0])*Ld[2]);

        K[0] = (EpsDif[0]*EpsPlus[0] + SigDif[0]*SigPlus[0])/(EpsPlus[0]*EpsPlus[0] + SigPlus[0]*SigPlus[0]);   /* CMx */
        K[1] = (EpsDif[1]*EpsPlus[1] + SigDif[1]*SigPlus[1])/(EpsPlus[1]*EpsPlus[1] + SigPlus[1]*SigPlus[1]);   /* CMy */
        K[2] = (EpsDif[2]*EpsPlus[2] + SigDif[2]*SigPlus[2])/(EpsPlus[2]*EpsPlus[2] + SigPlus[2]*SigPlus[2]);   /* CMz */
    }

    /* Generalized Clausius-Mosotti factors for single-shell oblate ellipsoid */
    else if(ANALYSIS == 11){
        vMatParam[1][1] = vMatParam[1][1]*C;
        EpsX = vMatParam[1][1]/Ld[0];
        EpsY = vMatParam[1][1]/Ld[1];
        EpsZ = vMatParam[1][1]/Ld[2];
        vMatParam[1][0] = vMatParam[1][0]*C;
        SigX = vMatParam[1][0]/Ld[0];
        SigY = vMatParam[1][0]/Ld[1];
        SigZ = vMatParam[1][0]/Ld[2];
    
        EpsDif[0] = w*(EpsX - vMatParam[0][1]);
        EpsDif[1] = w*(EpsY - vMatParam[0][1]);
        EpsDif[2] = w*(EpsZ - vMatParam[0][1]);
        SigDif[0] = SigX - vMatParam[0][0];
        SigDif[1] = SigY - vMatParam[0][0];
        SigDif[2] = SigZ - vMatParam[0][0];
        EpsPlus[0] = w*(vMatParam[0][1] + (EpsX - vMatParam[0][1])*Ld[0]);
        EpsPlus[1] = w*(vMatParam[0][1] + (EpsY - vMatParam[0][1])*Ld[1]);
        EpsPlus[2] = w*(vMatParam[0][1] + (EpsZ - vMatParam[0][1])*Ld[2]);
        SigPlus[0] = (vMatParam[0][0] + (SigX - vMatParam[0][0])*Ld[0]);
        SigPlus[1] = (vMatParam[0][0] + (SigY - vMatParam[0][0])*Ld[1]);
        SigPlus[2] = (vMatParam[0][0] + (SigZ - vMatParam[0][0])*Ld[2]);

        K[0] = (EpsDif[0]*EpsPlus[0] + SigDif[0]*SigPlus[0])/(EpsPlus[0]*EpsPlus[0] + SigPlus[0]*SigPlus[0]);   /* CMx */
        K[1] = (EpsDif[1]*EpsPlus[1] + SigDif[1]*SigPlus[1])/(EpsPlus[1]*EpsPlus[1] + SigPlus[1]*SigPlus[1]);   /* CMy */
        K[2] = (EpsDif[2]*EpsPlus[2] + SigDif[2]*SigPlus[2])/(EpsPlus[2]*EpsPlus[2] + SigPlus[2]*SigPlus[2]);   /* CMz */
    }

    /* Force using dipolar approximation */
    F[0] = A*(K[0]*E[1][0][0]*E[2][0][0] + K[1]*E[0][1][0]*E[1][1][0] + K[2]*E[0][0][1]*E[1][0][1]);    /* Fx */
    F[1] = A*(K[0]*E[1][0][0]*E[1][1][0] + K[1]*E[0][1][0]*E[0][2][0] + K[2]*E[0][0][1]*E[0][1][1]);    /* Fy */
    F[2] = A*(K[0]*E[1][0][0]*E[1][0][1] + K[1]*E[0][1][0]*E[0][1][1] + K[2]*E[0][0][1]*E[0][0][2]);    /* Fz */

    return 0;
}
