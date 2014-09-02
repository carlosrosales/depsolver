/******************************************************************************
* File     : force.c                                                          *
* Author   : Carlos Rosales Fernandez (carlos.rosales.fernandez(at)gmail.com) *
* Revision : 2.1 (2008-11-10)                                                 *
******************************************************************************/
/**
 * @brief Dielectrophoretic force calculation functions
 *
 * @file
 * Definitions of the six functions responsible for the calculation of
 * dielectrophoretic (DEP) force in depSolver: forceEllipsoid_tria3,
 * forceEllipsoid_tria6, forceMST_tria3, forceMST_tria3, forceMultipole_tria3,
 * forceMultipole_tria6.
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
#include "depolarization.h"
#include "field.h"
#include "force.h"
#include "forceIntegrals.h"
#include "globalToLocal.h"
#include "memHandler.h"
#include "shapeFunctions.h"

extern unsigned int nElems, nNodes;


/**
 * Calculates the DEP force at dielectric interfaces using a dipolar
 * approximation for an ellipsoid of semiaxis given by input array \a axis.
 * The force is returned in array \a F.
 *
 * Works for linear interpolation in triangular elements (3-noded triangles).
 *
 * @param ANALYSIS   : [ Input ] Analysis type, homogeneous (10) or shelled ellipsoid (11)
 * @param nOrder     : [ Input ] Multipole approximation order
 * @param axis       : [ Input ] Ellipsoid semi-axis lengths
 * @param Xeval      : [ Input ] Coordinates of evaluation point (ellipsoid center)
 * @param mNodes     : [ Input ] Coordinates of all nodes in the geometry
 * @param mElems     : [ Input ] Connectivity of all nodes in the geometry
 * @param vB         : [ Input ] Solution vector of the electrostatic problem
 * @param vMatParam  : [ Input ] Material properties
 * @param vProbParam : [ Input ] Frequency of the applied electric field
 * @param F          : [ Output ] Dielectrophoretic (DEP) force vector
 */
int forceEllipsoid_tria3(unsigned int ANALYSIS,
                         unsigned int nOrder,
                         double *axis,
                         double *Xeval,
                         double **mNodes,
                         unsigned int **mElems,
                         double *vB,
                         double **vMatParam,
                         double *vProbParam,
                         double *F)
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
    K[0] = K[1] = K[2] = 0.0;
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


/**
 * Calculates the DEP force at dielectric interfaces using a dipolar
 * approximation for an ellipsoid of semiaxis given by input array \a axis.
 * The force is returned in array \a F.
 *
 * Works for quadratic interpolation in triangular elements (6-noded triangles).
 *
 * @param ANALYSIS   : [ Input ] Analysis type, homogeneous (10) or shelled ellipsoid (11)
 * @param nOrder     : [ Input ] Multipole approximation order
 * @param axis       : [ Input ] Ellipsoid semi-axis lengths
 * @param Xeval      : [ Input ] Coordinates of evaluation point (ellipsoid center)
 * @param mNodes     : [ Input ] Coordinates of all nodes in the geometry
 * @param mElems     : [ Input ] Connectivity of all nodes in the geometry
 * @param vB         : [ Input ] Solution vector of the electrostatic problem
 * @param vMatParam  : [ Input ] Material properties
 * @param vProbParam : [ Input ] Frequency of the applied electric field
 * @param F          : [ Output ] Dielectrophoretic (DEP) force vector
 */
int forceEllipsoid_tria6(unsigned int ANALYSIS,
                         unsigned int nOrder,
                         double *axis,
                         double *Xeval,
                         double **mNodes,
                         unsigned int **mElems,
                         double *vB,
                         double **vMatParam,
                         double *vProbParam,
                         double *F)
{
    const unsigned int ELEMS = nElems, NODES_IN_ELEM = 6;
    unsigned int currentNode, i, j, k, p, q, r;
    double A, B, C, Eps, EpsX, EpsY, EpsZ, SigX, SigY, SigZ, w;
    double K[3], EpsDif[3], EpsPlus[3], SigDif[3], SigPlus[3], mE[6], Ld[3];
    double X[6][3];
    double E[7][7][7];
    double intE[6][7][7][7];

    /* Initialize */
    w = vProbParam[0]*pi2*eps0;
    Eps = vMatParam[0][1]*eps0;
    A = pi2*Eps*axis[0]*axis[1]*axis[2]/3.0;
    B = 1.0/(pi4*eps0);
    C = axis[0]*axis[2]/(2.0*axis[2] + axis[0]);
    F[0] = F[1] = F[2] = 0.0;
    K[0] = K[1] = K[2] = 0.0;
    for(i = 0; i < 7; i++)
        for(j = 0; j < 7 ; j++)
            for(k = 0; k < 7; k++) E[i][j][k] = 0.0;

    /* Electric field value at the required position */
    field_tria6(Xeval,mNodes,mElems,vB,mE);

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
        
        intDE_tria6(nOrder,X,Xeval,intE);
        for(j = 0; j < NODES_IN_ELEM; j++)
            for(p = 0; p < 7; p++)
                for(q = 0; q < 7 ; q++)
                    for(r = 0; r < 7; r++)
                        E[p][q][r] += B*intE[j][p][q][r]*vB[mElems[i][j] - 1];
    }

    /****************************************************/
    /*** CALCULATE FORCE USING DIPOLAR APPROXIMATION  ***/
    /****************************************************/

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
    F[0] = A*(K[0]*E[1][0][0]*E[2][0][0]+K[1]*E[0][1][0]*E[1][1][0]+K[2]*E[0][0][1]*E[1][0][1]);    /*  Fx */
    F[1] = A*(K[0]*E[1][0][0]*E[1][1][0]+K[1]*E[0][1][0]*E[0][2][0]+K[2]*E[0][0][1]*E[0][1][1]);    /* Fy */
    F[2] = A*(K[0]*E[1][0][0]*E[1][0][1]+K[1]*E[0][1][0]*E[0][1][1]+K[2]*E[0][0][1]*E[0][0][2]);    /* Fz */

    return 0;
}


/**
 * Calculates the DEP force at dielectric interfaces using the Maxwell Stress
 * Tensor method. Points where the force needs to be calculated are identified
 * by a boundary condition type 6. Returns the total force on the particle and
 * and application point (the particle center of mass) in arrays \a Fcm and
 * \a Xcm.
 *
 * Works for linear interpolation in triangular elements (3-noded triangles).
 *
 * @param mNodes  : [ Input ] Coordinates of all nodes in the geometry
 * @param mElems  : [ Input ] Connectivity of all nodes in the geometry
 * @param vBCType : [ Input ] Boundary condition type for each node in the computational domain
 * @param vB      : [ Input ] Solution vector of the electrostatic problem
 * @param eps     : [ Input ] Electric permittivity of the fluid
 * @param Fce     : [ Output ] DEP force vector at each of the elements in the particle surface
 * @param Xce     : [ Output ] Geometrical center of each of the elements in the surface of the particle
 * @param Fcm     : [ Output ] Total dielectrophoretic (DEP) force vector
 * @param Xcm     : [ Output ] Coordinates of the center of mass of the particle
 *
 * @warning
 * Notice that it is assumed that only one particle is considered in the
 * calculation, and that we assume that its mass density is uniform.
 */
int forceMST_tria3(double **mNodes,
                   unsigned int **mElems,
                   unsigned int *vBCType,
                   double *vB,
                   double Eps,
                   double **Fce,
                   double **Xce,
                   double *Fcm,
                   double *Xcm)
{
    const unsigned int ELEMS = nElems, NODES = nNodes, NODES_IN_ELEM = 3;
    unsigned int currentNode, i, j, N;
    double A, J;
    double L[2], M[3], mF[3], mE[6], Xin[3];
    double X[3][3], E[3][3];
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

            field_tria3(Xin,mNodes,mElems,vB,mE);

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
            intF_tria3(X,E,mF);

            /* Element centre */
            shape_tria3(L,M);
            J = X2L_tria3(X,Xin,L,M);

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


/**
 * Calculates the DEP force at dielectric interfaces using the Maxwell Stress
 * Tensor method. Points where the force needs to be calculated are identified
 * by a boundary condition type 6. Returns the total force on the particle and
 * and application point (the particle center of mass) in arrays \a F[] and
 * \a Xcm[].
 *
 * Works for quadratic interpolation in triangular elements (6-noded triangles).
 *
 * @param mNodes  : [ Input ] Coordinates of all nodes in the geometry
 * @param mElems  : [ Input ] Connectivity of all nodes in the geometry
 * @param vBCType : [ Input ] Boundary condition type for each node in the computational domain
 * @param vB      : [ Input ] Solution vector of the electrostatic problem
 * @param eps     : [ Input ] Electric permittivity of the fluid
 * @param Fce     : [ Output ] DEP force vector at each of the elements in the particle surface
 * @param Xce     : [ Output ] Geometrical center of each of the elements in the surface of the particle
 * @param Fcm     : [ Output ] Total dielectrophoretic (DEP) force vector
 * @param Xcm     : [ Output ] Coordinates of the center of mass of the particle
 *
 * @warning
 * Notice that it is assumed that only one particle is considered in the
 * calculation, and that we assume that its mass density is uniform.
 */
int forceMST_tria6(double **mNodes,
                   unsigned int **mElems,
                   unsigned int *vBCType,
                   double *vB,
                   double Eps,
                   double **Fce,
                   double **Xce,
                   double *Fcm,
                   double *Xcm)
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


/**
 * Calculates the DEP force at dielectric interfaces using a multipolar
 * approximation of order \a nOrder for a sphere of radius \a R. The force is
 * returned in array \a F.
 *
 * Works for linear interpolation in triangular elements (3-noded triangles).
 *
 * @param nOrder     : [ Input ] Multipole approximation order
 * @param R          : [ Input ] Radius of the spherical particle
 * @param Xeval      : [ Input ] Coordinates of evaluation point (sphere center)
 * @param mNodes     : [ Input ] Coordinates of all nodes in the geometry
 * @param mElems     : [ Input ] Connectivity of all nodes in the geometry
 * @param vB         : [ Input ] Solution vector of the electrostatic problem
 * @param vMatParam  : [ Input ] Material properties
 * @param vProbParam : [ Input ] Frequency of the applied electric field
 * @param F          : [ Output ] Dielectrophoretic (DEP) force vector
 */
int forceMultipole_tria3(unsigned int nOrder,
                         double R,
                         double *Xeval,
                         double **mNodes,
                         unsigned int **mElems,
                         double *vB,
                         double **vMatParam,
                         double *vProbParam,
                         double *F)
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
    for(i = 0; i < 6; i++) C[i] = 0.0;
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


/**
 * Calculates the DEP force at dielectric interfaces using a multipolar
 * approximation of order \a nOrder for a sphere of radius \a R. The force is
 *  returned in array \a F.
 *
 * Works for quadratic interpolation in triangular elements (6-noded triangles).
 *
 *
 * @param nOrder     : [ Input ] Multipole approximation order
 * @param R          : [ Input ] Radius of the spherical particle
 * @param Xeval      : [ Input ] Coordinates of evaluation point (sphere center)
 * @param mNodes     : [ Input ] Coordinates of all nodes in the geometry
 * @param mElems     : [ Input ] Connectivity of all nodes in the geometry
 * @param vB         : [ Input ] Solution vector of the electrostatic problem
 * @param vMatParam  : [ Input ] Material properties
 * @param vProbParam : [ Input ] Frequency of the applied electric field
 * @param F          : [ Output ] Dielectrophoretic (DEP) force vector
 */
int forceMultipole_tria6(unsigned int nOrder,
                         double R,
                         double *Xeval,
                         double **mNodes,
                         unsigned int **mElems,
                         double *vB,
                         double **vMatParam,
                         double *vProbParam,
                         double *F)
{
    const unsigned int ELEMS = nElems, NODES_IN_ELEM = 6;
    unsigned int currentNode, i, j, p, q, r;
    double B,CM1, CM2, CM3, CM4, CM5, Eps, EpsDif, EpsPlus, SigDif, SigPlus, w;
    double C[6], mE[6];
    double X[6][3];
    double E[7][7][7];
    double intE[6][7][7][7];
    
    /* Initialize */
    w = vProbParam[0]*pi2*eps0;
    Eps = vMatParam[0][1]*eps0;
    B = 1.0/(pi4*eps0);
    F[0] = F[1] = F[2] = 0.0;
    for(i = 0; i < 6; i++) C[i] = 0.0;
    for(p = 0; p < 7; p++)
        for(q = 0; q < 7 ; q++)
            for(r = 0; r < 7; r++) E[p][q][r] = 0.0;

    /* Electric field value at the required position */
    field_tria6(Xeval,mNodes,mElems,vB,mE);

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
        
        intDE_tria6(nOrder,X,Xeval,intE);
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
        C[2] = C[1]*R*R*0.333333333333333333;
        
        if(nOrder > 2){
            EpsPlus = w*(3.0*vMatParam[1][1] + 4.0*vMatParam[0][1]);
            SigPlus = 3.0*vMatParam[1][0] + 4.0*vMatParam[0][0];
            CM3 = (EpsDif*EpsPlus + SigDif*SigPlus)/(EpsPlus*EpsPlus + SigPlus*SigPlus);
            C[3] = C[2]*R*R*0.1;
            
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

