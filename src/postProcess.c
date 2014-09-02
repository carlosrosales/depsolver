/******************************************************************************
* File     : postProcess.c                                                    *
* Author   : Carlos Rosales Fernandez (carlos.rosales.fernandez(at)gmail.com) *
* Revision : 2.0 (2008-06-16)                                                 *
******************************************************************************/
/**
 * @brief Post-processing functions
 *
 * @file
 * This file contains the definitions of the seven post-processing
 * driver functions: headerSetup, postProcessSetup, multipoleSetup,
 * postProcessCol_tria3, postProcessCol_tria6, postProcess_tria3,
 * postProcess_tria6
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "constants.h"
#include "memHandler.h"
#include "field.h"
#include "force.h"
#include "postProcess.h"
#include "potential.h"

extern FILE *file_log;
extern unsigned int nElems, nInternalPoints, nFPoints, nCols, nColType,
                    NX, NY, NZ;
extern double deltaX, deltaY, deltaZ;


/**
 * Create headers for the output files depending on the values of the
 * \a ANALYSIS and \a cOutputType variables.
 *
 * @param ANALYSIS    : [ Input ] Type of post-processing to be done
 * @param cOutputType : [ Input ] Output format type
 * @param fp          : [ Output ] File pointer to post-processing files
 * @param Xinner      : [ Input ] Internal nodes where post-processing is to be done
 */
int headerSetup(unsigned int ANALYSIS,
                char *cOutputType,
                FILE **fp,
                double **Xinner)
{

    if( strcmp(cOutputType,"STD") == 0 ){
        switch(ANALYSIS){

            case(0):
            fp[0] = fopen("potential.dat","w");
            fprintf(fp[0],"#x\t\ty\t\tz\t\tRe[Pot(x,y,z)]\tIm[Pot(x,y,z)]\n");
            break;

            case(1):
            fp[1] = fopen("field.dat","w");
            fprintf(fp[1],"#x\t\ty\t\tz\t\tRe[Ex(x,y,z)]\tIm[Ex(x,y,z)]\tRe[Ey(x,y,z)]");
            fprintf(fp[1],"\tIm[Ey(x,y,z)]\tRe[Ez(x,y,z)]\tIm[Ez(x,y,z)]\n");
            break;

            default:
            fp[0] = fopen("potential.dat","w");
            fp[1] = fopen("field.dat","w");
            fprintf(fp[0],"#x\t\ty\t\tz\t\tRe[Pot(x,y,z)]\tIm[Pot(x,y,z)]\n");
            fprintf(fp[1],"#x\t\ty\t\tz\t\tRe[Ex(x,y,z)] \tIm[Ex(x,y,z)]\tRe[Ey(x,y,z)]");
            fprintf(fp[1],"\tIm[Ey(x,y,z)] \tRe[Ez(x,y,z)]\tIm[Ez(x,y,z)]\n");
            break;
        }
    }
    else if( strcmp(cOutputType,"VTK") == 0 ){
        fp[0] = fopen("results.vtk","w");
        fprintf(fp[0],"# vtk DataFile Version 2.0\n");
        fprintf(fp[0],"depSolver v2.0\n");
        fprintf(fp[0],"ASCII\n");
        fprintf(fp[0],"\n");
        fprintf(fp[0],"DATASET STRUCTURED_POINTS\n");
        fprintf(fp[0],"DIMENSIONS %d %d %d\n",NX,NY,NZ);
        fprintf(fp[0],"ORIGIN %e %e %e\n",Xinner[0][0],Xinner[0][1],Xinner[0][2]);
        fprintf(fp[0],"SPACING %e %e %e\n",deltaX,deltaY,deltaZ);
        fprintf(fp[0],"\n");
        fprintf(fp[0],"POINT_DATA %d\n",NX*NY*NZ);
    }

    if(ANALYSIS == 3 || ANALYSIS == 4){
        fp[2] = fopen("force-mst.dat","w");
        fprintf(fp[2],"#Xcm\t\tYcm\t\tZcm\t\tFXcm\t\tFYcm\t\tFZcm\n");
    }
    else if(ANALYSIS > 4){
        fp[2] = fopen("force-mp.dat","w");
        fprintf(fp[2],"#Xcm\t\tYcm\t\tZcm\t\tFXcm\t\tFYcm\t\tFZcm\n");
    }

    return 0;
}


/**
 * Filter input data to eliminate points inside columns from the post-processing
 * list of points if necessary (STD output only).
 *
 * @param cOutputType : [ Input ] Output format type
 * @param xCols       : [ Input ] Coordinates and radius of columns
 * @param XinerTemp   : [ Input ] Original internal nodes from input files
 * @param Xinner      : [ Output ] Internal nodes where post-processing is to be done
 */
int postProcessSetup(char *cOutputType,
                     double **xCols,
                     double **XinnerTemp,
                     double **Xinner)
{
    unsigned int i, j, nInside, count;
    double dx, dy, x, y, Rsq;

    /* Filter data to avoid calculation inside columns */
    if( nInternalPoints > 0 ){
        if( nCols > 0 && strcmp(cOutputType,"STD") == 0 ){

            nInside = 0;
            for(i = 0; i < nInternalPoints; i++){
                x = XinnerTemp[i][0];
                y = XinnerTemp[i][1];

                count = 0;
                for( j = 0; j < nCols; j++ ){
                    dx  = fabs( x - xCols[j][0] );
                    dy  = fabs( y - xCols[j][1] );
                    Rsq = xCols[j][2]*xCols[j][2];

                    if( nColType == 1 ){
                        if( ( dx*dx +dy*dy ) > Rsq ) count++;
                    }
                    else{
                        if( (dx > xCols[j][2]) && (dy > xCols[j][2]) ) count++;
                    }
                }
                if(count == nCols){
                    XinnerTemp[i][3] = 1.0;
                    nInside++;
                }
                else XinnerTemp[i][3] = -1.0;
            }

            count = 0;
            for( i = 0; i< nInternalPoints; i ++){
                if( XinnerTemp[i][3] > 0.0 ){
                    Xinner[count][0] = XinnerTemp[i][0];
                    Xinner[count][1] = XinnerTemp[i][1];
                    Xinner[count][2] = XinnerTemp[i][2];
                    count++;
                }
            }

            /* Redefine the number of evaluation points */
            nInternalPoints = nInside;
        }
        else{
            for( i = 0; i< nInternalPoints; i ++){
                Xinner[i][0] = XinnerTemp[i][0];
                Xinner[i][1] = XinnerTemp[i][1];
                Xinner[i][2] = XinnerTemp[i][2];
            }
        }
    }

    return 0;
}

/**
 * Set up the type of multipolar approximation (order and shape) required.
 *
 * @param ANALYSIS : [ Input ] Type of post-processing to be done
 * @param axis     : [ Input ] Semi-axis of ellipsoidal particle
 * @param nShape   : [ Output ] Sphere (0) or ellipsoid (1)
 * @param nOrder   : [ Output ] Multipolar approximation order
 * @param R        : [ Output ] Radius of spherical particle
 */
int multipoleSetup(unsigned int ANALYSIS,
                   double *axis,
                   unsigned int *nShape,
                   unsigned int *nOrder,
                   double *R)
{
    /* Default to sphere with dipolar approximation */
    (*nOrder) = 1;
    (*nShape) = 0;
    (*R)      = axis[0];

    switch(ANALYSIS){
        case 5:         /* Dipolar Approximation (Order 1) */
            (*nOrder) = 1;
            break;
        case 6:         /* Quadrupolar Approximation (Order 2) */
            (*nOrder) = 2;
            break;
        case 7:         /* Octupolar Approximation (Order 3) */
            (*nOrder) = 3;
            break;
        case 8:         /* Multipolar Approx Order 4 */
            (*nOrder) = 4;
            break;
        case 9:         /* Multipolar Approx Order 5 */
            (*nOrder) = 5;
            break;
        case 10:        /* Dipolar Approx for Homogeneous Ellipsoid */
            (*nOrder) = 1;
            (*nShape) = 1;
            break;
        case 11:        /* Dipolar Approx for Single-shelled Oblate Ellipsoid */
            (*nOrder) = 1;
            (*nShape) = 1;
            break;
    }

    return 0;
}

/**
 * Calculates potential and electric field at the required points, and stores
 * them in files 'potential.dat' and 'field.dat' using the standard format.
 *
 * It also calculates the force on the dielectric interface if required and
 * stores it in file 'force-mst.dat' or 'force-mp.dat' for MST and multipolar
 * approximations respectively.
 *
 * Works for linear interpolation in triangular elements (3-noded triangles).
 *
 * @param ANALYSIS    : [ Input ] Type of post-processing to be done
 * @param cOutputType : [ Input ] Output format type
 * @param axis        : [ Input ] Semi-axis of ellipsoidal particle
 * @param XF          : [ Input ] Points where the force should be calculated
 * @param Xinner      : [ Input ] Nodes where potential and/or field should be calculated
 * @param mNodes  : [ Input ] Coordinates of all nodes in the computational domain
 * @param mElems  : [ Input ] Connectivity of all nodes in the computational domain
 * @param vMatParam  : [ Input ] Electric properties of dielectric materials
 * @param vProbParam : [ Input ] Frequency of the applied electric field
 * @param vBCType : [ Input ] Boundary condition type for each node in the domain
 * @param vB      : [ Input ] Solution vector for the electrostatic problem
 */
int postProcess_tria3(unsigned int ANALYSIS,
                      char *cOutputType,
                      double *axis,
                      double **XF,
                      double **Xinner,
                      double **mNodes,
                      unsigned int **mElems,
                      double **vMatParam,
                      double *vProbParam,
                      unsigned int *vBCType,
                      double *vB)
{
    FILE *fp[3];
    unsigned int i, nOrder, nShape;
    double Eps, R;
    double Fcm[3], mPot[2], mE[6], Xin[3], Xcm[3], Xeval[3];
    double **Fce, **Xce;

    /* Setup the header of the output files */
    headerSetup(ANALYSIS,cOutputType,fp,Xinner);

    /* Electric field and potential values at the required points */
    if( nInternalPoints > 0 ){
        mPot[0] = mPot[1] = 0.0;
        for(i = 0; i < 6; i++) mE[i] = 0.0;

        for(i = 0; i < nInternalPoints; i++){
            Xin[0] = Xinner[i][0];
            Xin[1] = Xinner[i][1];
            Xin[2] = Xinner[i][2];

            if(ANALYSIS == 0){
                potential_tria3(Xin,mNodes,mElems,vB,mPot);
                fprintf(fp[0],"%e\t%e\t%e\t",Xin[0],Xin[1],Xin[2]);
                fprintf(fp[0],"%e\t%e\n",mPot[0],mPot[1]);
            }
            else if(ANALYSIS == 1){
                field_tria3(Xin,mNodes,mElems,vB,mE);
                fprintf(fp[1],"%e\t%e\t%e\t",Xin[0],Xin[1],Xin[2]);
                fprintf(fp[1],"%e\t%e\t%e\t",mE[0],mE[1],mE[2]);
                fprintf(fp[1],"%e\t%e\t%e\n",mE[3],mE[4],mE[5]);
            }
            else if(ANALYSIS != 4){
                potential_tria3(Xin,mNodes,mElems,vB,mPot);
                field_tria3(Xin,mNodes,mElems,vB,mE);
                fprintf(fp[0],"%e\t%e\t%e\t",Xin[0],Xin[1],Xin[2]);
                fprintf(fp[0],"%e\t%e\n",mPot[0],mPot[1]);
                fprintf(fp[1],"%e\t%e\t%e\t",Xin[0],Xin[1],Xin[2]);
                fprintf(fp[1],"%e\t%e\t%e\t",mE[0],mE[1],mE[2]);
                fprintf(fp[1],"%e\t%e\t%e\n",mE[3],mE[4],mE[5]);
            }
        }
        if(ANALYSIS != 0) fclose(fp[1]);
        if(ANALYSIS != 1) fclose(fp[0]);
    }

    /* Calculate Fdep at dielectric interfaces if required */
    if(ANALYSIS == 3 || ANALYSIS == 4){
        Fce = doubleMatrix(nElems,3,1);
        Xce = doubleMatrix(nElems,3,1);
        Eps = vMatParam[0][1]*eps0;
        forceMST_tria3(mNodes,mElems,vBCType,vB,Eps,Fce,Xce,Fcm,Xcm);
        fprintf(fp[2],"%e\t%e\t%e\t%e\t%e\t%e\n",Xcm[0],Xcm[1],Xcm[2],Fcm[0],Fcm[1],Fcm[2]);
        fprintf(fp[2],"#Xce\t\tYce\t\tZce\t\tFXce\t\tFYce\t\tFZce\n");
        for(i = 0; i < nElems; i++){
            if(vBCType[mElems[i][0]-1] == 6){
                fprintf(fp[2],"%e\t%e\t%e\t",Xce[i][0],Xce[i][1],Xce[i][2]);
                fprintf(fp[2],"%e\t%e\t%e\n",Fce[i][0],Fce[i][1],Fce[i][2]);
            }
        }
        fclose(fp[2]);
        freeDoubleMatrix(Fce,nElems);
        freeDoubleMatrix(Xce,nElems);
    }

    /* Calculate Fdep at particle centre using multipolar method if required */
    else if(ANALYSIS > 4){
        multipoleSetup(ANALYSIS,axis,&nShape,&nOrder,&R);
        for(i = 0; i < nFPoints; i++){
            Xeval[0] = XF[i][0];
            Xeval[1] = XF[i][1];
            Xeval[2] = XF[i][2];
            if(nShape == 0) forceMultipole_tria3(nOrder,R,Xeval,mNodes,mElems,vB,vMatParam,vProbParam,Fcm);
            else forceEllipsoid_tria3(ANALYSIS,nOrder,axis,Xeval,mNodes,mElems,vB,vMatParam,vProbParam,Fcm);
            fprintf(fp[2],"%e\t%e\t%e\t",Xeval[0],Xeval[1],Xeval[2]);
            fprintf(fp[2],"%e\t%e\t%e\n",Fcm[0],Fcm[1],Fcm[2]);
        }
        fclose(fp[2]);
    }

    return 0;
}


/**
 * Calculates potential and electric field at the required points, and stores
 * them in files 'potential.dat' and 'field.dat' using the standard format.
 *
 * It also calculates the force on the dielectric interface if required and
 * stores it in file 'force-mst.dat' or 'force-mp.dat' for MST and multipolar
 * approximations respectively.
 *
 * Works for quadratic interpolation in triangular elements (6-noded triangles).
 *
 * @param ANALYSIS    : [ Input ] Type of post-processing to be done
 * @param cOutputType : [ Input ] Output format type
 * @param axis        : [ Input ] Semi-axis of ellipsoidal particle
 * @param XF          : [ Input ] Points where the force should be calculated
 * @param Xinner      : [ Input ] Nodes where potential and/or field should be calculated
 * @param mNodes  : [ Input ] Coordinates of all nodes in the computational domain
 * @param mElems  : [ Input ] Connectivity of all nodes in the computational domain
 * @param vMatParam  : [ Input ] Electric properties of dielectric materials
 * @param vProbParam : [ Input ] Frequency of the applied electric field
 * @param vBCType : [ Input ] Boundary condition type for each node in the domain
 * @param vB      : [ Input ] Solution vector for the electrostatic problem
 */
int postProcess_tria6(unsigned int ANALYSIS,
                      char *cOutputType,
                      double *axis,
                      double **XF,
                      double **Xinner,
                      double **mNodes,
                      unsigned int **mElems,
                      double **vMatParam,
                      double *vProbParam,
                      unsigned int *vBCType,
                      double *vB)
{
    FILE *fp[3];
    unsigned int i, nOrder, nShape;
    double Eps, R;
    double Fcm[3], mPot[2], mE[6], Xin[3], Xcm[3], Xeval[3];
    double **Fce, **Xce;

    /* Setup the header of the output files */
    headerSetup(ANALYSIS,cOutputType,fp,Xinner);

    /* Electric field and potential values at the required points */
    if( nInternalPoints > 0 ){
        mPot[0] = mPot[1] = 0.0;
        for(i = 0; i < 6; i++) mE[i] = 0.0;

        for(i = 0; i < nInternalPoints; i++){
            Xin[0] = Xinner[i][0];
            Xin[1] = Xinner[i][1];
            Xin[2] = Xinner[i][2];

            if(ANALYSIS == 0){
                potential_tria6(Xin,mNodes,mElems,vB,mPot);
                fprintf(fp[0],"%e\t%e\t%e\t",Xin[0],Xin[1],Xin[2]);
                fprintf(fp[0],"%e\t%e\n",mPot[0],mPot[1]);
            }
            else if(ANALYSIS == 1){
                field_tria6(Xin,mNodes,mElems,vB,mE);
                fprintf(fp[1],"%e\t%e\t%e\t",Xin[0],Xin[1],Xin[2]);
                fprintf(fp[1],"%e\t%e\t%e\t",mE[0],mE[1],mE[2]);
                fprintf(fp[1],"%e\t%e\t%e\n",mE[3],mE[4],mE[5]);
            }
            else{
                potential_tria6(Xin,mNodes,mElems,vB,mPot);
                field_tria6(Xin,mNodes,mElems,vB,mE);
                fprintf(fp[0],"%e\t%e\t%e\t",Xin[0],Xin[1],Xin[2]);
                fprintf(fp[0],"%e\t%e\n",mPot[0],mPot[1]);
                fprintf(fp[1],"%e\t%e\t%e\t",Xin[0],Xin[1],Xin[2]);
                fprintf(fp[1],"%e\t%e\t%e\t",mE[0],mE[1],mE[2]);
                fprintf(fp[1],"%e\t%e\t%e\n",mE[3],mE[4],mE[5]);
            }
        }
        if(ANALYSIS != 0) fclose(fp[1]);
        if(ANALYSIS != 1) fclose(fp[0]);
    }

    /* Calculate Fdep at dielectric interfaces if required */
    if(ANALYSIS == 3 || ANALYSIS == 4){
        Fce = doubleMatrix(nElems,3,1);
        Xce = doubleMatrix(nElems,3,1);
        Eps = vMatParam[0][1]*eps0;
        forceMST_tria6(mNodes,mElems,vBCType,vB,Eps,Fce,Xce,Fcm,Xcm);
        fprintf(fp[2],"%e\t%e\t%e\t%e\t%e\t%e\n",Xcm[0],Xcm[1],Xcm[2],Fcm[0],Fcm[1],Fcm[2]);
        fprintf(fp[2],"#Xce\t\tYce\t\tZce\t\tFXce\t\tFYce\t\tFZce\n");
        for(i = 0; i < nElems; i++){
            if(vBCType[mElems[i][0]-1] == 6){
                fprintf(fp[2],"%e\t%e\t%e\t",Xce[i][0],Xce[i][1],Xce[i][2]);
                fprintf(fp[2],"%e\t%e\t%e\n",Fce[i][0],Fce[i][1],Fce[i][2]);
            }
        }
        fclose(fp[2]);
        freeDoubleMatrix(Fce,nElems);
        freeDoubleMatrix(Xce,nElems);
    }

    /* Calculate Fdep at particle centre using multipolar method if required */
    else if(ANALYSIS > 4){
        multipoleSetup(ANALYSIS,axis,&nShape,&nOrder,&R);
        for(i = 0; i < nFPoints; i++){
            Xeval[0] = XF[i][0];
            Xeval[1] = XF[i][1];
            Xeval[2] = XF[i][2];
            if(nShape == 0) forceMultipole_tria6(nOrder,R,Xeval,mNodes,mElems,vB,vMatParam,vProbParam,Fcm);
            else forceEllipsoid_tria6(ANALYSIS,nOrder,axis,Xeval,mNodes,mElems,vB,vMatParam,vProbParam,Fcm);
            fprintf(fp[2],"%e\t%e\t%e\t",Xeval[0],Xeval[1],Xeval[2]);
            fprintf(fp[2],"%e\t%e\t%e\n",Fcm[0],Fcm[1],Fcm[2]);
        }
        fclose(fp[2]);
    }

    return 0;
}


/**
 * Calculates potential and electric field at the required points, and stores
 * them in file 'results.vtk' using the vtk format.
 *
 * It also calculates the force on the dielectric interface if required and
 * stores it in file 'force-mst.dat' or 'force-mp.dat' for MST and multipolar
 * approximations respectively.
 *
 * Works for linear interpolation in triangular elements (3-noded triangles).
 *
 * @param ANALYSIS    : [ Input ] Type of post-processing to be done
 * @param cOutputType : [ Input ] Output format type
 * @param axis        : [ Input ] Semi-axis of ellipsoidal particle
 * @param XF          : [ Input ] Points where the force should be calculated
 * @param Xinner      : [ Input ] Nodes where potential and/or field should be calculated
 * @param mNodes  : [ Input ] Coordinates of all nodes in the computational domain
 * @param mElems  : [ Input ] Connectivity of all nodes in the computational domain
 * @param vMatParam  : [ Input ] Electric properties of dielectric materials
 * @param vProbParam : [ Input ] Frequency of the applied electric field
 * @param vBCType : [ Input ] Boundary condition type for each node in the domain
 * @param vB      : [ Input ] Solution vector for the electrostatic problem
 */
int vtkPostProcess_tria3(unsigned int ANALYSIS,
                         char *cOutputType,
                         double *axis,
                         double **XF,
                         double **Xinner,
                         double **mNodes,
                         unsigned int **mElems,
                         double **vMatParam,
                         double *vProbParam,
                         unsigned int *vBCType,
                         double *vB)
{
    FILE *fp[3];
    unsigned int i, nOrder, nShape;
    double Eps, R;
    double Fcm[3], mE[6], mPot[2], Xin[3], Xcm[3], Xeval[3];
    double **Fce, **Xce;

    /* Setup the header of the output files */
    headerSetup(ANALYSIS,cOutputType,fp,Xinner);

    /* Electric potential (real part) at the required points */
    if( nInternalPoints > 0 ){
        mPot[0] = mPot[1] = 0.0;
        for(i = 0; i < 6; i++) mE[i] = 0.0;

        fprintf(fp[0],"\n");
        fprintf(fp[0],"SCALARS Re[V] double\n");
        fprintf(fp[0],"LOOKUP_TABLE default\n");
        for(i = 0; i < nInternalPoints; i++){
            Xin[0] = Xinner[i][0];
            Xin[1] = Xinner[i][1];
            Xin[2] = Xinner[i][2];

            potential_tria3(Xin,mNodes,mElems,vB,mPot);
            fprintf(fp[0],"%e\n",mPot[0]);
        }

        /* Electric potential (imaginary part) at the required points */
        fprintf(fp[0],"\n");
        fprintf(fp[0],"SCALARS Im[V] double\n");
        fprintf(fp[0],"LOOKUP_TABLE default\n");
        for(i = 0; i < nInternalPoints; i++){
            Xin[0] = Xinner[i][0];
            Xin[1] = Xinner[i][1];
            Xin[2] = Xinner[i][2];

            potential_tria3(Xin,mNodes,mElems,vB,mPot);
            fprintf(fp[0],"%e\n",mPot[1]);
        }

        /* Electric field (real part) at the required points */
        fprintf(fp[0],"\n");
        fprintf(fp[0],"VECTORS Re[E] double\n");
        for(i = 0; i < nInternalPoints; i++){
            Xin[0] = Xinner[i][0];
            Xin[1] = Xinner[i][1];
            Xin[2] = Xinner[i][2];

            field_tria3(Xin,mNodes,mElems,vB,mE);
            fprintf(fp[0],"%e\t%e\t%e\n",mE[0],mE[2],mE[4]);
        }

        /* Electric field (imaginary part) at the required points */
        fprintf(fp[0],"\n");
        fprintf(fp[0],"VECTORS Im[E] double\n");
        for(i = 0; i < nInternalPoints; i++){
            Xin[0] = Xinner[i][0];
            Xin[1] = Xinner[i][1];
            Xin[2] = Xinner[i][2];

            field_tria3(Xin,mNodes,mElems,vB,mE);
            fprintf(fp[0],"%e\t%e\t%e\n",mE[1],mE[3],mE[5]);
        }
        fclose(fp[0]);
    }


    /* Calculate Fdep at dielectric interfaces if required */
    if(ANALYSIS == 3 || ANALYSIS == 4){
        Fce = doubleMatrix(nElems,3,1);
        Xce = doubleMatrix(nElems,3,1);
        Eps = vMatParam[0][1]*eps0;
        forceMST_tria3(mNodes,mElems,vBCType,vB,Eps,Fce,Xce,Fcm,Xcm);
        fprintf(fp[2],"%e\t%e\t%e\t%e\t%e\t%e\n",Xcm[0],Xcm[1],Xcm[2],Fcm[0],Fcm[1],Fcm[2]);
        fprintf(fp[2],"#Xce\t\tYce\t\tZce\t\tFXce\t\tFYce\t\tFZce\n");

        for(i = 0; i < nElems; i++){
            if(vBCType[mElems[i][0]-1] == 6){
                fprintf(fp[2],"%e\t%e\t%e\t",Xce[i][0],Xce[i][1],Xce[i][2]);
                fprintf(fp[2],"%e\t%e\t%e\n",Fce[i][0],Fce[i][1],Fce[i][2]);
            }
        }
        fclose(fp[2]);
        freeDoubleMatrix(Fce,nElems);
        freeDoubleMatrix(Xce,nElems);
    }

    /* Calculate Fdep at particle centre using multipolar method if required */
    else if(ANALYSIS > 4){
        multipoleSetup(ANALYSIS,axis,&nShape,&nOrder,&R);
        for(i = 0; i < nFPoints; i++){
            Xeval[0] = XF[i][0];
            Xeval[1] = XF[i][1];
            Xeval[2] = XF[i][2];
            if(nShape == 0) forceMultipole_tria3(nOrder,R,Xeval,mNodes,mElems,vB,vMatParam,vProbParam,Fcm);
            else forceEllipsoid_tria3(ANALYSIS,nOrder,axis,Xeval,mNodes,mElems,vB,vMatParam,vProbParam,Fcm);
            fprintf(fp[2],"%e\t%e\t%e\t",Xeval[0],Xeval[1],Xeval[2]);
            fprintf(fp[2],"%e\t%e\t%e\n",Fcm[0],Fcm[1],Fcm[2]);
        }
        fclose(fp[2]);
    }

    return 0;
}

/**
 * Calculates potential and electric field at the required points, and stores
 * them in file 'results.vtk' using the vtk format.
 *
 * It also calculates the force on the dielectric interface if required and
 * stores it in file 'force-mst.dat' or 'force-mp.dat' for MST and multipolar
 * approximations respectively.
 *
 * Works for quadratic interpolation in triangular elements (6-noded triangles).
 *
 * @param ANALYSIS    : [ Input ] Type of post-processing to be done
 * @param cOutputType : [ Input ] Output format type
 * @param axis        : [ Input ] Semi-axis of ellipsoidal particle
 * @param XF          : [ Input ] Points where the force should be calculated
 * @param Xinner      : [ Input ] Nodes where potential and/or field should be calculated
 * @param mNodes  : [ Input ] Coordinates of all nodes in the computational domain
 * @param mElems  : [ Input ] Connectivity of all nodes in the computational domain
 * @param vMatParam  : [ Input ] Electric properties of dielectric materials
 * @param vProbParam : [ Input ] Frequency of the applied electric field
 * @param vBCType : [ Input ] Boundary condition type for each node in the domain
 * @param vB      : [ Input ] Solution vector for the electrostatic problem
 */
int vtkPostProcess_tria6(unsigned int ANALYSIS,
                         char *cOutputType,
                         double *axis,
                         double **XF,
                         double **Xinner,
                         double **mNodes,
                         unsigned int **mElems,
                         double **vMatParam,
                         double *vProbParam,
                         unsigned int *vBCType,
                         double *vB)
{
    FILE *fp[3];
    unsigned int i, nOrder, nShape;
    double Eps, R;
    double Fcm[3], mE[6], mPot[2], Xin[3], Xcm[3], Xeval[3];
    double **Fce, **Xce;

    /* Setup the header of the output files */
    headerSetup(ANALYSIS,cOutputType,fp,Xinner);

    /* Electric potential at the required points */
    if( nInternalPoints > 0 ){
        mPot[0] = mPot[1] = 0.0;
        for(i = 0; i < 6; i++) mE[i] = 0.0;

        fprintf(fp[0],"\n");
        fprintf(fp[0],"SCALARS Re[V] double\n");
        fprintf(fp[0],"LOOKUP_TABLE default\n");
        for(i = 0; i < nInternalPoints; i++){
            Xin[0] = Xinner[i][0];
            Xin[1] = Xinner[i][1];
            Xin[2] = Xinner[i][2];

            potential_tria6(Xin,mNodes,mElems,vB,mPot);
            fprintf(fp[0],"%e\n",mPot[0]);
        }

        /* Electric potential at the required points */
        fprintf(fp[0],"\n");
        fprintf(fp[0],"SCALARS Im[V] double\n");
        fprintf(fp[0],"LOOKUP_TABLE default\n");
        for(i = 0; i < nInternalPoints; i++){
            Xin[0] = Xinner[i][0];
            Xin[1] = Xinner[i][1];
            Xin[2] = Xinner[i][2];

            potential_tria6(Xin,mNodes,mElems,vB,mPot);
            fprintf(fp[0],"%e\n",mPot[1]);
        }

        /* Electric field at the required points */
        fprintf(fp[0],"\n");
        fprintf(fp[0],"VECTORS Re[E] double\n");
        for(i = 0; i < nInternalPoints; i++){
            Xin[0] = Xinner[i][0];
            Xin[1] = Xinner[i][1];
            Xin[2] = Xinner[i][2];

            field_tria6(Xin,mNodes,mElems,vB,mE);
            fprintf(fp[0],"%e\t%e\t%e\n",mE[0],mE[2],mE[4]);
        }

        /* Electric field at the required points */
        fprintf(fp[0],"\n");
        fprintf(fp[0],"VECTORS Im[E] double\n");
        for(i = 0; i < nInternalPoints; i++){
            Xin[0] = Xinner[i][0];
            Xin[1] = Xinner[i][1];
            Xin[2] = Xinner[i][2];

            field_tria6(Xin,mNodes,mElems,vB,mE);
            fprintf(fp[0],"%e\t%e\t%e\n",mE[1],mE[3],mE[5]);
        }
        fclose(fp[0]);
    }


    /* Calculate Fdep at dielectric interfaces if required */
    if(ANALYSIS == 3 || ANALYSIS == 4){
        Fce = doubleMatrix(nElems,3,1);
        Xce = doubleMatrix(nElems,3,1);
        Eps = vMatParam[0][1]*eps0;
        forceMST_tria6(mNodes,mElems,vBCType,vB,Eps,Fce,Xce,Fcm,Xcm);
        fprintf(fp[2],"%e\t%e\t%e\t%e\t%e\t%e\n",Xcm[0],Xcm[1],Xcm[2],Fcm[0],Fcm[1],Fcm[2]);
        for(i = 0; i < nElems; i++){
            if(vBCType[mElems[i][0]-1] == 6){
                fprintf(fp[2],"%e\t%e\t%e\t",Xce[i][0],Xce[i][1],Xce[i][2]);
                fprintf(fp[2],"%e\t%e\t%e\n",Fce[i][0],Fce[i][1],Fce[i][2]);
            }
        }
        fclose(fp[2]);
        freeDoubleMatrix(Fce,nElems);
        freeDoubleMatrix(Xce,nElems);
    }

    /* Calculate Fdep at particle centre using multipolar method if required */
    else if(ANALYSIS > 4){
        multipoleSetup(ANALYSIS,axis,&nShape,&nOrder,&R);
        for(i = 0; i < nFPoints; i++){
            Xeval[0] = XF[i][0];
            Xeval[1] = XF[i][1];
            Xeval[2] = XF[i][2];
            if(nShape == 0) forceMultipole_tria6(nOrder,R,Xeval,mNodes,mElems,vB,vMatParam,vProbParam,Fcm);
            else forceEllipsoid_tria6(ANALYSIS,nOrder,axis,Xeval,mNodes,mElems,vB,vMatParam,vProbParam,Fcm);
            fprintf(fp[2],"%e\t%e\t%e\t",Xeval[0],Xeval[1],Xeval[2]);
            fprintf(fp[2],"%e\t%e\t%e\n",Fcm[0],Fcm[1],Fcm[2]);
        }
        fclose(fp[2]);
    }

    return 0;
}



