/******************************************************************************
* File     : depSolver.c                                                      *
* Author   : Carlos Rosales Fernandez (carlos.rosales.fernandez(at)gmail.com) *
* Revision : 2.0 (2008-06-16)                                                 *
******************************************************************************/
/**
 * @brief Driver for depSolver.
 *
 * @file
 * Main function for Laplace AC 3D multidomain electrical problem. The Gauss
 * solver requires formation of the full matrix, so the factor 2 memory saving
 * from the use of the symmetry in the matrix can't be used as in the GMRES.
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
#include <ctype.h>
#include <math.h>
#include <time.h>
#include "assembly_el.h"
#include "errorHandler.h"
#include "gaussBksb.h"
#include "gaussData.h"
#include "gmres_el.h"
#include "memHandler.h"
#include "postProcess.h"
#include "preProcess.h"

FILE *file_log;
char cElemType[6];
unsigned int nDim, nElems, nElemType, nFPoints, nInternalPoints, nInterfaces, 
             nMats, nNodes, nNodesAux, nNodesInElem, nCols, nColType, NX, NY, NZ;
double       deltaX, deltaY, deltaZ;

int main()
{
    FILE *fIn, *fAux;
    time_t currentTime;
    struct tm *localTime;
    double t1,t2,t3,t4,t_tot,t1s,t2s,t3s,t4s,t_tots;
    double axis[3], dr[3];
    int t1m, t2m, t3m, t4m, t_totm;
    unsigned int ANALYSIS, i, j, nBuffer, nInit, pNodes, preCond, size, internalPtsTemp;
    char cFilename[]="input.bem", cBuffer[33], cBCFile[33], cElemFile[33],
         cNodeFile[33], cInternalPointsFile[33], cSection[33], cTmp[33],
         solver[11], cForcePoints[33], cOutputType[4]="NAN";
    unsigned int *vBCType;
    unsigned int **vInterfaces, **mElems;
    double *vProbParam, *vB;
    double **vMatParam, **mBC, **mNodes, **Xinner, **XinnerTemp, **mA, **XF, **xCols;

    /* Open and initialize log file */
    file_log = fopen("bem.log","a");
    fprintf(file_log,"********************************************");
    fprintf(file_log,"************************************\n");
    fprintf(file_log,"depSolver Version 2.0.\n\n");
    fprintf(file_log,"Copyright 2006, 2008 Carlos Rosales Fernandez and IHPC (A*STAR).\n");
    fprintf(file_log,"License: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>.");
    fprintf(file_log,"This is free software: you are free to change and redistribute it.\n");
    fprintf(file_log,"There is NO WARRANTY, to the extent permitted by law.\n");
    fprintf(file_log,"********************************************");
    fprintf(file_log,"************************************\n\n");
    fprintf(file_log,"depSolver(): started @ ");
    currentTime = time(NULL);
    localTime = localtime(&currentTime);
    fputs(asctime(localTime),file_log);
    fclose(file_log);
    file_log = fopen("bem.log","a");

    /* Initialize pointers to arrays that may not need to be allocated */
    vInterfaces = NULL;
    xCols       = NULL;
    XF          = NULL;
    Xinner      = NULL;
    XinnerTemp  = NULL;

    /****************************************************************/
    /***     READ MAIN INPUT FILE, input.bem                      ***/
    /****************************************************************/

    /* Clean and open input file */
    comFilter(cFilename);
    fIn=fopen(cFilename,"r");

    /* Read node section */
    fscanf(fIn,"%32s %u %32s",cSection,&nNodes,cNodeFile);

    /* Read element section */
    fscanf(fIn,"%32s %u %5s %32s",cSection,&nElems,cElemType,cElemFile);
    for(i = 0; i < 5; i++) cElemType[i] = tolower(cElemType[i]);

    /* Read materials section */
    fscanf(fIn,"%32s %u",cSection,&nMats);
    vMatParam = doubleMatrix(nMats,2,0);
    for(i = 0; i < nMats; i++)
        fscanf(fIn,"%u %le %le",&nBuffer,&vMatParam[i][0],&vMatParam[i][1]);

    /* Read interfaces section */
    fscanf(fIn,"%32s %u",cSection,&nInterfaces);
    if(nInterfaces > 0){
        vInterfaces = uintMatrix(nInterfaces,2,0);
        for(i = 0; i < nInterfaces; i++)
            fscanf(fIn,"%u %u %u",&nBuffer,&vInterfaces[i][0],&vInterfaces[i][1]);
    }

    /* Read problems section */
    fscanf(fIn,"%32s",cSection);
    vProbParam = doubleVector(1,1);     /* This could be a double scalar      */
    fscanf(fIn,"%le",&vProbParam[0]);
    fscanf(fIn,"%32s",cBCFile);

    /* Get nodes in rigid body and displacement */
    fscanf(fIn,"%s",cSection);
    fscanf(fIn,"%u",&pNodes);
    fscanf(fIn,"%le %le %le",&dr[0],&dr[1],&dr[2]);

    /* Read analysis section - default solver is Gauss with backsubstitution */
    fscanf(fIn,"%32s %10s ",cSection,solver);
    for(i = 0; i < strlen(solver); i++) solver[i] = toupper(solver[i]);
    if(strcmp(solver,"GMRES") == 0) /* Get initialization method for gmres */
        fscanf(fIn,"%u %u",&preCond,&nInit);
    fscanf(fIn," %u ",&ANALYSIS);

    /* Multipolar approximations for sphere */
    if((ANALYSIS > 4) && (ANALYSIS < 10)){
        fscanf(fIn,"%u %le %32s",&nFPoints,&axis[0],cForcePoints);
        fAux = fopen(cForcePoints,"r");
        XF = doubleMatrix(nFPoints,3,0);
        for(i = 0; i < nFPoints; i++)
            fscanf(fAux,"%u %le %le %le",&nBuffer,&XF[i][0],&XF[i][1],&XF[i][2]);
        fclose(fAux);
    }

    /* Dipolar approximation for ellipsoid */
    else if(ANALYSIS == 10 || ANALYSIS == 11){
        fscanf(fIn,"%u %le %le %le %32s",&nFPoints,&axis[0],&axis[1],&axis[2],cForcePoints);
        fAux = fopen(cForcePoints,"r");
        XF = doubleMatrix(nFPoints,3,0);
        for(i = 0; i < nFPoints; i++)
            fscanf(fAux,"%u %le %le %le",&nBuffer,&XF[i][0],&XF[i][1],&XF[i][2]);
        fclose(fAux);
    }

    /* Read internal points section (points where V and E are evaluated) */
    internalPtsTemp = 0;
    fscanf(fIn,"%32s %u",cSection,&nInternalPoints);
    if(nInternalPoints > 0){
        internalPtsTemp = nInternalPoints;
        fscanf(fIn,"%3s",cOutputType);
        for(i = 0; i < strlen(cOutputType); i++) cOutputType[i] = toupper(cOutputType[i]);
        if( strncmp(cOutputType,"STD",3) == 0 || strncmp(cOutputType,"VTK",3) == 0 ){
            fscanf(fIn,"%32s",cInternalPointsFile);
            XinnerTemp = doubleMatrix(nInternalPoints,4,0);
        }
        else errorHandler("\n\nError: Unknown output type option.\n\n");
    }

    /* Read columns positions */
    fscanf(fIn,"%32s %u",cSection,&nCols);
    if( nCols > 0 ){
        fscanf(fIn,"%u",&nColType);
        xCols = doubleMatrix(nCols,3,0);
        for(i = 0; i < nCols; i++)
            fscanf(fIn,"%le %le %le",&xCols[i][0],&xCols[i][1],&xCols[i][2]);
    }
    fclose(fIn);

    /* Set element properties (nDim, nNodesInElem, nElemType) */
    nElemType = elemType();
    if( nElemType != 5 && nElemType != 6 )
        errorHandler("Error: Unknown element type option");

    /* Define the size of the problem */
    size = 2*nNodes;

    /* Memory allocation for element connectivity, nodes cordinates and bcs */
    mNodes  = doubleMatrix(nNodes,nDim,0);
    mElems  = uintMatrix(nElems,nNodesInElem,0);
    mBC     = doubleMatrix(size,3,1);
    vBCType = uintVector(size,1);


    /****************************************************************/
    /*** READ (nodes, elements, bcs, post-processing) INPUT FILES ***/
    /****************************************************************/

    /* Clean and read node file */
    comFilter(cNodeFile);
    fAux = fopen(cNodeFile,"r");
    for(i = 0; i < nNodes; i++){
        fscanf(fAux,"%u",&nBuffer);         /* Discard index */
        for(j = 0; j < nDim; j++) fscanf(fAux,"%le",&mNodes[i][j]);
    }
    fclose(fAux);
    /* Modify Rigid Body Position as required */
    for(i = pNodes; i < nNodes; i++){
        mNodes[i][0] += dr[0];
        mNodes[i][1] += dr[1];
        mNodes[i][2] += dr[2];
    }

    /* Clean and read element file */
    comFilter(cElemFile);
    fAux = fopen(cElemFile,"r");
    for(i = 0; i < nElems; i++){
        fscanf(fAux,"%u",&nBuffer);         /* Discard index */
        for(j = 0; j < nNodesInElem; j++) fscanf(fAux,"%u",&mElems[i][j]);
    }
    fclose(fAux);

    /* Clean and read boundary conditions file (bc stored as integer value) */
    comFilter(cBCFile);
    fAux = fopen(cBCFile,"r");
    for(j = 0; j < 2*nNodes; j++){
        fscanf(fAux,"%u",&nBuffer);     /* Discard index */
        fscanf(fAux,"%32s",cBuffer);
        if(strncmp(cBuffer,"C",1) == 0 || strncmp(cBuffer,"c",1) == 0){
            strncpy(cTmp,&cBuffer[1],strlen(cBuffer)-1);
            nBuffer = atoi(cTmp) + 5;
        }
        else nBuffer = atoi(cBuffer);
        vBCType[j] = nBuffer;

        /* Decide how many values to read depending on bc type */
        if((nBuffer > 2) && (nBuffer < 6))
            fscanf(fAux,"%le %le %le",&mBC[j][0],&mBC[j][1],&mBC[j][2]);
        else if(nBuffer == 0 || nBuffer == 6)
            fscanf(fAux,"%le %le",&mBC[j][0],&mBC[j][1]);
        else fscanf(fAux,"%le",&mBC[j][0]);
    }
    fclose(fAux);

    /* Clean and read internal points file */
    if(nInternalPoints > 0){
        comFilter(cInternalPointsFile);
        fAux = fopen(cInternalPointsFile,"r");
        if( strncmp(cOutputType,"VTK",3) == 0 ){
            fscanf(fAux,"%u %u %u",&NX,&NY,&NZ);
            fscanf(fAux,"%le %le %le",&deltaX,&deltaY,&deltaZ);
        }
        for(i = 0; i < nInternalPoints; i++){
            fscanf(fAux,"%u",&nBuffer);         /* Discard index */
            for(j = 0; j < nDim; j++) fscanf(fAux,"%le",&XinnerTemp[i][j]);
        }
        fclose(fAux);
    }

    /* Remove temporary files from current directory */
    system("rm -f ./*.tmp");

    /* Produce warning for specially problematic inputs */
    if ( strcmp(cOutputType,"VTK") == 0 && nCols > 0)
        warningHandler("The output options nCol > 0 and VTK may produce unexpected results");

    /* Record program status at this stage in file log */
    fprintf(file_log,"\nUsing %u nodes and %u %s ",nNodes,nElems,cElemType);
    fprintf(file_log,"elements to solve AC electrostatic problem");
    fprintf(file_log,"\nNumber of interfaces: %u\n",nInterfaces);
    fprintf(file_log,"Number of different materials: %u\n",nMats);
    fclose(file_log);
    file_log = fopen("bem.log","a");
    t1 = (double)clock()/CLOCKS_PER_SEC;


    /****************************************************************/
    /***     INPUT DATA OBTAINED, SOLUTION OF PROBLEM FOLLOWS     ***/
    /****************************************************************/

    /* Memory allocation for coefficient matrix and right hand side vector */
    if(strcmp(solver,"GMRES") == 0) mA = doubleMatrix(nNodes,size,1);
    else mA = doubleMatrix(size,size,1);
    vB = doubleVector(size,1);

    /* Assemble coefficient matrix */
    fprintf(file_log,"\ndepSolver(): performing matrix assembly ...");
    fclose(file_log);
    file_log = fopen("bem.log","a");
    if(nElemType == 5){
        if(strcmp(solver,"GMRES") == 0)
            electricFormA_tria3(mNodes,mElems,vBCType,vInterfaces,vMatParam,vProbParam,mA,mBC,vB);
        else
            electricFormFullA_tria3(mNodes,mElems,vBCType,vInterfaces,vMatParam,vProbParam,mA,mBC,vB);
    }
    else{
        if(strcmp(solver,"GMRES") == 0)
            electricFormA_tria6(mNodes,mElems,vBCType,vInterfaces,vMatParam,vProbParam,mA,mBC,vB);
        else
            electricFormFullA_tria6(mNodes,mElems,vBCType,vInterfaces,vMatParam,vProbParam,mA,mBC,vB);
    }
    t2 = (double)clock()/CLOCKS_PER_SEC-t1;

    /* Solve linear system */
    if(strcmp(solver,"GMRES") == 0){
        fprintf(file_log,"\ndepSolver(): solving linear system with GMRES ...");
        fclose(file_log);
        file_log = fopen("bem.log","a");
        solverGMRES_el(preCond,nInit,mA,vB);
    }
    else{
        fprintf(file_log,"\ndepSolver(): solving linear system with gaussBksb ...");
        fclose(file_log);
        file_log = fopen("bem.log","a");
        gaussBksb(size,mA,vB);
    }
    t3 = (float)clock()/CLOCKS_PER_SEC-t2-t1;


    /**************************************************************/
    /***      SOLUTION OBTAINED, POST-PROCESSING FOLLOWS        ***/
    /**************************************************************/

    /* Potential, electric field and force at required points */
    fprintf(file_log,"\ndepSolver(): performing postProcessing tasks ...");
    fclose(file_log);
    file_log = fopen("bem.log","a");
    Xinner = doubleMatrix(nInternalPoints,3,0);
    postProcessSetup(cOutputType,xCols,XinnerTemp,Xinner);
    if( nInternalPoints > 0 ) freeDoubleMatrix(XinnerTemp,internalPtsTemp);
    if(nElemType == 5){
        if( strcmp(cOutputType,"STD") == 0 ){
            postProcess_tria3(ANALYSIS,cOutputType,axis,XF,Xinner,mNodes,
            mElems,vMatParam,vProbParam,vBCType,vB);
        }
        else{
            vtkPostProcess_tria3(ANALYSIS,cOutputType,axis,XF,Xinner,mNodes,
            mElems,vMatParam,vProbParam,vBCType,vB);
        }
    }
    else{
        if( strcmp(cOutputType,"STD") == 0 ){
            postProcess_tria6(ANALYSIS,cOutputType,axis,XF,Xinner,mNodes,
            mElems,vMatParam,vProbParam,vBCType,vB);
        }
        else{
            vtkPostProcess_tria6(ANALYSIS,cOutputType,axis,XF,Xinner,mNodes,
            mElems,vMatParam,vProbParam,vBCType,vB);
        }
    }

    /* Save solution to file */
    fAux = fopen("solution.dat","w");
    for(i = 0; i < nNodes; i++){
        fprintf(fAux,"%le\t%le\t%le\t",mNodes[i][0],mNodes[i][1],mNodes[i][2]);
        fprintf(fAux,"%le\t%le\n",vB[i],vB[i + nNodes]);
    }
    fclose(fAux);
    t4 = (double)clock()/CLOCKS_PER_SEC-t3-t2-t1;


    /*******************************************************************/
    /***  POST-PROCESSING FINISHED, ANALYSIS OF PERFORMANCE FOLLOWS  ***/
    /*******************************************************************/

    fprintf(file_log,"\n\n*** depSolver(): performance analysis ***");
    if(t1 > 60){
        t1m = (int)(t1/60.0);
        t1s = t1-t1m*60.0;
        fprintf(file_log,"\nTIME READING INPUT: \t\t%u m %2.1f s",t1m,t1s);
    }
    else fprintf(file_log,"\nTIME READING INPUT: \t\t%2.1f seconds",t1);
    
    if(t2 > 60){
        t2m = (int)(t2/60.0);
        t2s = t2-t2m*60;
        fprintf(file_log,"\nTIME IN ASSEMBLY: \t\t%u m %2.1f s",t2m,t2s);
    }
    else fprintf(file_log,"\nTIME IN ASSEMBLY: \t\t%2.1f seconds",t2);
    
    if(t3 > 60){
        t3m = (int)(t3/60.0);
        t3s = t3-t3m*60.0;
        fprintf(file_log,"\nTIME IN SOLVER: \t\t%u m %2.1f s",t3m,t3s);
    }
    else  fprintf(file_log,"\nTIME IN SOLVER: \t\t%2.1f seconds",t3);
    
    if(t4 > 60){
        t4m = (int)(t4/60.0);
        t4s = t4-t4m*60.0;
        fprintf(file_log,"\nTIME IN POST-PROCESSING: \t%u m %2.1f s",t4m,t4s);
    }
    else fprintf(file_log,"\nTIME IN POST-PROCESSING: \t%2.1f seconds",t4);
    
    t_tot = t1+t2+t3+t4;
    if(t_tot > 60){
        t_totm = (int)(t_tot/60.0);
        t_tots = t_tot-t_totm*60.0;
        fprintf(file_log,"\nTOTAL EXECUTION TIME: \t%u m %2.1f s",t_totm,t_tots);
    }
    else fprintf(file_log,"\nTOTAL EXECUTION TIME: \t\t%2.1f seconds",t_tot);

    /************************************************************************/
    /***  PERFORMANCE ANALYSIS FINISHED, CLOSE ALL FILES AND FREE MEMORY  ***/
    /************************************************************************/

    free(vB);
    free(vBCType);
    free(vProbParam);
    freeUintMatrix(mElems,nElems);
    freeDoubleMatrix(vMatParam,nMats);
    freeDoubleMatrix(mBC,size);
    freeDoubleMatrix(mNodes,nNodes);

    /* Free only if previously allocated */
    if( ANALYSIS > 4 ) freeDoubleMatrix(XF,nFPoints);
    if( nInterfaces > 0 ) freeUintMatrix(vInterfaces,nInterfaces);
    if( nInternalPoints > 0 ) freeDoubleMatrix(Xinner,internalPtsTemp);
    if( nCols > 0 ) freeDoubleMatrix(xCols,nCols);
    if( strcmp(solver,"GMRES") == 0 ) freeDoubleMatrix(mA,nNodes);
    else freeDoubleMatrix(mA,size);

    /* Program has finished */
    fprintf(file_log,"\n\ndepSolver(): finished @ ");
    currentTime = time(NULL);
    localTime = localtime(&currentTime);
    fputs(asctime(localTime),file_log);
    fprintf(file_log,"\n\n********************************************");
    fprintf(file_log,"************************************\n\n");
    fclose(file_log);

    return 0;
}

