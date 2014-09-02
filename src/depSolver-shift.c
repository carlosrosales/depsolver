/******************************************************************************
* File      : depSolver-shift.c                                               *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* Main function for Laplace AC 3D multidomain electrical problem. Single      *
* precission complex version when using GMRES. the Gauss solver requires      *
* formation of the full matrix, so the factor 2 memory saving from the use of *
* the symmetry in the matrix can't be used easily.                            *
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

#include "gaussData.h"
#include "depSolver.h"

int main()
{
    FILE *fIn, *fAux;
    time_t currentTime;
    struct tm *localTime;
    char cFilename[]="input.bem", cBuffer[32], cBCFile[32], cElemFile[32], 
          cNodeFile[32], cInternalPointsFile[32], cSection[20], cTmp[32], 
          solver[10];
    int t1m, t2m, t3m, t4m, t_totm;
    unsigned int ANALYSIS, i, j, nBuffer, nInit, pNodes, preCond, size;
    unsigned int *vBCType, *indx;
    unsigned int **vInterfaces, **mElems;
    double exchanges, t1,t2,t3,t4,t_tot,t1s,t2s,t3s,t4s,t_tots;
    double axis[3], dr[3];
    double *vProbParam, *vB;
    double **vMatParam, **mBC, **mNodes, **Xinner, **mA, **XF, **xCols;

    /* Open and initialize log file */
    file_log = fopen("bem.log","a");
    fprintf(file_log,"depSolver-sft(): started @ ");
    currentTime = time(NULL);
    localTime = localtime(&currentTime);
    fputs(asctime(localTime),file_log);

    /* Clean and open input file */
    comFilter(cFilename);
    fIn=fopen(cFilename,"r");

    /* Read node section */
    fscanf(fIn,"%s %d %s",cSection,&nNodes,cNodeFile);

    /* Read element section */
    fscanf(fIn,"%s %d %s %s",cSection,&nElems,cElemType,cElemFile);
    for(i = 0; i < 6; i++) cElemType[i] = tolower(cElemType[i]);

    /* Read materials section */
    fscanf(fIn,"%s %d",cSection,&nMats);
    vMatParam = doubleMatrix(nMats,2,0);
    for(i = 0; i < nMats; i++)
        fscanf(fIn,"%d %le %le",&nBuffer,&vMatParam[i][0],&vMatParam[i][1]);

    /* Read interfaces section */
    fscanf(fIn,"%s %d",cSection,&nInterfaces);
    if(nInterfaces > 0){
        vInterfaces = uintMatrix(nInterfaces,2,0);
        for(i = 0; i < nInterfaces; i++)
            fscanf(fIn,"%d %d %d",&nBuffer,&vInterfaces[i][0],&vInterfaces[i][1]);
    }

    /* Read problems section */
    fscanf(fIn,"%s",cSection);
    size = 2*nNodes;                    /* Characteristic size of the problem */
    vProbParam = doubleVector(1,1);     /* This could be a double scalar      */
    vBCType = uintVector(size,1);
    mBC = doubleMatrix(size,3,1);

    /* Read problem parameters */
    fscanf(fIn,"%le",&vProbParam[0]);
    fscanf(fIn,"%s",cBCFile);

    /* Clean and read boundary conditions file (bc stored as integer value) */
    comFilter(cBCFile);
    fAux = fopen(cBCFile,"r");
    for(j = 0; j < size; j++){
        fscanf(fAux,"%d",&nBuffer);     /* Discard index */
        fscanf(fAux,"%s",cBuffer);
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

    /* Set element properties (nDim, nNodesInElem, nElemType) */
    nElemType = elemType();
    if(nElemType != 5 && nElemType != 6)
        errorHandler("ERROR: Incorrect Element Type");
        
    /* Get nodes in rigid body and displacement */
    fscanf(fIn,"%s",cSection);
    fscanf(fIn,"%d",&pNodes);
    fscanf(fIn,"%le %le %le",&dr[0],&dr[1],&dr[2]);

    /* Read analysis section - default solver if Gauss-Jordan */
    fscanf(fIn,"%s %s ",cBuffer,solver);
    for(i = 0; i < strlen(solver); i++) solver[i] = toupper(solver[i]);
    if(strcmp(solver,"GMRES") == 0) /* Get initialization method for gmres */
        fscanf(fIn,"%d %d",&preCond,&nInit);
    fscanf(fIn," %d ",&ANALYSIS);
    
    /* Multipolar approximations for sphere */
    if((ANALYSIS > 4) && (ANALYSIS < 10)){
        fscanf(fIn,"%d %le",&nFPoints,&axis[0]);
        XF = doubleMatrix(nFPoints,3,0);
        for(i = 0; i < nFPoints; i++)
            fscanf(fIn,"%le %le %le",&XF[i][0],&XF[i][1],&XF[i][2]);
    }
    
    /* Dipolar approximation for ellipsoid */
    else if(ANALYSIS == 10 || ANALYSIS == 11){
        fscanf(fIn,"%d %le %le %le",&nFPoints,&axis[0],&axis[1],&axis[2]);
        XF = doubleMatrix(nFPoints,3,0);
        for(i = 0; i < nFPoints; i++)
            fscanf(fIn,"%le %le %le",&XF[i][0],&XF[i][1],&XF[i][2]);
    }

    /* Read internal points section (points where V and E are evaluated) */
    fscanf(fIn,"%s %d",cBuffer,&nInternalPoints);
    if(nInternalPoints > 0){
        fscanf(fIn,"%s",cInternalPointsFile);
        Xinner = doubleMatrix(nInternalPoints,nDim,0);
    }

    /* Read columns positions */
    fscanf(fIn,"%s %d",cBuffer,&nCols);
    if(nCols){
        fscanf(fIn,"%d",&nColType);
        xCols = doubleMatrix(nCols,3,0);
        for(i = 0; i < nCols; i++)
            fscanf(fIn,"%le %le %le",&xCols[i][0],&xCols[i][1],&xCols[i][2]);
    }
    fclose(fIn);    
    
    /* Memory allocation for element connectivity and nodes cordinates */
    mNodes = doubleMatrix(nNodes,nDim,0);
    mElems = uintMatrix(nElems,nNodesInElem,0);

    /* Clean and read node file */
    comFilter(cNodeFile);
    fAux = fopen(cNodeFile,"r");
    for(i = 0; i < nNodes; i++){
        fscanf(fAux,"%d",&nBuffer);         /* Discard index */
        for(j = 0; j < nDim; j++) fscanf(fAux,"%le",&mNodes[i][j]);
    }
    fclose(fAux);

    /* Clean and read element file */
    comFilter(cElemFile);
    fAux = fopen(cElemFile,"r");
    for(i = 0; i < nElems; i++){
        fscanf(fAux,"%d",&nBuffer);         /* Discard index */
        for(j = 0; j < nNodesInElem; j++) fscanf(fAux,"%d",&mElems[i][j]);
    }
    fclose(fAux);

    /* Clean and read internal points file */
    if(nInternalPoints > 0){
        comFilter(cInternalPointsFile);
        fAux = fopen(cInternalPointsFile,"r");
        for(i = 0; i < nInternalPoints; i++){
            fscanf(fAux,"%d",&nBuffer);        /* Discard index */
            for(j = 0; j < nDim; j++) fscanf(fAux,"%le",&Xinner[i][j]);
        }
        fclose(fAux);
    }

    /* Record program status at this stage in file log */
    fprintf(file_log,"\nUsing %d Nodes and %d %s ",nNodes,nElems,cElemType);
    fprintf(file_log,"Elems to solve AC Electrostatic Problem");
    fprintf(file_log,"\nNumber of interfaces: %d\n",nInterfaces);
    fprintf(file_log,"Number of different materials: %d\n",nMats);
    t1 = (double)clock()/CLOCKS_PER_SEC;

    /****************************************************************/
    /***     INPUT DATA OBTAINED, SOLUTION OF PROBLEM FOLLOWS     ***/
    /****************************************************************/
    
    /* Modify Rigid Body Position as required */
    for(i = pNodes; i < nNodes; i++){
        mNodes[i][0] += dr[0];
        mNodes[i][1] += dr[1];
        mNodes[i][2] += dr[2];
    }

    /* Memory allocation for coefficient matrix and right hand side vector */
    if(strcmp(solver,"GMRES") == 0) mA = doubleMatrix(nNodes,size,1);
    else mA = doubleMatrix(size,size,1);
    vB = doubleVector(size,1);

    /* Assemble coefficient matrix */
    if(nElemType == 5){
        fprintf(file_log,"\ndepSolver-sft(): calling formMatrix_tria3() ...");
        if(strcmp(solver,"GMRES") == 0)
            electricFormA_tria3(mNodes,mElems,vBCType,vInterfaces,vMatParam,vProbParam,mA,mBC,vB);
        else
            electricFormFullA_tria3(mNodes,mElems,vBCType,vInterfaces,vMatParam,vProbParam,mA,mBC,vB);
    }
    else{
        fprintf(file_log,"\ndepSolver-sft(): calling formMatrix_tria6() ...");
        if(strcmp(solver,"GMRES") == 0)
            electricFormA_tria6(mNodes,mElems,vBCType,vInterfaces,vMatParam,vProbParam,mA,mBC,vB);
        else
            electricFormFullA_tria6(mNodes,mElems,vBCType,vInterfaces,vMatParam,vProbParam,mA,mBC,vB);
    }
    t2 = (double)clock()/CLOCKS_PER_SEC-t1;

    /* Solve linear system */
    if(strcmp(solver,"GMRES") == 0){
        fprintf(file_log,"\ndepSolver-sft(): calling solverGMRES() ...");
        solverGMRES_el(preCond,nInit,mA,vB);
    }
    else{
        fprintf(file_log,"\ndepSolver-sft(): calling gaussBksb() ...");
        gaussBksb(size,mA,vB);
    }
    t3 = (float)clock()/CLOCKS_PER_SEC-t2-t1;

    /**************************************************************/
    /***      SOLUTION OBTAINED, POST-PROCESSING FOLLOWS        ***/
    /**************************************************************/

    /* Potential, electric field and force at required points */
    if(nElemType == 5){
        fprintf(file_log,"\ndepSolver-sft(): calling postProcess_tria3() ...");
        if(nCols > 0){
            postProcessCol_tria3(ANALYSIS,axis,XF,Xinner,mNodes,mElems,
            vMatParam,vProbParam,vBCType,xCols,vB);
        }
        else{
            postProcess_tria3(ANALYSIS,axis,XF,Xinner,mNodes,mElems,vMatParam,
            vProbParam,vBCType,vB);
        }
    }
    else{
        fprintf(file_log,"\ndepSolver-sft(): calling postProcess_tria6() ...");
        if(nCols > 0){
            postProcessCol_tria6(ANALYSIS,axis,XF,Xinner,mNodes,mElems,
            vMatParam,vProbParam,vBCType,xCols,vB);
        }
        else{
            postProcess_tria6(ANALYSIS,axis,XF,Xinner,mNodes,mElems,vMatParam,
            vProbParam,vBCType,vB);
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

    fprintf(file_log,"\n\n*** depSolver-sft(): performance analysis ***");
    if(t1 > 60){
        t1m = (int)(t1/60.0);
        t1s = t1-t1m*60.0;
        fprintf(file_log,"\nTIME READING INPUT: \t\t%d m %2.1f s",t1m,t1s);
    }
    else fprintf(file_log,"\nTIME READING INPUT: \t\t%2.1f seconds",t1);
    
    if(t2 > 60){
        t2m = (int)(t2/60.0);
        t2s = t2-t2m*60;
        fprintf(file_log,"\nTIME IN ASSEMBLY: \t\t%d m %2.1f s",t2m,t2s);
    }
    else fprintf(file_log,"\nTIME IN ASSEMBLY: \t\t%2.1f seconds",t2);
    
    if(t3 > 60){
        t3m = (int)(t3/60.0);
        t3s = t3-t3m*60.0;
        fprintf(file_log,"\nTIME IN SOLVER: \t\t%d m %2.1f s",t3m,t3s);
    }
    else  fprintf(file_log,"\nTIME IN SOLVER: \t\t%2.1f seconds",t3);
    
    if(t4 > 60){
        t4m = (int)(t4/60.0);
        t4s = t4-t4m*60.0;
        fprintf(file_log,"\nTIME IN POST-PROCESSING: \t%d m %2.1f s",t4m,t4s);
    }
    else fprintf(file_log,"\nTIME IN POST-PROCESSING: \t%2.1f seconds",t4);
    
    t_tot = t1+t2+t3+t4;
    if(t_tot > 60){
        t_totm = (int)(t_tot/60.0);
        t_tots = t_tot-t_totm*60.0;
        fprintf(file_log,"\nTOTAL EXECUTION TIME: \t%d m %2.1f s",t_totm,t_tots);
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
    if(ANALYSIS > 4) freeDoubleMatrix(XF,nFPoints);
    if(nInterfaces) freeUintMatrix(vInterfaces,nInterfaces);
    if(nCols) freeDoubleMatrix(xCols,nCols);
    if(strcmp(solver,"GMRES") == 0) freeDoubleMatrix(mA,nNodes);
    else freeDoubleMatrix(mA,size);

    /* Program has finished */
    fprintf(file_log,"\n\ndepSolver-sft(): finshed @ ");
    currentTime = time(NULL);
    localTime = localtime(&currentTime);
    fputs(asctime(localTime),file_log);
    fprintf(file_log,"\n\n********************************************");
    fprintf(file_log,"************************************\n\n");
    fclose(file_log);

    return 0;
}
