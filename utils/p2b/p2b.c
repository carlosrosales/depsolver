/******************************************************************************
* File      : p2b.c                                                           *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* Creates the *.bem fields from the *.out inputs exported from Patran using   *
* bem_save(). If a node has two different boundary conditions (with one bc    *
* for interface and one for source) only the bc for the source is kept.       *
* When a node is shared between different elements and appears several times  *
* in the node file, both the node file, the element connectivity and the bc   *
* file are updated so that linearly dependent equations do not result. The    *
* number of nodes, elements and boundary conditions is printed to the file    *
* "mesh-info.dat". This is the version for electrical problems.               *
*                                                                             *
* NODES       : Number of nodes in file "nodes.out"                           *
* ELEMS       : Number of elements in fiel "elems.out"                        *
* ELEM_TYPE   : Type of element (tria3, tria6, quad4 or quad8)                *
* SCALING     : Scaling factor by wich the coordinate values are multiplied   *
* BCS_NUMBER  : Number of lines in file "bcs.out"                             *
* REORDER     : Need to reorder nodes in the elements? [0 == no | 1 == yes]   *
******************************************************************************/

/******************************************************************************
* COPYRIGHT & LICENSE INFORMATION                                             *
*                                                                             *
* Copyright 2006 Carlos Rosales Fernandez and The Institute of High           *
* Performance Computing (A*STAR)                                              *
*                                                                             *
* This file is part of p2b.                                                   *
*                                                                             *
* p2b is free software; you can redistribute it and/or modify                 *
* it under the terms of the GNU General Public License as published by        *
* the Free Software Foundation; either version 2 of the License, or           *
* (at your option) any later version.                                         *
*                                                                             *
* p2b is distributed in the hope that it will be useful,                      *
* but WITHOUT ANY WARRANTY; without even the implied warranty of              *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               *
* GNU General Public License for more details.                                *
*                                                                             *
* You should have received a copy of the GNU General Public License           *
* along with p2b; if not, write to the Free Software                          *
* Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA  *
******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "errorHandler.c"
#include "doubleMatrix.c"
#include "freeDoubleMatrix.c"
#include "freeUintMatrix.c"
#include "uintMatrix.c"
#include "uintVector.c"

int main(int argc, char *argv[])
{
    FILE    *bcs_in, *bcs_out, *elems_in, *elems_out, *info_out, *nodes_in, 
            *nodes_out;
    unsigned int    BCS_NUMBER, count, COLS, ROWS, i, id, j, k, m, NODES, ELEMS, 
                    FLAG, NODES_IN_ELEM, REORDER;
    unsigned int    **elemMat, *indx, *newIndx;
    double          SCALING_FACTOR, x, y, z;
    double          a[4][4], bc[4], **M1, **M2, **nodeMat;

    /* Check number and validity of arguments */
    if(argc == 2 && strcmp(argv[1],"-h") == 0){
        printf("\nCall as:\n\tp2b NODES ELEMS ELEM_TYPE SCALING BCS_NUMBER ");
        printf("REORDER\n\nCreates the *.bem fields from the *.out inputs ");
        printf("exported from Patran using\nbem_save(). If a node is repeated");
        printf(" (with one bc for interface and one for\nsource) only the bc");
        printf(" for the source is kept. This is the version for the\n");
        printf("electric case.\n\n");
        printf("NODES\t: Number of nodes in file 'nodes.out'\n");
        printf("ELEMS\t: Number of elements in fiel 'elems.out'\n");
        printf("ELEM_TYPE\t: Type of element (tria3, tria6, quad4 or quad8)\n");
        printf("SCALING\t: Scaling factor by wich the coordinate values will ");
        printf("be multiplied\nBCS_NUMBER\t: Number of lines in file ");
        printf("'bcs.out'\nREORDER\t: Need to reorder nodes in the elements? ");
        printf("[0 == no | 1 == yes]\n");
        exit(0);
    }
    else if(argc != 7){
        printf("Correct syntaxis is:\n\n");
        printf("./p2b NODES ELEMS ELEM_TYPE SCALING_FACTOR BCS_NUMBER REORDER\n");
        errorHandler("Type './p2b -h' for help.\n\n");
    }
    else{
        NODES = atoi(argv[1]);
        ELEMS = atoi(argv[2]);
        SCALING_FACTOR = atof(argv[4]);
        BCS_NUMBER = atoi(argv[5]);
        REORDER = atoi(argv[6]);
    }
    if( (NODES < 0) || (ELEMS < 0) )
        errorHandler("Error: Number of nodes and elements must be positive!!");
        
    /* Check that files can be open */
    if((bcs_in = fopen("bcs.out","r")) == NULL)
        errorHandler("Error: Can't open input file bcs.out");
    if((elems_in = fopen("elems.out","r")) == NULL)
        errorHandler("Error: Can't open input file elems.out");
    if((nodes_in = fopen("nodes.out","r")) == NULL)
        errorHandler("Error: Can't open input file nodes.out");
    if((bcs_out = fopen("bcs.bem","w")) == NULL)
        errorHandler("Error: Can't open output file bcs.bem");
    if((elems_out = fopen("elems.bem","w")) == NULL)
        errorHandler("Error: Can't open output file elems.bem");
    if((nodes_out = fopen("nodes.bem","w")) == NULL)
        errorHandler("Error: Can't open output file nodes.bem");
    if((info_out = fopen("mesh-info.dat","w")) == NULL)
        errorHandler("Error: Can't open output file mesh-info.dat");

    /* Get number of nodes in element */ 
    if(strcmp(argv[3],"tria3") == 0) NODES_IN_ELEM = 3;
    else if(strcmp(argv[3],"tria6") == 0) NODES_IN_ELEM = 6;
    else if(strcmp(argv[3],"quad4") == 0) NODES_IN_ELEM = 4;
    else if(strcmp(argv[3],"quad8") == 0) NODES_IN_ELEM = 8;
    else {
            fclose(nodes_in);
            fclose(elems_in);
            fclose(bcs_in);
            fclose(nodes_out);
            fclose(elems_out);
            fclose(bcs_out);
            errorHandler("Error: Incorrect element type");
    }
        
    /* Memory for temporary storage */
    ROWS = 2*NODES;
    COLS = 4;
    M1 = doubleMatrix(BCS_NUMBER,COLS,1);
    M2 = doubleMatrix(ROWS,COLS,1);
    nodeMat = doubleMatrix(NODES,4,1);
    elemMat = uintMatrix(ELEMS,NODES_IN_ELEM,1);
    indx = uintVector(NODES,1);
    newIndx = uintVector(NODES,1);

    /***************************************************/
    /***         FIX REPEATED NODES PROBLEM          ***/
    /***************************************************/
    
    /* Read node and element data into temporary storage */
    for(i = 0; i < NODES; i++){
        fscanf(nodes_in,"%d %le ",&id,&nodeMat[i][0]);
        fscanf(nodes_in,"%le %le",&nodeMat[i][1],&nodeMat[i][2]);
    }
    for(i = 0; i < ELEMS; i++){
        fscanf(elems_in,"%d",&id);
        for(j = 0; j < NODES_IN_ELEM; j++) fscanf(elems_in,"%d",&elemMat[i][j]);
    }
    fclose(nodes_in);
    fclose(elems_in);
    
    /* Set remaining indexes */
    for(i = 0; i < NODES; i++) indx[i] = i;
       
    /* Loop through all nodes */
    count = 0;
    for(i = 0; i < NODES; i++){
        
        /* Set evaluation node */
        x = nodeMat[i][0];
        y = nodeMat[i][1];
        z = nodeMat[i][2];
        
        /* compare to all other nodes */
        for(j = (i + 1); j < NODES; j++){
            
            /* Check repeated appearance of this node */
            if( nodeMat[j][3] == 0){
                if( (x == nodeMat[j][0]) && (y == nodeMat[j][1]) && 
                    (z == nodeMat[j][2]) ){
                    
                    /* Eliminate the node from the matrix */
                    nodeMat[j][3] = 1;
                    count++;
                    indx[j] = i;
                }
            }
        }
    }
    
    /************************************************************/
    /***    WRITE SCALED NODE FILE WITHOUT REPEATED NODES     ***/
    /************************************************************/
    
    count = 0;
    for(i = 0; i < NODES; i++){
        if(nodeMat[i][3] == 0){
            x = nodeMat[i][0]*SCALING_FACTOR;
            y = nodeMat[i][1]*SCALING_FACTOR;
            z = nodeMat[i][2]*SCALING_FACTOR;
            fprintf(nodes_out,"%d\t%le\t%le\t%le\n",++count,x,y,z);        
            newIndx[i] = count;
        }
        else newIndx[i] = newIndx[indx[i]];
    }
    fclose(nodes_out);

    /************************************************************/
    /*** CONVERT CONNECTIVITY FILE TO ACCEPTABLE INPUT FORMAT ***/
    /************************************************************/
    
    /* Substitute repeated nodes in element connectivity file */
    for(i = 0; i < NODES; i++){
        for(j = 0; j < ELEMS; j++){
            for(k = 0; k < NODES_IN_ELEM; k++){
                if( (elemMat[j][k] - 1) == i ) elemMat[j][k] = newIndx[i];
            }
        }
    }

    if(REORDER == 1){
        for(i = 0; i < ELEMS; i++){
            if(NODES_IN_ELEM == 3){
                fprintf(elems_out,"%d\t%d\t",i+1,elemMat[i][1]);
                fprintf(elems_out,"%d\t%d\n",elemMat[i][2],elemMat[i][0]);
            }
            else if(NODES_IN_ELEM = 6){
                fprintf(elems_out,"%d\t%d\t",i+1,elemMat[i][1]);
                fprintf(elems_out,"%d\t%d\t",elemMat[i][4],elemMat[i][2]);
                fprintf(elems_out,"%d\t%d\t",elemMat[i][5],elemMat[i][0]);
                fprintf(elems_out,"%d\n",elemMat[i][3]);
            }
            else if(NODES_IN_ELEM == 4){
                fprintf(elems_out,"%d\t%d\t",i+1,elemMat[i][1]);
                fprintf(elems_out,"%d\t%d\t",elemMat[i][2],elemMat[i][3]);
                fprintf(elems_out,"%d\n",elemMat[i][0]);
            }
            else if(NODES_IN_ELEM == 8){
                fprintf(elems_out,"%d\t%d\t",i+1,elemMat[i][1]);
                fprintf(elems_out,"%d\t%d\t",elemMat[i][5],elemMat[i][2]);
                fprintf(elems_out,"%d\t%d\t",elemMat[i][6],elemMat[i][3]);
                fprintf(elems_out,"%d\t%d\t",elemMat[i][7],elemMat[i][0]);
                fprintf(elems_out,"%d\n",elemMat[i][4]);
            }
        }
    }
    else{
        for(i = 0; i < ELEMS; i++){
            fprintf(elems_out,"%d",i+1);
            for(j = 0; j < NODES_IN_ELEM; j++){
                fprintf(elems_out,"\t%d",elemMat[i][j]);
            }
            fprintf(elems_out,"\n");
        }
    }
    fclose(elems_out);
    
    /***************************************************/
    /*** CONVERT BCS FILE TO ACCEPTABLE INPUT FORMAT ***/
    /***************************************************/

    for(i = 0; i < BCS_NUMBER; i++){
        fscanf(bcs_in,"%le %le ",&M1[i][0],&M1[i][1]);
        if(M1[i][1] == 0.0) fscanf(bcs_in,"%le %le ",&M1[i][2],&M1[i][3]);
        else fscanf(bcs_in,"%le ",&M1[i][2]);
    }
    fclose(bcs_in);

    /* GET BCS FOR REAL PART */
    for(i = 0; i < NODES; i++){    
        FLAG = 0;
        for(j = 0; j < BCS_NUMBER; j++){
            if(M1[j][0] == i+1){
                if(M1[j][1] == 0.0){
                    bc[FLAG] = 0.0;
                    a[FLAG][0] = i+1;
                    a[FLAG][1] = M1[j][1];
                    a[FLAG][2] = M1[j][2];
                    a[FLAG][3] = M1[j][3];
                }else{
                    bc[FLAG] = 1.0;
                    a[FLAG][0] = i+1;
                    a[FLAG][1] = M1[j][1];
                    a[FLAG][2] = M1[j][2];
                }
                FLAG = FLAG + 1;
            }
        }
        if(FLAG == 2){
            M2[i][0] = a[0][0];
            M2[i][1] = a[0][1];
            M2[i][2] = a[0][2];

            k = i + NODES;
            M2[k][0] = a[1][0];
            M2[k][1] = a[1][1];
            M2[k][2] = a[1][2];
        
            if(bc[0] == 0.0){
                M2[i][3] = a[0][3];
                M2[k][3] = a[1][3];
            }
        }else if(FLAG == 4){
            if(bc[0] == 1.0){
                M2[i][0] = a[0][0];
                M2[i][1] = a[0][1];
                M2[i][2] = a[0][2];

                k = i + NODES;
                M2[k][0] = a[1][0];
                M2[k][1] = a[1][1];
                M2[k][2] = a[1][2];     
            }else{
                M2[i][0] = a[2][0];
                M2[i][1] = a[2][1];
                M2[i][2] = a[2][2];

                k = i + NODES;
                M2[k][0] = a[3][0];
                M2[k][1] = a[3][1];
                M2[k][2] = a[3][2];     
            }
        }
        else errorHandler("Incorrect FLAG value");
    }

    /* Write bcs output file (only for non-repeated nodes) */
    for(i = 0; i < NODES; i++){
        if(nodeMat[i][3] == 0){
            if(M2[i][1] == 0.0){
                fprintf(bcs_out,"%d\t%d \t",(int)M2[i][0],(int)M2[i][1]);
                fprintf(bcs_out,"%d \t%d\n",(int)M2[i][2],(int)M2[i][3]);
            }
            else{
                fprintf(bcs_out,"%d\t%d\t",(int)M2[i][0],(int)M2[i][1]);
                fprintf(bcs_out,"%le\n",M2[i][2]);
            }
        }
    }
    for(i = 0; i < NODES; i++){
        j = i + NODES;
        if(nodeMat[i][3] == 0){
            if(M2[j][1] == 0.0){
                fprintf(bcs_out,"%d\t%d \t",(int)M2[j][0],(int)M2[j][1]);
                fprintf(bcs_out,"%d \t%d\n",(int)M2[j][2],(int)M2[j][3]);
            }
            else{
                fprintf(bcs_out,"%d\t%d\t",(int)M2[j][0],(int)M2[j][1]);
                fprintf(bcs_out,"%le\n",M2[j][2]);
            }
        }
    }
    fclose(bcs_out);
    
    fprintf(info_out,"Nodes :\t%d\nElems :\t%d\n",count,ELEMS);
    fprintf(info_out,"Nodes in each element :\t%d\n",NODES_IN_ELEM);
    fprintf(info_out,"Repeated nodes eliminated :\t%d\n",NODES - count);
    fclose(info_out);
    
    freeDoubleMatrix(M1,ROWS);
    freeDoubleMatrix(M2,ROWS);
    freeDoubleMatrix(nodeMat,NODES);
    freeUintMatrix(elemMat,ELEMS);
    free(indx);
    free(newIndx);

    return 1;
}
