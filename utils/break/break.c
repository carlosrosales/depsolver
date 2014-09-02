/******************************************************************************
* File      : break.c                                                         *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* Breaks a data file into several. Useful when calculations for several planes*
* in Z are done and need to be stored in independent files with only the X    *
* and Y information.                                                          *
* Takes the filename to break as input 1 and the number of files to produce as*
* input 2, followed by the number of columns to read in the original file and *
* the number of rows in each of the output data files. The files are saved as *
* "data1, data2,..." and saves columns (col_1,col_2,...col_NCOLS).            *
******************************************************************************/

/******************************************************************************
* COPYRIGHT & LICENSE INFORMATION                                             *
*                                                                             *
* Copyright 2006 Carlos Rosales Fernandez and The Institute of High           *
* Performance Computing (A*STAR)                                              *
*                                                                             *
* This file is part of Break.                                                 *
*                                                                             *
* Break is free software; you can redistribute it and/or modify               *
* it under the terms of the GNU General Public License as published by        *
* the Free Software Foundation; either version 2 of the License, or           *
* (at your option) any later version.                                         *
*                                                                             *
* Break is distributed in the hope that it will be useful,                    *
* but WITHOUT ANY WARRANTY; without even the implied warranty of              *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               *
* GNU General Public License for more details.                                *
*                                                                             *
* You should have received a copy of the GNU General Public License           *
* along with Break; if not, write to the Free Software                        *
* Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA  *
******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "errorHandler.c"

int main(int argc, char *argv[])
{
FILE *fIn, *fOut, *fTmp;
char filename[32],buffer[33],name[]="data",title[33];
unsigned int i, j, k, NFILES, NCOLS, NROWS;
double data;

if(argc == 2 && strcmp(argv[1],"-h") == 0){
    printf("\n\nCall as:\n\n\t./break filename nFiles nCols nRows\n\n");
    printf("Breaks a data file into several. Useful when calculations for ");
    printf("several planes\nin Z are done and need to be stored in ");
    printf("independent files with only the X\nand Y information.\n\n");
    printf("Takes the filename to break as input 1 and the number of files ");
    printf("to produce as\ninput 2, followed by the number of columns to ");
    printf("read in the original file and\nthe number of rows in each of the ");
    printf("output data files. The files are saved as\n'data1, data2,...' and");
    printf(" saves columns (X,Y,Z,T).\n\n");
    exit(0);
}
else if(argc != 5){
    printf("\n\nCorrect symtax is: ./break filename nFiles nCols nRows\n");
    errorHandler("Type './break -h' for help.\n");
}

/* INITIALIZATION */
fIn = fopen(argv[1],"r");
NFILES = atoi(argv[2]);
NCOLS = atoi(argv[3]);
NROWS = atoi(argv[4]);

for(i = 1; i <= NFILES; i++)
{
    sprintf(buffer,"%d",i);
    strcpy(filename,name);
    strcat(filename,buffer);
    fOut = fopen(filename,"w");
    if(i>1) fTmp = fopen("data1","r");
    for(j = 0; j <= NROWS; j++){
        if(i == 1 && j == 0){
            for(k = 0; k < NCOLS; k++){
                fscanf(fIn,"%s",&title);
                fprintf(fOut,"%s    ",title);
            }
        }else if(!j){
            for(k = 0; k < NCOLS; k++){
                fscanf(fTmp,"%s",&title);
                fprintf(fOut,"%s    ",title);
            }
        }else{
            for(k = 0; k < NCOLS; k++){
                fscanf(fIn,"%le",&data);
                fprintf(fOut,"%le   ",data);
            }
        }
        fprintf(fOut,"\n");
    }
    fclose(fOut);
    if(i>1) fclose(fTmp);
}
fclose(fIn);

return 1;
}

