/******************************************************************************
* File      : comFilter.c                                                     *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* Receives the name of the input file and filters out all C++ style comments, *
* then creates uncommented file filename.tmp. Input file is not changed.      *
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define LINELENGTH 320

void comFilter(char *FileName)
{
    FILE *fIn, *fOut;
    char cTest[LINELENGTH];
    int  i;

    if((fIn = fopen(FileName,"r")) != NULL){   
        if((fOut = fopen(strcat(FileName,".tmp"),"w")) != NULL){
            while(fgets(cTest,LINELENGTH,fIn) != NULL){
                for(i = 0; i < strlen(cTest); i++)
                    if((cTest[i] == '/') && cTest[i+1] == '/'){
                        if(i ==0) cTest[i]='\0';
                        else{
                            cTest[i]='\n';
                            cTest[i+1]='\0';
                        }
                    break;
                    }
                if(strlen(cTest) > 0) fprintf(fOut,"%s",cTest);
            }
        }
        else{
            printf("ERROR: Can't create output file %s.\n", FileName);
            exit(1);
        }
    }
    else{
        printf("ERROR: Can't open input file %s.\n", FileName);
        exit(1);
    }

    fclose(fIn);
    fclose(fOut);
}

