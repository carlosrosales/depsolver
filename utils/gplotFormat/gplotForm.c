/******************************************************************************
* File      : gplotForm.c                                                     *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* Inserts a blank row every nrows for all ncols in filename for all N rows in *
* the file. Input file is not altered in any way. Output file is by default   *
* called "temp.gnu".                                                          *              
******************************************************************************/

/******************************************************************************
* COPYRIGHT & LICENSE INFORMATION                                             *
*                                                                             *
* Copyright 2006 Carlos Rosales Fernandez and The Institute of High           *
* Performance Computing (A*STAR)                                              *
*                                                                             *
* gplotForm is free software; you can redistribute it and/or modify           *
* it under the terms of the GNU General Public License as published by        *
* the Free Software Foundation; either version 2 of the License, or           *
* (at your option) any later version.                                         *
*                                                                             *
* gplotForm is distributed in the hope that it will be useful,                *
* but WITHOUT ANY WARRANTY; without even the implied warranty of              *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               *
* GNU General Public License for more details.                                *
*                                                                             *
* You should have received a copy of the GNU General Public License           *
* along with gplotForm; if not, write to the Free Software                    *
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
    char cBuffer[32];
    unsigned int count, i, j, k, N, nrows, ncols;
    double data;

    if(argc == 2 && strcmp(argv[1],"-h") == 0){
        printf("\n\nCall as:\n\n\t./gplotForm filename N nrows ncols.\n");
        printf("\nInserts a blank row every nrows for all ncols in filename");
        printf(" for all N rows in\nthe file. Input file is not altered in ");
        printf("any way. Output file is by default\ncalled 'temp.gnu'.\n\n"); 
        exit(0);
    }
    else if(argc != 5){
        printf("\n\nCorrect syntax is:\n\n\t./gplotForm filename N nrows ncols.\n");
        errorHandler("Type './sep -h' for help.\n");
    }

    /* Initialize */
    fIn = fopen(argv[1],"r");
    N = atoi(argv[2]);
    nrows = atoi(argv[3]);
    ncols = atoi(argv[4]);
    for(i = 0; i < ncols; i++) fscanf(fIn,"%s",cBuffer);

	fOut = fopen("temp.gnu","w");
    count = 0;
	for(j = 0; j < N; j++){
		for(k = 0; k < ncols; k++){
    		fscanf(fIn,"%le",&data);
			fprintf(fOut,"%le\t",data);
		}
		fprintf(fOut,"\n");
		count++;
		if(count == nrows){
		  count = 0;
		  fprintf(fOut,"\n");
		}	  
	}
    fclose(fOut);
    fclose(fIn);

return 0;
}

