/******************************************************************************
* File      : meshgen.c                                                       *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* Generates a regular 2d mesh with limits [(xmin,xmax),(ymin,ymax)] in the    *
* plane constantPlane = r, where constantPlane takes the values x, y or z.    *
*                                                                             *
* lineNumber == 0 : line number will not be included in the file              *
* lineNumber == 1 : line number will be recorded in the file                  *
*                                                                             *
* Example:    mesh 1 z 100 20 -50 50 -10 10 0.0 1                             *
* Produces a mesh in the plane z = 0.0 in the range x = (-50,50),             *
* y = (-10,10), with 100 points in the x direction and 20 in the y            *
* direction. The line number is saved to the ouput file starting with 1.      *
*                                                                             *
* The output file is called "mesh.dat" for convenience. It can be called      *
* sucesively in order to produce a single file with all the necessary points, *
* as long as "a" contains the correct number of the first element of the      *
* plane (1 in the first call, 442 in the second if n = 21, et...).            *              
******************************************************************************/

/******************************************************************************
* COPYRIGHT & LICENSE INFORMATION                                             *
*                                                                             *
* Copyright 2006 Carlos Rosales Fernandez and The Institute of High           *
* Performance Computing (A*STAR)                                              *
*                                                                             *
* meshgen is free software; you can redistribute it and/or modify             *
* it under the terms of the GNU General Public License as published by        *
* the Free Software Foundation; either version 2 of the License, or           *
* (at your option) any later version.                                         *
*                                                                             *
* meshgen is distributed in the hope that it will be useful,                  *
* but WITHOUT ANY WARRANTY; without even the implied warranty of              *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               *
* GNU General Public License for more details.                                *
*                                                                             *
* You should have received a copy of the GNU General Public License           *
* along with meshgen; if not, write to the Free Software                      *
* Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA  *
******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "errorHandler.c"

int main(int argc, char *argv[])
{
    FILE *fout;
    char constantPlane[2];
    int a, count, i, j, nx, ny, nz, n, lines;
    double dx, dy, r, x, xmin, xmax, y, ymin, ymax, z;

    /* Check for input error */
    if(argc == 2 && strcmp(argv[1],"-h") == 0){
        printf("\n\nCall as:\n\n\t./meshgen a constantPlane ");
        printf("nx ny xmin xmax ymin ymax r lineNumber\n");
        printf("\nGenerates a regular 2d mesh with limits [(xmin,xmax),");
        printf("(ymin,ymax)] in the\nplane constantPlane = r, where ");
        printf("constantPlane takes the values x, y or z.\n\n");
        printf("lineNumber == 0 : the line number is not recorded in the file\n");
        printf("lineNumber == 1 : the line number is recorded in the file\n\n");
        printf("Example:\n\n\tmesh 1 z 100 20 -50 50 -10 10 0.0 1\n\n");
        printf("Produces a mesh in the plane z = 0.0 in the range x = (-50,50),\n");
        printf("y = (-10,10), with 100 points in the x direction and 20 in the y\n");
        printf("direction. The line number is saved to the ouput file ");
        printf("starting with 1.\n\nThe output file is called 'mesh.dat' ");
        printf("for convenience. It can be called\nsucesively in order to ");
        printf("produce a single file with all the necessary points,\n");
        printf("as long as 'a' contains the correct number of the first "); 
        printf("element of the\nplane (1 in the first call, 442 in the second ");
        printf("if n = 21, et...).\n\n");
        exit(0);
    }
    else if(argc != 11){
        printf("\n\nCorrect syntaxis is:\nmeshgen a constantPlane nx ny xmin ");
        printf("xmax ymin ymax r lineNumber\n");
        errorHandler("Call 'mesh -h' for help.\n");
    }
    else{
        a = atoi(argv[1]);
        strcpy(constantPlane,argv[2]);
        nx = atoi(argv[3]);
        ny = atoi(argv[4]);
        xmin = atof(argv[5]);
        xmax= atof(argv[6]);
        ymin = atof(argv[7]);
        ymax= atof(argv[8]);
        r = atof(argv[9]);
        lines = atoi(argv[10]);
    }
    if((fout = fopen("mesh.dat","a")) == NULL)
        errorHandler("Error: Can't open output file mesh.dat");

    /* Initialize */
    constantPlane[0] = toupper(constantPlane[0]);
    dx = (xmax - xmin)/(nx - 1);
    dy = (ymax - ymin)/(ny - 1);

    count = a;
    for(i = 0; i < nx; i++){
        x = xmin + i*dx;
        for(j = 0; j < ny; j++){
            y = ymin + j*dy;
            if(lines){
                if(strcmp(constantPlane,"X") == 0) fprintf(fout,"%d %le %le %le\n",count,r,x,y);
                else if(strcmp(constantPlane,"Y") == 0) fprintf(fout,"%d    %le %le %le\n",count,x,r,y);
                else fprintf(fout,"%d   %le %le %le\n",count,x,y,r);
            }else{
                if(strcmp(constantPlane,"X") == 0) fprintf(fout,"%le    %le %le\n",r,x,y);
                else if(strcmp(constantPlane,"Y") == 0) fprintf(fout,"%le   %le %le\n",x,r,y);
                else fprintf(fout,"%le  %le %le\n",x,y,r);
            }
            count++;
        }
    }
    fclose(fout);

    return 0;
}

 
