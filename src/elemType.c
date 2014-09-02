/******************************************************************************
* File      : elemType.c                                                      *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* Returns an index between 1 and 9 indicating the element type. If 0 is       *
* returned the element type is invalid. It also provides the number of nodes  *
* in each element, nNodesInElem, and the dimension of the problem, nDim.      *
* Recognised element types are line1, line2, line3, tria1, tria3, tria6,      *
* quad1, quad4 and quad8.                                                     *
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

extern char cElemType[];
extern unsigned int nDim, nNodesInElem; 

int elemType()
{
    /* Transform to lowercase for comparison */
    cElemType[0] = tolower(cElemType[0]);
    cElemType[1] = tolower(cElemType[1]);
    cElemType[2] = tolower(cElemType[2]);
    cElemType[3] = tolower(cElemType[3]);
    cElemType[4] = tolower(cElemType[4]);
    cElemType[5] = tolower(cElemType[5]);

    /* Find element characteristic from input */
    if(strcmp(cElemType,"line1") == 0){
        nDim = 2;
        nNodesInElem = 1;
        return 1;
    }
    else if(strcmp(cElemType,"line2") == 0){
        nDim = 2;
        nNodesInElem = 2;
        return 2;
    }
    else if(strcmp(cElemType,"line3") == 0){
        nDim = 2;
        nNodesInElem = 3;
        return 3;
    }
    else if(strcmp(cElemType,"tria1") == 0){
        nDim = 3;
        nNodesInElem = 1;
        return 4;
    }
    else if(strcmp(cElemType,"tria3") == 0){
        nDim = 3;
        nNodesInElem = 3;
        return 5;
    }
    else if(strcmp(cElemType,"tria6") == 0){
        nDim = 3;
        nNodesInElem = 6;
        return 6;
    }
    else if(strcmp(cElemType,"quad1") == 0){
        nDim = 3;
        nNodesInElem = 1;
        return 7;
    }
    else if(strcmp(cElemType,"quad4") == 0){
        nDim = 3;
        nNodesInElem = 4;
        return 8;
    }
    else if(strcmp(cElemType,"quad8") == 0){
        nDim = 3;
        nNodesInElem = 8;
        return 9;
    }

    /* This indicates an invalid element type */
    else return 0;
}
