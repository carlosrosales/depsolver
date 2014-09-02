/******************************************************************************
* File      : uintVector.c                                                    *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* Allocates memory for an unsigned int array of size LENGTH.                  *
* INIT == 0  : malloc function is used for the allocation.                    *
* INIT != 0  : calloc function is used, all elements are initialized to zero. *
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

unsigned int *uintVector(unsigned int LENGTH, unsigned int INIT)
{
    unsigned int *V;
    
    if(!INIT){
        if((V = (unsigned int *)malloc(LENGTH*sizeof(unsigned int))) == NULL)
            errorHandler("Error - In uintVector(): Can't allocate memory");
    }
    else{
        if((V = (unsigned int *)calloc(LENGTH,sizeof(unsigned int))) == NULL)
            errorHandler("Error - In uintVector(): Can't allocate memory");
    }
    
    return V;       
}
