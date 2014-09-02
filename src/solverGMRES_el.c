/******************************************************************************
* File      : solverGMRES_el.c                                                *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* It performs the iterative solving of the linear system using GMRES.         *
* size        : d.o.f of the problem.                                         *
* MAX_ITERS   : max. iterations in the outer loop.                            *
* tolerance   : the convergence criteria.                                     *
* x           : initial guess and also the solution returned by GMRES.        *
* convergence : converged solution (0) or not converged solution (1)          *
* preCond     : use Jacobi preconditioner (1) or not (0)                      *
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

#include "solverGMRES_el.h"

void solverGMRES_el(unsigned int preCond, unsigned int nInit, double **A, 
                    double *B)
{
    FILE *fin;
    unsigned int i, iter1 = 0, iter2 = 0, j, convergence = 1, size;
    double tolerance = 1.0, buffer, diag;
    double *x;

    fg = fopen("gmres.log","w");
    size = 2*nNodes;
    x = doubleVector(size,1);
    
    if(nInit != 0){
        fin = fopen("solution.init","r");
        for(i = 0; i < nInit; i++) fscanf(fin,"%le",&x[i]);
        fclose(fin);
    }
    
    /* Use Jacobi Preconditioner */
    if(preCond == 1){
        for(i = 0; i < nNodes; i++){
            diag = 1.0/A[i][i];
            B[i] = B[i]*diag;
            for(j = 0; j < size; j++) A[i][j] = diag*A[i][j];
        }
    }
    
    /* Set the tolerance w.r.t. the magnitute of the B */
    tolerance = L2Norm(size,B)*TOLERANCE;

    /* Printing the convergence performance of the solution progress */
    fprintf(fg,"\n\nRELATIVE TOLERANCE = %10.4e\n", tolerance);
    fprintf(fg,"\n%10s%21s\n", "ITER #    ", "RELATIVE RESIDUAL");
    fprintf(fg,"%10s%21s\n", "***************", "*******************");
        
    /* Iterations starts */
    for(i = 1; i <= MAX_OUTER_LOOP; i++){
        iter1++;
        iter2 = 0;
        
        /* Solving the linear system here */
        convergence = iterGMRES_el(size,iter1,&iter2,tolerance,A,x,B);

        if(convergence == 0){
            fprintf(fg,"\n\n***PROBLEM CONVERGES IN %d ITERATIONS***\n\n", MAX_ITERS*(iter1-1)+iter2);
    
            /* Return solution in B */    
            for(i = 0; i < size; i++) B[i] = x[i];
            free(x);
            break;
        }
    }   
        
    if(convergence == 1){
        fprintf(fg,"Problem fails to vonverge to given tolerance in %d iterations\n\n",MAX_OUTER_LOOP*MAX_ITERS);
        errorHandler("Linear system may be ill-conditioned");
    }
    fclose(fg);    
}
