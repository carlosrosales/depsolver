/******************************************************************************
* File      : iterGMRES_el.c                                                  *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* Iterative solution of the linear system Ax=B using GMRES.                   *
* c : cosine in the Givens Rotation                                           *
* s : sine in the Givens Rotation                                             *
* h : Hessenberg matrix                                                       *
* converge: (0) problem converges, (1) problem fails to converge within       *
*           MAX_ITERS, (2) problem converging too slowly                      *
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

#include "iterGMRES_el.h"

int iterGMRES_el(unsigned int size, unsigned int iter1, unsigned int *iter2, 
                double tolerance, double **A, double *x, double *B)
{	
	int i, j, k, iters = 0, converge = 1;
	double norm, htmp, h1, h2; 
	double error[2];
	double *res, *c, *s;
	double **h, **krylov;

  	s = doubleVector(MAX_ITERS,1);
  	c = doubleVector(MAX_ITERS,1);
  	res = doubleVector(MAX_ITERS,1);
  	krylov = doublePointer(MAX_ITERS,1);
  	h = doubleMatrix(MAX_ITERS,MAX_ITERS,1);
  	    	
	/* Determine the first Krylov subspace vector */
	krylov[0] = doubleVector(size,0);              
    initRes_el(size,A,x,B,krylov[0]);	
    norm = L2Norm(size,krylov[0]);
    for(i = 0; i < size; i++) krylov[0][i] = krylov[0][i]/norm;	
  	res[0] = norm;
  	error[0] = norm;
    
    /* If initial guess is accurate exit */
  	if(error[0] < tolerance) return 0;
	
	/* GMRES iterations starts */	    
	for(i = 0; i < MAX_ITERS; i++){
		iters++;
		
		/* Determine Krylov[i+1] vector using Krylov[i] vector */
		krylov[i+1] = doubleVector(size,0);
		matVectProd_el(A,krylov[i],krylov[i+1]);
		norm = L2Norm(size,krylov[i+1]);

		/* orthogonalization of krylov[i+1] w.r.t. all previous krylov[i] */
		for(j = 0; j <= i; j++){
			h[j][i] = dotProd(size,krylov[j],krylov[i+1]);
            for(k = 0; k < size; k++)
                krylov[i+1][k] = krylov[i+1][k] - h[j][i]*krylov[j][k];
		}
	
		h[i+1][i] = L2Norm(size,krylov[i+1]);	
		if(h[i+1][i] != 0.0)
            for(j = 0; j < size; j++) krylov[i+1][j] = krylov[i+1][j]/h[i+1][i];

		/* Compute the Givens rotation that reduces the    */ 
    	/* Hessenberg matrix to an upper triangular matrix */
		if(i > 0){
            for(j = 0; j <= (i-1); j++){
                h1 = c[j]*h[j][i] + s[j]*h[j+1][i];
                h2 = -s[j]*h[j][i] + c[j]*h[j+1][i];
        
                h[j][i] = h1;
                h[j+1][i] = h2;
            }
        }
		
		/* Perform the current Givens rotation */
        htmp = h[i+1][i]/h[i][i];
		if(norm != 0.0){
			c[i] = 1.0/sqrt(1.0 + htmp*htmp);
			s[i] = c[i]*htmp;
	    		
			h[i][i] = c[i]*h[i][i] + s[i]*h[i+1][i];
			h[i+1][i] = 0.0;

			res[i+1] = -s[i]*res[i];
			res[i] = c[i]*res[i];
		}   
    
		/* Update the error and test for convergence */ 
		error[1] = fabs(res[i+1]);
		fprintf(fg,"%8d%25.4e\n", MAX_ITERS*(iter1-1)+iters, error[1]);
		
		if(error[1] < tolerance){
			converge = 0;
			break;
		}
		else if(fabs((error[1] - error[0])/error[0]) < 1.0e-6*tolerance)
			errorHandler("CONVERGENCE RATE TOO SLOW");
		else error[0] = error[1];	
	}
    
	/* Solve the triangle system of equations for the weights */
  	for(i = iters-1; i >= 0; i--){
		for(j = i+1; j <= iters; j++) res[i] -= h[i][j]*res[j];		
  		res[i] = res[i]/h[i][i];
	}

	/* Return the solution in x */
	for(i = 0; i < size; i++)
		for(j = 0; j < iters; j++) x[i] += res[j]*krylov[j][i];  
			
	*iter2 = iters;
    
	/* Free dynamically allocated memory */
	free(s);
  	free(c);
  	free(res);
  	freeDoubleMatrix(h,MAX_ITERS);
  	freeDoublePointer(krylov,iters + 1);
	
	return converge;
}
