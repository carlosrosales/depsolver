/******************************************************************************
* File      : gaussBksb.c                                                     *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* Gauss elimination with partial pivoting (rows only) and backsubstitution.   *
* This solver is even faster than LU decomposition, but less accurate.        *
* A : Coefficient matrix                                                      *
* B : (Input) Independent Coefficients    |    (Output) Solution Vector       *
* N : Number of unknowns                                                      *
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
#include <math.h>

int gaussBksb(unsigned int N, double **A, double *B)
{
    unsigned int i, j, k, N1, i1, k1;
    double TOL, C, D;

    TOL = 1.0E-9;
    N1 = N - 1;

    for(i = 0; i < N1; i++){
        i1 = i+1;
        C = A[i][i];
        if((fabs(C)-TOL) <= 0){
            for(j = i1; j < N; j++){
                
                /* Exchange rows to get non-zero diagonal element */
                if(fabs(A[j][i])-TOL>0){
                    for(k = i; k < N; k++){ /* for ll=K:NN */
                        C = A[i][k];
                        A[i][k] = A[j][k];
                        A[j][k] = C;
                    }
                    C = B[i];
                    B[i] = B[j];
                    B[j] = C;
                    C = A[i][i];
                }
                else{    
                    printf("\nSingularity in row %d in system matrix", i);
                    exit(1);
                }
            }
        }

        /* Divide row by diagonal coefficient */
        C = A[i][i];
        for(j = i1; j < N; j++) A[i][j] = A[i][j]/C;
        B[i] = B[i]/C;

        /* Eliminate unknown X(K) from row ii */
        for(j = i1; j < N; j++){
            C = A[j][i];
            for(k = i1; k < N; k++) A[j][k] = A[j][k] - C*A[i][k];
            B[j] = B[j] - C*B[i];
        }
    }

    /* Compute last unknown */
    if(fabs(A[N1][N1]) > TOL){
        B[N1] = B[N1]/A[N1][N1];
    
        /* Apply backsubstitution process to remaining unknowns */
        for(i = 0; i < N1; i++){
            k = N1 - i - 1;
            k1 = k + 1;
            for(j = k1; j < N; j++) B[k] = B[k] - A[k][j]*B[j];
        
            /* Compute value of determinant */
            D = 1.0;
            for(j = 0; j< N; j++) D = D*A[j][j];
        }
    }
    else{
        printf("Singularity in system matrix, A[%d][%d] = %le",N1,N1,A[N1][N1]);
        exit(1);
    }
    
    return 0;
}
