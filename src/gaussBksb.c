/******************************************************************************
* File     : gaussBksb.c                                                      *
* Author   : Carlos Rosales Fernandez (carlos.rosales.fernandez(at)gmail.com) *
* Revision : 2.0 (2008-06-16)                                                 *
******************************************************************************/
/**
 * @brief Gauss elimination solver
 *
 * @file
 * Gauss elimination with partial pivoting (rows only) and backsubstitution. The
 * input constant vector \a B is overwritten with the solution at the output.
 *
 * @param A : [ Input ] Coefficient matrix
 * @param B : [ Input ] Independent coefficients; [ Output ] Solution vector
 * @param N : [ Input ] Number of unknowns
 */

/*******************************************************************************
* Copyright 2006, 2008 Carlos Rosales Fernandez and IHPC (A*STAR).             *
*                                                                              *
* This file is part of depSolver.                                              *
*                                                                              *
* depSolver is free software: you can redistribute it and/or modify it under   *
* the terms of the GNU GPL version 3 or (at your option) any later version.    *
*                                                                              *
* depSolver is distributed in the hope that it will be useful, but WITHOUT ANY *
* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS    *
* FOR A PARTICULAR PURPOSE. See the GNU General Public License for details.    *
*                                                                              *
* You should have received a copy of the GNU General Public License along with *
* depSolver. If not, see <http://gnu.org/licenses/gpl.html>.                   *
*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gaussBksb.h"

int gaussBksb(unsigned int N,
              double **A,
              double *B)
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
        printf("Singularity in system matrix, A[%d][%d] = %e",N1,N1,A[N1][N1]);
        exit(1);
    }
    
    return 0;
}
