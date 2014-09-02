/******************************************************************************
* File     : potential.h                                                      *
* Author   : Carlos Rosales Fernandez (carlos.rosales.fernandez(at)gmail.com) *
* Revision : 2.0 (2008-06-16)                                                 *
******************************************************************************/
/**
 * @brief Prototype declarations of the potential calculation functions
 *
 * @file
 * Prototype declarations of the potential calculation functions
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

int potential_tria3(double *Xin,
                    double **mNodes,
                    unsigned int **mElems, 
                    double *vB,
                    double *mPot);

int potential_tria6(double *Xin,
                    double **mNodes,
                    unsigned int **mElems,
                    double *vB,
                    double *mPot);

