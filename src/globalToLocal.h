/******************************************************************************
* File     : globalToLocal.h                                                  *
* Author   : Carlos Rosales Fernandez (carlos.rosales.fernandez(at)gmail.com) *
* Revision : 2.0 (2008-06-16)                                                 *
******************************************************************************/
/**
 * @brief Global to local transformation function prototype declarations
 *
 * @file
 * Global to local transformation function prototype declarations.
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

double X2L_line2(double X[2][2],
                 double *XL,
                 double *N);

inline double X2L_tria3(double X[3][3],
                        double *XL,
                        double *L,
                        double *N);

inline double X2L_tria6(double X[6][3],
                        double *XL,
                        double *L,
                        double *N);
