/******************************************************************************
* File     : force.h                                                          *
* Author   : Carlos Rosales Fernandez (carlos.rosales.fernandez(at)gmail.com) *
* Revision : 2.0 (2008-06-16)                                                 *
******************************************************************************/
/**
 * @brief Dielectrophoretic (DEP) force calculation function prototypes
 *
 * @file
 * Prototype declarations of the six functions responsible for the calculation
 * of dielectrophoretic (DEP) forces in depSolver: forceEllipsoid_tria3,
 * forceEllipsoid_tria6, forceMST_tria3, forceMST_tria3, forceMultipole_tria3,
 * forceMultipole_tria6.
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

int forceEllipsoid_tria3(unsigned int ANALYSIS,
                         unsigned int nOrder,
                         double *axis,
                         double *Xeval,
                         double **mNodes,
                         unsigned int **mElems,
                         double *vB,
                         double **vMatParam,
                         double *vProbParam,
                         double *F);

int forceEllipsoid_tria6(unsigned int ANALYSIS,
                         unsigned int nOrder,
                         double *axis,
                         double *Xeval,
                         double **mNodes,
                         unsigned int **mElems,
                         double *vB,
                         double **vMatParam,
                         double *vProbParam,
                         double *F);

int forceMST_tria3(double **mNodes,
                   unsigned int **mElems,
                   unsigned int *vBCType,
                   double *vB,
                   double Eps,
                   double **Fce,
                   double **Xce,
                   double *Fcm,
                   double *Xcm);

int forceMST_tria6(double **mNodes,
                   unsigned int **mElems,
                   unsigned int *vBCType,
                   double *vB,
                   double Eps,
                   double **Fce,
                   double **Xce,
                   double *Fcm,
                   double *Xcm);

int forceMultipole_tria3(unsigned int nOrder,
                         double R,
                         double *Xeval,
                         double **mNodes,
                         unsigned int **mElems,
                         double *vB,
                         double **vMatParam,
                         double *vProbParam,
                         double *F);

int forceMultipole_tria6(unsigned int nOrder,
                         double R,
                         double *Xeval,
                         double **mNodes,
                         unsigned int **mElems,
                         double *vB,
                         double **vMatParam,
                         double *vProbParam,
                         double *F);

