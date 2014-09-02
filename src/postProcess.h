/******************************************************************************
* File     : postProcess.h                                                    *
* Author   : Carlos Rosales Fernandez (carlos.rosales.fernandez(at)gmail.com) *
* Revision : 2.0 (2008-06-16)                                                 *
******************************************************************************/
/**
 * @brief Prototype declarations of post-processing functions
 *
 * @file
 * Prototype declarations of post-processing functions
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

int headerSetup(unsigned int ANALYSIS,
                char *cOutputType,
                FILE **fp,
                double **Xinner);

int postProcessSetup(char *cOutputType,
                     double **xCols,
                     double **XinnerTemp,
                     double **Xinner);

int multipoleSetup(unsigned int ANALYSIS,
                   double *axis,
                   unsigned int *nShape,
                   unsigned int *nOrder,
                   double *R);

int postProcess_tria3(unsigned int ANALYSIS,
                      char *cOutputType,
                      double *axis,
                      double **XF,
                      double **Xinner,
                      double **mNodes,
                      unsigned int **mElems,
                      double **vMatParam,
                      double *vProbParam,
                      unsigned int *vBCType,
                      double *vB);

int postProcess_tria6(unsigned int ANALYSIS,
                      char *cOutputType,
                      double *axis,
                      double **XF,
                      double **Xinner,
                      double **mNodes,
                      unsigned int **mElems,
                      double **vMatParam,
                      double *vProbParam,
                      unsigned int *vBCType,
                      double *vB);

int vtkPostProcess_tria3(unsigned int ANALYSIS,
                         char *cOutputType,
                         double *axis,
                         double **XF,
                         double **Xinner,
                         double **mNodes,
                         unsigned int **mElems,
                         double **vMatParam,
                         double *vProbParam,
                         unsigned int *vBCType,
                         double *vB);

int vtkPostProcess_tria6(unsigned int ANALYSIS,
                         char *cOutputType,
                         double *axis,
                         double **XF,
                         double **Xinner,
                         double **mNodes,
                         unsigned int **mElems,
                         double **vMatParam,
                         double *vProbParam,
                         unsigned int *vBCType,
                         double *vB);

