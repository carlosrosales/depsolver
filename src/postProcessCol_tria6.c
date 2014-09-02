/******************************************************************************
* File      : postProcess_tria6.c                                             *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* Calculates potential and electric field at the required points, and stores  *
* them in files 'potential.dat' and 'field.dat'. It also calculates the force *
* on the dielectric interface if required and stores it in files              *
* "force-mst.dat" and "force-mp.dat" for MST and multipolar approximations    *
* respectively. Avoids calculations inside columns when present.              *
* Works for quadratic interpolation in triangular elements (6-noded triangles)*
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

#include "postProcessCol_tria6.h"

int postProcessCol_tria6(unsigned int ANALYSIS, double *axis, double **XF, 
                         double **Xinner, double **mNodes, 
                         unsigned int **mElems, double **vMatParam, 
                         double *vProbParam, unsigned int *vBCType, 
                         double **xCols, double *vB)
{
    FILE *fE,*fF,*fP;
    unsigned int count, i, j, nOrder, nPointInside, nShape;
    double Eps, R;
    double Fcm[3], mPot[2], mE[6], Xin[3], Xcm[3], Xeval[3];
    double **Fce, **Xce;

    /* Initialize */
    Fce = doubleMatrix(nElems,3,1);
    Xce = doubleMatrix(nElems,3,1);
    switch(ANALYSIS){
        case(0):
            fprintf(file_log,"\n\tpostProcessCol_tria6(): ");
            fprintf(file_log,"doing potential calculation...");
            fP = fopen("potential.dat","w");
            fprintf(fP,"#x\t\ty\t\tz\t\tRe[Pot(x,y,z)]\tIm[Pot(x,y,z)]\n");
            break;
        case(1):
            fprintf(file_log,"\n\tpostProcessCol_tria6(): ");
            fprintf(file_log,"doing field calculation...");
            fE = fopen("field.dat","w");
            fprintf(fE,"#x\t\ty\t\tz\t\tRe[Ex(x,y,z)]\tIm[Ex(x,y,z)]\tRe[Ey(x,y,z)]");
            fprintf(fE,"\tIm[Ey(x,y,z)]\tRe[Ez(x,y,z)]\tIm[Ez(x,y,z)]\n");
            break;
        default:
            fprintf(file_log,"\n\tpostProcessCol_tria6(): ");
            fprintf(file_log,"doing potential and field calculation...");
            fE = fopen("field.dat","w");
            fP = fopen("potential.dat","w");
            fprintf(fP,"#x\t\ty\t\tz\t\tRe[Pot(x,y,z)]\tIm[Pot(x,y,z)]\n");
            fprintf(fE,"#x\t\ty\t\tz\t\tRe[Ex(x,y,z)] \tIm[Ex(x,y,z)]\tRe[Ey(x,y,z)]");
            fprintf(fE,"\tIm[Ey(x,y,z)] \tRe[Ez(x,y,z)]\tIm[Ez(x,y,z)]\n");
            break;
    }

    /* Electric field and potential values at the required points */
    for(i = 0; i < nInternalPoints; i++){
        Xin[0] = Xinner[i][0];
        Xin[1] = Xinner[i][1];
        Xin[2] = Xinner[i][2];
        mPot[0] = mPot[1] = 0.0;
        mE[0] = mE[1] = mE[2] = mE[3] = mE[4] = mE[5] = 0.0;
        
        /* Exclude points within the defined columns */
        count = 0;
        for(j = 0; j < nCols; j++){
            if(nColType == 1){
                if(sqrt(pow(Xin[0] - xCols[j][0],2.0) + pow(Xin[1] - 
                xCols[j][1],2.0)) > xCols[j][2]) count++;
            }
            else{
                if((fabs(Xin[0] - xCols[j][0])> xCols[j][2]) && 
                (fabs(Xin[1] - xCols[j][1]) > xCols[j][2])) count++;
            }
        }
        if(count == nCols) nPointInside = 1;
        else nPointInside = 0;

        if(ANALYSIS == 0){
            if(nPointInside) potential_tria6(Xin,mNodes,mElems,vB,mPot);
            fprintf(fP,"%le\t%le\t%le\t",Xin[0],Xin[1],Xin[2]);
            fprintf(fP,"%le\t%le\n",mPot[0],mPot[1]);
        }
        else if(ANALYSIS == 1){
            if(nPointInside) field_tria6(Xin,mNodes,mElems,vB,mE);
            fprintf(fE,"%le\t%le\t%le\t",Xin[0],Xin[1],Xin[2]);
            fprintf(fE,"%le\t%le\t%le\t",mE[0],mE[1],mE[2]);
            fprintf(fE,"%le\t%le\t%le\n",mE[3],mE[4],mE[5]);
        }
        else{
            if(nPointInside){
                potential_tria6(Xin,mNodes,mElems,vB,mPot);
                field_tria6(Xin,mNodes,mElems,vB,mE);
            }
            fprintf(fP,"%le\t%le\t%le\t",Xin[0],Xin[1],Xin[2]);
            fprintf(fP,"%le\t%le\n",mPot[0],mPot[1]);
            fprintf(fE,"%le\t%le\t%le\t",Xin[0],Xin[1],Xin[2]);
            fprintf(fE,"%le\t%le\t%le\t",mE[0],mE[1],mE[2]);
            fprintf(fE,"%le\t%le\t%le\n",mE[3],mE[4],mE[5]);
        }
    }
    if(ANALYSIS) fclose(fE);
    if(ANALYSIS != 1) fclose(fP);

    /* Calculate Fdep at dielectric interfaces if required */
    if(ANALYSIS == 3 || ANALYSIS == 4){
        fprintf(file_log,"\n\tpostProcess_tria6(): calling forceMST_tria6() ...");
        fF = fopen("force-mst.dat","w");
        Eps = vMatParam[0][1]*eps0;
        
        forceMST_tria6(mNodes,mElems,vBCType,vB,Eps,Fce,Xce,Fcm,Xcm);
        fprintf(fF,"#Xcm\t\tYcm\t\tZcm\t\tFXcm\t\tFYcm\t\tFZcm\n");
        fprintf(fF,"%le\t%le\t%le\t%le\t%le\t%le\n",Xcm[0],Xcm[1],Xcm[2],Fcm[0],Fcm[1],Fcm[2]);
        fprintf(fF,"#Xce\t\tYce\t\tZce\t\tFXce\t\tFYce\t\tFZce\n");
        
        for(i = 0; i < nElems; i++){
            if(vBCType[mElems[i][0]-1] == 6){
                fprintf(fF,"%le\t%le\t%le\t",Xce[i][0],Xce[i][1],Xce[i][2]);
                fprintf(fF,"%le\t%le\t%le\n",Fce[i][0],Fce[i][1],Fce[i][2]);
            }
        }
        fclose(fF);
    }
    
    /* Calculate Fdep at particle centre using multipolar method if required */
    else if(ANALYSIS > 4){
        fprintf(file_log,"\n\tpostProcessCol_tria6(): calling ");
        fprintf(file_log,"forceMultipole_tria6() ...");
        nShape = 0;
        switch(ANALYSIS){
            case 5:         /* Dipolar Approximation (Order 1) */
                nOrder = 1;
                R = axis[0];
                break;
            case 6:         /* Quadrupolar Approximation (Order 2) */
                nOrder = 2;
                R = axis[0];
                break;
            case 7:         /* Octupolar Approximation (Order 3) */
                nOrder = 3;
                R = axis[0];
                break;
            case 8:         /* Multipolar Approx Order 4 */
                nOrder = 4;
                R = axis[0];
                break;
            case 9:         /* Multipolar Approx Order 5 */
                nOrder = 5;
                R = axis[0];
                break;
            case 10:        /* Dipolar Approx for Homogeneous Ellipsoid */
                nOrder = 1;
                nShape = 1;
                break;
            case 11:            /* Dipolar Approx for Single-shelled Oblate Ellipsoid */
                nOrder = 1;
                nShape = 1;
                break;
        }
        fF = fopen("force-mp.dat","w");
        fprintf(fF,"#Xcm\t\tYcm\t\tZcm\t\tFXcm\t\tFYcm\t\tFZcm\n");
        for(i = 0; i < nFPoints; i++)
        {
            Xeval[0] = XF[i][0];
            Xeval[1] = XF[i][1];
            Xeval[2] = XF[i][2];
            if(!nShape) forceMultipole_tria6(nOrder,R,Xeval,mNodes,mElems,vB,vMatParam,vProbParam,Fcm);
            else forceEllipsoid_tria6(ANALYSIS,nOrder,axis,Xeval,mNodes,mElems,vB,vMatParam,vProbParam,Fcm);
            fprintf(fF,"%le\t%le\t%le\t",Xeval[0],Xeval[1],Xeval[2]);
            fprintf(fF,"%le\t%le\t%le\n",Fcm[0],Fcm[1],Fcm[2]);
        }
        fclose(fF);
    }

    freeDoubleMatrix(Fce,nElems);
    freeDoubleMatrix(Xce,nElems);

    return 0;
}
