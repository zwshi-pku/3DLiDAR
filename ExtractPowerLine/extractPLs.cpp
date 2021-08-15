/****************************************************************************
   Copyright 2020 Institute of Remote Sensing & GIS, 
   School of Earth and Space Sciences, Peking University

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*****************************************************************************/

/****************************************************************************
 * Please cite the following paper, If you use this code in your work.
 *
 * Shi, Zhenwei, Yi Lin, and Hui Li. "Extraction of urban power lines and 
 * potential hazard analysis from mobile laser scanning point clouds." 
 * International Journal of Remote Sensing 41, no. 9 (2020): 3411-3428.
 *
 * The paper can be downloaded from
 * https://www.tandfonline.com/doi/full/10.1080/01431161.2019.1701726
 * Copyright
 * Institute of Remote Sensing & GIS, 
 * School of Earth and Space Sciences, Peking University (www.pku.edu.cn)
 * Zhenwei Shi; Yi Lin; Hui Li
 * contact us: zwshi@pku.edu.cn; lihui@pku.edu.cn
*****************************************************************************/

#include "detectpowerline.h"
#include <mex.h>
#include <string>
using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *points = mxGetPr(prhs[0]);
	double radius   = mxGetScalar(prhs[1]);
    double angle   = mxGetScalar(prhs[2]);
    double Ls   = mxGetScalar(prhs[3]);	
	int rows = mxGetM(prhs[0]);

	list<Point> pointcloud;
	Point p;
	for (int i = 0; i < rows; i++) {
		p.x = points[i + 0 * rows];
		p.y = points[i + 1 * rows];
		p.z = points[i + 2 * rows];
		pointcloud.push_back(p);
	}
	
	list<int> isPLIndex;
	DetectPowerline dpl;
    dpl.radius_ = radius;
    dpl.angleThr_ = angle;
    dpl.LThr_ = Ls;    
	dpl.extractPL(pointcloud, isPLIndex);
    
    plhs[0] = mxCreateNumericMatrix(isPLIndex.size(),1, mxINT32_CLASS, mxREAL);  
	int* outputIsPLMatrix = (int *)mxGetData(plhs[0]);

	vector<int> isPLIndexV;
	isPLIndexV.assign(isPLIndex.begin(), isPLIndex.end());

	for (int i = 0; i < isPLIndexV.size(); i++)
		outputIsPLMatrix[i] = isPLIndexV[i];
}
