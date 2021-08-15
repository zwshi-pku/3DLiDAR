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

#include "pca.cpp"
#include "detectpowerline.h"
#include "kdtree.h"

#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
using namespace std;

void DetectPowerline::extractPL(list<Point> &_points, list<int> &_isPLIndex)
{
	KdTree tree(_points);
	tree.setNeighboursRadius(radius_);
	
	int id(0);
	for (list<Point>::iterator it = _points.begin(); it != _points.end(); ++it)
	{
		tree.kNearestNeighbor(it->x, it->y, it->z);
		list<Point> neigPoints = tree.getNearestNeighbor();
		if (neigPoints.size() < 3)
		{
			_isPLIndex.push_back(0);
			continue;
		}

		Matrix3d D, V;
		VectorXi indexs;
		PCA pca;
		pca.eigenDV(neigPoints, D, V, indexs);
		
		double eValue1(D(indexs(0), indexs(0))), eValue2(D(indexs(1), indexs(1))), eValue3(D(indexs(2), indexs(2)));
		double a(V(0, indexs(0))), b(V(1, indexs(0))), c(V(2, indexs(0)));
		double angle = acos(c / sqrt(a*a + b*b + c*c))*180.0 / 3.14159265;
				
		double L = (eValue1 - eValue2) / eValue1;
		double P = (eValue2 - eValue3) / eValue1;
		double S = eValue3 / eValue1;
		
		if (abs(angle - 90) < angleThr_ && L>LThr_)
		{
			_isPLIndex.push_back(1);
		}
		else
		{
			_isPLIndex.push_back(0);
		}
	}
}
