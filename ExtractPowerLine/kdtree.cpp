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

//#pragma once
#include "kdtree.h"
int KdTree::kNearestNeighbor(double &_x,double &_y,double &_z)
{
	pnts_.clear();
	nNearestNodes nNN(numberOfNeighbours_);
	nNN.setDistandIndx(radius_);
	nNN.setSearchPnt(_x,_y,_z);
	t_.locateNodes(&nNN,1);
	//cout<<"找到的点数："<<nNN.found<<endl;
	Point p;
	for(int i=1;i<=nNN.found;++i)
	{
		p.x=nNN.index[i]->pos[0];
		p.y=nNN.index[i]->pos[1];
		p.z=nNN.index[i]->pos[2];
		pnts_.push_back(p);
	}
	return pnts_.size();	//返回最临近点个数
}

list<Point> KdTree::kNearestNeighbor( Point &_point )
{
	pnts_.clear();
	nNearestNodes nNN(numberOfNeighbours_);
	nNN.setDistandIndx(radius_);
	nNN.setSearchPnt(_point.x,_point.y,_point.z);
	t_.locateNodes(&nNN,1);
	//cout<<"找到的点数："<<nNN.found<<endl;
	Point p;
	for(int i=1;i<=nNN.found;++i)
	{
		p.x=nNN.index[i]->pos[0];
		p.y=nNN.index[i]->pos[1];
		p.z=nNN.index[i]->pos[2];
		pnts_.push_back(p);
	}
	return pnts_;		//返回最临近点个数
}
