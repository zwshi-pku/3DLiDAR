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

#ifndef KDTREE_H
#define KDTREE_H
#include <list>
#include <math.h>
#include "3dKDtree.h"
#include "point.h"
using namespace std;
class KdTree
{
	kdTree t_;
	int numberOfNeighbours_;						
	double radius_;									
	list<Point> pnts_;								
	int id;											
public:
	KdTree(list<Point> &_p):t_(_p.size())			
	{
		list<Point>::iterator it;
		
		id=0;
		for(it=_p.begin();it!=_p.end();++it)
		{
			t_.store(it->x,it->y,it->z,id);
			++id;
		}
		t_.treeBalance();							
	}
	KdTree(vector<Point> &_p):t_(_p.size())			
	{
		vector<Point>::iterator it;

		id=0;
		for(it=_p.begin();it!=_p.end();++it)
		{
			t_.store(it->x,it->y,it->z,id);
			++id;
		}
		t_.treeBalance();							
	}

	KdTree(list<Point> &_p,string c):t_(_p.size())			
	{
		list<Point>::iterator it;

		id=0;
		
		if (c=="x")
		{
			for(it=_p.begin();it!=_p.end();++it)
			{
				t_.store(0,it->y,it->z,id);
				++id;
			}
			t_.treeBalance();							
		}else if (c=="y")
		{
			for(it=_p.begin();it!=_p.end();++it)
			{
				t_.store(it->x,0,it->z,id);
				++id;
			}
			t_.treeBalance();							
		}else if (c=="z")
		{
			for(it=_p.begin();it!=_p.end();++it)
			{
				t_.store(it->x,it->y,0,id);
				++id;
			}
			t_.treeBalance();							
		}else
		{
			cout<<"build tree error!";
		}

	}
	KdTree(list<Point> &_p,int _size):t_(_size)		
	{
		list<Point>::iterator it;

		id=0;
		for(it=_p.begin();it!=_p.end();++it)
		{
			t_.store(it->x,it->y,it->z,id);
			++id;
		}
		t_.treeBalance();							
	}

	void setNumberOfNeighbours(int _num)			
	{
		numberOfNeighbours_=_num;
		radius_=10.0;								
	}
	void setNeighboursRadius(double _radius)		
	{
		radius_=_radius;
		numberOfNeighbours_=5000;					
	}
	int kNearestNeighbor(double &_x,double &_y,double &_z);	
	list<Point> kNearestNeighbor(Point &_point);

	void pushBack(double &_x,double &_y,double &_z)
	{
		t_.store(_x,_y,_z,id);
		++id;
		t_.treeBalance();							
	}

	list<Point>& getNearestNeighbor()				
	{
		return pnts_;
	}
};

#endif