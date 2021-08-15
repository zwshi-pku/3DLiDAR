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

#ifndef POINT_H
#define POINT_H
#include <fstream>
#include <list>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <time.h>
#include <iomanip>
#include <math.h>
using namespace std;
struct Point_RGB
{
	int r,g,b;
	bool operator==(const Point_RGB &_point_RGB)
	{
		if(r==_point_RGB.r)
			if(g==_point_RGB.g)
				if(b==_point_RGB.b)
					return true;
		return false;
	}
};
class Point															//点类用于存放点云数据
{
public:
	double x,y,z;
	double angle, time_;

	friend bool operator<(const Point &_pnta,const Point &_pntb)	//可以用来sort和作为map的key
	{
		if(_pnta.x<_pntb.x){
			return true;
		}else if(_pnta.x>_pntb.x){
			return false;
		}else if(_pnta.y<_pntb.y){
			return true;
		}else if(_pnta.y>_pntb.y){
			return false;
		}else if(_pnta.z<_pntb.z){
			return true;
		}else if(_pnta.z>_pntb.z){
			return false;
		}else{
			return false;
		}
	}
	bool operator==(const Point &_pnt)								//可以用来unique、remove等，排序后删除相邻相同值
	{
		if(x==_pnt.x)
			if(y==_pnt.y)
				if(z==_pnt.z)
					return true;
		return false;
	}
	void operator=(const Point &_pnt)
	{
		x=_pnt.x;
		y=_pnt.y;
		z=_pnt.z;
		angle = _pnt.angle;
		time_ = _pnt.time_;
	}
	Point operator+(const Point &_point)
	{
		Point point;
		point.x=x+_point.x;
		point.y=y+_point.y;
		point.z=z+_point.z;
		return point;
	}
	Point operator-(const Point &_point)
	{
		Point point;
		point.x=x-_point.x;
		point.y=y-_point.y;
		point.z=z-_point.z;
		return point;
	}
	Point operator/(const int n)
	{
		Point point;
		point.x=x/(n*1.0);
		point.y=y/(n*1.0);
		point.z=z/(n*1.0);
		return point;
	}
	
	Point()															//构造函数
	{
		x = y = z = angle = time_ = 0;;
	}

	void input1(list<Point> &_points, string _ifilePath)								//数据读入到list
	{
		ifstream ifile;
		ifile.open(_ifilePath);
		if(!ifile)
		{
			cout<<"file data do not open ！"<<endl;
			return ;
		}		

		string tempStr;
		stringstream strings;
		strings.str(tempStr);
		Point tem;
		while (getline(ifile, tempStr))
		{
			strings.str(tempStr);
			strings >> tem.x >> tem.y >> tem.z;
			_points.push_back(tem);
			strings.clear();
		}
		ifile.close();
	}

	double point2point(Point &pnt1, Point &pnt2)
	{
		return sqrt(pow(pnt1.x - pnt2.x, 2) + pow(pnt1.y - pnt2.y, 2) + pow(pnt1.z - pnt2.z, 2));
	}
};

typedef list<Point> Points;

class PointRGB: public Point
{
public:
	int r_, g_, b_;
	PointRGB() :Point()
	{
		r_ = g_ = b_ = 0;
	}
	void output1(list<PointRGB> &_points)
	{
		string ofilePath;
		cout << "Write data:\n";
		cin >> ofilePath;

		ofstream ofile;
		ofile.open(ofilePath);
		if (!ofile)
		{
			cout << "file data do not open ！" << endl;
			return;
		}
		list<PointRGB>::iterator it;

		for (it = _points.begin(); it != _points.end(); ++it)
		{
			ofile << it->x << " " << it->y << " " << it->z << " " << it->r_ << " " << it->g_ << " " << it->b_ << endl;
		}
		ofile.close();
	}

};


#endif