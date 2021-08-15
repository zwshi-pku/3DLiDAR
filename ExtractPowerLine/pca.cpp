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

#pragma once
#include "pca.h"
//#include "../EigenMatrix/eigenmatrix.h"

PCA::PCA()
{

}
double PCA::planarity( list<Point> &_points )
{
	Matrix3d coVariaceMatrix;											//covariance_matrix of xyz
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			coVariaceMatrix(i, j) = 0;
		}
	}

	double meanx_, meany_, meanz_;										//相应维数的均值
	meanx_ = meany_ = meanz_ = 0.0;
	list<Point>::iterator it;
	for (it = _points.begin(); it != _points.end(); ++it)
	{
		meanx_ += it->x;
		meany_ += it->y;
		meanz_ += it->z;
	}
	meanx_ /= double(_points.size());
	meany_ /= double(_points.size());
	meanz_ /= double(_points.size());


	for (list<Point>::iterator it = _points.begin(); it != _points.end(); ++it)
	{
		coVariaceMatrix(0, 0) += (it->x - meanx_)*(it->x - meanx_);		//cov(x,x)
		coVariaceMatrix(0, 1) += (it->x - meanx_)*(it->y - meany_);		//cov(x,y)
		coVariaceMatrix(0, 2) += (it->x - meanx_)*(it->z - meanz_);		//cov(x,z)

		coVariaceMatrix(1, 0) = coVariaceMatrix(0, 1);					//cov(y,x)
		coVariaceMatrix(1, 1) += (it->y - meany_)*(it->y - meany_);		//cov(y,y)
		coVariaceMatrix(1, 2) += (it->y - meany_)*(it->z - meanz_);		//cov(y,z)

		coVariaceMatrix(2, 0) = coVariaceMatrix(0, 2);					//cov(x,z)
		coVariaceMatrix(2, 1) = coVariaceMatrix(1, 2);					//cov(y,z)
		coVariaceMatrix(2, 2) += (it->z - meanz_)*(it->z - meanz_);		//cov(z,z)
	}

	coVariaceMatrix(0, 0) /= double((_points.size() - 1));						//cov(x,x)
	coVariaceMatrix(0, 1) /= double((_points.size() - 1));						//cov(x,y)
	coVariaceMatrix(0, 2) /= double((_points.size() - 1));						//cov(x,z)

	coVariaceMatrix(1, 0) /= double((_points.size() - 1));						//cov(y,x)
	coVariaceMatrix(1, 1) /= double((_points.size() - 1));						//cov(y,y)
	coVariaceMatrix(1, 2) /= double((_points.size() - 1));						//cov(y,z)

	coVariaceMatrix(2, 0) /= double((_points.size() - 1));						//cov(x,z)
	coVariaceMatrix(2, 1) /= double((_points.size() - 1));						//cov(y,z)
	coVariaceMatrix(2, 2) /= double((_points.size() - 1));						//cov(z,z)
	
	MatrixXd eigvalue(1, 3);
	sortEigenvalue(coVariaceMatrix, eigvalue);
	
	return (eigvalue(0, 1) - eigvalue(0, 0)) / eigvalue(0, 2);
}

void PCA::egenValue(list<Point> &_points, double &_eValue1, double &_eValue2, double &_eValue3)
{
	Matrix3d coVariaceMatrix;											//covariance_matrix of xyz
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			coVariaceMatrix(i, j) = 0;
		}
	}

	double meanx_, meany_, meanz_;										//相应维数的均值
	meanx_ = meany_ = meanz_ = 0.0;
	list<Point>::iterator it;
	for (it = _points.begin(); it != _points.end(); ++it)
	{
		meanx_ += it->x;
		meany_ += it->y;
		meanz_ += it->z;
	}
	meanx_ /= double(_points.size());
	meany_ /= double(_points.size());
	meanz_ /= double(_points.size());
	
	for (list<Point>::iterator it = _points.begin(); it != _points.end(); ++it)
	{
		coVariaceMatrix(0, 0) += (it->x - meanx_)*(it->x - meanx_);		//cov(x,x)
		coVariaceMatrix(0, 1) += (it->x - meanx_)*(it->y - meany_);		//cov(x,y)
		coVariaceMatrix(0, 2) += (it->x - meanx_)*(it->z - meanz_);		//cov(x,z)

		coVariaceMatrix(1, 0) = coVariaceMatrix(0, 1);					//cov(y,x)
		coVariaceMatrix(1, 1) += (it->y - meany_)*(it->y - meany_);		//cov(y,y)
		coVariaceMatrix(1, 2) += (it->y - meany_)*(it->z - meanz_);		//cov(y,z)

		coVariaceMatrix(2, 0) = coVariaceMatrix(0, 2);					//cov(x,z)
		coVariaceMatrix(2, 1) = coVariaceMatrix(1, 2);					//cov(y,z)
		coVariaceMatrix(2, 2) += (it->z - meanz_)*(it->z - meanz_);		//cov(z,z)
	}

	coVariaceMatrix(0, 0) /= double((_points.size() - 1));						//cov(x,x)
	coVariaceMatrix(0, 1) /= double((_points.size() - 1));						//cov(x,y)
	coVariaceMatrix(0, 2) /= double((_points.size() - 1));						//cov(x,z)

	coVariaceMatrix(1, 0) /= double((_points.size() - 1));						//cov(y,x)
	coVariaceMatrix(1, 1) /= double((_points.size() - 1));						//cov(y,y)
	coVariaceMatrix(1, 2) /= double((_points.size() - 1));						//cov(y,z)

	coVariaceMatrix(2, 0) /= double((_points.size() - 1));						//cov(x,z)
	coVariaceMatrix(2, 1) /= double((_points.size() - 1));						//cov(y,z)
	coVariaceMatrix(2, 2) /= double((_points.size() - 1));						//cov(z,z)

	MatrixXd eigvalue(1, 3);
	sortEigenvalue(coVariaceMatrix, eigvalue);
	_eValue1 = eigvalue(0, 0);
	_eValue2 = eigvalue(0, 1);
	_eValue3 = eigvalue(0, 2);
}

/**对向量进行排序，从大到小
* vec: 待排序的向量
* sorted_vec: 排序的结果
* ind: 排序结果中各个元素在原始向量的位置
*/
void sort_vec(const VectorXd& vec, VectorXd& sorted_vec, VectorXi& ind){
	ind = VectorXi::LinSpaced(vec.size(), 0, vec.size() - 1);//[0 1 2 3 ... N-1]
	auto rule = [vec](int i, int j)->bool{
		return vec(i)>vec(j);
	};
	std::sort(ind.data(), ind.data() + ind.size(), rule);
	sorted_vec.resize(vec.size());
	for (int i = 0; i<vec.size(); i++){
		sorted_vec(i) = vec(ind(i));
	}
}

////测试
//int main(){
//	VectorXd x(5);
//	x << 3, 4, 1, 5, 6;
//	VectorXi ind;
//	VectorXd sorted_vec;
//	sort_vec(x, sorted_vec, ind);
//	cout << "原始向量:\n";
//	cout << x << endl << endl;
//	cout << "排序后:\n";
//	cout << sorted_vec << endl << endl;
//	cout << "排序后向量各元素对应的原始向量中的位置" << endl;
//	cout << ind << endl;
//
//	return 0;
//}
//-------------------- -
//作者：X_And_Y
//来源：CSDN
//原文：https ://blog.csdn.net/X_And_Y/article/details/83383520 
//版权声明：本文为博主原创文章，转载请附上博文链接！

void PCA::eigenDV(list<Point> &_points, Matrix3d &_D, Matrix3d &_V, VectorXi &_indexs)
{
	Matrix3d coVariaceMatrix;											//covariance_matrix of xyz
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			coVariaceMatrix(i, j) = 0;
		}
	}

	double meanx_, meany_, meanz_;										//相应维数的均值
	meanx_ = meany_ = meanz_ = 0.0;
	list<Point>::iterator it;
	for (it = _points.begin(); it != _points.end(); ++it)
	{
		meanx_ += it->x;
		meany_ += it->y;
		meanz_ += it->z;
	}
	meanx_ /= double(_points.size());
	meany_ /= double(_points.size());
	meanz_ /= double(_points.size());

	for (list<Point>::iterator it = _points.begin(); it != _points.end(); ++it)
	{
		coVariaceMatrix(0, 0) += (it->x - meanx_)*(it->x - meanx_);		//cov(x,x)
		coVariaceMatrix(0, 1) += (it->x - meanx_)*(it->y - meany_);		//cov(x,y)
		coVariaceMatrix(0, 2) += (it->x - meanx_)*(it->z - meanz_);		//cov(x,z)

		coVariaceMatrix(1, 0) = coVariaceMatrix(0, 1);					//cov(y,x)
		coVariaceMatrix(1, 1) += (it->y - meany_)*(it->y - meany_);		//cov(y,y)
		coVariaceMatrix(1, 2) += (it->y - meany_)*(it->z - meanz_);		//cov(y,z)

		coVariaceMatrix(2, 0) = coVariaceMatrix(0, 2);					//cov(x,z)
		coVariaceMatrix(2, 1) = coVariaceMatrix(1, 2);					//cov(y,z)
		coVariaceMatrix(2, 2) += (it->z - meanz_)*(it->z - meanz_);		//cov(z,z)
	}

	coVariaceMatrix(0, 0) /= double((_points.size() - 1));						//cov(x,x)
	coVariaceMatrix(0, 1) /= double((_points.size() - 1));						//cov(x,y)
	coVariaceMatrix(0, 2) /= double((_points.size() - 1));						//cov(x,z)

	coVariaceMatrix(1, 0) /= double((_points.size() - 1));						//cov(y,x)
	coVariaceMatrix(1, 1) /= double((_points.size() - 1));						//cov(y,y)
	coVariaceMatrix(1, 2) /= double((_points.size() - 1));						//cov(y,z)

	coVariaceMatrix(2, 0) /= double((_points.size() - 1));						//cov(x,z)
	coVariaceMatrix(2, 1) /= double((_points.size() - 1));						//cov(y,z)
	coVariaceMatrix(2, 2) /= double((_points.size() - 1));						//cov(z,z)

	EigenSolver<Matrix3d> es(coVariaceMatrix);
	_D = es.pseudoEigenvalueMatrix();
	_V = es.pseudoEigenvectors();

	//cout << "The pseudo-eigenvalue matrix D is:" << endl << _D << endl;
	//cout << "The pseudo-eigenvector matrix V is:" << endl << _V << endl;
		
	//_D.maxCoeff(&_row, &_col);

	//cout << _row << " " << _col << endl;

	VectorXd Ds(3);
	Ds << _D(0, 0), _D(1, 1), _D(2, 2);
	VectorXd sorted_Ds;
	sort_vec(Ds, sorted_Ds, _indexs);

	//cout << "原始向量:\n";
	//cout << Ds << endl << endl;
	//cout << "排序后:\n";
	//cout << sorted_Ds << endl << endl;
	//cout << "排序后向量各元素对应的原始向量中的位置" << endl;
	//cout << index << endl;	
}
