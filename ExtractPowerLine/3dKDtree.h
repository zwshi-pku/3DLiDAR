/****************************************************************************
   Copyright 2020 Original developer

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

#ifndef _3DKDTREE_H_
#define _3DKDTREE_H_

#define TYPE double

#include <iostream>
using namespace std;

//-------------------------------
//     KD树节点
//-------------------------------
typedef struct kdNode
{
	TYPE pos[3];  //节点坐标
	short plane;  //节点的分割平面
	int id;       //从UCD文件格式读取的节点的ID
}kdNode;

//-------------------------------------------
//    KD树
//-------------------------------------------
class nNearestNodes;
class kdTree
{
public:
	//
	//  构造函数
	//  nodes:树的最大节点数目，用来分配内存空间
	//
	kdTree(const int nodes);
	kdTree(){}	//无参构造函数
	//
	//  析构函数
	//
	~kdTree();
	//
	//   在树中存储一个节点,按读取顺序存储的
	//
	inline void store(const TYPE xCoord,const TYPE yCoord,const TYPE zCoord,const int pId);
	inline void store(const TYPE coords[3],const int pId);

	//
	//  平衡树
	//
	void treeBalance();

	//
	//  在index开始的子树中寻找邻近点
	//
	void locateNodes(nNearestNodes * const nNN,const int index)const;

	//
	//  
	//
	void balancePartition(kdNode** pBalanced,kdNode** pOriginal,const int index,const int start,const int end);

	//
	//
	//
	void medianPartition(kdNode** pOrig,const int start,const int end,const int median,const int axis);
protected:
	kdNode* kdNodes;

	int storedKDNodes;
	int halfStoredKDNodes;  //非叶子节点数
	int maxNumOfNodes;

	double boundrayMin[3];  //三个方向的最小值
	double boundrayMax[3];  //三个方向的最大值，这两个值相当于所有三维点的包围框
};

//----------------------------------------
//   最近临节点
//----------------------------------------
class nNearestNodes
{
public:
	int max;       //最近临节点的数目
	int found;     //搜索过程中找到的nNN的数目
	int got_Heap;  //寻找nNN时是否已经建堆  0表示没有建堆
	double pos[3]; //指定点的坐标
	double* dist2; //找到最近临点到指定点的距离的平方
	const kdNode** index;  //找到的最近临点

	//
	//   构造函数
	//
	nNearestNodes(const int maximum)
	{
		max = maximum;
		found = got_Heap = 0;
	}
	//
	//    析构函数
	//
	~nNearestNodes()
	{
		delete[]dist2;
		delete[]index;
	}

	//dis ：搜索距离的阈值
	void setDistandIndx(double dis = 1.0)
	{
		dist2 = new double[max + 1];
		dist2[0] = dis * dis;   //最小搜索距离的平方,小于这个距离的点认为是要搜索的
		//dist2[0] = dis;
		index = new const kdNode*[max + 1];
	}

	void setSearchPnt(const double x,const double y,const double z)
	{
		pos[0] = x;
		pos[1] = y;
		pos[2] = z;
	}
};

//-----------------------------------------------------------------
void kdTree::store(const TYPE xCoord,const TYPE yCoord,const TYPE zCoord,const int pId)
{
	const TYPE coords[3] = {xCoord,yCoord,zCoord};
	store(coords,pId);
}

void kdTree::store(const TYPE coords[3],const int pId)
{
	if(storedKDNodes > maxNumOfNodes)
	{
		cout << "No room for more nodes..." << endl;
		return;
	}

	storedKDNodes++;
	kdNode *const node = &kdNodes[storedKDNodes];
	for(int i = 0; i < 3; i++)
	{
		node->pos[i] = coords[i];

		if(node->pos[i] < boundrayMin[i])
			boundrayMin[i] = node->pos[i];
		if(node->pos[i] > boundrayMax[i])
			boundrayMax[i] = node->pos[i];
	}
	node->id = pId;
}

#endif