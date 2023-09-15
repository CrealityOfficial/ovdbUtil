#ifndef ADJACENTOCTREE_H
#define ADJACENTOCTREE_H

#include "trimesh2/TriMesh.h"

#include <string>
#include <map>
#include <list>
#include <set>


// 八叉树是基于模型网格的相对坐标建立的
// 八叉树的分裂，基于右手系，以 x 轴正方向为前，y 轴正方向为左，z 轴正方向为上
class MAOTreeNode
{
public:
	MAOTreeNode(std::string uid, trimesh::TriMesh* pModelMesh, int depth, MAOTreeNode* parent);
	~MAOTreeNode();

	void setExtent(double xMin, double xMax, double yMin, double yMax, double zMin, double zMax);
	int getDepth() { return m_depth; }
	int getSubTreeDepth() { return m_subTreeDepth; }

	bool isSpliting();
	void split();

	void updateAdjacency(int adjacencyIndex, MAOTreeNode* oldAdjacency, std::vector<MAOTreeNode*> newAdjacency);

	int addData(int faceID);
	void clearData();
	bool isContained(int faceID);

	bool hasChild();
	MAOTreeNode* getChild(bool up, bool left, bool front);
	void clearAll();

	void getAdjacentFaceID(int targetID, std::vector<int>& adjacencyIDArray);
	bool lineCollide(trimesh::dvec3 linePos, trimesh::dvec3 lineDir, std::vector<std::pair<int, trimesh::dvec3>>& intersectedFaceIDNPosArray);

private:
	inline int getChildIndex(bool up, bool left, bool front);

private:
	std::string m_uid;
	int m_depth;
	// 子树深度，包含自身
	int m_subTreeDepth;
	int m_splitThreshold;

	double m_xMin;
	double m_xMax;
	double m_yMin;
	double m_yMax;
	double m_zMin;
	double m_zMax;

	trimesh::dvec3 m_origin;
	trimesh::TriMesh* m_pModelMesh;

	std::set<int> m_faceIDSet;
	std::map<trimesh::dvec3, std::vector<int>> m_adjacentFaceMap;

	MAOTreeNode* m_pParent;
	std::vector<MAOTreeNode*> m_childArray;

	// std::vector(up, down, left, right, front, back)（暂时不使用邻接功能）
	//std::vector<std::list<MAOTreeNode*>> m_adjacencyArray;
};


class ModelAdjacentOctree
{
public:
	ModelAdjacentOctree(trimesh::TriMesh* modelMesh = nullptr);
	~ModelAdjacentOctree();

	void setModelMesh(trimesh::TriMesh* modelMesh);

	int getTreeDepth();
	void getAdjacentFaceID(int targetID, std::vector<int>& adjacencyIDArray);
	bool lineCollide(trimesh::dvec3 linePos, trimesh::dvec3 lineDir, std::vector<std::pair<int, trimesh::dvec3>>& intersectedFaceIDNPosArray);

private:
	void buildTree();

private:
	trimesh::TriMesh* m_pModelMesh;

	MAOTreeNode* m_pRootNode;
};

#endif // ADJACENTOCTREE_H