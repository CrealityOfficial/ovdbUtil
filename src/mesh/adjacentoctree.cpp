#include "adjacentoctree.h"

#include <algorithm>
#include <functional>


MAOTreeNode::MAOTreeNode(std::string uid, trimesh::TriMesh* pModelMesh, int depth, MAOTreeNode* parent)
	: m_uid(uid)
	, m_pModelMesh(pModelMesh)
	, m_pParent(parent)
	, m_depth(depth)
	, m_subTreeDepth(1)
	, m_splitThreshold(100)
{
	m_origin = m_pModelMesh->vertices[0];

	m_childArray.resize(8);
	//m_adjacencyArray.resize(6);
}

MAOTreeNode::~MAOTreeNode()
{
	m_pModelMesh = nullptr;
	m_pParent = nullptr;

	clearAll();
}

void MAOTreeNode::setExtent(double xMin, double xMax, double yMin, double yMax, double zMin, double zMax)
{
	m_xMin = xMin;
	m_xMax = xMax;
	m_yMin = yMin;
	m_yMax = yMax;
	m_zMin = zMin;
	m_zMax = zMax;
}

bool MAOTreeNode::isSpliting()
{
	return m_faceIDSet.size() > m_splitThreshold;
}

void MAOTreeNode::split()
{
	double halfXLength = (m_xMax - m_xMin) / 2.0;
	double halfYLength = (m_yMax - m_yMin) / 2.0;
	double halfZLength = (m_zMax - m_zMin) / 2.0;

	// 分裂子节点
	for (int up = 0; up < 2; ++up)
	{
		for (int left = 0; left < 2; ++left)
		{
			for (int front = 0; front < 2; ++front)
			{
				std::string childUid = m_uid + std::to_string(up) + std::to_string(left) + std::to_string(front);
				MAOTreeNode* childNode = new MAOTreeNode(childUid, m_pModelMesh, m_depth + 1, this);
				childNode->setExtent(
					m_xMin + front * halfXLength,
					m_xMax - (1 - front) * halfXLength,
					m_yMin + left * halfYLength,
					m_yMax - (1 - left) * halfYLength,
					m_zMin + up * halfZLength,
					m_zMax - (1 - up) * halfZLength
				);

				m_childArray[getChildIndex(up, left, front)] = childNode;
			}
		}
	}

	// 分配邻接节点（暂不使用）
	{
		//for (int i = 0; i < m_childArray.size(); ++i)
		//{
		//	MAOTreeNode* childNode = m_childArray[i];

		//	// 分配父节点的相邻节点
		//	int tempi = i;
		//	for (int j = 0; j < 3; ++j)
		//	{
		//		int adjacencyIndex = tempi % 2 + 2 * j;
		//		
		//		std::vector<MAOTreeNode*> childAdjacencyArray;
		//		for (int k = 0; k < m_adjacencyArray[adjacencyIndex].size(); ++k)
		//		{
		//			auto itr = m_adjacencyArray[adjacencyIndex].begin();
		//			std::advance(itr, k);
		//			MAOTreeNode* adjacentNode = *itr;
		//			if (adjacentNode->getDepth() >= m_depth)
		//			{
		//				childAdjacencyArray.push_back(adjacentNode);
		//			}
		//			else
		//			{
		//				bool adjacented =
		//					(j == 0 && adjacentNode->m_xMin >= childNode->m_xMin && adjacentNode->m_xMax <= childNode->m_xMax
		//						&& adjacentNode->m_yMin >= childNode->m_yMin && adjacentNode->m_yMax <= childNode->m_yMax) ||
		//					(j == 1 && adjacentNode->m_xMin >= childNode->m_xMin && adjacentNode->m_xMax <= childNode->m_xMax
		//						&& adjacentNode->m_zMin >= childNode->m_zMin && adjacentNode->m_zMax <= childNode->m_zMax) ||
		//					(j == 2 && adjacentNode->m_yMin >= childNode->m_yMin && adjacentNode->m_yMax <= childNode->m_yMax
		//						&& adjacentNode->m_zMin >= childNode->m_zMin && adjacentNode->m_zMax <= childNode->m_zMax);

		//				if (adjacented)
		//					childAdjacencyArray.push_back(adjacentNode);
		//			}
		//		}

		//		childNode->updateAdjacency(adjacencyIndex, nullptr, childAdjacencyArray);

		//		tempi = tempi >> 1;
		//	}

		//	// 分配兄弟节点里的相邻节点
		//	bool up = ((i >> 1) >> 1) % 2;
		//	bool left = (i >> 1) % 2;
		//	bool front = i % 2;
		//	tempi = i;
		//	for (int j = 0; j < 3; ++j)
		//	{
		//		int adjacencyIndex = (1 - tempi % 2) + 2 * j;
		//		int adjacentNodeIndex = getChildIndex((j == 0 ? -1 : 1) * up, (j == 1 ? -1 : 1) * left, (j == 2 ? -1 : 1) * front);
		//		childNode->updateAdjacency(adjacencyIndex, nullptr, { m_childArray[adjacentNodeIndex] });
		//	}
		//}
	}

	// 分配所包含的面
	int maxSubTreeDepth = 1;
	std::set<int> faceIDSetTemp = m_faceIDSet;
	for (auto itr = faceIDSetTemp.begin(); itr != faceIDSetTemp.end(); ++itr)
	{
		int faceID = *itr;
		for (int i = 0; i < m_childArray.size(); ++i)
		{
			int depth = m_childArray[i]->addData(faceID);
			if (depth > 0)
			{
				if (maxSubTreeDepth < depth)
					maxSubTreeDepth = depth;

				m_faceIDSet.erase(faceID);

				break;
			}
		}
	}

	m_subTreeDepth = maxSubTreeDepth + 1;
}

void MAOTreeNode::updateAdjacency(int adjacencyIndex, MAOTreeNode* oldAdjacency, std::vector<MAOTreeNode*> newAdjacency)
{
	//if (oldAdjacency)
	//	m_adjacencyArray[adjacencyIndex].remove(oldAdjacency);

	//for (int i = 0; i < newAdjacency.size(); ++i)
	//{
	//	m_adjacencyArray[adjacencyIndex].push_back(newAdjacency[i]);
	//}
}

int MAOTreeNode::addData(int faceID)
{
	// 判断是否包含该面
	if (!isContained(faceID))
		return -1;

	bool notExisted = m_faceIDSet.insert(faceID).second;
	// faceID 是否已存在
	if (notExisted)
	{
		// 是否有子节点，即是否已分裂过了
		if (hasChild())
		{
			int depth = m_subTreeDepth - 1;
			for (int i = 0; i < m_childArray.size(); ++i)
			{
				int newDepth = m_childArray[i]->addData(faceID);

				if (newDepth != -1)
				{
					m_faceIDSet.erase(faceID);

					if (depth < newDepth)
						depth = newDepth;

					break;
				}
			}

			m_subTreeDepth = depth + 1;
		}
		// 是否需要进行分裂
		else if (isSpliting())
		{
			split();

			m_adjacentFaceMap.clear();
		}
		// 更新局部邻接表（暂不使用）
		else
		{
			//trimesh::TriMesh::Face newFace = m_pModelMesh->faces[faceID];
			//for (int i = 0; i < 3; ++i)
			//{
			//	trimesh::point vertex = m_pModelMesh->vertices[newFace[i]];
			//	m_adjacentFaceMap[trimesh::dvec3(vertex.x, vertex.y, vertex.z)].push_back(faceID);
			//}

			m_subTreeDepth = 1;
		}
	}

	return m_subTreeDepth;
}

void MAOTreeNode::clearData()
{
	m_faceIDSet.clear();
	m_adjacentFaceMap.clear();
}

bool MAOTreeNode::isContained(int faceID)
{
	trimesh::TriMesh::Face newFace = m_pModelMesh->faces[faceID];
	for (int i = 0; i < 3; ++i)
	{
		trimesh::point vertex = m_pModelMesh->vertices[newFace[i]];
		trimesh::point relaVertex = vertex - trimesh::point(m_origin.x, m_origin.y, m_origin.z);

		if (relaVertex.x < m_xMin || relaVertex.x > m_xMax ||
			relaVertex.y < m_yMin || relaVertex.y > m_yMax ||
			relaVertex.z < m_zMin || relaVertex.z > m_zMax)
			return false;
	}

	return true;
}

bool MAOTreeNode::hasChild()
{
	return m_childArray[0] != nullptr;
}

MAOTreeNode* MAOTreeNode::getChild(bool up, bool left, bool front)
{
	return m_childArray[getChildIndex(up, left, front)];
}

void MAOTreeNode::clearAll()
{
	clearData();

	for (int i = 0; i < m_childArray.size(); ++i)
	{
		delete m_childArray[i];
		m_childArray[i] = nullptr;
	}
	
	//for (int i = 0; i < m_adjacencyArray.size(); i++)
	//{
	//	m_adjacencyArray[i].clear();
	//}
}

void MAOTreeNode::getAdjacentFaceID(int targetID, std::vector<int>& adjacencyIDArray)
{
}

bool MAOTreeNode::lineCollide(trimesh::dvec3 linePos, trimesh::dvec3 lineDir, std::vector<std::pair<int, trimesh::dvec3>>& intersectedFaceIDNPosArray)
{
	// 判断是否与自身相交
	// 边缘点投影到碰撞线上
	bool intersected = false;
	//double offset = 0.00001;
	trimesh::dvec3 linePosRela = linePos - m_origin;
	for (int i = 0; i < 8; i++)
	{
		double tempX = i % 2 ? m_xMin : m_xMax;
		double tempY = (i >> 1) % 2 ? m_yMin : m_yMax;
		double tempZ = ((i >> 1) >> 1) % 2 ? m_zMin : m_zMax;
		trimesh::dvec3 tempPos(tempX, tempY, tempZ);

		trimesh::dvec3 projPos = linePosRela + ((tempPos - linePosRela) DOT lineDir ) * lineDir;
		if (projPos.x >= m_xMin && projPos.x <= m_xMax &&
			projPos.y >= m_yMin && projPos.y <= m_yMax &&
			projPos.z >= m_zMin && projPos.z <= m_zMax )
		{
			intersected = true;
			break;
		}
	}

	if (!intersected)
		return false;

	// 与子节点求交点
	if (hasChild())
	{
		for (int i = 0; i < m_childArray.size(); i++)
		{
			m_childArray[i]->lineCollide(linePos, lineDir, intersectedFaceIDNPosArray);
		}
	}

	// 与该层三角面求交点
	for (auto itr = m_faceIDSet.begin(); itr != m_faceIDSet.end(); itr++)
	{
		int faceID = *itr;

		trimesh::TriMesh::Face face = m_pModelMesh->faces[faceID];
		trimesh::dvec3 A(m_pModelMesh->vertices[face[0]]);
		trimesh::dvec3 B(m_pModelMesh->vertices[face[1]]);
		trimesh::dvec3 C(m_pModelMesh->vertices[face[2]]);
		trimesh::dvec3 faceNormal = trimesh::normalized((B - A) TRICROSS (C - A));
		
		double temp = lineDir.x * faceNormal.x + lineDir.y * faceNormal.y + lineDir.z * faceNormal.z;
		if (temp == 0)
			continue;

		double t = ((A.x - linePos.x) * faceNormal.x + (A.y - linePos.y) * faceNormal.y + (A.z - linePos.z) * faceNormal.z) / temp;

		trimesh::dvec3 planeInterPos(lineDir.x * t + linePos.x, lineDir.y * t + linePos.y, lineDir.z * t + linePos.z);

		trimesh::dvec3 v0 = C - A;
		trimesh::dvec3 v1 = B - A;
		trimesh::dvec3 v2 = planeInterPos - A;
		double dot00 = v0 DOT v0;
		double dot01 = v0 DOT v1;
		double dot02 = v0 DOT v2;
		double dot11 = v1 DOT v1;
		double dot12 = v1 DOT v2;

		double bottom = 1.0 / (dot00 * dot11 - dot01 * dot01);
		double u = (dot11 * dot02 - dot01 * dot12) * bottom;
		double v = (dot00 * dot12 - dot01 * dot02) * bottom;
		if (u >= 0 && v >= 0 && u + v <= 1)
			intersectedFaceIDNPosArray.push_back({faceID, planeInterPos});
	}

	return intersectedFaceIDNPosArray.size() > 0;
}

int MAOTreeNode::getChildIndex(bool up, bool left, bool front)
{
	return up * 4 + left * 2 + front;
}

ModelAdjacentOctree::ModelAdjacentOctree(trimesh::TriMesh* modelMesh)
	: m_pRootNode(nullptr)
{
	setModelMesh(modelMesh);
}

ModelAdjacentOctree::~ModelAdjacentOctree()
{
	if (m_pRootNode)
	{
		delete m_pRootNode;
		m_pRootNode = nullptr;
	}

	m_pModelMesh = nullptr;
}

void ModelAdjacentOctree::setModelMesh(trimesh::TriMesh* modelMesh)
{
	if (modelMesh && m_pModelMesh != modelMesh)
	{
		m_pModelMesh = modelMesh;
		buildTree();
	}
}

int ModelAdjacentOctree::getTreeDepth()
{
	return m_pRootNode ? m_pRootNode->getSubTreeDepth() : 0;
}

void ModelAdjacentOctree::getAdjacentFaceID(int targetID, std::vector<int>& adjacencyIDArray)
{
	if (m_pRootNode)
		m_pRootNode->getAdjacentFaceID(targetID, adjacencyIDArray);
}

bool ModelAdjacentOctree::lineCollide(trimesh::dvec3 linePos, trimesh::dvec3 lineDir, std::vector<std::pair<int, trimesh::dvec3>>& intersectedFaceIDNPosArray)
{
	bool intersected = false;
	if (m_pRootNode)
	{
		intersected = m_pRootNode->lineCollide(linePos, lineDir, intersectedFaceIDNPosArray);
		if (intersected)
		{
			// 按距离排序
			std::map<trimesh::dvec3, double> distance2Map;
			std::function<bool(std::pair<int, trimesh::dvec3>, std::pair<int, trimesh::dvec3>)> compareFunc =
				[&linePos, &distance2Map](std::pair<int, trimesh::dvec3> A, std::pair<int, trimesh::dvec3> B)->bool 
			{
				double dis2A;
				auto itrA = distance2Map.find(A.second);
				if (itrA != distance2Map.end())
				{
					dis2A = itrA->second;
				}
				else
				{
					dis2A = trimesh::dist2(A.second, linePos);
					distance2Map[A.second] = dis2A;
				}

				double dis2B;
				auto itrB = distance2Map.find(B.second);
				if (itrB != distance2Map.end())
				{
					dis2B = itrB->second;
				}
				else
				{
					dis2B = trimesh::dist2(B.second, linePos);
					distance2Map[B.second] = dis2B;
				}

				return dis2A < dis2B;
			};

			std::sort(intersectedFaceIDNPosArray.begin(), intersectedFaceIDNPosArray.end(), compareFunc);

		}
	}

	return intersected;
}

void ModelAdjacentOctree::buildTree()
{
	if (m_pRootNode)
	{
		m_pRootNode->clearAll();
	}
	else
	{
		m_pRootNode = new MAOTreeNode("_", m_pModelMesh, 1, nullptr);
	}

	m_pModelMesh->need_bsphere();
	double diameter = m_pModelMesh->bsphere.r * 2;
	m_pRootNode->setExtent(-diameter, diameter, -diameter, diameter, -diameter, diameter);

	for (int i = 0; i < m_pModelMesh->faces.size(); ++i)
	{
		m_pRootNode->addData(i);
	}
}
