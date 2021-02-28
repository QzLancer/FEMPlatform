#include "FEM2DStaticSolver.h"

void FEM2DStaticSolver::solve()
{
	//测试数据是否正确传输到求解器中
	//cout << "Number of Nodes: " << m_num_nodes << endl;
	//cout << "Number of VtxElement: " << m_num_vtxele << endl;
	//cout << "Number of EdgElement: " << m_num_edgele << endl;
	//cout << "Number of TriElement: " << m_num_triele << endl;
	//for (int i = 0; i < m_num_nodes; ++i) {
	//	cout << "Nodes: " << mp_node[i].x << " , " << mp_node[i].y << endl;
	//}
	//for (int i = 0; i < m_num_vtxele; ++i) {
	//	cout << "vtxelement: " << mp_vtxele[i].n << ", domain: " << mp_vtxele[i].domain << endl;
	//}
	//for (int i = 0; i < m_num_edgele; ++i) {
	//	cout << "edgelement: " << mp_edgele[i].n[0] << ", " << mp_edgele[i].n[1] << ", domain: " << mp_edgele[i].domain << endl;
	//}
	//for (int i = 0; i < m_num_triele; ++i) {
	//	cout << "trielement: " << mp_triele[i].n[0] << ", " << mp_triele[i].n[1] << ", " << mp_triele[i].n[2] << ", domain: " << mp_triele[i].domain << endl;
	//}
	//cout << "number of materialmap: " << materialmap.size() << endl;
	//cout << "number of loadmap: " << loadmap.size() << endl;
	//cout << "number of boundarymap : " << boundarymap.size() << endl;

	//计算三角形单元几何部分
	for (int i = 0; i < m_num_triele; ++i) {
		makeTrangle(i);
	}

	//计算自由度
	//这里的boundarymap只用了第一类边界条件，只考虑线单元域
	//1.提取线所在全部域，2.提取域对应的单元编号，3.提取节点
	//4.去重 5.给node设置边界条件 6.将边界点排序到末尾
	int num_dof = m_num_nodes;
	int num_bdr;
	std::vector<int> boundarynodes;
	//原有版本的边界检索
	for (auto iter = boundarymap.begin(); iter != boundarymap.end(); ++iter) {
		if (iter->second->getBoundaryType() == 1) {
			int boundarynodedomain = iter->first;
			for (int i = 0; i < m_num_edgele; ++i) {
				if (mp_edgele[i].domain == boundarynodedomain) {
					boundarynodes.push_back(mp_edgele[i].n[0]);
					boundarynodes.push_back(mp_edgele[i].n[1]);
				}
			}
		}
		else {
			cout << "Invalid boundarytype!\n";
			exit(0);
		}
	}
	std::sort(boundarynodes.begin(), boundarynodes.end());
	boundarynodes.erase(unique(boundarynodes.begin(), boundarynodes.end()), boundarynodes.end());
	//直接通过几何检索边界
	//for (int i = 0; i < m_num_nodes; ++i) {
	//	double x = mp_node[i].x;
	//	double y = mp_node[i].y;
	//	if (mp_node[i].x <= 1e-8 || sqrt(x * x + y * y) > 0.05 - 1e-3) {
	//		boundarynodes.push_back(i);
	//	}
	//}

	num_bdr = boundarynodes.size();
	num_dof -= num_bdr;
	for (auto a : boundarynodes) {
		mp_node[a].bdr = 1;
	}
	cout << "num_bdr: " << num_bdr << ", num_dof: " << num_dof << endl;
	//将边界点排序到末尾
	std::vector<int> node_reorder(m_num_nodes);	//前num_dof个元素对应非边界节点，之后的元素对应第一类边界条件
	std::vector<int> node_pos(m_num_nodes);	//原节点编号对应的reorder后的节点编号
	int node_bdr = 0;
	for (int ibdr = 0; ibdr < m_num_nodes; ++ibdr) {
		if (mp_node[ibdr].bdr == 1) {
			++node_bdr;
			node_reorder[m_num_nodes - node_bdr] = ibdr;
			node_pos[ibdr] = m_num_nodes - node_bdr;
			mp_node[ibdr].A = 0;
		}
		else
		{
			node_reorder[ibdr - node_bdr] = ibdr;
			node_pos[ibdr] = ibdr - node_bdr;
		}
	}


	//设置单元材料
	for (int i_tri = 0; i_tri < m_num_triele; ++i_tri) {
		int domain = mp_triele[i_tri].domain;
		mp_triele[i_tri].material = materialmap[domain];
	}

	//设置单元负载
	for (int i_tri = 0; i_tri < m_num_triele; ++i_tri) {
		int domain = mp_triele[i_tri].domain;
		if (loadmap.find(domain) != loadmap.end()) {
			mp_triele[i_tri].J = loadmap[domain];
		}
	}
	//线性三角形一阶单元分析、装配
	////不考虑第一类边界条件的装配
	//std::vector<std::vector<int>> locs(2, std::vector<int>(9 * (size_t)(m_num_triele)));
	//std::vector<double> vals(9 * (size_t)(m_num_triele));
	//std::vector<double> F(m_num_nodes);
	//int pos = 0;
	//for (int i_tri = 0; i_tri < m_num_triele; ++i_tri) {
	//	CTriElement triele = mp_triele[i_tri];
	//	for (int i = 0; i < 3; ++i) {
	//		for (int j = 0; j < 3; ++j) {
	//			if (triele.material->getLinearFlag() == true) {
	//				double mu = triele.material->getMu();
	//				double Se = triele.C[i][j] / mu;
	//				locs[0][pos] = triele.n[i];
	//				locs[1][pos] = triele.n[j];
	//				vals[pos] = Se;
	//				++pos;
	//			}
	//		}
	//		double Fe = PI * triele.J * triele.area * (mp_node[triele.n[i]].x + 3 * triele.ydot) / 6;
	//		F[triele.n[i]] += Fe;
	//	}
	//}
	//考虑第一类边界条件的装配
	std::vector<std::vector<int>> locs(2, std::vector<int>(9 * (size_t)(m_num_triele)));
	std::vector<double> vals(9 * (size_t)(m_num_triele));
	std::vector<double> F(num_dof);
	int pos = 0;
	for (int i_tri = 0; i_tri < m_num_triele; ++i_tri) {
		CTriElement triele = mp_triele[i_tri];
		for (int i = 0; i < 3; ++i) {
			int n1 = triele.n[i];
			for (int j = 0; j < 3; ++j) {
				int n2 = triele.n[j];
				if (triele.material->getLinearFlag() == true) {
					if (mp_node[n1].bdr != 1 && mp_node[n2].bdr != 1) {
						double mu = triele.material->getMu();	//存在mu=0的情况
						double Se = triele.C[i][j] / mu;
						locs[0][pos] = node_pos[n1];
						locs[1][pos] = node_pos[n2];
						vals[pos] = Se;
						++pos;
					}
				}
			}
			if (mp_node[n1].bdr != 1) {
				double Fe = PI * triele.J * triele.area * (mp_node[triele.n[i]].x + 3 * triele.ydot) / 6;
				F[node_pos[n1]] += Fe;
			}
		}
	}

	locs[0].resize(pos);
	locs[1].resize(pos);
	cout << "pos: " << pos << endl;


	////求解器测试
	//vector<vector<int>> locs1(2,vector<int>(12));
	//vector<double> vals1(12);
	//vector<double> F1(5);
	//double s = 19, u = 21, p = 16, e = 5, r = 18, l = 12;
	//locs1[0][0] = 0; locs1[1][0] = 0; vals1[0] = s;
	//locs1[0][1] = 0; locs1[1][1] = 2; vals1[1] = u;
	//locs1[0][2] = 0; locs1[1][2] = 3; vals1[2] = u;
	//locs1[0][3] = 1; locs1[1][3] = 0; vals1[3] = l;
	//locs1[0][4] = 1; locs1[1][4] = 1; vals1[4] = u;
	//locs1[0][5] = 2; locs1[1][5] = 1; vals1[5] = l;
	//locs1[0][6] = 2; locs1[1][6] = 2; vals1[6] = p;
	//locs1[0][7] = 3; locs1[1][7] = 3; vals1[7] = e;
	//locs1[0][8] = 3; locs1[1][8] = 4; vals1[8] = u;
	//locs1[0][9] = 4; locs1[1][9] = 0; vals1[9] = l;
	//locs1[0][10] = 4; locs1[1][10] = 1; vals1[10] = l;
	//locs1[0][11] = 4; locs1[1][11] = 4; vals1[11] = r;
	//for (int i = 0; i < 5; ++i) F1[i] = 1.0;
	//double* res = matsolver->solveMatrix(locs1, vals1, F1, 12, 5);
	//for (int i = 0; i < 5; ++i) cout << res[i] << endl;

	//求解
	double* res1 = matsolver->solveMatrix(locs, vals, F, pos, num_dof);
	A.resize(m_num_nodes);
	for (int i = 0; i < num_dof; ++i) {
		int index = node_reorder[i];
		A[index] = res1[i];
	}

}
