#include "FEM2DStaticSolver.h"

void FEM2DStaticSolver::solve()
{
	////测试数据是否正确传输到求解器中
	//cout << "Number of Nodes: " << m_num_nodes << endl;
	//cout << "Number of VtxElement: " << m_num_vtxele << endl;
	//cout << "Number of EdgElement: " << m_num_edgele << endl;
	//cout << "Number of TriElement: " << m_num_triele << endl;
	//for (int i = 0; i < m_num_nodes; ++i) {
	//	cout << "Nodes: " << mp_node[i].x << " , " << mp_node[i].y << endl;
	//}
	//for (int i = 0; i < m_num_vtxele; ++i) {
	//	cout << "VtxElement: " << mp_vtxele[i].n << ", Domain: " << mp_vtxele[i].domain << endl;
	//}
	//for (int i = 0; i < m_num_edgele; ++i) {
	//	cout << "EdgElement: " << mp_edgele[i].n[0] << ", " << mp_edgele[i].n[1] << ", Domain: " << mp_edgele[i].domain << endl;
	//}
	//for (int i = 0; i < m_num_triele; ++i) {
	//	cout << "TriElement: " << mp_triele[i].n[0] << ", " << mp_triele[i].n[1] << ", " << mp_triele[i].n[2] << ", Domain: " << mp_triele[i].domain << endl;
	//}
	//cout << "Number of Materialmap: " << materialmap.size() << endl;
	//cout << "Number of Loadmap: " << loadmap.size() << endl;
	//cout << "Number of Boundarymap : " << boundarymap.size() << endl;

	//计算三角形单元几何部分
	for (int i = 0; i < m_num_triele; ++i) {
		makeTrangle(i);
	}

	//计算自由度
	//这里的boundarymap只用了第一类边界条件，只考虑线单元域
	//1.提取线所在全部域，2.提取域对应的单元编号，3.提取节点，4.去重
	int num_dof = m_num_nodes;
	std::vector<int> boundarynodes;
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
	num_dof -= boundarynodes.size();

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

	//定义求解过程中的变量
	std::vector<std::vector<int>> locs(2, std::vector<int>(9*(size_t)(m_num_triele)));
	std::vector<double> vals(9*(size_t)(m_num_triele));
	std::vector<double> F(m_num_nodes);
	int pos = 0;

	//线性三角形一阶单元装配
	for (int i_tri = 0; i_tri < m_num_triele; ++i_tri) {
		CTriElement triele = mp_triele[i_tri];
		for (int i = 0; i < 3; ++i) {
			for (int j = 0; j < 3; ++j) {
				if (triele.material->getLinearFlag() == true) {
					double mu = triele.material->getMu();
					double Se = triele.C[i][j] / mu;
					locs[0][pos] = triele.n[i];
					locs[1][pos] = triele.n[j];
					vals[pos] = Se;
					++pos;
				}
			}
			double fe = PI * triele.J * triele.area * (mp_node[triele.n[i]].x + 3 * triele.ydot) / 6;
			F[triele.n[i]] += fe;
		}
	}

	cout << pos << endl;
	//非线性迭代
	for (int iter = 0; iter < maxitersteps; ++iter) {
		//非线性单元装配过程
		int pos = 0;
		

	}
}
