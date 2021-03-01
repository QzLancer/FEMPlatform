#include "FEM2DSolver.h"

void FEM2DSolver::makeTrangle(int index)
{
	int k, m, n;
	double q0, q1, q2, r0, r1, r2, area;
	k = mp_triele[index].n[0];
	m = mp_triele[index].n[1];
	n = mp_triele[index].n[2];

	q0 = mp_node[m].y - mp_node[n].y;
	mp_triele[index].Q[0] = q0;
	q1 = mp_node[n].y - mp_node[k].y;
	mp_triele[index].Q[1] = q1;
	q2 = mp_node[k].y - mp_node[m].y;
	mp_triele[index].Q[2] = q2;

	r0 = mp_node[n].x - mp_node[m].x;
	mp_triele[index].R[0] = r0;
	r1 = mp_node[k].x - mp_node[n].x;
	mp_triele[index].R[1] = r1;
	r2 = mp_node[m].x - mp_node[k].x;
	mp_triele[index].R[2] = r2;

	area = 0.5 * std::abs(q1 * r2 - r1 * q2);
	mp_triele[index].area = area;

	mp_triele[index].rc = (mp_node[k].x +
		mp_node[m].x +
		mp_node[n].x) / 3;
	mp_triele[index].zc = (mp_node[k].y +
		mp_node[m].y +
		mp_node[n].y) / 3;

	int flag = 0;
	for (int f = 0; f < 3; f++) {
		if (mp_node[mp_triele[index].n[f]].x < 1e-7) {
			flag++;
		}
	}

	//�������������İ뾶
	if (flag == 2) {
		mp_triele[index].ydot = mp_triele[index].rc;
	}
	else {
		mp_triele[index].ydot = 1 / (mp_node[k].x + mp_node[m].x);
		mp_triele[index].ydot += 1 / (mp_node[k].x + mp_node[n].x);
		mp_triele[index].ydot += 1 / (mp_node[m].x + mp_node[n].x);
		mp_triele[index].ydot = 1.5 / mp_triele[index].ydot;
	}

	//����һ����������ԳƵ�Ԫϵ������
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			mp_triele[index].C[i][j] = (PI * mp_triele[index].ydot * 
				(mp_triele[index].R[i] * mp_triele[index].R[j] + mp_triele[index].Q[i] * mp_triele[index].Q[j])) / (2 * mp_triele[index].area);
		}
	}
}

void FEM2DSolver::processBoundaryCondition()
{
	//�����boundarymapֻ���˵�һ��߽�������ֻ�����ߵ�Ԫ��
	//1.��ȡ������ȫ����2.��ȡ���Ӧ�ĵ�Ԫ��ţ�3.��ȡ�ڵ�
	//4.ȥ�� 5.��node���ñ߽����� 6.���߽������ĩβ
	num_freenodes = m_num_nodes;
	std::vector<int> boundarynodes;
	//ԭ�а汾�ı߽����
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

	num_freenodes -= boundarynodes.size();
	for (auto a : boundarynodes) {
		mp_node[a].bdr = 1;
	}
	cout << "num_bdr: " << boundarynodes.size() << ", num_freenodes: " << num_freenodes << endl;
	//���߽������ĩβ

	node_reorder.resize(m_num_nodes);
	node_pos.resize(m_num_nodes);
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
}

void FEM2DSolver::processMaterial()
{
	for (int i_tri = 0; i_tri < m_num_triele; ++i_tri) {
		int domain = mp_triele[i_tri].domain;
		mp_triele[i_tri].material = materialmap[domain];
	}
}

void FEM2DSolver::processLoad()
{
	//���õ�Ԫ����
	for (int i_tri = 0; i_tri < m_num_triele; ++i_tri) {
		int domain = mp_triele[i_tri].domain;
		if (loadmap.find(domain) != loadmap.end()) {
			mp_triele[i_tri].J = loadmap[domain];
		}
	}
}

void FEM2DSolver::updateB()
{
	for( int i_tri = 0; i_tri < m_num_triele; ++i_tri ) {
		double bx = 0, by = 0;
		for (int i = 0; i < 3; ++i) {
			int n = mp_triele[i_tri].n[i];
			bx += mp_triele[i_tri].R[i] * A[n];
			by += mp_triele[i_tri].Q[i] * A[n];
		}
		bx = bx / 2 / mp_triele[i_tri].area;
		Bx[i_tri] = bx;
		by = -by / 2 / mp_triele[i_tri].area;
		By[i_tri] = by;
		B[i_tri] = sqrt(bx * bx + by * by);
	}
}

